#!/usr/bin/env Rscript

## GENE.INFO.FILE <- "result/step1/gene.info.gz"
## GENE           <- 2
## DATA.DIR       <- "result/step6/eqtl/data/"
## PLINK.DIR      <- "data/rosmap_geno/"
## TEMP.DIR       <- "temp"
## OUT.HDR        <- "temp"
DO.PERMUTE <- FALSE

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 6) q()

GENE.INFO.FILE <- argv[1]
GENE           <- as.integer(argv[2])
DATA.DIR       <- argv[3]
PLINK.DIR      <- argv[4]
TEMP.DIR       <- argv[5]
OUT.HDR        <- argv[6]

if(length(argv) > 6){
    DO.PERMUTE <- as.logical(argv[7])
}

################################################################

source("Util-geno.R")
library(tidyverse)
library(data.table)

out.files <-
    c(OUT.HDR %&% "_stat.bed.gz",
      OUT.HDR %&% "_celltype.bed.gz")

if(all(file.exists(out.files))) {
    log.msg("files exist --> done")
    q()
}

################################################################
CIS.DIST <- 7e5 # cis distance

dir.create(dirname(OUT.HDR), recursive = TRUE, showWarnings = FALSE)
TEMP.DIR <- TEMP.DIR %&% OUT.HDR
dir.create(TEMP.DIR, recursive = TRUE, showWarnings = FALSE)

.fread <- function(file.name, ...) {
    .s <- str_length(file.name)
    if(str_sub(file.name, .s - 2, .s) == ".gz") {
        return(fread("gzip -cd " %&% file.name, ...))
    }
    return(fread(file.name, ...))
}

.bed.write <- function(out.dt, out.file) {
    out.tsv.file <- str_remove(out.file, ".gz$")

    if(file.exists(out.file)) {
        unlink(out.file)
    }

    if(file.exists(out.tsv.file)) {
        unlink(out.tsv.file)
    }

    fwrite(out.dt,
           file = out.tsv.file,
           sep = "\t",
           eol = "\n",
           row.names = FALSE,
           col.names = TRUE,
           na = "NA",
           quote = FALSE)

    Rsamtools::bgzip(out.tsv.file, dest=out.file)
    Rsamtools::indexTabix(out.file, format="bed")
    unlink(out.tsv.file)
}


GENE.INFO <- .fread(GENE.INFO.FILE, header=TRUE) %>%
    as.data.frame

read.gene.data <- function(g,
                           data.dir = DATA.DIR,
                           gene.info = GENE.INFO,
                           temp.dir = TEMP.DIR) {

    expr.files <- list.files(data.dir,
                             full.names = TRUE,
                             pattern="bed.gz.tbi") %>%
        str_remove_all(".tbi$")

    if(length(expr.files) < 1) return(NULL)

    .info <- gene.info[g, ]
    .chr <- .info["chr"]
    .tss <- .info["transcript_start.hg19"] %>% as.integer
    .tes <- .info["transcript_end.hg19"] %>% as.integer
    .qq <- .chr %&% ":" %&% .tss %&% "-" %&% .tes
    .ensg <- .info["ensembl_gene_id"] %>% as.character

    if(!(.chr %in% as.character(1:22))) return(NULL) # sorry only autosomal

    .read.tab <- function(.file){
        .cmd <- "tabix -h " %&% .file %&% " " %&% .qq
        .ret <-
            suppressMessages(fread(.cmd, header=TRUE)) %>%
            mutate(`#chr` = as.character(`#chr`)) %>%
            mutate(transcript_start = as.integer(transcript_start)) %>%
            mutate(transcript_end = as.integer(transcript_end)) %>%
            mutate(ensembl_gene_id = as.character(ensembl_gene_id)) %>%
            mutate(hgnc_symbol = as.character(hgnc_symbol)) %>%
            mutate(celltype = as.character(celltype)) %>%
            as.data.table

        .id.vars <- colnames(.ret)[1:6]
        .ret <-
            melt(.ret, id.vars = .id.vars) %>%
            mutate(variable = as.character(variable)) %>%
            as.data.table

        .id.list <- c(.id.vars, "variable")
        .ret <- .ret[, .(value = mean(value)), by = .id.list]
    }

    expr.dt <-
        lapply(expr.files, .read.tab) %>%
        bind_rows %>%
        filter(ensembl_gene_id == .ensg) %>%
        as_tibble %>%
        spread(key = variable, value = value) %>%
        as.data.frame

    if(is.null(expr.dt) || nrow(expr.dt) < 1) return(NULL)

    .yy <- expr.dt[, -(1:6), drop = FALSE]

    .chr.num <- str_remove(.chr, "chr")
    .lb <- pmax(.tss - CIS.DIST, 0)
    .ub <- (.tes + CIS.DIST)

    .plink <- subset.plink(PLINK.DIR %&% "/chr" %&% .chr.num,
                           .chr.num, .lb, .ub, TEMP.DIR)

    if(is.null(.plink)) return(NULL)

    x.pos.df <-
        .plink$FAM %>%
        select(iid) %>%
        mutate(iid = str_remove_all(iid, "MAP")) %>%
        mutate(iid = str_remove_all(iid, "ROS")) %>%
        mutate(x.pos = 1:n()) %>%
        group_by(iid) %>%
        slice(which.min(x.pos)) %>%
        ungroup

    y.pos.df <- tibble(iid = names(.yy)) %>%
        mutate(iid = str_remove_all(iid, "MAP")) %>%
        mutate(iid = str_remove_all(iid, "ROS")) %>%
        mutate(y.pos = 1:n())

    pos.df <-
        left_join(y.pos.df, x.pos.df, by = "iid") %>%
        na.omit %>%
        as.data.table

    .qnorm <- function(x) {
        .loc <- is.finite(x)
        x.safe <- x[.loc]
        nn <- length(x.safe)
        qq <- qnorm((1:nn)/(nn + 1))
        x.safe[order(x.safe)] <- qq
        ret <- x
        ret[.loc] <- x.safe
        return(ret)
    }

    xx <- .plink$BED[pos.df$x.pos, , drop = FALSE]
    yy <- t(as.matrix(.yy[, pos.df$y.pos, drop = FALSE]))

    ## Deal with an empty gene
    n.valid <- apply(is.finite(yy), 2, sum)

    ## convert back to pseudo-counting data
    yy <- scale(yy) %>% exp() %>% as.matrix
    yy[, n.valid < 10] <- NA

    list(x = xx,
         y = yy,
         snp.info = .plink$BIM,
         gene.info = expr.dt[, 1:6],
         pos = pos.df)
}

.data <- read.gene.data(GENE)

if(is.null(.data)) {
    write_tsv(data.frame(), OUT.HDR %&% "_stat.bed.gz")
    write_tsv(data.frame(), OUT.HDR %&% "_celltype.bed.gz")
    log.msg("nothing to do")
    q()
}

opts <- list(do.hyper = TRUE,
             svd.init = TRUE,
             jitter = 1e-4,
             pi.ub = 0,
             pi.lb = log(.1/.9),
             gammax = 1e3,
             vbiter = 3500,
             print.interv = 100,
             out.residual = FALSE,
             rate = 1e-2,
             decay = 0,
             tol = 1e-6)

y <- .data$y %>% as.matrix

svd.knockoff <- function(x){
    .svd <- zqtl::take.ld.svd(x, eigen.tol = 0, eigen.reg = 0)    

    uu <- .svd$U
    .offset <- sample(ncol(uu), 1)
    .col <- seq(.offset, ncol(uu) + .offset - 1) %% ncol(uu) + 1
    uu.ko <- uu[, .col, drop = FALSE]

    xx <- (sweep(uu, 2, .svd$D, `*`) %*% .svd$V.t) %>%
        scale

    xx.ko <- (sweep(uu.ko, 2, .svd$D, `*`) %*% .svd$V.t) %>%
        scale

    list(x = xx, ko = xx.ko)
}

if(max(apply(y, 2, sd, na.rm=TRUE), na.rm=TRUE) < 1e-4) {
    write_tsv(data.frame(), OUT.HDR %&% "_stat.bed.gz")
    write_tsv(data.frame(), OUT.HDR %&% "_celltype.bed.gz")
    log.msg("Too small variance")
    q()
}

.ko <- svd.knockoff(.data$x)
x <- .ko$x
x.ko <- .ko$ko

if(DO.PERMUTE){
    y <- apply(y, 2, sample) %>% as.matrix
}

.fqtl <- fqtl::fit.fqtl(y=y,
                        x.mean = cbind(x, x.ko),
                        factored = TRUE,
                        k=min(ncol(y), 2 * ncol(x)),
                        model = "nb",
                        options = opts)

.snp.info <- .data$snp.info %>%
    mutate(x.col = 1:n()) %>%
    mutate(ko = 0)

ptot <- nrow(.snp.info)

.snp.info.ko <- .snp.info %>%
    mutate(x.col = x.col + ptot) %>%
    mutate(ko = 1)

.snp.info <- rbind(.snp.info, .snp.info.ko)

.gene.info <- .data$gene.info %>%
    mutate(y.col = 1:n()) %>%
    mutate(`#chr` = as.integer(`#chr`))

.ensg <- unique(.gene.info$ensembl_gene_id)
.hgnc <- unique(.gene.info$hgnc_symbol)

.dt <- melt.fqtl.effect(.fqtl$mean.left) %>% 
    dplyr::rename(x.col = .row, k.col = .col) %>%
    (function(x) left_join(.snp.info, x)) %>%
    as.data.table

.dt[order(.dt$lodds, decreasing = TRUE),
    efdr := pmin(cumsum(ko) / seq(1, .N), 1),
    by = .(k.col)]

.left.dt <- .dt %>%
    mutate(theta.sd = sqrt(theta.var)) %>% 
    dplyr::select(-theta.var) %>% 
    dplyr::rename(`#chr` = `chr`) %>% 
    mutate(`#chr` = as.integer(`#chr`)) %>% 
    mutate(`start` = snp.loc - 1) %>% 
    mutate(`stop` = snp.loc) %>% 
    mutate(ensembl_gene_id = .ensg, hgnc_symbol = .hgnc) %>% 
    arrange(`#chr`, `stop`, `k.col`) %>% 
    dplyr::rename(k = `k.col`) %>% 
    dplyr::select(`#chr`, `start`, `stop`,
                  starts_with("plink"), `k`,
                  `ensembl_gene_id`, `hgnc_symbol`, `ko`,
                  theta, theta.sd, lodds, efdr) %>% 
    as.data.table

.right.dt <- melt.fqtl.effect(.fqtl$mean.right) %>%
    mutate(theta.sd = sqrt(theta.var)) %>% 
    dplyr::select(-theta.var) %>% 
    dplyr::rename(y.col = .row, k.col = .col) %>%
    (function(x) left_join(.gene.info, x)) %>%
    dplyr::rename(k = `k.col`) %>% 
    dplyr::select(- y.col) %>% 
    as.data.table

.bed.write(.left.dt, OUT.HDR %&% "_stat.bed.gz")

.bed.write(.right.dt, OUT.HDR %&% "_celltype.bed.gz")

log.msg("done")
