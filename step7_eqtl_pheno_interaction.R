#!/usr/bin/env Rscript

## GENE.INFO.FILE <- "result/step1/gene.info.gz"
## GENE           <- 67
## DATA.DIR       <- "result/step6/eqtl/data/"
## PLINK.DIR      <- "data/rosmap_geno/"
## TEMP.DIR       <- "temp/"
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
      OUT.HDR %&% "_geno.bed.gz",
      OUT.HDR %&% "_pheno.bed.gz")

if(all(file.exists(out.files))) {
    log.msg("files exist --> done")
    q()
}

################################################################
CIS.DIST <- 1e6 # cis distance
PHENO.FILE <- "data/rosmap_phenotype.csv.gz"

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
                           pheno.file = PHENO.FILE,
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

    .scale <- function(x, p.min = .1, p.max = .9){
        .d <- p.max - p.min
        ret <- (x - min(x,na.rm=TRUE))/diff(range(x,na.rm=TRUE))
        p.min + ret * .d
    }

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

    .pheno <- .fread(pheno.file, na.strings = "-9")

    .pheno[, age.death := exp(.qnorm(age_death))]
    .pheno[, np.sqrt := exp(.qnorm(np_sqrt))]
    .pheno[, nft.sqrt := exp(.qnorm(nft_sqrt))]
    .pheno[, cogdx := exp(.qnorm(cogn_ep_random_slope))]
    .pheno[, apoe.e4 := .pheno$apoe4n]

    .pheno <- .pheno %>% 
        mutate(iid = as.character(projid)) %>%
        select(iid, pathoAD, msex, np.sqrt, nft.sqrt, age.death, apoe.e4, cogdx)

    pos.df <- pos.df[as.character(iid) %in% as.character(.pheno$iid)]

    phi <-
        left_join(pos.df[, .(iid)], .pheno) %>% 
        select(-iid) %>%
        as.matrix

    xx <- scale(.plink$BED[pos.df$x.pos, , drop = FALSE])
    yy <- t(as.matrix(.yy[, pos.df$y.pos, drop = FALSE]))

    ## Deal with an empty gene
    n.valid <- apply(is.finite(yy), 2, sum)

    ## convert back to pseudo-counting data
    yy <- scale(yy)
    yy[yy > 3] <- 3
    yy[yy < -3] <- 3
    yy <- exp(yy) %>% as.matrix
    yy[, n.valid < 10] <- 0

    list(x = xx,
         y = yy,
         phi = phi,
         phenotypes = colnames(phi),
         snp.info = .plink$BIM,
         gene.info = expr.dt[, 1:6],
         pos = pos.df)
}

.data <- read.gene.data(GENE)

if(is.null(.data)) {
    write_tsv(data.frame(), OUT.HDR %&% "_stat.bed.gz")
    write_tsv(data.frame(), OUT.HDR %&% "_geno.bed.gz")
    write_tsv(data.frame(), OUT.HDR %&% "_pheno.bed.gz")
    log.msg("nothing to do")
    q()
}

opts <- list(do.hyper=FALSE,
             pi=-0,
             tau=-4,
             svd.init=TRUE,
             jitter=1e-4,
             gammax=1e3,
             vbiter=7500,
             print.interv=100,
             out.residual=FALSE,
             rate = 1e-2,
             tol = 1e-6)

y <- .data$y %>% as.matrix
x <- .data$x %>% as.matrix
phi <- .data$phi %>% as.matrix

if(max(apply(y, 2, sd, na.rm=TRUE)) < 1e-4) {
    write_tsv(data.frame(), OUT.HDR %&% "_stat.bed.gz")
    write_tsv(data.frame(), OUT.HDR %&% "_geno.bed.gz")
    write_tsv(data.frame(), OUT.HDR %&% "_pheno.bed.gz")
    log.msg("Too small variance")
    q()
}

## Interaction QTL analysis

if(DO.PERMUTE){

    y.perm <- apply(y, 2, sample) %>% as.matrix
    phi.perm <- apply(phi, 2, sample) %>% as.matrix
    log.msg("Permuted samples in the phenotypes")

    .fqtl <- fqtl::fit.fqtl(y=y.perm,
                            x.mean = x,
                            factored = TRUE,
                            weight.nk = phi.perm,
                            model = "nb",
                            c.mean = cbind(x, 1),
                            options = opts)

} else {

    .fqtl <- fqtl::fit.fqtl(y=y, x.mean=x,
                            factored = TRUE,
                            weight.nk = phi,
                            model = "nb",
                            c.mean = cbind(x, 1),
                            options = opts)

}

.snp.info <- .data$snp.info %>%
    mutate(x.col = 1:n())

.gene.info <- .data$gene.info %>%
    mutate(y.col = 1:n())

.pheno.info <- data.table(pheno = .data$phenotypes) %>% 
    mutate(k.col = 1:n())

.left.dt <- melt.fqtl.effect(.fqtl$mean.left) %>%
    mutate(theta.sd = sqrt(theta.var)) %>% 
    dplyr::select(-theta.var) %>% 
    dplyr::rename(x.col = .row, k.col = .col) %>%
    (function(x) left_join(.snp.info, x)) %>% 
    dplyr::rename(`#chr` = `chr`) %>% 
    left_join(.pheno.info) %>%
    mutate(start = snp.loc - 1) %>% 
    mutate(stop = snp.loc) %>% 
    arrange(`#chr`, `stop`, `pheno`) %>% 
    dplyr::select(`#chr`, `start`, `stop`,
                  starts_with("plink"), `pheno`,
                  pheno, theta, theta.sd, lodds) %>% 
    as.data.table

.right.dt <- melt.fqtl.effect(.fqtl$mean.right) %>% 
    mutate(theta.sd = sqrt(theta.var)) %>% 
    dplyr::select(-theta.var) %>% 
    dplyr::rename(y.col = .row, k.col = .col) %>% 
    (function(x) left_join(.gene.info, x)) %>%
    left_join(.pheno.info) %>%
    dplyr::select(-ends_with(".col")) %>% 
    as.data.table

.cov.dt <-
    melt.fqtl.effect(.fqtl$mean.cov) %>% 
    rename(x.col = .row, y.col = .col) %>%
    mutate(theta.sd = sqrt(theta.var)) %>% 
    left_join(.snp.info) %>% 
    na.omit %>% 
    left_join(.gene.info) %>% 
    mutate(start = `snp.loc` - 1) %>% 
    mutate(stop = `snp.loc`) %>% 
    arrange(`#chr`, `stop`) %>% 
    select(`#chr`, `start`, `stop`, starts_with("plink"), ensembl_gene_id,
           hgnc_symbol, celltype, theta, theta.sd, lodds) %>% 
    as.data.table

.bed.write(.left.dt, OUT.HDR %&% "_stat.bed.gz")
.bed.write(.cov.dt, OUT.HDR %&% "_geno.bed.gz")
.bed.write(.right.dt, OUT.HDR %&% "_pheno.bed.gz")

log.msg("done")
