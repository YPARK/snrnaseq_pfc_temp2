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
CIS.DIST <- 1e6 # cis distance

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

    ## quantlie normalization
    yy <- apply(yy, 2, .qnorm) %>% as.matrix
    yy[, n.valid < 10] <- 0

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

opts <- list(do.hyper=FALSE, pi=0, gammax=1e3, vbiter=5000,
             print.interv=100, out.residual=FALSE, tol = 1e-6)

y <- .data$y %>% as.matrix
x <- .data$x %>% as.matrix

if(max(apply(y, 2, sd, na.rm=TRUE)) < 1e-4) {
    write_tsv(data.frame(), OUT.HDR %&% "_stat.bed.gz")
    write_tsv(data.frame(), OUT.HDR %&% "_celltype.bed.gz")
    log.msg("Too small variance")
    q()
}

if(DO.PERMUTE){
    y <- apply(y, 2, sample) %>% as.matrix
    x <- x[sample(nrow(x)), , drop = FALSE]
}

.fqtl <- fqtl::fit.fqtl(y=y, x.mean = x, factored = TRUE,
                        k=min(ncol(y), ncol(x)),
                        options = opts)

.snp.info <- .data$snp.info %>%
    mutate(x.col = 1:n())

.gene.info <- .data$gene.info %>%
    mutate(y.col = 1:n())

.left.dt <- melt.fqtl.effect(.fqtl$mean.left) %>%
    dplyr::rename(x.col = .row, k.col = .col) %>%
    left_join(.snp.info) %>%
    as.data.table

.right.dt <- melt.fqtl.effect(.fqtl$mean.right) %>%
    dplyr::rename(y.col = .row, k.col = .col) %>%
    (function(x) left_join(.gene.info, x)) %>%
    as.data.table

out.stat <-
    left_join(.left.dt, .right.dt,
              suffix = c(".geno", ".celltype"),
              by = c("k.col")) %>%
    mutate(start = snp.loc - 1) %>%
    mutate(stop = snp.loc) %>%
    arrange(`#chr`, `stop`, `celltype`) %>%
    select(`#chr`, `start`, `stop`, starts_with("plink"),
           ensembl_gene_id, hgnc_symbol,
           celltype,
           starts_with("theta"),
           starts_with("lodds")) %>%
    as.data.table

.bed.write(out.stat, OUT.HDR %&% "_stat.bed.gz")

.bed.write(.right.dt, OUT.HDR %&% "_celltype.bed.gz")

log.msg("done")
