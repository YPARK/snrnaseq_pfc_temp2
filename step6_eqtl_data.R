#!/usr/bin/env Rscript

## data.file <- "result/step5/pseudobulk/In-SST.ln_mu.gz"
## col.file <- "result/step5/pseudobulk/In-SST.mu_cols.gz"
## row.file <- "result/step5/filtered/In-SST.rows.gz"
## info.file <- "result/step1/gene.info.gz"
## n.pc <- 10
## geno.dir <- "data/rosmap_geno/"
## TEMP.DIR <- "temp"
## gene.idx <- "1,10"

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) != 9) {
    q()
}

data.file <- argv[1]
col.file <- argv[2]
row.file <- argv[3]
info.file <- argv[4]
n.pc <- as.integer(argv[5])
geno.dir <- argv[6]
TEMP.DIR <- argv[7]
gene.idx <- argv[8]
out.file <- argv[9]

if(file.exists(out.file)) {
    q()
}

################################################################

library(tidyverse)
library(data.table)
source("Util-geno.R")

`%&%` <- function(a,b) paste0(a, b)

.fread <- function(file.name, ...) {
    .s <- str_length(file.name)
    if(str_sub(file.name, .s - 2, .s) == ".gz") {
        return(fread("gzip -cd " %&% file.name, ...))
    }
    return(fread(file.name, ...))
}

gene.idx <- unlist(str_split(gene.idx, pattern=",")) %>% as.integer

celltype <-
    .fread(col.file, header=FALSE, nrow=1) %>%
    unlist %>%
    str_split("_") %>%
    (function(x) x[[1]][2])

info.tab <- .fread(info.file, header=TRUE, sep = "\t", fill=TRUE)

read.data <- function(dat.file, col.file, rows){
    .cols <- .fread(col.file, header=FALSE, col.names="sample")
    .cols[, c("iid", "celltype") := tstrsplit(sample, split="_")]
    ret <- .fread(dat.file, header=FALSE, col.names=.cols$iid)
    rownames(ret) <- rows$gene
    return(ret)
}

.rows <- .fread(row.file, header=FALSE, col.names="gene")
.data <- read.data(data.file, col.file, rows = .rows)

genes <- .rows %>%
    mutate(j = 1:n()) %>%
    left_join(info.tab) %>%
    mutate(chr = as.integer(chr)) %>%
    na.omit %>%
    arrange(chr, transcript_start.hg19) %>%
    as.data.table

Y <- t(as.matrix(.data[genes$j, ]))

ngene <- nrow(genes)

if(!any(gene.idx %in% 1:ngene)) {
    log.msg("Nothing to do")
    q()
}

.remove.outlier <- function(y) { # deal with zero count spots
    return(log(1 + exp(y)))
}

Y <- apply(Y, 2, .remove.outlier) %>% as.matrix

ind.names <- colnames(.data)

stopifnot(all(ind.names == rownames(Y)))

match.ind <- function(.plink, ind.names) {
    x.pos.df <- .plink$FAM %>%
        select(iid) %>%
        mutate(iid = str_remove_all(iid, "MAP")) %>%
        mutate(iid = str_remove_all(iid, "ROS")) %>%
        mutate(x.pos = 1:n())

    y.pos.df <- tibble(iid = ind.names) %>%
        mutate(iid = str_remove_all(iid, "MAP")) %>%
        mutate(iid = str_remove_all(iid, "ROS")) %>%
        mutate(y.pos = 1:n())

    pos.df <- left_join(y.pos.df, x.pos.df, by = "iid") %>%
        na.omit
}

adjust.confounder <- function(g, max.rank=30, ntop=100, cis.dist=1e6, p.cutoff=5e-4) {

    .symb <- unlist(genes[g, .(hgnc_symbol)])
    .chr <- unlist(genes[g, .(chr)])

    temp.dir <- TEMP.DIR %&% "/" %&% celltype %&% "/gene_g" %&% g %&% "_" %&% .symb
    .mkdir(temp.dir)

    tss <- genes[g, .(transcript_start.hg19)] %>% unlist
    tes <- genes[g, .(transcript_end.hg19)] %>% unlist

    .plink <- subset.plink(plink.hdr = geno.dir %&% "/chr" %&% .chr,
                           chr = .chr,
                           plink.lb = max(tss - cis.dist, 0),
                           plink.ub = tes + cis.dist,
                           temp.dir = temp.dir)

    if(is.null(.plink)) return(NULL)

    unlink(temp.dir, recursive = TRUE)

    pos.df <- match.ind(.plink, ind.names)

    gg <-
        .plink$BED[pos.df$x.pos, , drop = FALSE] %>%
        scale

    yy <- Y[pos.df$y.pos, g, drop = FALSE]
    y <- Y[, g, drop = FALSE]

    ## 1. (G) Read geno type matrix
    ## 2. (Y0) Find most correlated genes
    ## 3. Test if a control gene can be used as a confounder variable, or not

    y0.idx <- genes$chr != .chr
    Y0 <- Y[, y0.idx, drop = FALSE]

    y0.stat <-
        zqtl::calc.qtl.stat(Y0, y) %>%
        arrange(p.val) %>%
        head(ntop) %>%
        filter(p.val < p.cutoff)

    if(nrow(y0.stat) < 1) {
        return(y)
    }

    valid <- c()

    for(j in y0.stat$x.col) {

        y0 <- Y0[pos.df$y.pos, j, drop = FALSE]

        ## (case 1) If "Y0" is a true confounder,
        ## i.e., G -> Y and Y0 -> Y
        ## this will artificially introduce a collider bias
        ## (case 2) If "Y0" is at the downstream,
        ## i.e., G -> Y -> Y0
        ## this will regress out all the genetic effects &
        ## there will be no association with genetic variants

        r0 <- lm(y0 ~ yy) %>%
            residuals %>%
            matrix(ncol=1)

        rx.stat <- zqtl::calc.qtl.stat(gg, r0) %>%
            filter(p.val < p.cutoff)

        if(nrow(rx.stat) > 0) {
            log.msg("Found %d, %.2e", j, min(rx.stat$p.val))
            valid <- c(valid, j)
        }
    }

    log.msg("Found %d control genes", length(valid))

    if(length(valid) < 1) {
        return(y)
    }

    y0.valid <- scale(Y0[, valid, drop = FALSE])
    y0.valid[is.na(y0.valid)] <- 0

    if(length(valid) > max.rank) {
        .svd <- rsvd::rsvd(y0.valid, k=max.rank)
        y0.valid <- .svd$u
    }

    rr <- lm(y ~ y0.valid) %>% residuals

    return(rr)
}

residual.mat <- matrix(NA, nrow=length(gene.idx), ncol=length(ind.names))

for(gi in 1:length(gene.idx)) {
    g <- gene.idx[gi]
    .r <- adjust.confounder(g, max.rank = n.pc, ntop = 50)
    if(!is.null(.r))
        residual.mat[gi, ] <- .r
}

out.bed <- genes[gene.idx] %>%
    rename(transcript_start = transcript_start.hg19,
           transcript_end = transcript_end.hg19,
           `#chr` = chr) %>%
    select(`#chr`,
           transcript_start,
           transcript_end,
           ensembl_gene_id,
           hgnc_symbol) %>%
    mutate(celltype)

out.dt <- data.table(residual.mat)
colnames(out.dt) <- ind.names

out.dt <- cbind(out.bed, out.dt) %>%
    as_tibble %>%
    arrange(`#chr`, transcript_start) %>%
    as.data.table

log.msg("out.dt ready")

out.tsv.file <- str_remove(out.file, ".gz$")

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
log.msg("Done")
