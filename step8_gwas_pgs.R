#!/usr/bin/env Rscript

## LD.FILE <- "LD.info.txt"
## LD.IDX <- 1
## GWAS.FILE <- "data/GWAS/pgc_scz.bed.gz"
## GENO.DIR <- "data/rosmap_geno/"
## TEMP.DIR <- "temp"
## OUT.FILE <- "out.bed.gz"

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) != 6) q()

LD.FILE   <- argv[1]
LD.IDX    <- as.integer(argv[2])
GWAS.FILE <- argv[3]
GENO.DIR  <- argv[4]
TEMP.DIR  <- argv[5]
OUT.FILE  <- argv[6]

source("Util-geno.R")
library(data.table)
library(tidyverse)

if(file.exists(OUT.FILE)) {
    log.msg("File %s exists", OUT.FILE)
    q()
}

.bed.write <- function(out.dt, out.file) {
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
}

.mkdir(TEMP.DIR)

ld.tab <- fread(LD.FILE) %>%
    mutate(`#CHR`=str_remove(chr, "chr")) %>%
    mutate(`#CHR`=as.integer(`#CHR`)) %>%
    mutate(qq = `#CHR` %&% ":" %&% `start` %&% "-" %&% `stop`) %>%
    mutate(ld = 1:n()) %>%
    as.data.frame

svd.knockoff <- function(.svd){
    uu <- .svd$U
    .offset <- sample(ncol(uu), 1)
    .col <- seq(.offset, ncol(uu) + .offset - 1) %% ncol(uu) + 1
    uu.ko <- uu[, .col, drop = FALSE]
    ret <- .svd
    ret$U <- uu.ko
    return(ret)
}

fit.pgs.ld <- function(ld.idx) {

    log.msg("\n\nFitting LD %d\n\n", ld.idx)

    qq <- ld.tab[ld.idx, "qq"]
    .chr <- ld.tab[ld.idx, "#CHR"]
    .lb <- ld.tab[ld.idx, "start"]
    .ub <- ld.tab[ld.idx, "stop"]

    .gwas.tab <- rbind(fread("tabix -h " %&% GWAS.FILE %&% " " %&% qq),
                       fread("tabix -h " %&% GWAS.FILE %&% " chr" %&% qq))

    .plink <- subset.plink(GENO.DIR %&% "/chr" %&% .chr, .chr, .lb, .ub, TEMP.DIR)

    if(nrow(.gwas.tab) < 1) { return(data.table()) }

    if(is.null(.plink)) { return(data.table()) }

    if("lodds" %in% names(.gwas.tab) && !("beta" %in% names(.gwas.tab))){
        .gwas.tab <- .gwas.tab %>%
            rename(beta = lodds) %>%
            as.data.table
    }

    .gwas.tab <- .gwas.tab %>%
        mutate(snp.loc = `stop`) %>%
        as.data.table

    valid.snps <- .plink$BIM %>%
        select(-missing, -rs) %>%
        mutate(x.pos = 1:n()) %>%
        (function(x) left_join(.gwas.tab, x, by = 'snp.loc')) %>%
        filter(!is.na(x.pos)) %>%
        filter((plink.a1 == a1 & plink.a2 == a2) | (plink.a2 == a1 & plink.a1 == a2)) %>%
        mutate(gwas.flip = if_else(plink.a1 == a1, 1.0, -1.0)) %>%
        mutate(beta = gwas.flip * beta) %>%
        select(-gwas.flip) %>%
        as.data.table

    if(nrow(valid.snps) < 1) {
        log.msg("skip this LD block")
        return(data.table())
    }

    valid.snps <-
        valid.snps[order(abs(valid.snps$beta), decreasing=TRUE),
                   head(.SD, 1),
                   by = .(rs)] %>%
        arrange(snp.loc) %>%
        as.data.table

    X <- .plink$BED %c% valid.snps$x.pos %>% scale

    zz <- valid.snps %>%
        mutate(z = beta/se) %>%
        select(z) %>%
        as.matrix()

    .svd <- zqtl::take.ld.svd(X, eigen.tol = .01, eigen.reg = 0)
    zz <- zqtl::scale.zscore(zz, .svd$V.t, .svd$D)

    if(is.null(.svd)) return(data.table())

    .svd.ko <- svd.knockoff(.svd)

    y.obs <- pred.prs(.svd, zz, 0.01) %>% as.numeric
    y.ko <- pred.prs(.svd.ko, zz, 0.01) %>% as.numeric



    uu <- .svd$U

    pv <-
        sapply(1:ncol(uu), function(k) {
            y.up <- y.tot[uu[, k] > 0]
            y.down <- y.tot[uu[, k] <= 0]
            wilcox.test(y.down, y.up)$p.value
        })

    valid.comp <- which(p.adjust(pv) < pv.cutoff)

    if(length(valid.comp) < 1) {

        y.hat <- y.tot * 0

    } else {

        .svd.valid <- select.svd(.svd, valid.comp)

        y.hat <- pred.prs(.svd.valid, zz, 0) %>% as.numeric

    }

    .dt <- data.table(iid = .plink$FAM$iid, y.hat = y.hat, y.tot = y.tot, ld = ld.idx)

    return(.dt)
}

n.ld <- nrow(ld.tab)
out.dt <- fit.pgs.ld(LD.IDX, pv.cutoff = 0.01 / n.ld)

.mkdir(dirname(OUT.FILE))
if(nrow(out.dt) > 1) {
    out.dt <- ld.tab[LD.IDX, ] %>% left_join(out.dt) %>%
        select(`#CHR`, `start`, `stop`, ld, iid, y.hat, y.tot) %>% 
        as.data.table
    .bed.write(out.dt, OUT.FILE)
} else {
    write_tsv(out.dt, OUT.FILE)
}

unlink(TEMP.DIR)

log.msg("done")
