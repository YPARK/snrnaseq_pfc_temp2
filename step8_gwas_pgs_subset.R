#!/usr/bin/env Rscript

## LD.FILE <- "LD.info.txt"
## LD.IDX <- 1012
## GWAS.FILE <- "data/GWAS/ctg_ad.bed.gz"
## GENO.DIR <- "data/rosmap_geno/"
## EQTL.DIR <- "result/step7"
## TEMP.DIR <- "temp"
## OUT.FILE <- "TEMP.bed.gz"

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) != 7) q()

LD.FILE   <- argv[1]
LD.IDX    <- as.integer(argv[2])
GWAS.FILE <- argv[3]
GENO.DIR  <- argv[4]
EQTL.DIR  <- argv[5]
TEMP.DIR  <- argv[6]
OUT.FILE  <- argv[7]

## Take PGS vectors induced by a subset of eQTL loci

source("Util-geno.R")

.mkdir(dirname(OUT.FILE))
.mkdir(TEMP.DIR)

library(data.table)
library(tidyverse)

ld.tab <- fread(LD.FILE) %>%
    mutate(`#CHR`=str_remove(chr, "chr")) %>%
    mutate(`#CHR`=as.integer(`#CHR`)) %>%
    mutate(qq = `#CHR` %&% ":" %&% `start` %&% "-" %&% `stop`) %>%
    mutate(ld = 1:n()) %>%
    as.data.frame

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

take.prs.celltype <- function(ld.idx, window.size, lodds.cutoff) {

    qq <- ld.tab[ld.idx, "qq"]
    .chr <- ld.tab[ld.idx, "#CHR"]
    .lb <- ld.tab[ld.idx, "start"]
    .ub <- ld.tab[ld.idx, "stop"]

    .eqtl.file <- EQTL.DIR %&% "/chr" %&% .chr %&% "_stat.bed.gz"

    .gwas.dt <- rbind(fread("tabix -h " %&% GWAS.FILE %&% " " %&% qq),
                      fread("tabix -h " %&% GWAS.FILE %&% " chr" %&% qq))

    if(nrow(.gwas.dt) < 1) { return(NULL) }

    if("lodds" %in% names(.gwas.dt) && !("beta" %in% names(.gwas.dt))){
        .gwas.dt <- .gwas.dt %>%
            rename(beta = lodds) %>%
            as.data.table
    }

    .gwas.dt <- .gwas.dt %>%
        mutate(snp.loc = `stop`) %>%
        as.data.table

    .plink <- subset.plink(GENO.DIR %&% "/chr" %&% .chr,
                           .chr, .lb, .ub, TEMP.DIR)

    if(is.null(.plink))  { return(NULL) }

    valid.snps <- .plink$BIM %>%
        select(-missing, -rs) %>%
        mutate(x.pos = 1:n()) %>%
        (function(x) left_join(.gwas.dt, x, by = 'snp.loc')) %>%
        filter(!is.na(x.pos)) %>%
        filter((plink.a1 == a1 & plink.a2 == a2) | (plink.a2 == a1 & plink.a1 == a2)) %>%
        mutate(gwas.flip = if_else(plink.a1 == a1, 1.0, -1.0)) %>%
        mutate(beta = gwas.flip * beta) %>%
        select(-gwas.flip) %>%
        as.data.table

    if(nrow(valid.snps) < 1) {
        log.msg("skip this LD block... no valid SNPs found")
        return(NULL)
    }

    valid.snps <-
        valid.snps[order(abs(valid.snps$beta), decreasing=TRUE),
                   head(.SD, 1),
                   by = .(rs)] %>%
        arrange(snp.loc) %>%
        select(beta, se, x.pos, start, stop, plink.a1, plink.a2) %>% 
        as.data.table

    .svd <- zqtl::take.ld.svd(.plink$BED, eigen.tol = .01, eigen.reg = .01)
    iid <- .plink$FAM$iid

    .prs.by.subset <- function(.dt, .svd){

        zz <- .dt %>%
            mutate(z = beta/se) %>%
            select(z) %>%
            as.matrix

        .x.pos <- .dt$x.pos
        .svd.subset <- .svd
        .svd.subset$V.t <- .svd.subset$V.t[, .x.pos, drop=FALSE]

        zz <- zqtl::scale.zscore(zz, .svd.subset$V.t, .svd.subset$D)

        pred.prs(.svd.subset, zz, 0.01) %>% as.numeric
    }

    ytot <- .prs.by.subset(valid.snps, .svd)

    ret <- data.table(iid, y = ytot, celltype = "tot")

    ## Select GWAS SNPs by +/- window from eQTLs
    .eqtl.dt <- fread("tabix -h " %&% .eqtl.file %&% " " %&% qq) %>%
        filter(lodds > lodds.cutoff) %>%
        as.data.table

    .celltypes <- unique(.eqtl.dt$celltype)

    for(.ct in .celltypes) {

        .bed <- .eqtl.dt %>% filter(celltype == .ct) %>%
            mutate(lb = `start` - window.size,
                   ub = `stop` + window.size) %>% 
            select(lb, ub)

        .dt.ct <- data.table()

        for(b in 1:nrow(.bed)) {
            .lb <- .bed[b, ]$lb
            .ub <- .bed[b, ]$ub
            .temp <- valid.snps[`start` >= .lb & `stop` <= .ub]
            .dt.ct <- rbind(.dt.ct, .temp)
        }

        .dt.ct <- .dt.ct[, head(.SD,1), by=.(x.pos)]

        if(nrow(.dt.ct) < 1) next

        y.ct <- .prs.by.subset(.dt.ct, .svd)

        ret <- rbind(ret, data.table(iid, y = y.ct, celltype = .ct))
    }

    ret <- ret %>% mutate(ld = ld.idx)

    return(as.data.table(ret))
}

out.dt <- data.table()

for(w in c(1000, 1e5)) {
    for(lo in c(0, log(.9/.1))){
        .dt <- take.prs.celltype(LD.IDX, w = 100, lo = 0)
        if(!is.null(.dt)){
            .dt[, l10.window := log10(w)]
            .dt[, lodds.cutoff := lo]
            out.dt <- rbind(out.dt, .dt)
        }
    }
}

.mkdir(dirname(OUT.FILE))
if(nrow(out.dt) > 1) {
    out.dt <-
        ld.tab[LD.IDX, ] %>%
        left_join(out.dt) %>%
        select(`#CHR`, `start`, `stop`, ld, iid, y, celltype,
               l10.window, lodds.cutoff) %>% 
        as.data.table
    .bed.write(out.dt, OUT.FILE)
} else {
    write_tsv(out.dt, OUT.FILE)
}

unlink(TEMP.DIR)

log.msg("done")
