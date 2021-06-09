#!/usr/bin/env Rscript

## ln.resid.file   <- "result/step5/cfa/Microglia-broad.ln_resid_mu.gz"
## ln.obs.file     <- "result/step5/cfa/Microglia-broad.ln_obs_mu.gz"
## col.file     <- "result/step5/cfa/Microglia-broad.mu_cols.gz"
## row.file     <- "result/step5/filtered/Microglia-broad.rows.gz"
## pheno.file   <- "data/rosmap_phenotype.csv.gz"

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 6) {
    q()
}

ln.resid.file <- argv[1]
ln.obs.file    <- argv[2]
col.file      <- argv[3]
row.file      <- argv[4]
pheno.file    <- argv[5]
out.file      <- argv[6]

if(!all(file.exists(argv[-6]))) {
    q()
}

if(file.exists(argv[6])) {
    q()
}

################################################################

library(tidyverse)
library(data.table)
source("Util.R")

################################################################

.fread.melt <- function(x, .name, rows, cols = NULL, col.file = NULL) {

    if(is.null(cols)){
        cols <- fread(col.file, header = FALSE, col.names = "sample")
        cols[, c("projid", "celltype") := tstrsplit(sample, "_")]
        cols <- cols$projid %>% unlist
    }

    .ret <- fread(x, header = FALSE, col.names = as.character(cols))
    .ret[, gene := rows]
    ret <- melt.data.table(.ret, id.vars = "gene",
                           variable.name = "projid",
                           value.name = .name)
    ret[, projid := as.character(projid)]
    ret[, projid := as.integer(projid)]
}

rows <- fread(row.file, header = FALSE, col.names = "gene")
rows[, c("ensg", "hgnc") := tstrsplit(gene, "_")]

.pheno <- fread(pheno.file, header=TRUE, na.strings = "-9")

.celltype <- fread(col.file, nrow=1, header=F) %>%
    unlist %>%
    (function(x) str_split(x, pattern="[_]")[[1]][2])

## Reading the aggregate results
obs.dt <- .fread.melt(ln.obs.file, "obs", rows$gene, col.file = col.file)

## Reading the counterfactually-adjusted results
cfa.dt <- .fread.melt(ln.resid.file, "resid", rows$gene, col.file = col.file)

## filter out samples with too much or little adjusted
sample.stat <- cfa.dt[, .(m = mean(resid), s = sd(resid)), by = .(projid)]
.invalid.samples <- sample.stat[abs(m/s) > 3, .(projid)] %>% unlist

################################################################

combined.dt <- merge(obs.dt, cfa.dt) %>%
    filter(!(projid %in% .invalid.samples)) %>%
    merge(rows, by = "gene") %>%
    left_join(.pheno[, .(projid, pathoAD)], by = "projid")

## treat the spots with too little values as missing
combined.dt <- combined.dt[obs > -4 & resid > -4]

################################################################

.scan.test <- function(gg, xx, yy) {

    x.obs <- apply(t(is.finite(xx)), 2, mean)
    y.obs <- apply(t(is.finite(yy)), 2, mean)

    .valid <- which(x.obs > .25 & y.obs > .25)

    .t.test <- lapply(.valid, function(j) t.test(yy[j, ], xx[j, ]))
    .w.test <- lapply(.valid, function(j) wilcox.test(yy[j, ], xx[j, ]))

    pv.t <- sapply(.t.test, function(x) x$p.value)
    tt <- sapply(.t.test, function(x) x$statistic)
    se <- sapply(.t.test, function(x) x$stderr)
    pv <- sapply(.w.test, function(x) x$p.value)

    gg[.valid, ] %>%
        mutate(pv = pv, pv.t =pv.t, t = tt, se = se) %>%
        as.data.table
}

run.test <- function(.dt) {
    gg <- .dt %>% select(ensg, hgnc) %>% mutate(celltype = .celltype) %>% as_tibble

    xx <- .dt %>% select(starts_with("0_")) %>% as.matrix
    yy <- .dt %>% select(starts_with("1_")) %>% as.matrix

    return(.scan.test(gg, xx, yy))
}

resid.stat <- dcast(combined.dt, hgnc + ensg ~ pathoAD + projid, value.var = "resid", fill=NA) %>%
    run.test %>%
    mutate(method = "resid")

obs.stat <- dcast(combined.dt, hgnc + ensg ~ pathoAD + projid, value.var = "obs", fill=NA) %>%
    run.test %>%
    mutate(method = "obs")

out.dt <- rbind(resid.stat, obs.stat)

fwrite(out.dt, file = out.file, sep = "\t", row.names=FALSE)
