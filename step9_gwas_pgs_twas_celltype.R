#!/usr/bin/env Rscript

## LD.FILE <- "LD.info.txt"
## LD.IDX <- 1
## GWAS.NAME <- "ctg_ad"
## GWAS.DIR <- "result/step8/subset"
## EQTL.DIR <- "result/step7"
## OUT.HDR <- "TEMP"
## DO.PERMUTE <- FALSE

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 6) q()

LD.FILE <- argv[1]
LD.IDX <- as.integer(argv[2])
GWAS.NAME <- argv[3]
GWAS.DIR <- argv[4]
EQTL.DIR <- argv[5]
OUT.HDR <- argv[6]

if(length(argv) > 6){
    DO.PERMUTE <- as.logical(argv[7])
}

WINDOW <- 1e4
CUTOFF <- 0

source("Util-geno.R")

out.stat.file <- OUT.HDR %&% ".stat.gz"
out.pve.file <- OUT.HDR %&% ".pve.gz"

.mkdir(dirname(OUT.HDR))

if(all(file.exists(c(out.stat.file, out.pve.file)))){
    log.msg("Files exist...")
    q()
}

library(data.table)
library(tidyverse)

ld.tab <- fread(LD.FILE) %>%
    mutate(`#CHR`=str_remove(chr, "chr")) %>%
    mutate(`#CHR`=as.integer(`#CHR`)) %>%
    mutate(qq = `#CHR` %&% ":" %&% `start` %&% "-" %&% `stop`) %>%
    mutate(ld = 1:n()) %>%
    as.data.frame

take.twas.data <- function(ld.idx) {
    qq <- ld.tab[ld.idx, "qq"]
    .chr <- ld.tab[ld.idx, "#CHR"]
    .lb <- ld.tab[ld.idx, "start"]
    .ub <- ld.tab[ld.idx, "stop"]

    .gwas.file <- sprintf("%s/%s/%04d.bed.gz", GWAS.DIR, GWAS.NAME, ld.idx)
    .eqtl.file <- EQTL.DIR %&% "/chr" %&% .chr %&% "_poly.bed.gz"
    .gwas.dt <- fread("tabix -h " %&% .gwas.file %&% " " %&% qq)
    .eqtl.dt <- fread("tabix -h " %&% .eqtl.file %&% " " %&% qq)

    .eqtl.dt <- .eqtl.dt[, head(.SD, 1), by = .(iid, hgnc_symbol, celltype)]

    ## Matching the same cell type
    .tot <-
        .gwas.dt[celltype != "tot"] %>%
        rename(gwas.tot = y) %>% 
        filter(l10.window == log10(WINDOW)) %>% 
        filter(lodds.cutoff == CUTOFF) %>% 
        left_join(.eqtl.dt) %>% 
        as.data.table
    
    if(nrow(.tot) < 1) return(NULL)
    
    return(.tot)
}

.data <- take.twas.data(LD.IDX)

if(is.null(.data)) {
    log.msg("Empty data")
    write_tsv(data.frame(), out.stat.file)
    write_tsv(data.frame(), out.pve.file)
    q()
}

.var <- function(...) var(unlist(...), na.rm=TRUE)

.pve <- function(.dt, x, y) {
    .x <- x[, .dt$x.col, drop = FALSE]
    .y <- y[, unique(.dt$y.col), drop = FALSE]
    .hat <- .x %*% .dt$theta
    list(.var(.hat) / .var(.y))
}

.pve.uni <- function(.dt, x, y) {
    .x <- x[, .dt$x.col, drop = FALSE]
    .y <- y[, unique(.dt$y.col), drop = FALSE]
    .x[!is.finite(.x)] <- 0
    .hat <- lm(.y ~ .x) %>% predict
    list(.var(.hat) / .var(.y))
}

out.dt <- data.table()
pve.df <- data.frame()

for(.ct in unique(.data$celltype)){

    .data.ct <- 
        dcast(.data[celltype == .ct], iid + gwas.tot ~ hgnc_symbol,
              value.var = "yhat", fun.aggregate = mean) %>%
        select(-iid) %>%
        as.data.table

    y <- .data.ct %>% select(gwas.tot) %>% as.matrix %>% scale
    x <- .data.ct %>% select(-gwas.tot) %>% as.matrix %>% scale

    if(DO.PERMUTE){
        y <- apply(y, 2, function(.y) sample(.y))
        x <- apply(x, 2, function(.x) sample(.x))
    }

    y.dt <- data.table(gwas = GWAS.NAME, y.col = 1)

    x.dt <- data.table(hgnc_symbol = unlist(colnames(.data.ct)[-1])) %>%
        mutate(x.col = 1:n()) %>%
        as.data.table

    .stat <- zqtl::calc.qtl.stat(x, y)

    #######################
    ## multivariate stat ##
    #######################

    opts <- list(do.hyper=FALSE, pi=0, gammax=1e3, vbiter=5000,
                 print.interv=100, out.residual=FALSE, tol = 1e-6)

    .fit <- fqtl::fit.fqtl(y, x, options = opts)

    .mult.stat <-
        melt.fqtl.effect(.fit$mean) %>%
        rename(x.col = .row, y.col = .col)

    .dt <-
        left_join(.mult.stat, .stat) %>% 
        left_join(x.dt) %>%
        left_join(y.dt) %>%
        mutate(celltype = .ct) %>% 
        as.data.table

    out.dt <- rbind(out.dt, .dt)

    ##############
    ## variance ##
    ##############

    vtot <- .var(.data$`gwas.tot`)

    .df <-
        .dt[,
            .(pve = .pve(.SD, x, y),
              pve.uni = .pve.uni(.SD, x, y),
              ngene = .N),
            by = .(celltype)] %>%
        mutate(vtot) %>%
        mutate(pve = unlist(pve)) %>% 
        mutate(pve.uni = unlist(pve.uni)) %>% 
        as_tibble
    
    pve.df <- rbind(pve.df, .df)
}

out.df <-
    out.dt[, .(gwas, celltype, hgnc_symbol,
               theta, theta.var, lodds,
               beta, se, n, p.val)] %>%
    as_tibble

write_tsv(out.df, out.stat.file)
write_tsv(pve.df, out.pve.file)
