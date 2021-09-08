#!/usr/bin/env Rscript

## LD.FILE <- "LD.info.txt"
## LD.IDX <- 1
## GWAS.FILE <- "result/step8/ctg_ad.bed.gz"
## EQTL.DIR <- "result/step7"
## OUT.HDR <- "TEMP"

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) != 5) q()

LD.FILE <- argv[1]
LD.IDX <- as.integer(argv[2])
GWAS.FILE <- argv[3]
EQTL.DIR <- argv[4]
OUT.HDR <- argv[5]

source("Util-geno.R")

out.stat.file <- OUT.HDR %&% ".stat.gz"
out.pve.file <- OUT.HDR %&% ".pve.gz"
lodds.cutoff <- log(.9/.1)

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

take.ld.data <- function(ld.idx) {
    qq <- ld.tab[ld.idx, "qq"]
    .chr <- ld.tab[ld.idx, "#CHR"]
    .lb <- ld.tab[ld.idx, "start"]
    .ub <- ld.tab[ld.idx, "stop"]

    .eqtl.file <- EQTL.DIR %&% "/chr" %&% .chr %&% "_poly.bed.gz"

    .gwas.dt <- fread("tabix -h " %&% GWAS.FILE %&% " " %&% qq) %>%
        rename(gwas.hat = y.hat, gwas.tot = y.tot)

    .eqtl.dt <- fread("tabix -h " %&% .eqtl.file %&% " " %&% qq)

    .dt <- left_join(.eqtl.dt, .gwas.dt) %>% 
        as.data.table

    if(nrow(.dt) < 1) return(NULL)

    dcast(.dt, iid + gwas.tot + gwas.hat ~ hgnc_symbol + celltype,
          value.var = "yhat",
          fun.aggregate = mean) %>%
        select(-iid) %>%
        as.data.table
}

.data <- take.ld.data(LD.IDX)

if(is.null(.data)) {
    log.msg("Empty data")
    write_tsv(data.frame(), out.stat.file)
    write_tsv(data.frame(), out.pve.file)
    q()
}

y <- .data %>% select(gwas.hat) %>% as.matrix %>% scale
x <- .data %>% select(-gwas.tot, -gwas.hat) %>% as.matrix %>% scale

gwas.name <- str_remove(basename(GWAS.FILE), ".bed.gz")

y.dt <- data.table(gwas = gwas.name, y.col = 1)

x.dt <- data.table(xx = unlist(colnames(.data)[-(1:2)])) %>%
    mutate(x.col = 1:n()) %>%
    separate(xx, c("hgnc_symbol", "celltype"), sep="[_]") %>% 
    as.data.table

.stat <- zqtl::calc.qtl.stat(x, y)

opts <- list(do.hyper=FALSE, pi=0, gammax=1e3, vbiter=5000,
             print.interv=100, out.residual=FALSE, tol = 1e-6)

.fit <- fqtl::fit.fqtl(y, x, options = opts)

.mult.stat <-
    melt.fqtl.effect(.fit$mean) %>%
    rename(x.col = .row, y.col = .col)

out.dt <-
    left_join(.mult.stat, .stat) %>% 
    left_join(x.dt) %>%
    left_join(y.dt) %>%
    as.data.table

.var <- function(...) var(unlist(...), na.rm=TRUE)

vtot <- .var(.data$`gwas.tot`)
vtot.hat <- .var(.data$`gwas.hat`)

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

if(nrow(out.dt[lodds > lodds.cutoff]) > 0) {

    pve.df <-
        out.dt[lodds > lodds.cutoff,
               .(pve = .pve(.SD, x, y),
                 pve.uni = .pve.uni(.SD, x, y),
                 ngene = .N),
               by = .(celltype)] %>%
        mutate(vtot, vtot.hat) %>%
        mutate(pve = unlist(pve)) %>% 
        mutate(pve.uni = unlist(pve.uni)) %>% 
        as_tibble

} else {
    pve.df <- data.frame()
}

out.df <-
    out.dt[, .(gwas, celltype, hgnc_symbol,
               theta, theta.var, lodds,
               beta, se, n, p.val)] %>%
    as_tibble

write_tsv(out.df, out.stat.file)
write_tsv(pve.df, out.pve.file)
