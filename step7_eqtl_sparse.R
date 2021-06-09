#!/usr/bin/env Rscript

## GENE.INFO.FILE <- "result/step1/gene.info.gz"
## GENE           <- 2
## DATA.DIR       <- "result/step6/eqtl/data/"
## PLINK.DIR      <- "data/rosmap_geno/"
## TEMP.DIR       <- "temp"
## OUT.HDR        <- "temp"

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) != 6) q()

GENE.INFO.FILE <- argv[1]
GENE           <- as.integer(argv[2])
DATA.DIR       <- argv[3]
PLINK.DIR      <- argv[4]
TEMP.DIR       <- argv[5]
OUT.HDR        <- argv[6]

################################################################
source("Util-geno.R")
library(tidyverse)
library(data.table)

out.files <- 
    c(OUT.HDR %&% "_stat.bed.gz",
      OUT.HDR %&% "_poly.bed.gz")

if(all(file.exists(out.files))) {
    log.msg("files exist --> done")
    q()
}

################################################################
CIS.DIST <- 1e6 # cis distance
NCV      <- 7 # number of cross validation

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

#' @param k.select
#' @param .svd
#' @param z.dt
.pgs.rank <- function(k.select, .svd, z.dt) {
    zz <- z.dt %>% select(-x.col) %>% as.matrix
    y.col <- colnames(zz)
    if(length(k.select) < 1) {
        y.hat <- matrix(0, nrow=nrow(.svd$U), ncol=1)
        return(y.hat)
    }

    dd <- .svd$D[k.select]
    y.hat <- 
        sweep(.svd$U[, k.select, drop = FALSE], MARGIN=2, STATS=dd, FUN=`/`) %*%
        (.svd$V.t[k.select, z.dt$x.col, drop = FALSE]) %*%
        zz
    colnames(y.hat) <- y.col
    return(y.hat)
}

#' @param .svd
#' @param max.K
run.varimax <- function(.svd, max.K) {
    V <- t(.svd$V.t)

    max.K <- min(ncol(V), max.K)

    if(max.K < 1) return(NULL)

    dd <- .svd$D[1:max.K, , drop = FALSE]
    V <- V[, 1:max.K, drop = FALSE]
    .vm <- varimax(V, eps = 1e-4, normalize=FALSE)

    R <- .vm$rotmat
    U <- .svd$U[, 1:max.K, drop = FALSE] %*% R
    V <- V %*% R

    list(U = U, D = dd, V = V, V.t = t(V))
}

#' @param yy
#' @param ncv
build.cv.idx <- function(yy, ncv = NCV) {

    ntot <- nrow(yy)
    ntest <- ceiling(ntot/NCV)
    rand.idx <- sample(ntot)

    .fun <- function(j) {

        .test <- rand.idx[seq(j * ntest + 1, min((j + 1) * ntest, ntot))]
        .train <- setdiff(1:ntot, .test)

        list(test = .test, train = .train)
    }

    lapply(seq(0, (ncv - 1)), .fun)
}

GENE.INFO <- .fread(GENE.INFO.FILE, header=TRUE) %>%
    as.data.frame

read.gene.data <- function(g, data.dir = DATA.DIR, gene.info = GENE.INFO, temp.dir = TEMP.DIR) {

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
        na.omit

    xx <- .plink$BED[pos.df$x.pos, , drop = FALSE]

    yy <- t(as.matrix(.yy[, pos.df$y.pos, drop = FALSE]))

    ## Total SVD for all the individuals
    svd.tot <-
        zqtl::take.ld.svd(.plink$BED, eigen.reg = 1e-2) %>%
        run.varimax(max.K = 300)

    ## SVD for those who have gene expressions
    svd.obj <- svd.tot
    iid.tot <- x.pos.df$iid
    svd.tot$U <- svd.tot$U[x.pos.df$x.pos, , drop = FALSE]

    iid <- pos.df$iid
    svd.obj$U <- svd.obj$U[pos.df$x.pos, , drop = FALSE]

    ## Deal with an empty gene
    yy <-  scale(yy)
    n.valid <- apply(is.finite(yy), 2, sum)
    yy[, n.valid < 10] <- 0

    list(x = xx,
         y = yy,
         iid = iid,
         iid.tot = iid.tot,
         svd = svd.obj,
         svd.tot = svd.tot,
         snp.info = .plink$BIM,
         gene.info = expr.dt[, 1:6],
         pos = pos.df)
}

.data <- read.gene.data(GENE)

if(is.null(.data)) {
    write_tsv(data.frame(), OUT.HDR %&% "_stat.bed.gz")
    write_tsv(data.frame(), OUT.HDR %&% "_poly.bed.gz")
    log.msg("nothing to do")
    q()
}

U <- scale(.data$svd$U)
y <- .data$y
x <- .data$x

#################################
## select relevant SVD factors ##
#################################

opts <- list(do.hyper=FALSE, pi=0, gammax=1e3, vbiter=5000,
             print.interv=100, out.residual=FALSE, tol = 1e-6)

.fit <- fqtl::fit.fqtl(y, U, options = opts)
.fqtl.dt <- melt.fqtl.effect(.fit$mean)

.stat <-
    zqtl::calc.qtl.stat(.data$x, .data$y) %>% 
    mutate(z = beta / se) %>%
    na.omit

##################################
## polygenic prediction results ##
##################################

.gene.info <- .data$gene.info %>%
    mutate(y.col = 1:n())

svd.obj <- .data$svd
svd.tot <- .data$svd.tot

iid.tot <- .data$iid.tot
out.hat <- tibble()

.cutoff <- (log(.1) - log(.9))

.temp <- .fqtl.dt %>%
    group_by(.col) %>%
    slice(which.max(lodds)) %>%
    ungroup %>%
    filter(lodds > .cutoff)

for(g in .temp$.col) {

    k.g <-
        .fqtl.dt %>%
        filter(.col == g) %>%
        filter(lodds > .cutoff) %>%
        select(.row) %>%
        unlist

    y.g <- .data$y[, g]

    zz.g <- .stat %>% 
        filter(y.col == g) %>% 
        as.data.table %>% 
        dcast(x.col ~ y.col,
              value.var = "z",
              fun.aggregate = mean)

    y.hat.g <- .pgs.rank(k.g, svd.obj, zz.g)

    ## compute PVE
    .hat.g <- predict(lm(y.g ~ y.hat.g))
    pve.g <- var(.hat.g, na.rm=TRUE) / var(y.g, na.rm=TRUE)

    ## PVE bootstrap
    pve.boot.g <-
        sapply(1:200,
               function(b) {
                   .idx <- sample(length(y.g), length(y.g), replace=TRUE)
                   .hat.boot <- predict(lm(y.g[.idx] ~ y.hat.g[.idx]))
                   .boot <- var(.hat.boot, na.rm=TRUE) / var(y.g, na.rm=TRUE)
               })
    
    ## compute total PGS
    y.hat.g <- .pgs.rank(k.g, svd.tot, zz.g) %>%
        unlist %>%
        as.numeric

    out.g <- tibble(iid = iid.tot,
                    yhat = y.hat.g,
                    pve = pve.g,
                    pve.sd = sd(pve.boot.g),
                    y.col = g)

    out.hat <- rbind(out.hat, out.g)
}

out.hat <-
    left_join(.gene.info, out.hat) %>%
    select(-y.col) %>%
    as.data.table

.bed.write(out.hat, OUT.HDR %&% "_poly.bed.gz")

#########################################
## marginal and sparse eQTL statistics ##
#########################################

opts <- list(do.hyper=FALSE, pi=0, gammax=1e3, vbiter=5000,
             print.interv=100, out.residual=FALSE, tol = 1e-6)

.fit <- fqtl::fit.fqtl(y, x, options = opts)

.sparse.dt <- melt.fqtl.effect(.fit$mean) %>%
    rename(x.col = .row, y.col = .col) %>%
    mutate(theta.sd = sqrt(theta.var))

.gene.info <- .data$gene.info %>%
    mutate(y.col = 1:n())

.snp.info <- .data$snp.info %>%
    mutate(x.col = 1:n())

out.stat <-
    .stat %>%
    left_join(.snp.info) %>% 
    left_join(.gene.info) %>% 
    left_join(.sparse.dt) %>% 
    mutate(start = snp.loc - 1) %>% 
    mutate(stop = snp.loc) %>% 
    arrange(`#chr`, `stop`) %>% 
    select(`#chr`, `start`, `stop`, starts_with("plink"), ensembl_gene_id, hgnc_symbol, celltype, beta, se, p.val, n, theta, theta.sd, lodds) %>% 
    as.data.table

.bed.write(out.stat, OUT.HDR %&% "_stat.bed.gz")

log.msg("done")
