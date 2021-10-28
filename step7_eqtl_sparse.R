#!/usr/bin/env Rscript

## GENE.INFO.FILE <- "result/step1/gene.info.gz"
## GENE           <- 67
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

out.files <- OUT.HDR %&% "_stat.bed.gz"

if(all(file.exists(out.files))) {
    log.msg("files exist --> done")
    q()
}

################################################################
CIS.DIST <- 5e5 # cis distance

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

#' @param yy
#' @param ncv
build.cv.idx <- function(yy, ncv) {

    ntot <- nrow(yy)
    ntest <- ceiling(ntot/ncv)
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

    iid <- pos.df$iid

    ## Deal with an empty gene
    n.valid <- apply(is.finite(yy), 2, sum)

    ## convert back to pseudo-counting data
    yy <- scale(yy) %>% exp() %>% as.matrix
    yy[, n.valid < 10] <- 0

    list(x = xx,
         y = yy,
         iid = iid,
         snp.info = .plink$BIM,
         gene.info = expr.dt[, 1:6],
         pos = pos.df)
}

.data <- read.gene.data(GENE)

if(is.null(.data)) {
    write_tsv(data.frame(), OUT.HDR %&% "_stat.bed.gz")
    log.msg("nothing to do")
    q()
}

.maf <-
    data.table(maf.1 = apply(.data$x, 2, mean),
               maf.2 = apply(2 - .data$x, 2, mean)) %>%
    mutate(maf = pmin(maf.1, maf.2)/2) %>%
    mutate(x.col =1:n()) %>%
    select(maf, x.col)

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

y <- .data$y
.ko <- svd.knockoff(.data$x)
x <- .ko$x
x.ko <- .ko$ko

if(max(apply(y, 2, sd, na.rm=TRUE), na.rm=TRUE) < 1e-4) {
    write_tsv(data.frame(), OUT.HDR %&% "_stat.bed.gz")
    log.msg("nothing to do")
    q()
}

if(DO.PERMUTE){
    y <- apply(y, 2, sample)
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

y.qnorm <- apply(y, 2, .qnorm)

.stat <-
    zqtl::calc.qtl.stat(x, y.qnorm) %>% 
    mutate(z = beta / se) %>%
    na.omit

#########################################
## marginal and sparse eQTL statistics ##
#########################################

opts <- list(do.hyper = TRUE,
             svd.init = TRUE,
             jitter = 0,
             pi.ub = 0,
             pi.lb = log(.1/.9),
             gammax = 1e3,
             vbiter = 3500,
             print.interv = 100,
             out.residual = FALSE,
             rate = 1e-2,
             decay = 0,
             tol = 1e-6)

#######################
## cross-validataion ##
#######################

.cv.set <- build.cv.idx(y, 5)
cv.tab <- data.table()

for(.cv in .cv.set) {

    y.train <- y[.cv$train, , drop = FALSE]
    x.train <- x[.cv$train, , drop = FALSE]
    y.test <- y[.cv$test, , drop = FALSE]
    x.test <- x[.cv$test, , drop = FALSE]

    .train <- fqtl::fit.fqtl(y = y.train,
                             x.mean = x.train,
                             model = "nb",
                             options = opts)

    y.pred <- x.test %*% .train$mean$theta

    rr <-
        sapply(1:ncol(y.pred), function(j) {
            cor(y.pred[,j], y.test[,j],
                method = "spearman",
                use = "pairwise.complete.obs")
        })

    cv.tab <- rbind(cv.tab, data.table(rr, y.col = 1:length(rr)))
}

cv.dt <- cv.tab[, .(r = mean(rr)), by = .(y.col)]

###################
## full training ##
###################

.fit <- fqtl::fit.fqtl(y = y,
                       x.mean = x,
                       c.mean = x.ko,
                       model = "nb",
                       options = opts)

####################
## output results ##
####################

.sparse.dt <- melt.fqtl.effect(.fit$mean) %>%
    rename(x.col = .row, y.col = .col) %>%
    mutate(theta.sd = sqrt(theta.var))

.ko.dt <- melt.fqtl.effect(.fit$mean.cov) %>%
    rename(x.col = .row, y.col = .col) %>%
    mutate(theta.sd = sqrt(theta.var))

.sparse.dt <- .sparse.dt %>%
    left_join(.ko.dt,
              by = c("x.col", "y.col"),
              suffix = c("", ".ko"))

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
    left_join(.maf) %>% 
    left_join(cv.dt) %>% 
    arrange(`#chr`, `stop`) %>% 
    select(`#chr`, `start`, `stop`, starts_with("plink"),
           ensembl_gene_id, hgnc_symbol, celltype,
           beta, se, p.val, n, maf, r,
           theta, theta.sd, lodds,
           theta.ko, theta.sd.ko, lodds.ko) %>% 
    as.data.table

.bed.write(out.stat, OUT.HDR %&% "_stat.bed.gz")

log.msg("done")
