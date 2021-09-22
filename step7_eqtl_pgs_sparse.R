#!/usr/bin/env Rscript

## GENE.INFO.FILE <- "result/step1/gene.info.gz"
## GENE           <- 2
## DATA.DIR       <- "result/step6/eqtl/data/"
## PLINK.DIR      <- "data/rosmap_geno/"
## TEMP.DIR       <- "temp"
## OUT.HDR        <- "temp"

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) != 7) q()

GENE.INFO.FILE <- argv[1]
GENE           <- as.integer(argv[2])
DATA.DIR       <- argv[3]
PLINK.DIR      <- argv[4]
TEMP.DIR       <- argv[5]
LODDS.CUTOFF   <- argv[6]
OUT.HDR        <- argv[7]

source("Util-geno.R")
library(tidyverse)
library(data.table)

.stat.file <- OUT.HDR %&% "_stat.bed.gz"

if(!(file.exists(.stat.file))) {
    write_tsv(data.frame(), OUT.HDR %&% "_poly_sparse.bed.gz")
    q()
}

.stat <- .fread(.stat.file)
cts <- .stat[lodds > LODDS.CUTOFF, .(celltype)] %>% unlist

if(length(cts) < 1){
    write_tsv(data.frame(), OUT.HDR %&% "_poly_sparse.bed.gz")
    q()
}

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

read.gene.data <- function(g, data.dir = DATA.DIR,
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
        na.omit

    xx <- .plink$BED[pos.df$x.pos, , drop = FALSE]

    yy <- t(as.matrix(.yy[, pos.df$y.pos, drop = FALSE]))

    iid <- pos.df$iid


    ## Deal with an empty gene
    yy <-  scale(yy)
    n.valid <- apply(is.finite(yy), 2, sum)
    yy[, n.valid < 10] <- 0

    list(x = xx,
         y = yy,
         iid = iid,
         x.tot = scale(.plink$BED),
         iid.tot = .plink$FAM$iid,
         snp.info = .plink$BIM,
         gene.info = expr.dt[, 1:6],
         pos = pos.df)
}

.data <- read.gene.data(GENE)
x <- .data$x.tot
x[is.na(x)] <- 0
iid.tot <- .data$iid.tot
out.hat <- data.frame()

for(.ct in cts) {
    .stat.ct <- .stat[celltype == .ct & lodds > LODDS.CUTOFF]
    .pos.ct <- match(.stat.ct$`stop`, .data$snp.info$snp.loc)
    .x <- x[, .pos.ct, drop = FALSE]
    .y.ct <- .x %*% matrix(.stat.ct$theta, ncol = 1)

    out.g <- tibble(iid = iid.tot,
                    yhat = .y.ct,
                    celltype = .ct)

    out.hat <- rbind(out.hat, out.g)
}

out.hat <-
    .data$gene.info %>%
    left_join(out.hat) %>%
    na.omit

.bed.write(out.hat, OUT.HDR %&% "_poly_sparse.bed.gz")
