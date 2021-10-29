#!/usr/bin/env Rscript
argv <- commandArgs(trailingOnly = TRUE)

.chr <- argv[1]
.pheno <- argv[2]

library(tidyverse)
library(data.table)
source("Util-geno.R")

.list.files <- function(...) list.files(..., full.names=TRUE)

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

eFDR <- function (stat, stat0)
{
    m <- length(stat)
    n <- ncol(stat0)
    m0 <- length(stat0)
    v <- c(rep(TRUE, m), rep(FALSE, m0))
    v <- v[order(c(stat, stat0), decreasing = TRUE)]
    u <- 1:length(v)
    w <- 1:m
    p <- (u[v == TRUE] - w)/m0
    p <- p[rank(-stat)]
    p <- pmax(p, 1/m0)
}

gene.info <- .fread("result/step1/gene.info.gz") %>%
    mutate(g = 1:n())

.out.chr.file <- str_c("share/2021-10-29/interaction/", .pheno, "/chr", .chr, ".geno.bed.gz")

.mkdir(dirname(.out.chr.file))

.genes <-
    gene.info %>%
    filter(chr == .chr) %>%
    select(g) %>%
    mutate(g = sprintf("%05d",g)) %>%
    unlist

max.g <- gene.info %>%
    filter(chr == .chr) %>%
    select(g) %>%
    unlist %>%
    max

if(!file.exists(.out.chr.file)) {

    .files <- str_c("result/step7/interaction/", .genes, "_stat.bed.gz")

    .dt.chr <-
        lapply(.files, function(x){
            .sz <- file.info(x)$size
            if(.sz < 100) return(data.table())

            .dt <- .fread(x) %>%
                mutate(`start` = as.integer(`start`)) %>%
                mutate(`stop` = as.integer(`stop`)) %>%
                filter(pheno == .pheno) %>%
                as.data.table

            .out <- .dt[ko == 0]

            stat1 <- .out[, .(lodds)] %>% unlist
            stat0 <- .dt[ko == 1, .(lodds)] %>% unlist
            .efdr <- eFDR(stat1, stat0)

            .out[, efdr := .efdr]

            g <- basename(x) %>% str_remove("_stat.bed.gz$") %>% as.integer
            log.msg("[%05d / %05d]", g, max.g)

            return(.out)
        }) %>%
    do.call(what = rbind)

    .dt.chr <- .dt.chr[order(.dt.chr$`stop`), ]
    .bed.write(.dt.chr, .out.chr.file)
}

.out.chr.file <- str_c("share/2021-10-29/interaction/", .pheno, "/chr", .chr, ".celltype.bed.gz")
.mkdir(dirname(.out.chr.file))

if(!file.exists(.out.chr.file)){
    .files <- str_c("result/step7/interaction/", .genes, "_pheno.bed.gz")

    .dt.chr <-
        lapply(.files, function(x){
            .sz <- file.info(x)$size
            if(.sz < 100) return(data.table())

            .dt <- .fread(x) %>%
                mutate(`transcript_start` = as.integer(`transcript_start`)) %>%
                mutate(`transcript_end` = as.integer(`transcript_end`)) %>%
                filter(pheno == .pheno) %>%
                as.data.table

            g <- basename(x) %>% str_remove("_pheno.bed.gz$") %>% as.integer
            log.msg("[%05d / %05d]", g, max.g)

            return(.dt)
        }) %>%
    do.call(what = rbind)

    .dt.chr <- .dt.chr[order(.dt.chr$`transcript_start`), ]
    .bed.write(.dt.chr, .out.chr.file)
}
