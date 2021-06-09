#!/usr/bin/env/Rscript --vanilla
argv <- commandArgs(trailingOnly = TRUE)

row.file <- argv[1] # e.g., row.file <- "result/step1/merged.rows.gz"
out.file <- argv[2] # e.g., "result/step1/gene.info.gz"

library(dplyr)
library(data.table)

.rows <- fread(row.file, header=FALSE, col.names="gene")
.rows[, c("ensembl_gene_id", "hgnc_symbol") := tstrsplit(gene,split="_")]

read.mart <- function(.ensg, .ensembl.hs) {

    .attr <- c("ensembl_gene_id",
               "chromosome_name",
               "transcript_start",
               "transcript_end")

    .temp <- biomaRt::getBM(attributes=.attr,
                            filters="ensembl_gene_id",
                            values=.ensg,
                            mart=.ensembl.hs,
                            useCache = FALSE)

    .temp <- as.data.table(.temp)

    .temp <- .temp[, .(chromosome_name = unique(chromosome_name),
                       transcript_start = min(transcript_start),
                       transcript_end = max(transcript_end)),
                   by = .(ensembl_gene_id)]

    return(.temp)
}

ensembl.hg19 <- biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL",
                                 host="grch37.ensembl.org",
                                 path="/biomart/martservice",
                                 dataset="hsapiens_gene_ensembl")

ensembl.hs.hg19 = biomaRt::useDataset("hsapiens_gene_ensembl",
                                      mart=ensembl.hg19)

.info.hg19 <- read.mart(.rows$ensembl_gene_id, ensembl.hs.hg19)

ensembl <- biomaRt::useMart(biomart = "ensembl", 
                            dataset = "hsapiens_gene_ensembl")

ensembl.hs = biomaRt::useDataset("hsapiens_gene_ensembl",
                                 mart=ensembl)

.info.recent <- read.mart(.rows$ensembl_gene_id, ensembl.hs)

.info <- merge(.info.hg19, .info.recent,
               by="ensembl_gene_id",
               suffixes=c(".hg19", ".recent"))

.info[, chr := chromosome_name.hg19]
.info[!(chr %in% c(1:22, "X", "Y")), chr := chromosome_name.recent]

out <- .rows %>%
    left_join(.info[, .(ensembl_gene_id, chr,
                        transcript_start.hg19,
                        transcript_end.hg19,
                        transcript_start.recent,
                        transcript_end.recent)])
             
fwrite(out, out.file, row.names=FALSE, col.names=TRUE, sep="\t")

