argv <- commandArgs(trailingOnly = TRUE)

marker.file <- argv[1] # e.g., "data/PsychENCODE.marker"
row.file <- argv[2] # e.g., "result/step1/merged.rows.gz"
out.file <- argv[3] # e.g., "result/step3/marker_broad.txt.gz"

`%&%` <- function(a,b) paste0(a,b)

library(tidyverse)
library(data.table)

dir.create("result/step3", recursive = TRUE, showWarnings = FALSE)

genes.dt <- fread(row.file, header=FALSE, col.names = "gene")
genes.dt[, c("ensg","hgnc") := tstrsplit(gene, "_", fixed=TRUE)]

ct.dt <- fread(marker.file, col.names =  c("hgnc", "celltype"))

annot.dt <- left_join(ct.dt, genes.dt, by="hgnc") %>%
    na.omit

.dt <- annot.dt[, .(gene, celltype)]
fwrite(.dt, out.file, sep = " ", col.names = FALSE)

