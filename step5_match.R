#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

col.file <- argv[1]
pheno.file <- argv[2]
out.file <- argv[3]

library(tidyverse)
library(data.table)

col.dt <- fread(paste0("gzip -cd ", col.file), header=FALSE)
pheno.dt <- fread(paste0("gzip -cd ", pheno.file), header=FALSE) %>% distinct
out.dt <- left_join(col.dt, pheno.dt, by = "V1") %>% as.data.table
fwrite(out.dt[, .(V2)], file = out.file, col.names = FALSE, row.names = FALSE)

