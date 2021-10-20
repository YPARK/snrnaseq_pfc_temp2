batch.file  <- "result/step1/merged.columns.gz"
annot.file  <- "result/step3/bbknn.annot.gz"
factor.file <- "result/step3/bbknn.factors.gz"
out.file    <- "result/step3/bbknn.umap.gz"

argv <- commandArgs(trailingOnly = TRUE)

batch.file  <- argv[1] # e.g., "result/step1/merged.columns.gz"
annot.file  <- argv[2] # "result/step3/bbknn.annot.gz"
factor.file <- argv[3] # "result/step3/bbknn.factors.gz"
out.file    <- argv[4] # "result/step3/bbknn.umap.gz"

library(data.table)
library(tidyverse)

svd_d.file <- str_replace(factor.file, "factors", "svd_D")
D <- fread(svd_d.file, header=FALSE) %>% unlist %>% as.numeric

V <- fread(factor.file, header=FALSE) %>% as.matrix

batch <- fread(batch.file, col.names = c("Barcode","batch"), header=FALSE)

.cols <- c("Barcode", "celltype", "prob", "ln.prob")
annot <- fread(annot.file, col.names = .cols) %>%
    mutate(j = 1:n()) %>%
    left_join(batch) %>% 
    as.data.table

annot[, c("barcode","projid") := tstrsplit(Barcode,split="_")]

VD <- sweep(V, 2, D, `*`)

umap.mat <- uwot::umap(VD[annot$j, ],
                       metric = "cosine",
                       fast_sgd = TRUE,
                       verbose = TRUE,
                       spread = 3,
                       n_components = 2,
                       n_threads = 15)

colnames(umap.mat) <- paste0("UMAP", 1:ncol(umap.mat))

out.dt <- cbind(annot, as.data.table(umap.mat))
fwrite(out.dt, out.file, sep="\t")

