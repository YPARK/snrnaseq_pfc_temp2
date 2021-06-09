#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
`%&%` <- function(a,b) paste0(a,b)


.mu.file <- "result/step4/aggregate/subtype.mu.gz"
.sum.file <- "result/step4/aggregate/subtype.sum.gz"
.col.file <- "result/step4/aggregate/subtype.mu_cols.gz"

.col.dt <- fread(.col.file, header=FALSE, col.names = c("col"))
.col.dt[, c("projid", "celltype") := tstrsplit(col, "_")]

.mu.dt <- fread(.mu.file, header=FALSE, col.names = .col.dt$col)
.sum.dt <- fread(.sum.file, header=FALSE, col.names = .col.dt$col)

## remove samples with too few observations
nobs <- apply(as.matrix(.sum.dt), 2, sum)
.valid <- which(nobs >= 100)
.mu <- as.matrix(.mu.dt)[, .valid]
.col <- .col.dt[.valid, ]

X <- scale(.mu, center=FALSE, scale=TRUE) * 1e4
log.X <- scale(t(log(X + 1)))

celltype <- .col$celltype

.umap.dt <- uwot::umap(log.X,
                       fast_sgd = TRUE,
                       n_neighbors = 50,
                       metric = "cosine",
                       n_threads = 10,
                       spread = 5,
                       verbose = TRUE) %>%
    as.data.table



colnames(.umap.dt) <- "umap" %&% 1:ncol(.umap.dt)
.umap.dt[, celltype := celltype]

.lab.dt <- .umap.dt[, .(umap1=median(umap1), umap2=median(umap2)), by = .(celltype)]

plt <-
    ggplot(.umap.dt, aes(x=umap1, y=umap2, fill = celltype)) +
    theme_classic() +
    geom_point(pch=21, stroke=.1) +
    geom_text(aes(label=celltype), data=.lab.dt, size=4)

.file <- "result/step4/Fig_neuron_celltype_layers.pdf"
ggsave(.file, plot=plt, units="in", width=5, height=5)
