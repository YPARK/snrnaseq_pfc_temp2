.row.file <- "result/step1/merged.rows.gz"
.col.file <- "result/step3/aggregate/broad.mu_cols.gz"
.mu.file <- "result/step3/aggregate/broad.mu.gz"
.marker.file <- "result/step3/marker_broad.txt.gz"

argv <- commandArgs(trailingOnly = TRUE)

.row.file <- argv[1] # e.g., "result/step1/merged.rows.gz"
.col.file <- argv[2] # e.g., "result/step3/aggregate/broad.mu_cols.gz"
.mu.file <- argv[3] # e.g., "result/step3/aggregate/broad.mu.gz"
.marker.file <- argv[4] # e.g., "result/step3/marker_broad.txt.gz"

library(tidyverse)
library(data.table)
library(patchwork)
source("Util.R")

.num.sci <- function(x) format(x, digits=2, scientific=TRUE)
.num.int <- function(x) format(as.integer(x), big.mark = ',')
.num.round <- function(x) round(x, digits=2)
`%&%` <- function(a,b) paste0(a,b)
.fread <- function(...) fread(..., header = FALSE)
setDTthreads(16)


.rows <- .fread(.row.file) %>% unlist
.cols <- .fread(.col.file) %>% unlist

.mu.dt <- .fread(.mu.file, colClasses="double")

ln.mat <- log(1 + t(as.matrix(.mu.dt)))

.svd <- rsvd::rsvd(ln.mat, k = 10)
colnames(.svd$u) <- "PC" %&% 1:ncol(.svd$u)

.svd.dt <- data.table(.svd$u, col = .cols)
.svd.dt[, c("projid", "celltype") := tstrsplit(col, "_", fixed=TRUE)]

plot.scatter <- function(gg) {
    gg +
        theme_classic() +
        theme(title = element_text(size = 8)) +
        theme(legend.text = element_text(size = 8)) +
        theme(legend.key.height = unit(.5, "lines")) +
        theme(legend.key.width = unit(.2, "lines")) +
        theme(legend.text = element_text(size=8)) +
        geom_point(size = .5) +
        scale_colour_brewer(palette = "Paired")
}

.tsne <- Rtsne::Rtsne(.svd$u, dims = 2, initial_dims = 10,
                      check_duplicates = FALSE, perplexity = 30,
                      num_threads = 16,
                      pca = FALSE, verbose = TRUE)

colnames(.tsne$Y) <- "tSNE" %&% 1:ncol(.tsne$Y)

.tsne.dt <- data.table(.tsne$Y, celltype = .svd.dt$celltype)

p1 <- ggplot(.svd.dt, aes(x = PC1, y = PC2, colour = celltype)) %>%
    plot.scatter

p2 <- ggplot(.svd.dt, aes(x = PC2, y = PC3, colour = celltype)) %>%
    plot.scatter

p3 <- ggplot(.tsne.dt, aes(x = tSNE1, y = tSNE2, colour = celltype)) %>%
    plot.scatter

plt <- p1 | p2 | p3

ggsave(file = "result/step3/celltype_aggregate.pdf", plot=plt,
       width = 9, height = 3)

marker.genes <- .fread(.marker.file)
celltypes <- sapply(.cols, function(x) tstrsplit(x, "_", fixed=TRUE)[[2]])

.gene.pos <- unique(match(marker.genes$V1, .rows)) %>% sort

.gene.names <- .rows[.gene.pos]

.row.order <- .cols[order(celltypes)] # now samples = rows

## visualize top marker genes

.dt <-
    .mu.dt[.gene.pos, ] %>% 
    (function(x) { colnames(x) <- .cols; x }) %>%
    mutate(col = .gene.names) %>%
    melt(id.vars = "col", variable.name = "row", value.name = "weight") %>% 
    mutate(weight = log(1 + weight)) %>% 
    na.omit %>%
    as.data.table

.dt[, c("projid","celltype") := tstrsplit(row, "_")]

## take top genes per each cell type

.sum <- .dt[, .(mu = mean(weight)), by = .(col, celltype)]

.top.genes <-
    .sum[order(.sum$mu), tail(.SD, 30), by = .(celltype)] %>%
    select(col) %>%
    unique %>%
    unlist

.dt <- .dt[col %in% .top.genes] %>% 
    col.order(.ro = .row.order, ret.tab = TRUE) %>%
    as.data.table

p1 <-
    ggplot(.dt, aes(y = col, x = row, fill = pmin(pmax(weight, 0), 4))) +
    theme_classic() +
    theme(axis.text = element_blank()) +
    theme(axis.ticks = element_blank()) +
    ylab(length(.top.genes) %&% " marker genes (top 30 per type)") +
    xlab("pseudo-bulk samples") +
    scale_fill_distiller(palette = "RdPu", direction = 1, guide=FALSE) +
    geom_tile()

.dt.col <- .dt %>%
    select(row, celltype) %>%
    distinct %>%
    mutate(y = 1)

p0 <-
    ggplot(.dt.col, aes(y=y, x=row, fill=celltype)) +
    theme_void() +
    geom_tile() +
    geom_text(aes(label=celltype), data=.dt.col[, head(.SD, 1), by = .(celltype)],
              size=4, hjust=0) +
    scale_fill_brewer(palette = "Paired", guide=FALSE)

plt <- (p0 / p1) + plot_layout(heights=c(.05, .95))

ggsave(file = "result/step3/celltype_pseudo_heatmap.pdf", plot=plt,
       width = 8, height = 7)
