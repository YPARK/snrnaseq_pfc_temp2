
library(tidyverse)
library(data.table)
library(patchwork)

.num.sci <- function(x) format(x, digits=2, scientific=TRUE)

.files <- list.files("result/step1/", pattern = "*.score.gz", full.names = TRUE)

.fread <- function(xx) {
    batch <- basename(xx) %>% str_remove(".score.gz")
    ret <- fread(xx, header = FALSE)
    colnames(ret) <- c("gene","nnz", "nmax",".mean",".sd",".cv")

    ret[, batch := batch]
    ret[, c("ensg", "hgnc") := tstrsplit(gene, "_", fixed=TRUE)]

    return(ret)
}

.plot.hist <- function(.dt) {

    .dt <- .dt[order(.dt$p), ]

    ggplot(.dt, aes(x = p, y = log(n), group=batch)) +
        theme_classic() +
        geom_line(colour = "gray40", size=.2) +
        theme(legend.title = element_text(size=6)) +
        theme(legend.text = element_text(size=5)) +
        theme(legend.key.height = unit(.15, "lines")) +
        scale_y_continuous(labels = function(x) .num.sci(exp(x))) +
        scale_x_continuous(labels = function(x) .num.sci(exp(x)))
}

.plot.all <- function(scores.dt) {

    hist.dt <- copy(scores.dt)[, p := round(5*log((nnz + 1)/(nmax + 1)))/5]
    hist.dt <- hist.dt[, .(n = .N), by = .(batch, p)]

    p1 <- .plot.hist(hist.dt) +
        geom_vline(xintercept = log(1e-4), colour="red", lty = 2) +
        xlab("Frequency of the detected within a batch") + ylab("Frequency (# genes)")

    hist.dt <- copy(scores.dt)[, p := round(5 * log(.mean))/5]
    hist.dt <- hist.dt[, .(n = .N), by = .(batch, p)]

    p2 <- .plot.hist(hist.dt) +
        geom_vline(xintercept = log(1e-4), colour="red", lty = 2) +
        xlab("Average expression") + ylab("Frequency (# genes)")

    hist.dt <- copy(scores.dt)[, p := round(5 * log(.cv))/5]
    hist.dt <- hist.dt[, .(n = .N), by = .(batch, p)]

    p3 <- .plot.hist(hist.dt) +
        geom_vline(xintercept = log(1.1), colour="red", lty = 2) +
        xlab("Coefficient of Variation") + ylab("Frequency (# genes)")

    hist.dt <- copy(scores.dt)[, p := round(5 * log(.sd))/5]
    hist.dt <- hist.dt[, .(n = .N), by = .(batch, p)]

    p4 <- .plot.hist(hist.dt) +
        geom_vline(xintercept = log(1e-2), colour="red", lty = 2) +
        xlab("Standard Deviation") + ylab("Frequency (# genes)")

    (p1 | p2) / (p3 | p4)
}

scores.dt <- lapply(.files, .fread) %>% bind_rows
plt <- .plot.all(scores.dt)
ggsave(filename = "result/step1/gene_stat.pdf", plot=plt, width=6, height=6)

valid.genes <-
    scores.dt[str_sub(hgnc, start=1, end=3) != "MT-" &
              nnz / nmax > 1e-4 &
              .cv > 1. &
              .sd > 1e-2,
              .(gene, batch)]

valid.genes <-
    valid.genes[, .(n = .N), by=.(gene)] %>%
    filter(n >= 14) %>% 
    select(gene) %>% 
    arrange %>%
    unique

fwrite(valid.genes, file="result/step1/features.tsv.gz", col.names = FALSE)
