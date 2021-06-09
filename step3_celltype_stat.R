argv <- commandArgs(trailingOnly = TRUE)

annot.file <- argv[1]  # "result/step3/bbknn.annot.gz" 
score.file <- argv[2]  # "result/step2/merged.score.gz"

library(tidyverse)
library(data.table)
library(patchwork)

.num.sci <- function(x) format(x, digits=2, scientific=TRUE)
.num.int <- function(x) format(as.integer(x), big.mark = ',')
.num.round <- function(x) round(x, digits=2)
`%&%` <- function(a,b) paste0(a,b)
.fread <- function(...) fread(..., header = FALSE)

.plot.hist <- function(.dt, log.x = FALSE, density.y = FALSE) {


    .dt <- .dt[order(.dt$p), ]
    .denom <- .dt[, .(ntot = sum(n)), by = .(celltype)]
    .dt <- merge(.dt, .denom, by = "celltype")
    .dt <- .dt[order(.dt$ntot), ]

    if(density.y){
        .aes <- aes(x = p, y = n/ntot, colour=celltype)
    } else {
        .aes <- aes(x = p, y = n, colour=celltype)
    }

    ret <- ggplot(.dt, .aes) +
        theme_classic() +
        theme(axis.title = element_text(size=6)) +
        theme(axis.text = element_text(size=5)) +
        theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
        theme(legend.title = element_text(size=6)) +
        theme(legend.text = element_text(size=5)) +
        theme(legend.key.height = unit(.15, "lines")) +
        geom_line(size=.5)

    if(log.x){
        .by <- (max(.dt$p) - min(.dt$p))/10
        ret <- ret +
            scale_x_continuous(limits = range(.dt$p),
                               breaks = seq(min(.dt$p), max(.dt$p), by = .by),
                               labels = function(x) .num.round(exp(x)))
    }

    if(density.y){
        ret <- ret + ylab("density")
        ret <- ret + scale_y_continuous(labels = function(x) .num.sci(x))
    } else {
        ret <- ret + ylab("count")        
        ret <- ret + scale_y_continuous(trans = "log10", labels = function(x) .num.sci(x))
    }

    ret <- ret + scale_colour_brewer(palette = "Paired")

    return(ret)
}

.plot.all <- function(annot.info.dt) {

    nnz.stat <- annot.info.dt[, .(n = .N), by = .(p = round(log(nnz)*5)/5, celltype)]

    p1 <- .plot.hist(nnz.stat, log.x=TRUE) +
        xlab("# Detected genes")

    mean.stat <- annot.info.dt[, .(n = .N), by = .(p = round(log(.mean)*5)/5, celltype)]
    
    p2 <- .plot.hist(mean.stat, log.x=TRUE) +
        xlab("Average expression")

    cv.stat <- annot.info.dt[, .(n = .N), by = .(p = round(log(.cv)*10)/10, celltype)]

    p3 <- .plot.hist(cv.stat, log.x=TRUE) +
        xlab("Coefficient of variation")

    sd.stat <- annot.info.dt[, .(n = .N), by = .(p = round(log(.sd)*10)/10, celltype)]

    p4 <- .plot.hist(sd.stat, log.x=TRUE) +
        xlab("Standard Deviation")

    (p1 | p2) / (p3 | p4)
}

dir.create("result/step3", recursive = TRUE, showWarnings = FALSE)

annot.dt <- .fread(annot.file)
colnames(annot.dt) <- c("cell", "celltype", "prob", "ln.prob")

info.dt <- .fread(score.file)
colnames(info.dt) <- c("cell","nnz", "nmax",".mean",".sd",".cv")

annot.info.dt <- merge(annot.dt, info.dt, by = "cell")

sum.dt <- annot.info.dt[, .(
    ncell = .num.int(.N),
    nnz = .num.int(mean(nnz)) %&% " (" %&% .num.int(sd(nnz)) %&% ")",
    mean = .num.round(mean(.mean)) %&% " (" %&% .num.round(sd(.mean)) %&% ")",
    cv = .num.round(mean(.cv)) %&% " (" %&% .num.round(sd(.cv)) %&% ")"),
    by = .(celltype)]

plt <- .plot.all(annot.info.dt)

ggsave(filename = "result/step3/celltype_stat.pdf", plot=plt, width=7, height=6)

write_tsv(sum.dt, file = "result/step3/celltype_stat.txt")
