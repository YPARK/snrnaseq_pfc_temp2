library(tidyverse)
library(data.table)
library(patchwork)

.num.sci <- function(x) format(x, digits=2, scientific=TRUE)
.num.int <- function(x) format(x, big.mark = ',')
.num.round <- function(x) round(x, digits=2)

.fread <- function(xx) {
    batch <- basename(xx) %>% str_remove(".score.gz")
    ret <- fread(xx, header = FALSE)
    colnames(ret) <- c("cell","nnz", "nmax",".mean",".sd",".cv")
    return(ret)
}

.plot.hist <- function(.dt, log.x = FALSE, log.y = FALSE) {

    .dt <- .dt[order(.dt$p), ]

    ret <- ggplot(.dt, aes(x = p, y = n)) +
        theme_classic() +
        theme(axis.title = element_text(size=6)) +
        theme(axis.text = element_text(size=5)) +
        theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
        geom_bar(stat="identity", fill="gray60") +
        geom_line(colour = "gray20", size=.2)

    if(log.x){
        .by <- (max(.dt$p) - min(.dt$p))/10
        ret <- ret +
            scale_x_continuous(limits = range(.dt$p),
                               breaks = seq(min(.dt$p), max(.dt$p), by = .by),
                               labels = function(x) .num.round(exp(x)))
    }

    if(log.y) {
        ret <- ret + scale_y_continuous(labels = function(x) .num.sci(exp(x)))
    } else {
        ret <- ret + scale_y_continuous(labels = function(x) .num.int(as.integer(x)))
    }

    return(ret)
}

.plot.all <- function(qc.dt) {

    mt.stat <- qc.dt %>%
        mutate(p = 100 * (1 + tot.mt)/(1 + tot.mt + tot.auto))

    mt.stat <- mt.stat[, .(n = .N), by = .(p = round(10 * log(p))/10)]

    p1 <- .plot.hist(mt.stat, log.x=TRUE) +
        xlab("Mitochondrial activity (%)") +
        ylab("# cells")

    nnz.stat <- qc.dt[, .(n = .N), by = .(p = round(log(nnz.auto)*10)/10)]

    p2 <- .plot.hist(nnz.stat, log.x=TRUE) +
        xlab("# Detected genes (non-MT)") +
        ylab("# cells")

    mean.stat <- qc.dt[, .(n = .N), by = .(p = round(log(.mean.auto)*10)/10)]

    p3 <- .plot.hist(mean.stat, log.x=TRUE) +
        xlab("Average expression (non-MT)") +
        ylab("# cells")

    cv.stat <- qc.dt[, .(n = .N), by = .(p = round(log(.cv.auto)*50)/50)]

    p4 <- .plot.hist(cv.stat, log.x=TRUE) +
        xlab("Coefficient of variation (non-MT)") +
        ylab("# cells")

    sd.stat <- qc.dt[, .(n = .N), by = .(p = round(log(.sd.auto)*50)/50)]

    p5 <- .plot.hist(sd.stat, log.x=TRUE) +
        xlab("Standard Deviation (non-MT)") +
        ylab("# cells")

    p6 <-  ggplot() + theme_void() + geom_blank()

    (p1 | p2 | p3)/ (p4 | p5 | p6)

}

.mt <- .fread("result/step2/mt_genes.score.gz")
.auto <- .fread("result/step2/merged.score.gz")

qc.dt <- left_join(.auto, .mt, by = "cell", suffix=c(".auto",".mt")) %>%
    na.omit %>%
    mutate(tot.mt = nmax.mt * .mean.mt) %>%
    mutate(tot.auto = nmax.auto * .mean.auto) %>%
    mutate(p.mt = tot.mt/(tot.mt + tot.auto) * 100) %>%
    as.data.table

plt <- .plot.all(qc.dt)
ggsave(filename = "result/step2/cell_stat.pdf", plot=plt, width=9, height=6)
