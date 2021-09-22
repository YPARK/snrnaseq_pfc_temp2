.pheno.file <- "data/ROSMAP_clinical.csv" 
.ct.file <- "result/step3/bbknn.annot.gz" 

argv <- commandArgs(trailingOnly = TRUE)

.pheno.file <- argv[1]  # e.g., .pheno.file <- "data/ROSMAP_clinical.csv" 
.ct.file <- argv[2]  # e.g., .ct.file <- "result/step3/bbknn.annot.gz" 
.out.file <- argv[3] # e.g., .out.file <- 

dir.create("result/step3", recursive = TRUE, showWarnings = FALSE)

library(tidyverse)
library(data.table)
library(patchwork)
source("Util.R")

read.celltype.order <- function(.color.file = "data/brain_colors.txt") {

    .dt <- fread(.color.file)

    .order <- c("Ex",
                "In",
                "OPC", "Oligo", "Microglia", "Astro",
                "Vasc")         

    data.table(celltype = .order) %>%
        left_join(.dt)
}

.num.sci <- function(x) format(x, digits=2, scientific=TRUE)
.num.int <- function(x) format(as.integer(x), big.mark = ',')
.num.round <- function(x) round(x, digits=2)

`%&%` <- function(a,b) paste0(a,b)
.fread <- function(...) fread(..., header = FALSE)
.scale <- function(x) { ret <- scale(x); ret[is.na(ret)] <- 0; ret }
   
.age.code <- function(x) {
    if(nchar(x) == 0) return(NA)
    if(x == "90+") return(2);
    if(as.numeric(x) > 80) return(1);
    return(0)
}

.apoe.code <- function(x) {
    if(is.na(x)) return(NA) 
    if(x == 44) return(2)
    if(x %in% c(24, 34)) return(1)
    return(0) 
}

.ad.dx <- function(b) {
    ret <- rep(NA, length(b))
    ret[!is.na(b)] <- 0
    ret[b > 3] <- 1
    return(ret)
}

.pheno <- fread(.pheno.file, header=TRUE)

.pheno[, age.death := sapply(age_death, .age.code)]
.pheno[, age.dx := sapply(age_first_ad_dx, .age.code)]
.pheno[, apoe.e4 := sapply(apoe_genotype, .apoe.code)]
.pheno[, pathoAD := sapply(braaksc, .ad.dx)]

.pheno <- .pheno[, .(projid, msex, educ, pmi, pathoAD, age.dx, age.death, apoe.e4, cogdx)]

.ct.tab <- fread(.ct.file, header=FALSE,
                 col.names = c("cell", "celltype", "prob", "ln.prob"))

.ct.tab[, c("barcode", "projid") := tstrsplit(cell, "_", fixed=TRUE)]
.ct.tab[, projid := as.integer(projid)]

ct.prop.tab <- .ct.tab[, .(nct = .N), by = .(celltype, projid)]
ct.prop.tab[, ntot := sum(nct), by = .(projid)]
ct.prop.tab[, prop := nct/ntot, by = .(projid)]

.df <- ct.prop.tab %>%
    dcast(projid ~ celltype, value.var = "prop", fill = 0) %>%
    left_join(.pheno, by = "projid")

xx <- .df %>% select(-sort(unique(ct.prop.tab$celltype)), -projid)
yy <- .df %>% select(sort(unique(ct.prop.tab$celltype)))

x.tib <- tibble(pheno = colnames(xx)) %>% mutate(x.col = 1:n())
y.tib <- tibble(celltype = colnames(yy)) %>% mutate(y.col = 1:n())

K <- ncol(yy)

xy.stat <-
    zqtl::calc.qtl.stat(xx, log(yy + 1/K)) %>% 
    left_join(y.tib, by = "y.col") %>%
    left_join(x.tib, by = "x.col") %>%
    select(-x.col, -y.col)    

.dt <- xy.stat %>%
    mutate(row = pheno, col = celltype, weight = -log10(p.val)) %>%
    order.pair(ret.tab = TRUE)

.dt.lab <- .dt %>% filter(p.val < 0.05)

plt <-
    ggplot(.dt, aes(y = row, x = col, fill = beta, size = weight)) +
    theme_bw() +
    theme(title = element_text(size=8)) +
    theme(axis.title = element_blank()) +
    theme(axis.text.x = element_text(angle=90, vjust=0, hjust=0)) +
    theme(legend.key.height = unit(.5, "lines")) +
    theme(legend.key.width = unit(.2, "lines")) +
    theme(legend.text = element_text(size=8)) +
    geom_point(pch = 22, stroke = .1) +
    geom_text(data = .dt.lab, size = 3, label = "*", col = "white") +
    scale_fill_gradient2(low="blue", high="red") +
    scale_x_discrete(position = "top") +
    scale_size_continuous(range = c(0,5),
                          breaks = seq(0,2,.3),
                          labels = function(x) .num.sci(10^(-x)))

ggsave(file = "result/step3/celltype_pheno.pdf", plot=plt,
       width = 4, height = 3)

fwrite(xy.stat, file = .out.file, sep = " ")

################################################################

.ggplot <- function(...) {
    ggplot(...) +
    theme_classic() +
    theme(title = element_text(size = 8)) +
    theme(legend.text = element_text(size = 8)) +
    theme(axis.text.x = element_blank()) +
    theme(axis.ticks.x = element_blank())
}

.ct <- read.celltype.order()

.co <-
    ct.prop.tab[celltype == "Ex", ] %>% 
    (function(x) { x[order(x$prop), .(projid)] }) %>%
    unlist

.dt <- ct.prop.tab %>%
    mutate(col = factor(projid, .co)) %>%
    mutate(celltype = factor(celltype, .ct$celltype))

p1 <-
    .ggplot(.dt, aes(x = col, y = prop, fill = celltype)) +
    ylab("Proportion of cell types") +
    xlab(length(.co) %&% " individuals") +
    geom_bar(stat="identity") +
    scale_fill_manual(values = .ct$hex)

.dt.tot <- .dt[, .(col, ntot)] %>% distinct

ntot.med <- median(.dt.tot$ntot)

.brk <- quantile(.dt.tot$ntot, c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 1))

p0 <-
    .ggplot(.dt.tot, aes(x = col, y = ntot)) +
    ylab("# cells") +
    xlab(length(.co) %&% " individuals") +
    geom_hline(yintercept = ntot.med, col = 2) +
    geom_bar(stat="identity") +
    scale_y_continuous(breaks=.brk, labels = .num.int)

plt <- (p0/p1) + plot_layout(heights=c(2,3))

ggsave(file = "result/step3/celltype_fraction.pdf", plot=plt,
       width = 9, height = 5)

