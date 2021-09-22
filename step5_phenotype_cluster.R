pheno.file <- "data/rosmap_phenotype.csv.gz"
col.file <- "result/step1/merged.cols.gz"

argv <- commandArgs(trailingOnly = TRUE)

pheno.file <- argv[1]
col.file <- argv[2]

dir.create("result/step5", recursive = TRUE, showWarnings = FALSE)

library(tidyverse)
library(data.table)
source("Util.R")

`%&%` <- function(a,b) paste0(a,b)
.fread <- function(...) fread(..., header = FALSE)
.scale <- function(x) { ret <- scale(x); ret[is.na(ret)] <- 0; ret }

.apoe.code <- function(x) {
    if(is.na(x)) return(NA) 
    if(x %in% c(24, 44, 34)) return(1)
    return(0) 
}

.pheno <- fread(pheno.file, header=TRUE, na.strings = "-9")
.pheno[.pheno == -9] <- NA
.pheno[.pheno <= -9] <- NA

.pheno[!is.na(apoe_genotype), apoe.e4 := sapply(apoe_genotype, .apoe.code)]
.pheno <- .pheno[, .(projid, msex, educ, pathoAD, age_death, apoe.e4, np_sqrt, nft_sqrt)]

.cols <- .fread(col.file)
.cols[, c("barcode", "projid") := tstrsplit(V1, "_", fixed=TRUE)]
.projid <- data.table(projid = as.integer(unique(.cols$projid)))

.pheno.proj <- left_join(.projid, .pheno, by = "projid")
.pheno.proj <- .pheno.proj %>% filter(!is.na(pathoAD))

.dat <- .pheno.proj %>% select(-projid) %>% as.matrix %>% .scale

.kmeans <- kmeans(.dat, centers=7, nstart = 500)

.melt <- .pheno.proj %>%
    mutate(k = .kmeans$cluster) %>% 
    melt.data.table(id.vars=c("projid", "k"))

.melt[, x := .scale(value), by = .(variable)]

.order <- .melt[, .(x = mean(x)), by = .(k, variable)] %>%
    rename(row = k, col = variable, weight = x) %>%
    order.pair()

.df <- .melt %>%
    mutate(pheno = factor(variable, .order$cols)) %>% 
    mutate(k = factor(k, .order$rows))

pheno.kmean.df <- .df

.aes <- aes(x = as.factor(projid), y = pheno, fill = pmin(pmax(x, -2), 2))

plt <- p1 <-
    ggplot(.df, .aes) +
    theme_classic() +
    facet_grid(.~ k, space="free", scales="free") +
    xlab(length(unique(.df$projid)) %&% " samples") +
    ylab("phenotypes/covariates") +
    theme(axis.text.x = element_blank()) +
    theme(axis.ticks.x = element_blank()) +
    theme(legend.key.size = unit(.2,"lines")) +
    theme(legend.title = element_text(size=6)) +
    theme(legend.text = element_text(size=6)) +
    geom_tile() +
    scale_fill_distiller("normalized\nweights",
                         palette = "RdBu", direction = -1)

ggsave(plot=plt, filename="result/step5/phenotype_cluster.pdf", width = 9, height = 2)

################################################################

out.dt <- .pheno.proj[, .(projid, pathoAD)]
out.dt[, AD := factor(pathoAD, c(0,1), c("Con", "AD"))]

col.dt <- copy(.cols)[, projid := as.integer(projid)]
col.dt <- merge(col.dt[projid %in% out.dt$projid], out.dt) %>%
    select(V1, AD)

fwrite(col.dt, file="result/step5/phenotype_cluster.txt.gz", sep=" ", col.names = FALSE)


################################################################

read.celltype.order <- function(.color.file = "data/brain_colors.txt") {

    .dt <- fread(.color.file)

    .order <- c("Ex-L2or3", "Ex-L4", "Ex-L5or6", "Ex-L5or6-CC",
                "In-PV", "In-SST", "In-SV2C", "In-VIP",
                "OPC", "Oligo", "Microglia", "Astro",
                "Endo", "Per", "Fib", "SMC")         

    data.table(celltype = .order) %>%
        left_join(.dt)
}

.annot.file <- "result/step4/subtype.annot.gz"

.ct.tab <- fread(.ct.file, header=FALSE, col.names = c("cell", "celltype"))

.ct.tab[, c("barcode", "projid") := tstrsplit(cell, "_", fixed=TRUE)]
.ct.tab[, projid := as.integer(projid)]

ct.prop.tab <- .ct.tab[, .(nct = .N), by = .(celltype, projid)]
ct.prop.tab[, ntot := sum(nct), by = .(projid)]
ct.prop.tab[, prop := nct/ntot, by = .(projid)]

.df <-
    pheno.kmean.df[, .(projid, k)] %>%
    unique %>% 
    left_join(ct.prop.tab) %>%
    na.omit %>% 
    mutate(celltype = factor(celltype, .ct.order$celltype))

plt <- p2 <-
    ggplot(.df, aes(y = prop, x = as.factor(projid), fill = celltype)) +
    theme_classic() +
    facet_grid(.~ k, space="free", scales="free") +
    ylab("Proportion of cell types") +
    xlab("407 samples") +
    geom_bar(stat="identity") +
    scale_fill_manual(values = .ct.order$hex) +
    theme(axis.text.x = element_blank()) +
    theme(axis.ticks.x = element_blank()) +
    theme(legend.key.size = unit(.2,"lines")) +
    theme(legend.title = element_text(size=6)) +
    theme(legend.text = element_text(size=6))

ggsave(plot=plt, filename="result/step5/phenotype_cluster_cf.pdf", width = 9, height = 2)

plt <- p1/p2

ggsave(plot=plt, filename="result/step5/phenotype_cluster_combined.pdf", width = 9, height = 4)

