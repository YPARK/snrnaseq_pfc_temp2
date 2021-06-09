library(data.table)
library(tidyverse)

assoc.file <- "result/step7/phenotype_association.txt.gz"

pheno.file <- "data/rosmap_phenotype.csv.gz"

pheno.dt <-
    fread(pheno.file, na.strings = "-9") %>%
    mutate(iid = as.integer(projid)) %>% 
    filter(!is.na(iid)) %>% 
    select(iid, pathoAD, apoe_genotype, np_sqrt, nft_sqrt, cogn_ep_random_slope, educ, msex, age_death)

pheno.melt <- melt(pheno.dt, id.vars = "iid")

if(!file.exists(assoc.file)) {

    read.assoc <- function(poly.file) {

        poly.dt <- fread(poly.file)
        
        .xx <- poly.dt %>%
            dcast(iid ~ hgnc_symbol + ensembl_gene_id + celltype, value.var = "yhat") %>%
            mutate(iid = as.integer(iid)) %>%
            as.data.table
        
        .yy <- .xx[, .(iid)] %>%
            left_join(pheno.dt, by = "iid") %>%
            mutate(apoe.e4 = if_else(str_ends(apoe_genotype, "4"), 1, 0)) %>% 
            select(-apoe_genotype) %>% 
            as.data.table
        
        x.mat <- as.matrix(as.data.frame(.xx)[, -1]) %>% scale
        y.mat <- as.matrix(as.data.frame(.yy)[, -1]) %>% scale
        
        x.tib <- tibble(x = colnames(.xx)[-1]) %>%
            mutate(x.col = 1:n()) %>%
            separate(x, c("hgnc_symbol", "ensGene", "celltype"), sep="[_]")
        
        y.tib <- tibble(pheno = colnames(.yy)[-1]) %>%
            mutate(y.col = 1:n())
        
        zqtl::calc.qtl.stat(x.mat, y.mat) %>%
            left_join(x.tib) %>%
            left_join(y.tib) %>%
            as.data.table
    }

    assoc.dt <-
        str_c("result/step7/chr", 1:22, "_poly.bed.gz") %>%
        lapply(FUN = read.assoc) %>%
        do.call(what=rbind) %>%
        select(-x.col, -y.col)

    fit.ashr <- function(.dt) {
        .beta <- .dt$beta
        .se <- .dt$se
        .ash <- ashr::ash(.beta, .se)
        list(.ash$result$lfsr, .ash$result$svalue,
             .ash$result$lfdr, .ash$result$qvalue)
    }

    assoc.dt[,
             c("lfsr", "svalue", "lfdr", "qvalue") := fit.ashr(.SD),
             by = .(pheno)]

    assoc.dt[, q.bh := p.adjust(p.val, "fdr"), by = .(pheno)]
    
    fwrite(assoc.dt, assoc.file, sep="\t", quote=FALSE)
}

gene.info <- fread("result/step1/gene.info.gz")
col.info <- fread("data/brain_colors.txt")

assoc.dt <- fread(assoc.file, sep="\t", header=TRUE) %>%
    left_join(gene.info)

.pheno <- "nft_sqrt"

.dt <- assoc.dt[pheno == .pheno]

  

.df.show <-
    .dt[order(.dt$p.val),
        head(.SD, 5),
        by = .(celltype)] %>%
    filter(qvalue < .2) %>%
    as.data.frame

.col <- .df.show %>%
    left_join(col.info) %>%
    select(celltype, hex) %>%
    unique %>% 
    arrange(celltype)

plt <-
    ggplot(.dt, aes(x=transcript_start.hg19, y=-log10(p.val))) +
    theme_bw() +
    theme(axis.text.x = element_blank()) +
    theme(legend.position = "none") +
    facet_grid(.~ as.integer(chr), space="free", scales="free") +
    geom_point(aes(colour = as.factor(as.integer(chr) %% 2)), stroke = 0, size = .5) +
    scale_colour_brewer(guide=FALSE) +
    geom_point(aes(fill=celltype), data = .df.show, stroke = .5, pch = 21) +
    scale_fill_manual(values = .col$hex) +
    ggrepel::geom_text_repel(aes(label = hgnc_symbol), data = .df.show, size=2)

ggsave("temp.pdf", plt, width = 10, height = 2)

g <- 9

.chr <- .df.show[g, c("chr")]
.tss <- .df.show[g, c("transcript_start.hg19")]
.tes <- .df.show[g, c("transcript_end.hg19")]

.gene <- .df.show[g, "hgnc_symbol"]
.qq <- str_c(.chr, ":", .tss, "-", .tes)

.file <- str_c("result/step7/chr", .chr, "_poly.bed.gz")

.poly <- fread(str_c("tabix -h ", .file, " ", .qq))

.poly.pheno.dt <- .poly %>%
    left_join(pheno.melt[variable==.pheno]) %>%
    group_by(celltype) %>%
    mutate(yhat = scale(yhat)) %>% 
    mutate(yobs = scale(value)) %>% 
    ungroup

ggplot(.poly.pheno.dt, aes(x = yhat, y = yobs)) +
    xlab(.gene) +
    ylab(.pheno) +
    theme_bw() +
    facet_wrap(~celltype, scales="free") +
    geom_point(stroke = 0, colour = "gray30") +
    geom_density2d(colour="green", size=.2) +
    geom_smooth(method = "lm", se=FALSE, size = .5, colour = "red")  +
    ggpubr::stat_cor()

