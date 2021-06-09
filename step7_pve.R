
assoc.file <- "result/step7/gene_pve_stat.txt.gz"

if(!file.exists(assoc.file)) {

    pheno.file <- "data/rosmap_phenotype.csv.gz"

    pheno.dt <-
        fread(pheno.file, na.strings = "-9") %>%
        mutate(iid = as.integer(projid)) %>% 
        filter(!is.na(iid)) %>% 
        select(iid, pathoAD, apoe_genotype, np_sqrt, nft_sqrt, cogn_ep_random_slope, educ, msex, age_death)

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

}
