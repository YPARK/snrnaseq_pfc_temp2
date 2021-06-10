library(data.table)
library(tidyverse)

read.assoc <- function(poly.file, .glm.fun, .term) {

    require(data.table)
    require(tidyverse)

    pheno.file <- "data/rosmap_phenotype.csv.gz"

    pheno.dt <-
        fread(pheno.file, na.strings = "-9") %>%
        mutate(iid = as.integer(projid)) %>%
        filter(!is.na(iid)) %>%
        mutate(apoe.e4 = if_else(str_ends(apoe_genotype, "4"), 1, 0)) %>%
        select(-apoe_genotype) %>%
        select(iid, pathoAD, apoe.e4, np_sqrt, nft_sqrt,
               cogn_ep_random_slope, educ, msex, age_death)

    print(poly.file)

    poly.dt <- fread(poly.file)

    .post.inter <- function(.glm) {
        ret <-
            summary(.glm) %>%
            coefficients %>%
            as.data.frame

        ret <- ret %>%
            mutate(term = rownames(ret)) %>%
            filter(term == .term) %>%
            select(Estimate, `Std. Error`) %>%
            as.list
    }

    .fun <- function(.dt){
        .lm.dt <- .dt[, .(iid, yhat)] %>%
            left_join(pheno.dt, by = "iid")
        .glm.fun(.lm.dt) %>% .post.inter
    }

    .ret <- poly.dt[!is.na(yhat), .fun(.SD),
                    by = .(hgnc_symbol, celltype)]
}

fit.ashr <- function(.dt) {
    .beta <- .dt$`Estimate`
    .se <- .dt$`Std. Error`
    .ash <- ashr::ash(.beta, .se)
    list(.ash$result$lfsr, .ash$result$svalue,
         .ash$result$lfdr, .ash$result$qvalue)
}

###################################################
## testing singleton terms without APOE genotype ##
###################################################

.np.fun0 <- function(.lm.dt) {
    lm(np_sqrt ~ yhat, data = .lm.dt)
}

assoc.file <- "result/step7/phenotype_np_twas_unc.txt.gz"

if(!file.exists(assoc.file)) {

    clust <- parallel::makeCluster(12)

    .files <- str_c("result/step7/chr", 1:22, "_poly.bed.gz")
    .out <- parallel::parLapply(clust, .files,
                                fun = read.assoc,
                                .glm.fun = .np.fun0,
                                .term = "yhat")

    out.dt <- do.call(rbind, .out)

    out.dt[,
           c("lfsr", "svalue", "lfdr", "qvalue") := fit.ashr(.SD),
           by = .(celltype)]

    fwrite(out.dt, assoc.file, sep="\t", quote=FALSE)
}

.nft.fun0 <- function(.lm.dt) {
    lm(nft_sqrt ~ yhat, data = .lm.dt)
}

assoc.file <- "result/step7/phenotype_nft_twas_unc.txt.gz"

if(!file.exists(assoc.file)) {

    clust <- parallel::makeCluster(12)

    .files <- str_c("result/step7/chr", 1:22, "_poly.bed.gz")
    .out <- parallel::parLapply(clust, .files,
                                fun = read.assoc,
                                .glm.fun = .nft.fun0,
                                .term = "yhat")

    out.dt <- do.call(rbind, .out)

    out.dt[,
           c("lfsr", "svalue", "lfdr", "qvalue") := fit.ashr(.SD),
           by = .(celltype)]

    fwrite(out.dt, assoc.file, sep="\t", quote=FALSE)
}

.ad.fun0 <- function(.lm.dt) {
    glm(pathoAD ~ yhat, family = "binomial", data = .lm.dt)
}

assoc.file <- "result/step7/phenotype_ad_twas_unc.txt.gz"

if(!file.exists(assoc.file)) {

    clust <- parallel::makeCluster(12)

    .files <- str_c("result/step7/chr", 1:22, "_poly.bed.gz")
    .out <- parallel::parLapply(clust, .files,
                                fun = read.assoc,
                                .glm.fun = .ad.fun0,
                                .term = "yhat")

    out.dt <- do.call(rbind, .out)

    out.dt[,
           c("lfsr", "svalue", "lfdr", "qvalue") := fit.ashr(.SD),
           by = .(celltype)]

    fwrite(out.dt, assoc.file, sep="\t", quote=FALSE)
}

#############################
## testing singleton terms ##
#############################

.np.fun <- function(.lm.dt) {
    lm(np_sqrt ~ apoe.e4 + yhat, data = .lm.dt)
}

assoc.file <- "result/step7/phenotype_np_twas.txt.gz"

if(!file.exists(assoc.file)) {

    clust <- parallel::makeCluster(12)

    .files <- str_c("result/step7/chr", 1:22, "_poly.bed.gz")
    .out <- parallel::parLapply(clust, .files,
                                fun = read.assoc,
                                .glm.fun = .np.fun,
                                .term = "yhat")

    out.dt <- do.call(rbind, .out)

    out.dt[,
           c("lfsr", "svalue", "lfdr", "qvalue") := fit.ashr(.SD),
           by = .(celltype)]

    fwrite(out.dt, assoc.file, sep="\t", quote=FALSE)
}

.nft.fun <- function(.lm.dt) {
    lm(nft_sqrt ~ apoe.e4 + yhat, data = .lm.dt)
}

assoc.file <- "result/step7/phenotype_nft_twas.txt.gz"

if(!file.exists(assoc.file)) {

    clust <- parallel::makeCluster(12)

    .files <- str_c("result/step7/chr", 1:22, "_poly.bed.gz")
    .out <- parallel::parLapply(clust, .files,
                                fun = read.assoc,
                                .glm.fun = .nft.fun,
                                .term = "yhat")

    out.dt <- do.call(rbind, .out)

    out.dt[,
           c("lfsr", "svalue", "lfdr", "qvalue") := fit.ashr(.SD),
           by = .(celltype)]

    fwrite(out.dt, assoc.file, sep="\t", quote=FALSE)
}

.ad.fun <- function(.lm.dt) {
    glm(pathoAD ~ apoe.e4 + yhat, family = "binomial", data = .lm.dt)
}

assoc.file <- "result/step7/phenotype_ad_twas.txt.gz"

if(!file.exists(assoc.file)) {

    clust <- parallel::makeCluster(12)

    .files <- str_c("result/step7/chr", 1:22, "_poly.bed.gz")
    .out <- parallel::parLapply(clust, .files,
                                fun = read.assoc,
                                .glm.fun = .ad.fun,
                                .term = "yhat")

    out.dt <- do.call(rbind, .out)

    out.dt[,
           c("lfsr", "svalue", "lfdr", "qvalue") := fit.ashr(.SD),
           by = .(celltype)]

    fwrite(out.dt, assoc.file, sep="\t", quote=FALSE)
}

###############################
## testing interaction terms ##
###############################

.np.inter.fun <- function(.lm.dt) {
    lm(np_sqrt ~ apoe.e4 + apoe.e4:yhat + yhat, data = .lm.dt)
}

assoc.file <- "result/step7/phenotype_np_interaction.txt.gz"

if(!file.exists(assoc.file)) {

    clust <- parallel::makeCluster(12)

    .files <- str_c("result/step7/chr", 1:22, "_poly.bed.gz")
    .out <- parallel::parLapply(clust, .files,
                                fun = read.assoc,
                                .glm.fun = .np.inter.fun,
                                .term = "apoe.e4:yhat")

    out.dt <- do.call(rbind, .out)

    out.dt[,
           c("lfsr", "svalue", "lfdr", "qvalue") := fit.ashr(.SD),
           by = .(celltype)]

    fwrite(out.dt, assoc.file, sep="\t", quote=FALSE)
}

.nft.inter.fun <- function(.lm.dt) {
    lm(nft_sqrt ~ apoe.e4 + apoe.e4:yhat + yhat, data = .lm.dt)
}

assoc.file <- "result/step7/phenotype_nft_interaction.txt.gz"

if(!file.exists(assoc.file)) {

    clust <- parallel::makeCluster(12)

    .files <- str_c("result/step7/chr", 1:22, "_poly.bed.gz")
    .out <- parallel::parLapply(clust, .files,
                                fun = read.assoc,
                                .glm.fun = .nft.inter.fun,
                                .term = "apoe.e4:yhat")

    out.dt <- do.call(rbind, .out)

    out.dt[,
           c("lfsr", "svalue", "lfdr", "qvalue") := fit.ashr(.SD),
           by = .(celltype)]

    fwrite(out.dt, assoc.file, sep="\t", quote=FALSE)
}

.ad.inter.fun <- function(.lm.dt) {
    glm(pathoAD ~ apoe.e4 + apoe.e4:yhat + yhat,
        family = "binomial", data = .lm.dt)
}

assoc.file <- "result/step7/phenotype_ad_interaction.txt.gz"

if(!file.exists(assoc.file)) {

    clust <- parallel::makeCluster(12)

    .files <- str_c("result/step7/chr", 1:22, "_poly.bed.gz")
    .out <- parallel::parLapply(clust, .files,
                                fun = read.assoc,
                                .glm.fun = .ad.inter.fun,
                                .term = "apoe.e4:yhat")

    out.dt <- do.call(rbind, .out)

    out.dt[,
           c("lfsr", "svalue", "lfdr", "qvalue") := fit.ashr(.SD),
           by = .(celltype)]

    fwrite(out.dt, assoc.file, sep="\t", quote=FALSE)
}
