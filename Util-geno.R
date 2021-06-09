options(stringsAsFactors = FALSE)

`%c%` <- function(mat, cols) mat[, cols, drop = FALSE]
`%r%` <- function(mat, rows) mat[rows, , drop = FALSE]
`%&%` <- function(a,b) paste(a, b, sep = '')

.unlist <- function(...) unlist(..., use.names = FALSE)
.zeros <- function(n1, n2) matrix(0, n1, n2)

sigmoid <- function(x) 1/(1+exp(-x))
logit <- function(x) log(x) - log(1- x)

.eval <- function(str) eval(parse(text = str))

log.msg <- function(...) {
    ss = as.character(date())
    cat(sprintf('[%s] ', ss), sprintf(...), '\n', file = stderr(), sep = '')
    flush(stderr())
}

load.data <- function(fileName){
    load(fileName)
    mget(ls()[ls() != "fileName"])
}

.mkdir <- function(...) {
    dir.create(..., recursive=TRUE, showWarnings=FALSE) 
}

RNORM <- function(d1, d2) {
    matrix(rnorm(d1 * d2), nrow = d1, ncol = d2)
}

#' @param plink.hdr
#' @param chr
#' @param plink.lb
#' @param plink.ub
#' @param temp.dir
subset.plink <- function(plink.hdr, chr, plink.lb, plink.ub, temp.dir) {

    require(dplyr)

    .error <- function(e) {
        print(e)
        log.msg('No QTL here!\n')
        return(NULL)
    }

    chr = unlist(chr)
    plink.lb = unlist(plink.lb) %>% as.integer()
    plink.ub = unlist(plink.ub) %>% as.integer()

    dir.create(temp.dir, recursive=TRUE, showWarnings=FALSE)

    .subset <- function(plink.hdr, chr, plink.lb, plink.ub, temp.dir) {

        chr.num = gsub(pattern = 'chr', replacement = '', chr) %>%
            as.integer()

        out.hdr = temp.dir %&% '/plink'

        plink.cmd = sprintf('rm -f %s.*; ./bin/plink --bfile %s --make-bed --geno 0.05 --maf 0.05 --chr %d --from-bp %d --to-bp %d --out %s', out.hdr, plink.hdr, chr.num, plink.lb, plink.ub, out.hdr)

        ## print(plink.cmd)

        system(plink.cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)

        plink = zqtl::read.plink(out.hdr)
        colnames(plink$BIM) = c('chr', 'rs', 'missing', 'snp.loc', 'plink.a1', 'plink.a2')
        colnames(plink$FAM) = c('fam', 'iid', 'father', 'mother', 'sex.code', '.pheno')

        if(any(is.logical(plink$BIM$plink.a1))) {
            plink$BIM$plink.a1 = 'T'
        }

        if(any(is.logical(plink$BIM$plink.a2))) {
            plink$BIM$plink.a2 = 'T'
        }
        return(plink)
    }

    plink = tryCatch(.subset(plink.hdr, chr, plink.lb, plink.ub, temp.dir),
                      error = .error)
    return(plink)
}

# Match the right-hand-side PLINK to the left-hand-side
#' @param .plink.lhs
#' @param .plink.rhs
match.plink <- function(.plink.lhs, .plink.rhs) {

    if(is.null(.plink.lhs) && !is.null(.plink.rhs)) { return(.plink.rhs) }
    if(!is.null(.plink.rhs) && is.null(.plink.rhs)) { return(.plink.lhs) }
    if(is.null(.plink.lhs) && is.null(.plink.rhs)) { return(NULL) }

    lhs = .plink.lhs$BIM %>%
        as_tibble() %>%
        mutate(lhs.pos = 1:n()) %>%
        rename(lhs.plink.a1 = plink.a1, lhs.plink.a2 = plink.a2) %>%
        select(-missing, -rs)

    rhs = .plink.rhs$BIM %>%
        as_tibble() %>%
        mutate(rhs.pos = 1:n()) %>%
        rename(rhs.plink.a1 = plink.a1, rhs.plink.a2 = plink.a2) %>%
        select(-missing, -rs)

    .matched = left_join(lhs, rhs, by = c("chr", "snp.loc")) %>% 
        na.omit()

    if(nrow(.matched) < 1) return(NULL)

    .matched = .matched %>%
        dplyr::filter(((lhs.plink.a1 == rhs.plink.a1) & (lhs.plink.a2 == rhs.plink.a2)) |
                      ((lhs.plink.a2 == rhs.plink.a1) & (lhs.plink.a1 == rhs.plink.a2))) %>%
        arrange(chr, snp.loc)

    if(nrow(.matched) < 1) return(NULL)

    ret.lhs = .plink.lhs
    ret.rhs = .plink.rhs

    ret.lhs$BIM = ret.lhs$BIM[.matched$lhs.pos, , drop = FALSE]
    ret.lhs$BED = ret.lhs$BED[ , .matched$lhs.pos, drop = FALSE]

    ret.rhs$BIM = ret.rhs$BIM[.matched$rhs.pos, , drop = FALSE]
    ret.rhs$BED = ret.rhs$BED[ , .matched$rhs.pos, drop = FALSE]

    flip.tab = ret.lhs$BIM %>%
        mutate(lhs.pos = 1:n()) %>%
        left_join(ret.rhs$BIM %>% mutate(rhs.pos = 1:n()),
                  by = c("chr", "snp.loc"),
                  suffix = c(".lhs", ".rhs")) %>%
        filter(plink.a1.lhs != plink.a1.rhs)

    ret.rhs$BIM[flip.tab$rhs.pos, ] <- ret.lhs$BIM[flip.tab$lhs.pos, ]

    flip.bed = ret.rhs$BED[, flip.tab$rhs.pos]
    zero.idx = flip.bed <= 0.5
    two.idx = flip.bed >= 1.5
    flip.bed[two.idx] = 0
    flip.bed[zero.idx] = 2
    ret.rhs$BED[, flip.tab$rhs.pos] = flip.bed

    list(lhs = ret.lhs, rhs = ret.rhs)
}

#' @param .svd SVD result
#' @param tau regularization parameter
take.DinvVt <- function(.svd, tau) {
    D = .svd$D
    U = .svd$U
    V.t = .svd$V.t
    .tau = tau / sqrt(nrow(U))
    sweep(V.t, 1, D + .tau, `/`)
}

#' @param .svd SVD result
#' @param zz z-score vector/matrix
#' @param tau regularization parameter
pred.coeff <- function(.svd, zz, tau) {
    W.t = take.DinvVt(.svd, tau)
    .beta = t(W.t) %*% (W.t %*% zz)
}

#' @param .svd SVD result
#' @param zz z-score vector/matrix
#' @param tau regularization parameter
pred.prs <- function(.svd, zz, tau) {
    if(is.null(.svd)) return(NULL)
    W.t = take.DinvVt(.svd, tau)
    .svd$U %*% (W.t %*% zz)
}

#' @param .svd SVD result
select.svd <- function(.svd, .valid){
    .svd$U <- .svd$U[, .valid, drop = FALSE]
    .svd$D <- .svd$D[.valid]
    .svd$V.t <- .svd$V.t[.valid, , drop = FALSE]
    .svd$V <- t(.svd$V.t)
    return(.svd)
}

#' @param .effect
melt.fqtl.effect <- function(.effect) {
    .names <- names(.effect)

    .take.j <- function(j){
        .mat <- .effect[[j]]
        colnames(.mat) <- 1:ncol(.mat)
        ret <- as.data.table(.mat) %>%
            mutate(.row = 1:n()) %>%
            melt(id.vars = ".row", variable.name = ".col", value.name = .names[j]) %>%
            mutate(.col = as.integer(.col))
        return(ret)
    }

    .list <- lapply(1:length(.effect), .take.j)

    if(length(.list) > 1) {
        ret <- Reduce(function(x, y) left_join(x, y), .list[-1], .list[[1]])
    } else {
        ret <- .list[[1]]
    }
    return(ret)
}

