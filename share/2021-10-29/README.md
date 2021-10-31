### column names for the eqtl/${celltype}/${chr}.bed.gz

* `chr`: chromosome
* `start`: SNP location -1
* `stop`: SNP location
* `plink.a1`: A1 in plink genotype
* `plink.a2`: A2 in plink genotype
* `ensembl_gene_id`: ENSEMBL Gene ID
* `hgnc_symbol`: HGNC Gene symbol
* `celltype`: Cell type name
* `beta`: Univariate eQTL effect size
* `se`: Univariate eQTL standard error
* `p.val`: Univariate eQTL p-value
* `n`: Effective sample size
* `maf`: Minor Allele Frequency
* `r`: 5-fold cross validation R
* `theta`: Multivariate effect (Negative Binomial)
* `theta.sd`: Multivariate effect standard deviation
* `lodds`: Multivariate effect log odds ratio (PIP)
* `efdr`: Multivariate effect empirical FDR (within each gene)

### column names for the interaction/${phenotype}/${chr}.{geno,celltype}.bed.gz

* `chr`: chromosome
* `start`: SNP location -1
* `stop`: SNP location
* `plink.a1`: A1 in plink genotype
* `plink.a2`: A2 in plink genotype
* `pheno`: Interacting phenotype name
* `celltype`: Cell type name
* `ko`: Set `1` if this covarite is a knock-off control variate
* `ensembl_gene_id`: ENSEMBL Gene ID
* `hgnc_symbol`: HGNC Gene symbol
* `theta`: Multivariate effect (Negative Binomial)
* `theta.sd`: Multivariate effect standard deviation
* `lodds`: Multivariate effect log odds ratio (PIP)
* `efdr`: Multivariate effect empirical FDR (within each gene)
