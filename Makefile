batch_names := $(shell ls -1 ./data/batches/*/filtered* | awk '!/Meta/' | awk -F'/' '{ print $$4 }')

all:

################################################################
# Step 1. Data integration
# Step 1a. Select features consistently expressed across all the batches
# Step 1b. Combine all the data

FEATURE_SCORES := $(foreach b, $(batch_names), result/step1/$(b).score.gz)
BATCH_ROWS := $(foreach b, $(batch_names), result/step1/$(b).rows.gz)
BATCH_COLS := $(foreach b, $(batch_names), result/step1/$(b).cols.gz)

STEP1 := $(FEATURE_SCORES) \
    result/step1/features.tsv.gz $(BATCH_COLS) result/step1/batches.txt \
    result/step1/merged.mtx.gz result/step1/merged.cols.gz result/step1/merged.rows.gz \
    result/step1/mt_genes.mtx.gz result/step1/mt_genes.tsv.gz result/step1/mt_genes.cols.gz \
    result/step1/gene.info.gz

step1: $(STEP1)

result/step1/gene.info.gz: step1_ensembl.R result/step1/merged.rows.gz
	Rscript --vanilla $^ $@

result/step1/merged.mtx.gz: result/step1/batches.txt result/step1/features.tsv.gz $(BATCH_COLS) $(BATCH_ROWS)
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	mmutil_merge_col result/step1/features.tsv.gz 500 result/step1/merged $$(cat $< | tr '\n' ' ')

result/step1/%.cols.gz: result/step1/%.columns.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	gzip -cd $^ | awk '{ print $$1 }' | gzip -c > $@

result/step1/merged.rows.gz: result/step1/merged.mtx.gz
	[ -f $@ ] || exit 1

result/step1/features.tsv.gz: step1_feature_scores.R $(FEATURE_SCORES)
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	Rscript --vanilla $<

result/step1/mt_genes.tsv.gz: $(BATCH_ROWS)
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	cat $^ | gzip -d | awk '/MT-/ { genes[$$1] ++ } END { for(g in genes) print g }' | gzip -c > $@

result/step1/mt_genes.mtx.gz: result/step1/batches.txt result/step1/mt_genes.tsv.gz $(BATCH_COLS) $(BATCH_ROWS)
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	mmutil_merge_col result/step1/mt_genes.tsv.gz 0 result/step1/mt_genes $$(cat $< | tr '\n' ' ')

# Calculate a feature score vector for each batch
result/step1/%.score.gz: data/batches/%/filtered_feature_bc_matrix/matrix.mtx.gz result/step1/%.rows.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	mmutil_row_score $^ $@

result/step1/%.rows.gz: data/batches/%/filtered_feature_bc_matrix/features.tsv.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	gzip -cd $^ | awk '{ print $$1 "_" $$2 }' | gzip -c > $@

result/step1/batches.txt: $(BATCH_COLS) $(BATCH_ROWS)
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	echo $(batch_names) | awk '{ for(j=1; j<=NF; ++j) { printf "data/batches/" $$j "/filtered_feature_bc_matrix/matrix.mtx.gz"; printf FS "result/step1/" $$j ".rows.gz"; print FS "result/step1/" $$j ".cols.gz"; } }' > $@

# Append ${project_id} to for each cell
result/step1/%.cols.gz: data/batches/%/filtered_feature_bc_matrix/barcodes.tsv.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	gzip -cd $< | awk -F '[-_]' -v M=data/batches/$*/filtered_feature_bc_matrix/projid.txt 'BEGIN{ while(("cat " M " " | getline line) > 0) { split(line,larr,"\t"); map[larr[1]] = larr[2] } }{ print $$1 "_" map[$$2] }' | gzip -c > $@


################################################################
# Step 2. Cell Q/C by mitochondrial mRNA content
# --> Found no need for Q/C

STEP2 := result/step2/merged.score.gz result/step2/mt_genes.score.gz

step2: $(STEP2) result/step2/cell_stat.pdf

result/step2/%.score.gz: result/step1/%.mtx.gz result/step1/%.cols.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	mmutil_col_score $^ $@

result/step2/cell_stat.pdf: step2_cell_scores.R $(STEP2)
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	Rscript --vanilla $<

################################################################
# Step 3. Cell type identification
# Step 3a. Generate feature scores

BROAD_DEF := In Oligo Astro Ex Vasc OPC Microglia

BBKNN_FILES := mtx factors svd_D svd_V svd_U cols

STEP3 := result/step3/marker_broad.txt.gz \
	$(foreach t, $(BBKNN_FILES), result/step3/bbknn.$(t).gz) \
	result/step3/bbknn.annot.gz \
	$(foreach ct, $(BROAD_DEF), $(foreach t, mtx cols, result/step3/sorted/$(ct).$(t).gz)) \
	result/step3/bbknn.umap.gz \
	result/step3/celltype_stat.txt \
	result/step3/celltype_stat.pdf \
	result/step3/celltype_pheno_stat.txt \
	result/step3/celltype_fraction.pdf \
	result/step3/celltype_pheno.pdf

step3: $(STEP3)

result/step3/bbknn.batch.gz: result/step1/merged.cols.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	gzip -cd $< | awk 'BEGIN { while("gzip -cd result/step1/merged.columns.gz" | getline line) { split(line, larr); batch[larr[1]] = larr[2]; } } { print batch[$$1] }' | gzip -c > $@

result/step3/bbknn.mtx.gz: result/step1/merged.mtx.gz result/step1/merged.cols.gz result/step3/bbknn.batch.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	mmutil_bbknn --mtx result/step1/merged.mtx.gz --col result/step1/merged.cols.gz --batch result/step3/bbknn.batch.gz --block_size 10000 --knn 150 --rank 50 --verbose --out result/step3/bbknn

result/step3/bbknn.factors.gz: result/step3/bbknn.mtx.gz
	[ -f $@ ] || exit 1
result/step3/bbknn.factors_doublet.gz: result/step3/bbknn.mtx.gz
	[ -f $@ ] || exit 1
result/step3/bbknn.cols.gz: result/step3/bbknn.mtx.gz
	[ -f $@ ] || exit 1
result/step3/bbknn.svd_D.gz: result/step3/bbknn.mtx.gz
	[ -f $@ ] || exit 1
result/step3/bbknn.svd_V.gz: result/step3/bbknn.mtx.gz
	[ -f $@ ] || exit 1
result/step3/bbknn.svd_U.gz: result/step3/bbknn.mtx.gz
	[ -f $@ ] || exit 1

result/step3/bbknn.annot.gz: result/step3/marker_broad.txt.gz $(foreach t, $(BBKNN_FILES), result/step3/bbknn.$(t).gz)
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	mmutil_annotate_col  --svd_u result/step3/bbknn.svd_U.gz --svd_d result/step3/bbknn.svd_D.gz --svd_v result/step3/bbknn.factors.gz --col result/step3/bbknn.cols.gz --row result/step1/merged.rows.gz --ann $< --em_iter 100 --batch_size 10000 --em_tol 1e-4 --verbose --out result/step3/bbknn

result/step3/marker_broad.txt.gz: step3_celltype_markers.R data/PsychENCODE.marker result/step1/merged.rows.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	Rscript --vanilla $^ result/step3/marker_broad-temp.txt
	cat result/step3/marker_broad-temp.txt | awk '/Endo/ || /Per/ { print $$1 "\tVasc" } !(/Endo/ || /Per/) { print $$0 }' | gzip -c > $@
	rm result/step3/marker_broad-temp.txt

#####################
# cell sorting part #
#####################

result/step3/sorted/%.select.gz: result/step3/bbknn.annot.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	gzip -cd $< | awk '$$2 == "$*"{ print $$1 }' | gzip -c > $@

result/step3/sorted/%.mtx.gz: result/step1/merged.mtx.gz result/step1/merged.cols.gz result/step3/sorted/%.select.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	if [ -f $@ ]; then touch $@; else mmutil_select_col $^ $(shell echo $@ | sed 's/.mtx.gz//g'); fi

result/step3/sorted/%.cols.gz: result/step3/sorted/%.mtx.gz
	[ -f $@ ] || exit 1

##############
# statistics #
##############

result/step3/bbknn.umap.gz: step3_celltype_umap.R result/step1/merged.columns.gz result/step3/bbknn.annot.gz result/step3/bbknn.factors.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	Rscript --vanilla $^ $@

result/step3/celltype_stat.txt: step3_celltype_stat.R result/step3/bbknn.annot.gz result/step2/merged.score.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	Rscript --vanilla $^

result/step3/celltype_stat.pdf: result/step3/celltype_stat.txt
	[ -f $@ ] || exit 1

result/step3/celltype_fraction.pdf: result/step3/celltype_pheno_stat.txt
	[ -f $@ ] || exit 1

result/step3/celltype_pheno_stat.txt: step3_celltype_phenotype.R data/ROSMAP_clinical.csv result/step3/bbknn.annot.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	Rscript --vanilla $^ $@

result/step3/celltype_pheno.pdf: result/step3/celltype_pheno_stat.txt
	[ -f $@ ] || exit 1

################################################################
# Step 4. Define subtypes and break down large cell types

STEP4_BBKNN := $(foreach t, Ex In Vasc, $(foreach x, mtx factors svd_D svd_V svd_U cols batch, result/step4/bbknn/$(t).$(x).gz))
STEP4_ANNOT := $(foreach t, In Ex Vasc, $(foreach s, bbknn, result/step4/subtype/$(s)_$(t).annot.gz))
STEP4_UMAP := $(foreach t, Ex In Vasc, result/step4/bbknn/$(t).umap.gz)

step4: $(STEP4_BBKNN) $(STEP4_ANNOT) $(STEP4_UMAP) result/step4/subtype.annot.gz

# % = {Ex, In}
result/step4/bbknn/%.batch.gz: result/step3/sorted/%.cols.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	gzip -cd $< | awk 'BEGIN{ while("gzip -cd result/step1/merged.columns.gz " | getline line) { split(line, larr); batch[larr[1]] = larr[2]; } } { print batch[$$1] }' | gzip -c > $@

result/step4/bbknn/%.mtx.gz: result/step3/sorted/%.mtx.gz result/step3/sorted/%.cols.gz result/step4/bbknn/%.batch.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	if [ -f $@ ]; then touch $@; else mmutil_bbknn --mtx result/step3/sorted/$*.mtx.gz  --col result/step3/sorted/$*.cols.gz --batch result/step4/bbknn/$*.batch.gz --knn 100 --rank 50 --out result/step4/bbknn/$*; fi

result/step4/bbknn/%.factors.gz: result/step4/bbknn/%.mtx.gz
	[ -f $@ ] || exit 1

result/step4/bbknn/%.factors_doublet.gz: result/step4/bbknn/%.mtx.gz
	[ -f $@ ] || exit 1

result/step4/bbknn/%.cols.gz: result/step4/bbknn/%.mtx.gz
	[ -f $@ ] || exit 1
result/step4/bbknn/%.svd_D.gz: result/step4/bbknn/%.mtx.gz
	[ -f $@ ] || exit 1
result/step4/bbknn/%.svd_V.gz: result/step4/bbknn/%.mtx.gz
	[ -f $@ ] || exit 1
result/step4/bbknn/%.svd_U.gz: result/step4/bbknn/%.mtx.gz
	[ -f $@ ] || exit 1

result/step4/subtype/Vasc.marker.gz: data/ROSMAP.VascularCells.celltypemMarkers.txt
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	cat $< | awk 'NR > 1 || $$3 > 1 { print $$1 "\t" $$7 }' > result/step4/subtype/Vasc-temp.txt
	Rscript --vanilla step4_celltype_markers.R result/step4/subtype/Vasc-temp.txt result/step1/merged.rows.gz $@
	rm result/step4/subtype/Vasc-temp.txt

result/step4/subtype/In.marker.gz: data/Velmeshev_etal_science_2019_markers.csv
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	cat $< | awk -F',' '$$1 ~ /^IN-/ && $$5 > 1 { print $$3 "\t" $$1 }' | awk '{ gsub("^IN","In",$$2); print }' > result/step4/subtype/In-temp.txt
	Rscript --vanilla step4_celltype_markers.R result/step4/subtype/In-temp.txt result/step1/merged.rows.gz $@
	rm result/step4/subtype/In-temp.txt

result/step4/subtype/Ex.marker.gz: data/Velmeshev_etal_science_2019_markers.csv
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	cat $< | awk -F',' '$$1 ~ /^L[1-9]+/ && $$5 > 1 { print $$3 "\tEx-" $$1 }' | awk '{ gsub("/","or"); print }' > result/step4/subtype/Ex-temp.txt
	Rscript --vanilla step4_celltype_markers.R result/step4/subtype/Ex-temp.txt result/step1/merged.rows.gz $@
	rm result/step4/subtype/Ex-temp.txt

result/step4/subtype/bbknn_%.annot.gz: result/step4/subtype/%.marker.gz $(foreach x, mtx factors svd_D svd_V svd_U cols, result/step4/bbknn/%.$(x).gz)
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	if [ -f $@ ]; then touch $@; else mmutil_annotate_col --svd_u result/step4/bbknn/$*.svd_U.gz --svd_d result/step4/bbknn/$*.svd_D.gz --svd_v result/step4/bbknn/$*.factors.gz --col result/step4/bbknn/$*.cols.gz --row result/step1/merged.rows.gz --ann $< --em_iter 100 --batch_size 5000 --em_tol 1e-8 --verbose --out $(shell echo $@ | sed 's/.annot.gz//'); fi

result/step4/subtype/%.label_names.gz: result/step4/subtype/%.annot.gz
	[ -f $@ ] || exit 1

result/step4/bbknn/%.umap.gz: step4_celltype_umap.R result/step1/merged.columns.gz result/step4/subtype/bbknn_%.annot.gz result/step4/bbknn/%.factors.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	if [ -f $@ ]; then touch $@; else Rscript --vanilla $^ $@; fi

############################
# final annotation results #
############################

CUTOFF := .9

result/step4/subtype.annot.gz: $(STEP4_ANNOT) 
	cat result/step4/subtype/bbknn_*.annot.gz result/step3/bbknn.annot.gz | gzip -d | awk '$$2 != "Ex"' | awk '$$2 != "In"' | awk '$$2 != "Vasc"' | awk '$$3 > $(CUTOFF){ print $$1 FS $$2 }' | gzip -c > $@

################################################################

SUBTYPES := Astro Endo Ex-L2or3 Ex-L4 Ex-L5or6 Ex-L5or6-CC Fib In-PV In-SST In-SV2C In-VIP Microglia OPC Oligo Per SMC

MTX_EXT := mtx rows cols
STEP5_DATA := $(foreach d, sorted filtered, $(foreach t, $(SUBTYPES), $(foreach ext, $(MTX_EXT),result/step5/$(d)/$(t).$(ext).gz)))
STEP5_EXT := ln_mu ln_mu_sd sum mean mu_cols
STEP5_PB := $(foreach t, $(SUBTYPES), $(foreach xx, $(STEP5_EXT), result/step5/pseudobulk/$(t).$(xx).gz))

step5: $(STEP5_DATA) $(STEP5_PB)

result/step5/sorted/%.select.gz: 
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	gzip -cd result/step4/subtype.annot.gz | awk '($$2 == "$*"){ print $$1 }' | gzip -c > $@

result/step5/sorted/%.mtx.gz: result/step5/sorted/%.select.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	mmutil_select_col result/step1/merged.mtx.gz result/step1/merged.cols.gz result/step5/sorted/$*.select.gz $(shell echo $@ | sed 's/.mtx.gz//g')

result/step5/sorted/%.rows.gz: result/step1/merged.rows.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	cp $< $@

result/step5/sorted/%.cols.gz: result/step5/sorted/%.mtx.gz
	[ -f $@ ] || exit 1

result/step5/filtered/%.mtx.gz: result/step5/sorted/%.mtx.gz result/step5/sorted/%.cols.gz result/step5/sorted/%.rows.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	mmutil_filter_row --mtx result/step5/sorted/$*.mtx.gz --row result/step5/sorted/$*.rows.gz --col result/step5/sorted/$*.cols.gz  --score MEAN --cutoff .05 --col_cutoff 500 --out result/step5/filtered/$*

result/step5/filtered/%.cols.gz: result/step5/filtered/%.mtx.gz
	[ -f $@ ] || exit 1

result/step5/filtered/%.rows.gz: result/step5/filtered/%.mtx.gz
	[ -f $@ ] || exit 1

############################################
# confounder identification and adjustment #
############################################

result/step5/pseudobulk/%.ln_mu.gz: result/step5/filtered/%.mtx.gz result/step5/filtered/%.cols.gz result/step5/pseudobulk/%.ind.gz result/step5/pseudobulk/%.lab.gz result/step5/pseudobulk/%.annot.gz
	rm -f result/step5/filtered/$*.mtx.gz.index
	mmutil_aggregate_col --mtx result/step5/filtered/$*.mtx.gz --col result/step5/filtered/$*.cols.gz --annot result/step5/pseudobulk/$*.annot.gz --lab result/step5/pseudobulk/$*.lab.gz --ind result/step5/pseudobulk/$*.ind.gz --verbose --out $(shell echo $@ | sed 's/.ln_mu.gz//g')

result/step5/pseudobulk/%.mu_cols.gz: result/step5/pseudobulk/%.ln_mu.gz
	[ -f $@ ] || exit 1

result/step5/pseudobulk/%.ln_mu_sd.gz: result/step5/pseudobulk/%.ln_mu.gz
	[ -f $@ ] || exit 1

result/step5/pseudobulk/%.sum.gz: result/step5/pseudobulk/%.ln_mu.gz
	[ -f $@ ] || exit 1

result/step5/pseudobulk/%.mean.gz: result/step5/pseudobulk/%.ln_mu.gz
	[ -f $@ ] || exit 1

result/step5/pseudobulk/%.ind.gz: result/step5/filtered/%.cols.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	gzip -cd $< | awk -F'_' '{ print $$2 }' | gzip -c > $@

result/step5/pseudobulk/%.annot.gz: result/step5/filtered/%.cols.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	gzip -cd $< | awk -vA=$* '{ print $$1 FS A }' | gzip -c > $@

result/step5/pseudobulk/%.lab.gz:  # just a single cell type
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	echo $* | gzip -c > $@

#####################
# prepare eQTL data #
#####################

BLOCK_SIZE=20
TEMPDIR=/broad/hptmp/ypp/snrnaseq_pfc/
NPC=15

step6: $(foreach ct, $(SUBTYPES), jobs/step6/eqtl_data_$(ct).jobs.gz)

step6_long: $(foreach ct, $(SUBTYPES), jobs/step6/eqtl_data_$(ct).long.gz)

step6_post:  $(foreach ct, $(SUBTYPES), $(foreach t, bed.gz bed.gz.tbi, result/step6/eqtl/data/$(ct).$(t)))

jobs/step6/eqtl_data_%.jobs.gz:
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	mkdir -p $(TEMPDIR)
	mkdir -p result/step6/eqtl/data/$*/
	awk -vB=$(BLOCK_SIZE) -vN=$(shell gzip -cd result/step5/filtered/$*.rows.gz | wc -l| sed 's/ //g') 'BEGIN{ M=int(N/B); if(M < N) M++; for(j=0; j<M; ++j) { printf (j*B + 1); for(b=(j*B + 1); b<(j+1)*B; ++b) printf "," (b + 1); printf "\n"; } }' | awk  -vEXE=step6_eqtl_data.R -vC=$* '{ printf "Rscript --vanilla %s result/step5/pseudobulk/%s.ln_mu.gz result/step5/pseudobulk/%s.mu_cols.gz result/step5/filtered/%s.rows.gz result/step1/gene.info.gz $(NPC) data/rosmap_geno/ $(TEMPDIR) %s result/step6/eqtl/data/$*/%04d.bed.gz\n", EXE, C, C, C, $$1, NR }' | gzip -c > $@
	[ $$(zless $@ | wc -l) -eq 0 ] || qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=4g -l h_rt=4:00:00 -b y -j y -N STEP6_$*_DATA -t 1-$$(zless $@ | wc -l) ./run_jobs.sh $@

jobs/step6/eqtl_data_%.long.gz: jobs/step6/eqtl_data_%.jobs.gz
	gzip -cd $< | awk 'system(" ! [ -f " $$NF " ]") == 0' | gzip > $@
	[ $$(zless $@ | wc -l) -eq 0 ] || qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=16g -l h_rt=24:00:00 -b y -j y -N STEP6_EQTL -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@


result/step6/eqtl/data/%.bed.gz:
	cat result/step6/eqtl/data/$*/*.bed.gz | bgzip -d | awk 'NR == 1 || $$1 != "#chr"' | awk '$$2 != "NA"' | sort -k1,1 -k2,2n | bgzip -c > $@

result/step6/eqtl/data/%.bed.gz.tbi: result/step6/eqtl/data/%.bed.gz
	tabix -p bed $< 


#######################
# run eQTL CV and PGS #
#######################
step7: jobs/step7/eqtl_sparse.jobs.gz

CHR := $(shell seq 1 22)
EXT := bed.gz bed.gz.tbi
STAT := stat poly

step7_post: $(foreach chr, $(CHR), $(foreach et, $(EXT), $(foreach tt, $(STAT), result/step7/chr$(chr)_$(tt).$(et)))) $(foreach et, $(EXT), result/step7/gene_pve.$(et))

jobs/step7/eqtl_sparse.jobs.gz:
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	mkdir -p $(TEMPDIR)
	mkdir -p result/step7/eqtl/
	awk -vDATA="result/step6/eqtl/data/" -vGENO="data/rosmap_geno/" -vTEMP=$(TEMPDIR) -vINFO="result/step1/gene.info.gz" -vN=$(shell gzip -cd result/step1/gene.info.gz | wc -l) -vEXE=step7_eqtl_sparse.R 'BEGIN{ for(j=1; j<=N; ++j) printf "Rscript --vanilla %s %s %d %s %s %s result/step7/eqtl/%05d\n", EXE, INFO, j, DATA, GENO, TEMP, j }' | gzip -c > $@
	[ $$(zless $@ | wc -l) -eq 0 ] || qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=2g -l h_rt=0:30:00 -b y -j y -N STEP7 -t 1-$$(zless $@ | wc -l) ./run_jobs.sh $@

jobs/step7/eqtl_sparse.long.gz: jobs/step7/eqtl_sparse.jobs.gz
	gzip -cd $< | awk 'system(" ! [ -f " $$NF "_stat.bed.gz ]") == 0' | gzip > $@
	[ $$(zless $@ | wc -l) -eq 0 ] || qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=16g -l h_rt=24:00:00 -b y -j y -N STEP7_EQTL -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

result/step7/chr%_stat.bed.gz:
	find result/step7/eqtl -name '*stat.bed.gz' -type f -exec gzip -cd {} + | awk 'NR == 1 { print $$0 } (NR > 1 && $$1 == $*) { printf("%d\t%d\t%d", $$1, $$2, $$3); for(j=4; j<=NF; ++j) printf "\t" $$j; printf "\n"; }' | sort -k1,1 -k2,2n --parallel=8 -S 8G | bgzip -c > $@

result/step7/%.bed.gz.tbi: result/step7/%.bed.gz
	tabix -p bed $< 

result/step7/chr%_poly.bed.gz:
	find result/step7/eqtl -name '*poly.bed.gz' -type f -exec gzip -cd {} + | awk 'NR == 1 { print $$0 } (NR > 1 && $$1 == $*) { printf("%d\t%d\t%d", $$1, $$2, $$3); for(j=4; j<=NF; ++j) printf "\t" $$j; printf "\n"; }' | sort -k1,1 -k2,2n --parallel=8 -S 8G | bgzip -c > $@

result/step7/gene_pve.bed.gz:
	find result/step7/eqtl -name '*poly.bed.gz' -type f -exec gzip -cd {} + | awk -F'\t' '{ k=$$1 FS $$2 FS $$3 FS $$5 FS $$6; if(!(k in data)) data[k] = $$9 FS $$11 FS $$12  } END { for(k in data) print k FS data[k] }' | sort -k1,1 -k2,2n --parallel=8 -S 8G | bgzip -c > $@

#########################
# GWAS polygenic scores #
#########################

GWAS := $(shell ls -1 data/GWAS/*.bed.gz | xargs -I file basename file .bed.gz)
LDFILE := LD.info.txt
NLD := $(shell tail -n+2 $(LDFILE) | wc -l | awk '{ print $$1 }')

step8: $(foreach gwas, $(GWAS), jobs/step8/gwas_$(gwas).jobs.gz)

step8_post: $(foreach gwas, $(GWAS), $(foreach et, $(EXT), result/step8/$(gwas).$(et)))

jobs/step8/gwas_%.jobs.gz: data/GWAS/%.bed.gz 
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	mkdir -p result/step8/$*/
	awk -vLDFILE=$(LDFILE) -vDATA="data/GWAS/$*.bed.gz" -vGENO="data/rosmap_geno/" -vTEMP="$(TEMPDIR)/gwas/$*" -vN=$(NLD) -vEXE=step8_gwas_pgs.R 'BEGIN{ for(j=1; j<=N; ++j) printf "Rscript --vanilla %s %s %d %s %s %s result/step8/$*/%04d.bed.gz\n", EXE, LDFILE, j, DATA, GENO, (TEMP "_" j), j }' | gzip -c > $@
	[ $$(zless $@ | wc -l) -eq 0 ] || qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=2g -l h_rt=0:30:00 -b y -j y -N $* -t 1-$$(zless $@ | wc -l) ./run_jobs.sh $@

result/step8/%.bed.gz:
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	ls -1 result/step8/$*/*.bed.gz | sort | xargs -I file gzip -cd file | awk 'NR == 1 || $$1 != "#CHR"' | bgzip -c > $@

result/step8/%.bed.gz.tbi: result/step8/%.bed.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	tabix -p bed $<

##########################
# Interpretation of GWAS #
##########################

GWAS := $(shell ls -1 result/step8/*.bed.gz | xargs -I file basename file .bed.gz)
TT := obs.stat null.stat obs.pve null.pve

step9: $(foreach gwas, $(GWAS), jobs/step9/interpretation_$(gwas).jobs.gz)
step9_long: $(foreach gwas, $(GWAS), jobs/step9/interpretation_$(gwas).long.gz)
step9_post: $(foreach gwas, $(GWAS), $(foreach tt, $(TT), result/step9/$(gwas).$(tt).gz))

jobs/step9/interpretation_%.jobs.gz: 
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	awk -vLDFILE=$(LDFILE) -vGWAS="result/step8/$*.bed.gz" -vEQTL="result/step7/" -vN=$(NLD) -vEXE=step9_gwas_pgs_twas.R 'BEGIN{ for(j=1; j<=N; ++j) printf "Rscript --vanilla %s %s %d %s %s result/step9/obs/$*/%04d\n", EXE, LDFILE, j, GWAS, EQTL, j }' | gzip -c > $@
	awk -vLDFILE=$(LDFILE) -vGWAS="result/step8/$*.bed.gz" -vEQTL="result/step7/" -vN=$(NLD) -vEXE=step9_gwas_pgs_twas_null.R 'BEGIN{ for(j=1; j<=N; ++j) printf "Rscript --vanilla %s %s %d %s %s result/step9/null/$*/%04d\n", EXE, LDFILE, j, GWAS, EQTL, j }' | gzip -c >> $@
	[ $$(zless $@ | wc -l) -eq 0 ] || qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=2g -l h_rt=0:30:00 -b y -j y -N $* -t 1-$$(zless $@ | wc -l) ./run_jobs.sh $@

jobs/step9/%.long.gz: jobs/step9/%.jobs.gz
	gzip -cd $< | awk 'system(" ! [ -f " $$NF ".stat.gz ]") == 0' | gzip > $@
	[ $$(zless $@ | wc -l) -eq 0 ] || qsub -P compbio_lab -o /dev/null -binding linear:1 -cwd -V -l h_vmem=8g -l h_rt=4:00:00 -b y -j y -N STEP9_$* -t 1-$$(zcat $@ | wc -l) ./run_jobs.sh $@

result/step9/%.obs.stat.gz:
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	gzip -cd result/step9/obs/$*/*.stat.gz | awk 'NR == 1 || $$1 != "gwas"' | gzip -c > $@

result/step9/%.null.stat.gz:
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	gzip -cd result/step9/null/$*/*.stat.gz | awk 'NR == 1 || $$1 != "gwas"' | gzip -c > $@

result/step9/%.obs.pve.gz:
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	gzip -cd result/step9/obs/$*/*.pve.gz | awk 'NR == 1 || $$1 != "celltype"' | gzip -c > $@

result/step9/%.null.pve.gz:
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	gzip -cd result/step9/null/$*/*.pve.gz | awk 'NR == 1 || $$1 != "celltype"' | gzip -c > $@
