
SCENIC_PATH=/kimdata/livinlrg/scAnalysis/WS290/Objects/scenic
SCRATCH_PATH=/scratch/livinlrg/WS290 # on kimclust52

# elegans test run
arboreto_with_multiprocessing.py \
	-m grnboost2 \
	--seed 777 \
	--num_workers 32 \
	-o $SCENIC_PATH/test/WS290_elegans_pyscenic_grn_tpm.grnboost2.csv \
	$SCENIC_PATH/cel_expression.csv \
	$SCENIC_PATH/cel_tfs.txt \
	-t

pyscenic ctx \
	--expression_mtx_fname $SCENIC_PATH/cel_expression.csv  \
	--mask_dropouts \
	--num_workers 32 \
	--annotations_fname $SCENIC_PATH/cel_chipseq_motif_table.tbl \
	-o $SCENIC_PATH/test/WS290_elegans_tpm_reg.csv \
	$SCENIC_PATH/test/WS290_elegans_pyscenic_grn_tpm.grnboost2.csv \
	$SCENIC_PATH/elegans_WS290_cistarget/elegans_WS290_500bp.genes_vs_motifs.rankings.feather \
	$SCENIC_PATH/elegans_WS290_cistarget/WS290_elegans_2000bp_gene_body.genes_vs_motifs.rankings.feather \
	$SCENIC_PATH/elegans_WS290_cistarget/elegans_chipseq_WS290.genes_vs_tracks.rankings.feather \
	-t
	
pyscenic ctx \
	--expression_mtx_fname $SCENIC_PATH/cel_expression.csv  \
	--mask_dropouts \
	--num_workers 32 \
	--annotations_fname $SCENIC_PATH/cel_chipseq_motif_table.tbl \
	-o $SCENIC_PATH/test/WS290_elegans_tpm_reg.csv \
	$SCENIC_PATH/test/WS290_elegans_pyscenic_grn_tpm.grnboost2.csv \
	$SCENIC_PATH/elegans_WS290_cistarget/elegans_WS290_500bp.genes_vs_motifs.rankings.feather \
	$SCENIC_PATH/elegans_WS290_cistarget/WS290_elegans_2000bp_gene_body.genes_vs_motifs.rankings.feather \
	$SCENIC_PATH/elegans_WS290_cistarget/elegans_chipseq_nonpeaked_WS290.genes_vs_tracks.rankings.feather \
	-t
	
# briggsae test run
arboreto_with_multiprocessing.py \
	-m grnboost2 \
	--seed 777 \
	--num_workers 32 \
	-o $SCENIC_PATH/Joint/WS260_briggsae_pyscenic_grn_tpm.grnboost2.csv \
	$SCENIC_PATH/cbr_expression.csv \
	$SCENIC_PATH/cbr_tfs.txt \
	-t

pyscenic ctx \
	--expression_mtx_fname $SCENIC_PATH/cel_expression.csv  \
	--mask_dropouts \
	--num_workers 32 \
	--annotations_fname $SCENIC_PATH/Cel_motif_names.tbl \
	-o $SCENIC_PATH/WS290_elegans_tpm_reg.csv \
	$SCENIC_PATH/WS290_elegans_pyscenic_grn_tpm.grnboost2.csv \
	$SCENIC_PATH/elegans_WS290_cistarget/briggsae_WS290_500bp.genes_vs_motifs.rankings.feather \
	$SCENIC_PATH/elegans_WS290_cistarget/WS290_briggsae_2000bp_gene_body.genes_vs_motifs.rankings.feather \
	-t

## run snakemake for cel and for cbr
# then aggregate the outputs
# and run aucell here

pyscenic aucell \
	$SCRATCH_PATH/WS290_elegans_expression.loom \
	$SCENIC_PATH/elegans_motif_chip.aggregate.csv \
	--output $SCENIC_PATH/WS290_elegans.cel_grn.pyscenic.csv \
	--num_workers 10
	
pyscenic aucell \
	$SCRATCH_PATH/WS290_elegans_expression.loom \
	$SCENIC_PATH/briggsae_motif.aggregate.csv \
	--output $SCENIC_PATH/WS290_elegans.cbr_grn.pyscenic.csv \
	--num_workers 10
	
pyscenic aucell \
	$SCRATCH_PATH/WS290_briggsae_expression.loom \
	$SCENIC_PATH/elegans_motif_chip.aggregate.csv \
	--output $SCENIC_PATH/WS290_briggsae.cel_grn.pyscenic.csv \
	--num_workers 10
	
pyscenic aucell \
	$SCRATCH_PATH/WS290_briggsae_expression.loom \
	$SCENIC_PATH/briggsae_motif.aggregate.csv \
	--output $SCENIC_PATH/WS290_briggsae.cbr_grn.pyscenic.csv \
	--num_workers 10
	
	
	
