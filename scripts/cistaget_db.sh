
WORK_DIR=/kimdata/livinlrg/scAnalysis/WS290/Objects/scenic
CISTARGET_DIR=~/Software/create_cisTarget_databases
cd $WORK_DIR

#WS290 genome
wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/releases/WS290/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS290.genomic.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/releases/WS290/species/c_briggsae/PRJNA10731/c_briggsae.PRJNA10731.WS290.genomic.fa.gz

gzip -d *.fa.gz

awk '/^>/{if (l!="") print l; print; l=0; next}{l+=length($0)}END{print l}' $WORK_DIR/c_elegans.PRJNA13758.WS290.genomic.fa > elegans_genome_length.txt

awk '/^>/{if (l!="") print l; print; l=0; next}{l+=length($0)}END{print l}' $WORK_DIR/c_briggsae.PRJNA10731.WS290.genomic.fa > briggsae_genome_length.txt

# C. elegans regions
bedtools getfasta -name -fi $WORK_DIR/c_elegans.PRJNA13758.WS290.genomic.fa -bed $WORK_DIR/WS290_elegans_500.bed -fo WS290_elegans_500_gene.fa
sed -i 's/^\([^:]*\)::.*/\1/' $WORK_DIR/WS290_elegans_500_gene.fa

bedtools getfasta -name -fi $WORK_DIR/c_elegans.PRJNA13758.WS290.genomic.fa -bed $WORK_DIR/WS290_elegans_2000_gene_body.bed -fo WS290_elegans_2000_gene_body.fa
sed -i 's/^\([^:]*\)::.*/\1/' $WORK_DIR/WS290_elegans_2000_gene_body.fa

# C. briggsae regions
bedtools getfasta -name -fi $WORK_DIR/c_briggsae.PRJNA10731.WS290.genomic.fa -bed $WORK_DIR/WS290_briggsae_500.bed -fo WS290_briggsae_500_gene.fa
sed -i 's/^\([^:]*\)::.*/\1/' $WORK_DIR/WS290_briggsae_500_gene.fa

bedtools getfasta -name -fi $WORK_DIR/c_briggsae.PRJNA10731.WS290.genomic.fa -bed $WORK_DIR/WS290_briggsae_2000_gene_body.bed -fo WS290_briggsae_2000_gene_body.fa
sed -i 's/^\([^:]*\)::.*/\1/' $WORK_DIR/WS290_briggsae_2000_gene_body.fa

# elegans
$CISTARGET_DIR/create_cistarget_motif_databases.py \
	-f $WORK_DIR/WS290_elegans_500_gene.fa \
	-g ^gene_id \
	-M $WORK_DIR/pwms_cel_motifs_cb/ \
	-m $WORK_DIR/joint_ele_motif_names.txt \
	-o $WORK_DIR/elegans_WS290_cistarget/elegans_WS290_500bp \
	-c ~/Software/cluster-buster/cbust \
	-t 20

$CISTARGET_DIR/create_cistarget_motif_databases.py \
	-f $WORK_DIR/WS290_elegans_2000_gene_body.fa \
	-g ^gene_id \
	-M $WORK_DIR/pwms_cel_motifs_cb/ \
	-m $WORK_DIR/joint_ele_motif_names.txt \
	-o $WORK_DIR/elegans_WS290_cistarget/WS290_elegans_2000bp_gene_body \
	-c ~/Software/cluster-buster/cbust \
	-t 20

# briggsae
$CISTARGET_DIR/create_cistarget_motif_databases.py \
        -g ^gene_id \
        -f $WORK_DIR/WS290_briggsae_500_gene.fa \
        -M $WORK_DIR/pwms_cbr_motifs_cb/ \
        -m $WORK_DIR/joint_cbr_motif_names.txt \
        -o $WORK_DIR/briggsae_WS290_cistarget/briggsae_WS290_500bp \
        -c ~/Software/cluster-buster/cbust \
	-t 20

$CISTARGET_DIR/create_cistarget_motif_databases.py \
        -g ^gene_id \
        -f $WORK_DIR/WS290_briggsae_2000_gene_body.fa \
        -M $WORK_DIR/pwms_cbr_motifs_cb/ \
        -m $WORK_DIR/joint_cbr_motif_names.txt \
        -o $WORK_DIR/briggsae_WS290_cistarget/WS290_briggsae_2000bp_gene_body \
        -c ~/Software/cluster-buster/cbust \
	-t 20

