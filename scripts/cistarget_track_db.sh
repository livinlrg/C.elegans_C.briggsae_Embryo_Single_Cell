
WORK_DIR=/kimdata/livinlrg/scAnalysis/WS290/Objects/scenic
CISTARGET_DIR=~/Software/create_cisTarget_databases
UCSC_DIR=~/Software/UCSC

cd $WORK_DIR

# WS290 genome acquired in cis_target_db.sh
# genome length in two column for bedgraph
samtools faidx $WORK_DIR/c_elegans.PRJNA13758.WS290.genomic.fa
cut -f1,2 $WORK_DIR/c_elegans.PRJNA13758.WS290.genomic.fa.fai > c_elegans.sizes.genome
sed -i -e 's/^/chr/' c_elegans.sizes.genome
# custom remove chr in front of mtdna and change to chrM

samtools faidx $WORK_DIR/c_briggsae.PRJNA10731.WS290.genomic.fa
cut -f1,2 $WORK_DIR/c_briggsae.PRJNA10731.WS290.genomic.fa.fai > c_briggsae.sizes.genome

# merge peaks that are overlapping
mkdir files_bedgraph_merge
for file in bedGraphFiles/*.bedgraph; do
	filename=$(basename $file)
	bedtools sort -i $file | bedtools merge -c 4 -d 0 -o max > ./files_bedgraph_merge/${filename%.*}.bedgraph
done

# bedgraph to bigWig
# bedGraphToBigWig myFile_sorted.bedgraph myChrom.sizes myBigWig.bw
mkdir bw_files
for file in files_bedgraph_merge/*; do
	filename=$(basename $file)
	
	$UCSC_DIR/bedGraphToBigWig $file c_elegans.sizes.genome ./bw_files/${filename%.*}.bw
done

# Create the region bed files using R
# Then use this script to create the database
$CISTARGET_DIR/create_cistarget_track_databases.py \
	-b $WORK_DIR/WS290_track_multipeak.bed \
	-T $WORK_DIR/bw_files \
	-d $WORK_DIR/track_file_names.txt \
	-o $WORK_DIR/elegans_WS290_cistarget/elegans_chipseq_WS290 \
	-a $UCSC_DIR/bigWigAverageOverBed \
	-t 20 \
	-g "_[0-9]+$"

# Create the region bed files using R
# Then use this script to create the database
$CISTARGET_DIR/create_cistarget_track_databases.py \
	-b $WORK_DIR/WS290_track.bed \
	-T $WORK_DIR/bw_files \
	-d $WORK_DIR/track_file_names.txt \
	-o $WORK_DIR/elegans_WS290_cistarget/elegans_chipseq_nonpeaked_WS290 \
	-a $UCSC_DIR/bigWigAverageOverBed \
	-t 20 \
	-g ^gene_id
