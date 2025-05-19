

# Using fastq-multx downloaded through https://github.com/ExpressionAnalysis/ea-utils
# the reads were re-demultiplexed based on read 1 downloaded from SRA.
# The different runs are then split into the different libraries for subsequent alignment.

while read lines; do
	echo $lines
	$SRATOOLS_DIR/fastq-dump \
	-O $FASTQ_DIR/large/ \
	--split-files \
	$lines
done < /kimdata/livinlrg/atlas/external/large_SraAccList.csv

gzip $FASTQ_DIR/large/*.fastq

while read lines; do
echo $lines
fastq-multx \
	-m 2 \
	-B barcodes/${lines}_barcodes.txt \
	$FASTQ_DIR/${lines}_1.fastq.gz \
	$FASTQ_DIR/${lines}_2.fastq.gz \
	$FASTQ_DIR/${lines}_3.fastq.gz \
	-o $FASTQ_DIR/${lines}_%_1.fastq.gz \
	-o $FASTQ_DIR/${lines}_%_2.fastq.gz \
	-o $FASTQ_DIR/${lines}_%_3.fastq.gz > ${lines}.log.txt
done < /kimdata/livinlrg/atlas/external/large_SraAccList.csv