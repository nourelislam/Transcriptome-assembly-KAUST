## Assembly the merged bam files ##
for file in *merged.sorted.bam; do sample=${file%-merged.sorted.bam}; stringtie -p 5 $file -o $sample.gtf -l $sample;done
