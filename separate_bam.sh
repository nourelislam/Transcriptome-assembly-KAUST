### merging all bam files ##
stringtie --merge -p 8 *gtf -o merged.gtf
###then stringtie assembly one more time using merged as the reference file ##
for file in *bam; do sample=${file%.sorted.bam}; stringtie -p 8 -G GTF/merged.gtf $file -o $sample.guided.gtf -l $sample;done
