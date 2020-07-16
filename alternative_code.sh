## Assembly the merged bam files ##
for file in *merged.sorted.bam; do sample=${file%-merged.sorted.bam}; stringtie -p 5 $file -o $sample.gtf -l $sample;done

## Guided assembly ##
for file in *merged.sorted.bam; do sample=${file%-merged.sorted.bam}; stringtie -G ../gtf/Oryza_sativa.IRGSP-1.0.44.gtf -p 5 $file -o $sample.gtf -l $sample;done
