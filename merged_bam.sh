## merged ##
stringtie --merge -p12 -G ~/Documents/DS/Whole_genome/Oryza_sativa.IRGSP-1.0.44.gtf *gtf -o merged.gtf

## compare #
gffcompare merged.gtf -r ~/Documents/DS/Whole_genome/Oryza_sativa.IRGSP-1.0.44.gtf -V

# feature count ##
featureCounts -T10 -a gffcmp.annotated.gtf -g transcript_id -o transcript_counts.txt ../*.bam  -G ~/Documents/DS/Whole_genome/whole_genome.fa
