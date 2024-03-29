for file in *.sam;do sample=${file%.sam}; samtools view -@ 7 -hbo $sample.bam $sample.sam; rm $sample.sam ; samtools sort -@ 7 $sample.bam -o $sample.sorted.bam; rm $sample.bam ;done
### merging all gtf files with reference gtf ##
stringtie --merge -p10 -G ~/Documents/DS/Whole_genome/Oryza_sativa.IRGSP-1.0.44.gtf *gtf -o merged.gtf

# compare #
gffcompare merged.gtf -r ~/Documents/DS/Whole_genome/Oryza_sativa.IRGSP-1.0.44.gtf -V


### featurecount ##
featureCounts -T10 -a gffcmp.annotated.gtf -g gene_id -o gene_counts.txt ../../*.sorted.bam  -G ~/Documents/DS/Whole_genome/whole_genome.fa

cat gene_counts.txt | cut -f 1,7-30 > simple_gene_counts.txt
cat transcript_counts.txt | cut -f 1,7-30 > simple_transcript_counts.txt
