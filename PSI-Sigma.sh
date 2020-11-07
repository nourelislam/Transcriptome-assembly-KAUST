### STAR aligner ##
for R1 in *R1*; do R2=${R1//R1_001.fastq.gz/R2_001.fastq.gz}; sample=${R1%_R1_001.fastq.gz}; STAR --runThreadN 8 --outSAMtype BAM SortedByCoordinate --outFilterIntronMotifs RemoveNoncanonical --genomeDir ../Whole_genome/GenomeDir --readFilesIn <(gunzip -c $R1 $R2) --outFileNamePrefix PSI-Sigma/$sample; done

## indexing ##
for file in *.bam;do samtools index -@8 $file ;done
# sorting #
(grep "^#" Oryza_sativa.IRGSP-1.0.48.gtf; grep -v "^#" Oryza_sativa.IRGSP-1.0.48.gtf | sort -k1,1 -k4,4n) > Oryza_sativa.IRGSP.sorted.gtf
