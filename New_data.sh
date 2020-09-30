## alignment ##
for R1 in *R1*; do R2=${R1//R1_001.fastq/R2_001.fastq}; sample=${R1%_R1_001.fastq}; hisat2 -x ../../Whole_genome/whole -q -1 $R1 -2 $R2 -S $sample.sam -t --dta --novel-splicesite-outfile $sample.tsv -p 8; rm $sample.sam; done
#### 2nd iteration ##
for R1 in *R1*; do R2=${R1//R1_001.fastq/R2_001.fastq}; sample=${R1%_R1_001.fastq}; hisat2 -x ../../Whole_genome/whole -q -1 $R1 -2 $R2 -S $sample.sam -t  --novel-splicesite-infile $sample.tsv -p 8; done
### denovo assembly ###
for file in *bam; do sample=${file%.sorted.bam}; stringtie -p 8 $file -o $sample.gtf -l $sample;done
### Assembly ##
stringtie --merge -p 8 M-19-4801_WT-S-6-3_D706-D505_L001.gtf M-19-4801_WT-S-6-3_D706-D505_L002.gtf M-19-4802_WT-S-6-4_D706-D506_L001.gtf M-19-4802_WT-S-6-4_D706-D506_L002.gtf -o WT-S-merged.gtf -l WT-S

### gff compare ##
gffcompare *merged.gtf -r ~/Documents/DS/Whole_genome/Oryza_sativa.IRGSP-1.0.44.gtf -R -V
