## alignment ##
for R1 in *R1*; do R2=${R1//R1_001.fastq/R2_001.fastq}; sample=${R1%_R1_001.fastq}; hisat2 -x ../../Whole_genome/whole -q -1 $R1 -2 $R2 -S $sample.sam -t --dta --novel-splicesite-outfile $sample.tsv -p 8; rm $sample.sam; done
#### 2nd iteration ##
for R1 in *R1*; do R2=${R1//R1_001.fastq/R2_001.fastq}; sample=${R1%_R1_001.fastq}; hisat2 -x ../../Whole_genome/whole -q -1 $R1 -2 $R2 -S $sample.sam -t  --novel-splicesite-infile $sample.tsv -p 8; done
