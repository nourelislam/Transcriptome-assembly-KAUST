#downloading genome chromosomes from Ensembl, then concatenate all of them into one file####

cat Oryza_sativa.IRGSP-1.0.dna.chromosome.1.fa Oryza_sativa.IRGSP-1.0.dna.chromosome.2.fa Oryza_sativa.IRGSP-1.0.dna.chromosome.3.fa Oryza_sativa.IRGSP-1.0.dna.chromosome.4.fa Oryza_sativa.IRGSP-1.0.dna.chromosome.5.fa Oryza_sativa.IRGSP-1.0.dna.chromosome.6.fa Oryza_sativa.IRGSP-1.0.dna.chromosome.7.fa Oryza_sativa.IRGSP-1.0.dna.chromosome.8.fa Oryza_sativa.IRGSP-1.0.dna.chromosome.9.fa Oryza_sativa.IRGSP-1.0.dna.chromosome.10.fa Oryza_sativa.IRGSP-1.0.dna.chromosome.11.fa Oryza_sativa.IRGSP-1.0.dna.chromosome.12.fa > Oryza_sativa.IRGSP-1.0.dna.fa
#template#
#hisat2 -x hisat2_idx/Oryza_sativa -q -1 M_19_0207_33-C-3-1_D701-D501_L001_R1_001.fastq M_19_0207_33-C-3-1_D701-D501_L001_R2_001.fastq -S M_19_0207_33-C-3-1_novel.sam -t --summary-file M_19_0207_33-C-3-1_with_novel.summary.txt --novel-splicesite-infile M_19_0207_33-C-3-1.tsv -p 4

# 1st iteration of alignment ###
for R1 in *R1*; do R2=${R1//R1_001.fastq/R2_001.fastq}; sample=${R1%_R1_001.fastq}; hisat2 -x hisat2_idx/Oryza_sativa -q -1 $R1 -2 $R2 -S $sample.sam -t --summary-file $sample.summary.txt --dta --novel-splicesite-outfile $sample.tsv -p 9; done
# --dta for long anchor alignment 
#--novel-splicesite-outfile to use it on the seceond run of alignment as an input of splice sites

### 2nd iteration alignment with the splicesites file in ########
for R1 in *R1*; do R2=${R1//R1_001.fastq/R2_001.fastq}; sample=${R1%_R1_001.fastq}; hisat2 -x hisat2_idx/Oryza_sativa -q -1 $R1 -2 $R2 -S $sample.sam -t --summary-file $sample-withNovel.summary.txt --novel-splicesite-infile $sample.tsv -p 8; samtools view -hbo $sample.bam $sample.sam; rm $sample.sam; done
for file in *.sam;do sample=${file%.sam}; samtools view -hbo $sample.bam $sample.sam;done

#stringtie assemly#
for file in *bam; do sample=${file%_sorted.bam}; stringtie -p 5 $file -o $sample.gtf -A $sample.tsv -l $sample;done

#stringtie merge#
stringtie --merge -p 4 M_19_0207_33-C-3-1_D701.gtf M_19_0208_33-C-3-2_D702.gtf M_19_0213_33-C-6-1_D708.gtf M_19_0214_33-C-6-2_D709-D501_L001.gtf -o 33-C-merged.gtf -l 33-C

#gff-compare#
gffcompare *merged.gtf -r Oryza_sativa.IRGSP-1.0.44.gtf -R -V

stringtie -p 4 -B -G $GTF -o ballgown/M_19_0207_33-C-3-1.gtf M_19_0207_33-C-3-1_D701_sorted.bam

##Adding the RG.ID 
picard=$CONDA_PREFIX/share/picard-2.20.8-0/picard.jar
java -jar $picard AddOrReplaceReadGroups I=M_19_0207_33-C-3-1_D701_sorted.bam O=M_19_0207_33-C-3-1_D701_rg_sorted.bam  RGID=@K00235:165:H2VVWBBXY:1 RGLB=M_19_0207_33-C-3-1 RGPL=illumina RGPU=@K00235:165:H2VVWBBXY:1.M_19_0207_33-C-3-1 RGSM=M_19_0207
