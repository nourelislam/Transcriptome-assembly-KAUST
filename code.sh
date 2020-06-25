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
##for file in *.sam;do sample=${file%.sam}; samtools view -hbo $sample.bam $sample.sam;done
## sorting ##
for file in *.bam;do sample=${file%.bam}; samtools sort $sample.bam -o $sample.sorted.bam;done

##Adding the RG.ID 
picard=$CONDA_PREFIX/share/picard-2.20.8-0/picard.jar
java -jar $picard AddOrReplaceReadGroups I=M_19_0207_33-C-3-1_D701_sorted.bam O=M_19_0207_33-C-3-1_D701_rg_sorted.bam  RGID=@K00235:165:H2VVWBBXY:1 RGLB=M_19_0207_33-C-3-1 RGPL=illumina RGPU=@K00235:165:H2VVWBBXY:1.M_19_0207_33-C-3-1 RGSM=M_19_0207
#merge bam files#
java -jar $picard MergeSamFiles I=M_19_0219_WT-C-3-1_D702-D502_L002_rg_sorted.bam I=M_19_0220_WT-C-3-2_D703-D502_L002_rg_sorted.bam I=M_19_0225_WT-C-6-1_D708-D502_L003_rg_sorted.bam I=M_19_0226_WT-C-6-2_D709-D502_L003_rg_sorted.bam O=WT-C-merged.bam

##sort merged files##
samtools sort 33-F-merged.bam -o 33-F-merged.sorted.bam # as an example 

#stringtie assemly#
for file in *bam; do sample=${file%_sorted.bam}; stringtie -p 5 $file -o $sample.gtf -A $sample.tsv -l $sample;done

#stringtie merge#
stringtie --merge -p 4 M_19_0207_33-C-3-1_D701.gtf M_19_0208_33-C-3-2_D702.gtf M_19_0213_33-C-6-1_D708.gtf M_19_0214_33-C-6-2_D709-D501_L001.gtf -o 33-C-merged.gtf -l 33-C

#gff-compare#
gffcompare *merged.gtf -r Oryza_sativa.IRGSP-1.0.44.gtf -R -V
-R If -r was specified, this option causes gffcompare to ignore reference transcripts that are not overlapped by any transcript in one of input1.gtf,…,inputN.gtf. Useful for ignoring annotated transcripts that are not present in your RNA-Seq samples and thus adjusting the “sensitivity” calculation in the accuracy report written in the file.

stringtie -p 4 -B -G $GTF -o ballgown/M_19_0207_33-C-3-1.gtf M_19_0207_33-C-3-1_D701_sorted.bam



###stringtie B ###
stringtie -p5 -B -G gffcmp.combined.gtf -o ballgown/RNAseq_33-F/33-F.gtf -A ballgown/RNAseq_33-F/33-F-gene_abundance.tsv /media/nourelislam/Personal/splice_site-project/sorted/merged/33-F-merged.sorted.bam -l 33-F

##### featurecounts Estimation##
featureCounts -T4 -a /media/nourelislam/Personal/splice_site-project/sorted/gtf/merged/gffcmp.combined.gtf -g gene_id -o gene_counts.txt *.bam -J -G $FASTA

###### Extracting the gene & ref annotation from combined.gtf #####
grep -e "TCONS_00007307" -e "TCONS_00001556" -e "TCONS_00001668" -e "TCONS_00004700" -e "TCONS_00004802" -e "TCONS_00004807" -e "TCONS_00012118" -e "TCONS_00011013" -e "TCONS_00015975" -e "TCONS_00019568" -e "TCONS_00023659" -e "TCONS_00027566" -e "TCONS_00027605" -e "TCONS_00029910" -e "TCONS_00030288" -e "TCONS_00034537" -e "TCONS_00036057" -e "TCONS_00033431" -e "TCONS_00037598" -e "TCONS_00038309" -e "TCONS_00038404" -e "TCONS_00039046" -e "TCONS_00040926" -e "TCONS_00043218" -e "TCONS_00040133" -e "TCONS_00040191" -e "TCONS_00043945" -e "TCONS_00044035" -e "TCONS_00045039" -e "TCONS_00045040" -e "TCONS_00047673" -e "TCONS_00049518" -e "TCONS_00049576" -e "TCONS_00050661" -e "TCONS_00052824" -e "TCONS_00052901" -e "TCONS_00054407" -e "TCONS_00054540" -e "TCONS_00057639" -e "TCONS_00059914" -e "TCONS_00060095" -e "TCONS_00060216" -e "TCONS_00065908" -e "TCONS_00065122" -e "TCONS_00065666" -e "TCONS_00064822" gffcmp.combined.gtf


############################################
cat path.txt | while read id; do grep $id gffcmp.combined.gtf >> frs_genes.gtf; done
##filteration #
grep -w "transcript" frs_genes.gtf | awk '{print $10, $20}' > filt.gtf
cat path.txt | while read id; do grep $id filt.gtf >> annotated.gtf; done
#####################

seqkit stats *gz -T | csvtk pretty -t > lane1.tsv












