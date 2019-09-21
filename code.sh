#downloading genome chromosomes from Ensembl, then concatenate all of them into one file####

cat Oryza_sativa.IRGSP-1.0.dna.chromosome.1.fa Oryza_sativa.IRGSP-1.0.dna.chromosome.2.fa Oryza_sativa.IRGSP-1.0.dna.chromosome.3.fa Oryza_sativa.IRGSP-1.0.dna.chromosome.4.fa Oryza_sativa.IRGSP-1.0.dna.chromosome.5.fa Oryza_sativa.IRGSP-1.0.dna.chromosome.6.fa Oryza_sativa.IRGSP-1.0.dna.chromosome.7.fa Oryza_sativa.IRGSP-1.0.dna.chromosome.8.fa Oryza_sativa.IRGSP-1.0.dna.chromosome.9.fa Oryza_sativa.IRGSP-1.0.dna.chromosome.10.fa Oryza_sativa.IRGSP-1.0.dna.chromosome.11.fa Oryza_sativa.IRGSP-1.0.dna.chromosome.12.fa > Oryza_sativa.IRGSP-1.0.dna.fa

#hisat2 -x hisat2_idx/Oryza_sativa -q -1 M_19_0207_33-C-3-1_D701-D501_L001_R1_001.fastq M_19_0207_33-C-3-1_D701-D501_L001_R2_001.fastq -S M_19_0207_33-C-3-1_novel.sam -t --summary-file M_19_0207_33-C-3-1_with_novel.summary.txt --novel-splicesite-infile M_19_0207_33-C-3-1.tsv -p 4

for R1 in *R1*; do R2=${R1//R1_001.fastq/R2_001.fastq}; sample=${R1%-D501_L001_R1_001.fastq}; hisat2 -x hisat2_idx/Oryza_sativa -q -1 $R1 -2 $R2 -S $sample.sam -t --summary-file $sample.summary.txt --dta --novel-splicesite-outfile $sample.tsv -p 5; done
# --dta for long anchor alignment 
#--novel-splicesite-outfile to use it on the seceond run of alignment as an input of splice sites

for R1 in *R1*; do R2=${R1//R1_001.fastq/R2_001.fastq}; sample=${R1%-D501_L001_R1_001.fastq}; hisat2 -x hisat2_idx/Oryza_sativa -q -1 $R1 -2 $R2 -S $sample.sam -t --summary-file $sample-withNovel.summary.txt --novel-splicesite-infile $sample.tsv -p 9; done
