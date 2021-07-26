# preparing the genome #
cat Oryza_sativa.IRGSP-1.0.dna.chromosome.1.fa Oryza_sativa.IRGSP-1.0.dna.chromosome.2.fa Oryza_sativa.IRGSP-1.0.dna.chromosome.3.fa Oryza_sativa.IRGSP-1.0.dna.chromosome.4.fa Oryza_sativa.IRGSP-1.0.dna.chromosome.5.fa Oryza_sativa.IRGSP-1.0.dna.chromosome.6.fa Oryza_sativa.IRGSP-1.0.dna.chromosome.7.fa Oryza_sativa.IRGSP-1.0.dna.chromosome.8.fa Oryza_sativa.IRGSP-1.0.dna.chromosome.9.fa Oryza_sativa.IRGSP-1.0.dna.chromosome.10.fa Oryza_sativa.IRGSP-1.0.dna.chromosome.11.fa Oryza_sativa.IRGSP-1.0.dna.chromosome.12.fa Oryza_sativa.IRGSP-1.0.dna.chromosome.Mt.fa Oryza_sativa.IRGSP-1.0.dna.chromosome.Pt.fa > Oryza_sativa.IRGSP.fa

## generate my own transcriptome reference sequence #
gffread -w Transcriptome.fa -g Oryza_sativa.IRGSP.fa Oryza_sativa.IRGSP-1.0.48.gtf

### STAR index #
STAR --runMode genomeGenerate --runThreadN 8 --genomeSAindexNbases 13 --genomeFastaFiles Oryza_sativa.IRGSP.fa --sjdbGTFfile Oryza_sativa.IRGSP-1.0.48.gtf

## STAR alignment ##
STAR --runThreadN 8 --genomeDir ../Whole_genome/GenomeDir --readFilesIn M_19_0211_33-F-3-1_D705-D501_L001_R1_001.fastq M_19_0211_33-F-3-1_D705-D501_L001_R2_001.fastq --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM --outFileNamePrefix STAR_BAM/M_19_0211
for R1 in *R1*; do R2=${R1//R1_001.fastq.gz/R2_001.fastq.gz}; sample=${R1%_R1_001.fastq.gz}; STAR --runThreadN 8 --genomeDir ../Whole_genome/GenomeDir --readFilesIn <(gunzip -c $R1 $R2) --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM --outFileNamePrefix STAR_align/$sample; done

# salmon quant #
salmon quant -t ../../Whole_genome/Transcriptome.fa -l A -a M_19_0211Aligned.toTranscriptome.out.bam -o salmon_quant -p 8
for file in *bam; do sample=${file%Aligned.toTranscriptome.out.bam}; salmon quant -t ../../Whole_genome/Transcriptome.fa -l A -a $file -o $sample -p 8;done
## indexing the col 1 and 4 with TSV##
awk '{print $1, $4}' OFS='\t' quant.sf > simple_quant.tsv
for file in *M_19_02*; do cd $file ;  awk '{print $1, $4}' OFS='\t' quant.sf > $file.tsv; cd ../ ; done
for file in *M_19_02*; do cd $file; rm *M_19_02*; done
for file in *M_19_02*; do cd $file; cat *tsv | cut -f2 > $file.tsv; cd ../; done
for file in *M_19_02*; do cd $file; mv *M_19_02* ../; cd ../; done
for file in *M_19_02*; do cd $file; rm *.tsv; cd ../; done
Paste *.tsv > expression_data.tsv  # edit the expression data to only have sample names col -1 col numbers by "nano" ###
### change col names #
sed -e '1s/TPM/33-C-3/' M_19_0207_33-C-3-1_D701-D501_L001 > 33-C-3.tsv

##generate events ##
suppa.py generateEvents -i ../Oryza_sativa.IRGSP-1.0.48.gtf -o Oryza_sativa.IRGSP -f ioe -e SE SS MX RI FL

##concatenate together ##
awk '
    FNR==1 && NR!=1 { while (/^<header>/) getline; }
    1 {print}
' *.ioe > ensembl_hg19.events.ioe


### after adding rownames on the first col of expression data by nano ###
awk '{print $1,$4,$5,$6,$7,$12,$13,$14,$15,$16,$17,$18,$19,$22,$23,$24,$25}' OFS='\t' expression_data.tsv > biotic_expression.tsv
awk '{print $1,$2,$3,$4,$5,$8,$9,$10, $11}' OFS='\t' biotic_expression.tsv > 33.EX.biotic.tsv
awk '{print $1,$6,$7,$12,$13,$14,$15,$16,$17}' OFS='\t' biotic_expression.tsv > WT.EX.biotic.tsv

## SUPPA2 event detection ##
suppa.py psiPerEvent --ioe-file ~/Documents/DS/Whole_genome/new_SUPPA/ensembl_hg19.events.ioe --expression-file expression_data.tsv -o events

##diffenretial experssion of local events ##
suppa.py diffSplice --method empirical --input ~/Documents/DS/Whole_genome/new_SUPPA/ensembl_hg19.events.ioe --psi 33.salt_event.psi WT.salt_event.psi --tpm 33.salt.tsv WT.salt.tsv --area 1000 --lower-bound 0.05 -gc -o salt_cond

#cold_Hours#
suppa.py diffSplice --method empirical --input ~/Documents/DS/Whole_genome/new_SUPPA/ensembl_hg19.events.ioe --psi 33.cold3_event.psi 33.cold6_event.psi --tpm 33.cold3.tsv 33.cold6.tsv --area 1000 --lower-bound 0.05 -gc -o cold_hours
suppa.py diffSplice --method empirical --input ~/Documents/DS/Whole_genome/new_SUPPA/ensembl_hg19.events.ioe --psi WT.cold3_event.psi WT.cold6_event.psi --tpm WT.cold3.tsv WT.cold6.tsv --area 1000 --lower-bound 0.05 -gc -o WT.cold_hours
##salt_hours ##
suppa.py diffSplice --method empirical --input ~/Documents/DS/Whole_genome/new_SUPPA/ensembl_hg19.events.ioe --psi 33.salt3_event.psi 33.salt6_event.psi --tpm 33.salt3.tsv 33.salt6.tsv --area 1000 --lower-bound 0.05 -gc -o 33.salt_hours
suppa.py diffSplice --method empirical --input ~/Documents/DS/Whole_genome/new_SUPPA/ensembl_hg19.events.ioe --psi WT.salt3_event.psi WT.salt6_event.psi --tpm WT.salt3.tsv WT.salt6.tsv --area 1000 --lower-bound 0.05 -gc -o WT.salt_hours
## control cond##
suppa.py diffSplice --method empirical --input ~/Documents/DS/Whole_genome/new_SUPPA/ensembl_hg19.events.ioe --psi WT.cont_event.psi 33.cont_event.psi --tpm WT.cont_expression.tsv 33.cont_expression.tsv --area 1000 --lower-bound 0.05 -gc -o cont_cond

#######

awk -F '[; : \t]' '{print $1,$2,$(NF-1), $NF}' OFS="\t" cold_cond.dpsi| head

grep "SE" simple_event.psi > SE_events.psi
awk -F '[; : \t]' '{print $1,$2,$(NF-11),$(NF-10),$(NF-9),$(NF-8),$(NF-7),$(NF-6),$(NF-5),$(NF-4),$(NF-3),$(NF-2),$(NF-1),$NF}' OFS="\t" SE_events.psi| head






