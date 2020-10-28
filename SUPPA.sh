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
for file in *M_19_02*; do cd $file ;awk '{print $1, $4}' OFS='\t' quant.sf > simple_quant.tsv; cd ../; done
for file in *M_19_02*; do cd $file ;  awk '{print $1, $4}' OFS='\t' quant.sf > $file.tsv; cd ../ ;done
### change col names #
sed -i -e '1s/TPM/WT-F-6/' simple_quant.tsv

##generate events ##
suppa.py generateEvents -i ../Oryza_sativa.IRGSP-1.0.48.gtf -o Oryza_sativa.IRGSP -f ioe -e SE SS MX RI FL

##concatenate together ##
awk '
    FNR==1 && NR!=1 { while (/^<header>/) getline; }
    1 {print}
' *.ioe > ensembl_hg19.events.ioe

## SUPPA2 event detection ##
suppa.py psiPerEvent --ioe-file ensembl_hg19.events.ioe --expression-file ~/Documents/DS/fastq/STAR_BAM/salmon_quant/simple_quant.tsv -o
