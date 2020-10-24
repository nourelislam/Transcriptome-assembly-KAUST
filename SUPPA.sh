# preparing the genome #
cat Oryza_sativa.IRGSP-1.0.dna.chromosome.1.fa Oryza_sativa.IRGSP-1.0.dna.chromosome.2.fa Oryza_sativa.IRGSP-1.0.dna.chromosome.3.fa Oryza_sativa.IRGSP-1.0.dna.chromosome.4.fa Oryza_sativa.IRGSP-1.0.dna.chromosome.5.fa Oryza_sativa.IRGSP-1.0.dna.chromosome.6.fa Oryza_sativa.IRGSP-1.0.dna.chromosome.7.fa Oryza_sativa.IRGSP-1.0.dna.chromosome.8.fa Oryza_sativa.IRGSP-1.0.dna.chromosome.9.fa Oryza_sativa.IRGSP-1.0.dna.chromosome.10.fa Oryza_sativa.IRGSP-1.0.dna.chromosome.11.fa Oryza_sativa.IRGSP-1.0.dna.chromosome.12.fa Oryza_sativa.IRGSP-1.0.dna.chromosome.Mt.fa Oryza_sativa.IRGSP-1.0.dna.chromosome.Pt.fa > Oryza_sativa.IRGSP.fa

## generate my own transcriptome reference sequence #
gffread -w Transcriptome.fa -g Oryza_sativa.IRGSP.fa Oryza_sativa.IRGSP-1.0.48.gtf

### STAR index #
STAR --runMode genomeGenerate --runThreadN 8 --genomeSAindexNbases 13 --genomeFastaFiles Oryza_sativa.IRGSP.fa --sjdbGTFfile Oryza_sativa.IRGSP-1.0.48.gtf

## STAR alignment ##
STAR --runThreadN 8 --genomeDir ../Whole_genome/GenomeDir --readFilesIn M_19_0211_33-F-3-1_D705-D501_L001_R1_001.fastq M_19_0211_33-F-3-1_D705-D501_L001_R2_001.fastq --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM --outFileNamePrefix STAR_BAM/M_19_0211

# salmon quant #
salmon quant -t ../../Whole_genome/Transcriptome.fa -l A -a M_19_0211Aligned.toTranscriptome.out.bam -o salmon_quant -p 8
