## merged ##
stringtie --merge -p12 -G ~/Documents/DS/Whole_genome/Oryza_sativa.IRGSP-1.0.44.gtf *gtf -o merged.gtf

## compare #
gffcompare merged.gtf -r ~/Documents/DS/Whole_genome/Oryza_sativa.IRGSP-1.0.44.gtf -V
