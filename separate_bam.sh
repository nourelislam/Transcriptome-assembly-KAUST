### merging all bam files ##
stringtie --merge -p 8 *gtf -o merged.gtf
###then stringtie assembly one more time using merged as the reference file ##
