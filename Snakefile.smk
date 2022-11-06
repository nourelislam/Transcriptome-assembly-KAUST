configfile: "config/config.yaml"

include: "workflow/rules/basics.smk"

# run rule all_gzip at the end of a run.
#
# Mappability takes a while...
rule all:
    input:
        #expand("01_rm/done__rm_{GE}", GE = config['genomes']),
        #expand("01_satDNA/done__satDNA_{GE}", GE = config['genomes']),
        #expand("01_busco/done__busco_{GE}", GE=config['genomes']),
        #expand("01_quast/{GE}", GE = config['genomes']),
        expand("01_EDTA/{GE}", GE = config['genomes']),
        #expand("01_assembly_stats/{GE}.stats", GE = config['genomes']),
        #expand("01_mappability/done__{GE}_bw", GE = config['genomes']),
        #expand("01_blastn_mito/done__{GE}_blastn_mito", GE = config['genomes'])
        #"01_naming/done__01_naming",
        #"02_general/done__02_general",
        #"03_completeness/done__03_completeness",
        #"04_annotation/done__04_annotation",
        ##"06_variation/done__06_variation",
        #"07_viz/done__07_viz",
        #"08_report/done__08_report",
        #"04_annotation/done__04_annotation"

rule hard_reset:
    shell:
        """
        rm logs/*
        rm -rdf 0*_*/*
        """
