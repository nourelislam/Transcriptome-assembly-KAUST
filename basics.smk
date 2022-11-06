
rule gunzip:
    input:
        config['genomeDIR'] + "{GE}" + config['genomeEXT'] + ".gz"
    output:
        config['genomeDIR'] + "{GE}" + config['genomeEXT']
    shell:
        """
        gunzip {input}
        """


rule gzip:
    input:
        config['genomeDIR'] + "{GE}" + config['genomeEXT']
    output:
        touch(config['genomeDIR'] + "done__gzip_{GE}")
    shell:
        """
        gzip {input}
        """

rule all_gzip:
    input:
        expand(config['genomeDIR'] + "{GE}" + config['genomeEXT'], GE=config['genomes'])

rule assembly_stats:
    input:
        fas = config['genomeDIR'] + "{GE}" + config['genomeEXT']
    output:
        "01_assembly_stats/{GE}.stats"
    params:
        "-t"
    conda:
        "../envs/assembly_stats.yaml"
    shell:
        """
        assembly-stats {params} {input} > {output}
        """


rule quast:
    input:
        fas = config['genomeDIR'] + "{GE}" + config['genomeEXT']
    output:
        directory("01_quast/{GE}")
    log:
        "01_quast/log__quast_{GE}"
    params:
        "01_quast"
    threads:
        config['threads_quast']
    conda:
        "../envs/quast.yaml"
    shell:
        """
        quast -o {output} -t {threads} {input} > {log} 2>&1
        """


rule busco:
    input:
        fas = config['genomeDIR'] + "{GE}" + config['genomeEXT']
    output:
        touch("01_busco/done__busco_{GE}")
    log:
        os.path.join(workflow.basedir, "01_busco/log__busco_{GE}")
    params:
        args = "-m genome -o busco -l metazoa",
        outdir = "01_busco"
    threads:
        config['threads_busco']
    conda:
        "../envs/busco.yaml"
    shell:
        """
        d={params.outdir}/{wildcards.GE}
        mkdir -p $d && 
        cd $d && 
        busco -i ../../{input.fas} -c {threads} {params.args} > {log} 2>&1
        """

rule EDTA:
    input:
        fasta = config['genomeDIR'] + '{GE}' + config['genomeEXT']
    output:
        "01_EDTA/{GE}"
    params:
        arg = "--species others --step all --sensitive 1 -anno 1 --evaluate 1"
    conda:
        "../envs/EDTA.yaml"
    threads:
        config['threads_EDTA']
    shell:
        """
        EDTA.pl --genome {input.fasta} -t {threads} {params.arg}
        """

rule rm:
    input:
        config['genomeDIR'] + "{GE}" + config['genomeEXT']
    output:
        touch("01_rm/done__rm_{GE}")
    params:
        outdir = "01_rm",
        lib = os.path.join(workflow.basedir, config['repeatmasker_lib']),
        args = "-a -xsmall -no_is -s -e rmblast -dir ./ -gff"
    log:
        "logs/01_rm_{GE}"
    threads: 
        config['threads_repeatmasker']
    conda:
        "../envs/repeatmasker.yaml"
    shell:
        """
        mkdir -p {params.outdir} && cd {params.outdir} && 
        RepeatMasker -pa {threads} {params.args} -lib {params.lib} ../{input} > ../{log} 2>&1
        """

# Mappability with gem
rule gem_index:
    input:
        config['genomeDIR'] + "{GE}" + config['genomeEXT']
    output:
        idx = directory("01_mappability/{GE}_gem_index")
    params:
        outdir = "01_mappability",
        idx = "01_mappability/gem_index",
        gem_bin = os.path.join(workflow.basedir, "workflow/scripts/bin")
    log:
        "logs/log__01_mappability_{GE}_gem_index"
    threads:
        config['threads_gem_index']
    shell:
        """
        export PATH={params.gem_bin}:$PATH
        mkdir -p {params.outdir} && 
        workflow/scripts/bin/gem-indexer -T {threads} -c dna -i {input} -o {output.idx} > {log} 2>&1
        """

rule gem:
    input:
        config['genomeDIR'] + "{GE}" + config['genomeEXT'],
        rules.gem_index.output
    output:
        touch("01_mappability/done__{GE}_gem")
    params:
        idx = "01_mappability/gem_index",
        out = "01_mappability/gem_out",
        gem_bin = os.path.join(workflow.basedir, "workflow/scripts/bin")
    log:
        "logs/log__01_mappability_{GE}_gem"
    threads:
        config['threads_gem']
    shell:
        """
        export PATH={params.gem_bin}:$PATH
        gem-mappability -T {threads} -I {params.idx}.gem -l 100 -o {params.out} > {log} 2>&1
        """
 
rule gem_to_bw:
    input:
        genome = config['genomeDIR'] + "{GE}" + config['genomeEXT'],
        gem = rules.gem.output
    output:
        touch("01_mappability/done__{GE}_bw")
    params:
        idx = "01_mappability/gem_index",
        gem_out = rules.gem.params.out,
        gem_bin = os.path.join(workflow.basedir, "workflow/scripts/bin")
    log:
        "logs/log__01_mappability_{GE}_bw"
    shell:
        """
        export PATH={params.gem_bin}:$PATH
        samtools faidx {input.genome} && 
        cat {input.genome}.fai | awk '{{print$1"\t"$2}}' > {input.genome}.chrsizes &&
        gem-2-wig -I {params.idx}.gem -i {params.gem_out}.mappability -o {params.gem_out} 2>> {log} && 
        cat {params.gem_out}.wig | awk -v OFS="\t" '{{print$1,$2,$5}}' > {params.gem_out}1.wig &&
        wigToBigWig {params.gem_out}1.wig {input.genome}.chrsizes {params.gem_out}.bw > {log} 2>&1
        """



###############################################################################

# BLAST SMALL LIST OF COMPLETE FLATWORM MITOCHONDRIAL GENOMES AGAINST ASSEMBLY

###############################################################################

rule make_genome_blastdb:
    input:
        genome = config['genomeDIR'] + "{GE}" + config['genomeEXT']
    output:
        touch("01_blastn_mito/done__{GE}_make_genome_blastdb")
    params:
        db_name = "{GE}",
        assembly_name = "{GE}" + config['genomeEXT'],
        db_prefix = "{GE}"
    log:
        os.path.join(workflow.basedir, "log/log__{GE}_make_genome_blastdb")
    conda:
        "../envs/blastn.yaml"
    shell:
        """
        mkdir -p tmp/blastDB_genome/
        cp {input} tmp/blastDB_genome/
        cd tmp/blastDB_genome/
        makeblastdb \
        -dbtype nucl \
        -input_type fasta \
        -out {params.db_prefix} \
        -title {params.db_name} \
        -in {params.assembly_name} > {log} 2>&1
        """


rule unpack_mito:
    input:
        mito_db = "workflow/resources/plan_mitogenomes_kraken_db.fa.gz"
    output:
        temp("tmp/mito_query.fa")
    shell:
        """
        zcat {input} > {output}
        """


rule blastn_mito_to_genome:
    input:
        rules.make_genome_blastdb.output,
        assembly = config['genomeDIR'] + "{GE}" + config['genomeEXT'],
        mito_db = rules.unpack_mito.output
    output:
        touch("01_blastn_mito/done__{GE}_blastn_mito"),
        tbl = report("01_blastn_mito/mito_to_{GE}_blast.out")
        #caption = "../report/04_06_blast_mito_to_genome.rst", category = "04_annotation", subcategory = "06_blastn_mito")
    params:
        db_name = "tmp/blastDB_genome/{GE}",
    threads:
        config['threads_blastn']
    conda:
        "../envs/blastn.yaml"
    shell:
        """
        blastn -db {params.db_name} \
        -num_threads {threads} \
        -query {input.mito_db} \
        -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qcovs sstart send evalue bitscore' > {output.tbl}
        """


rule satDNA:
    # This is simply a repeatmasker run with a sattelite repeat library
    input:
        config['genomeDIR'] + "{GE}" + config['genomeEXT']
    output:
        touch("01_satDNA/done__satDNA_{GE}")
    params:
        outdir = "01_satDNA",
        lib = os.path.join(workflow.basedir, config['satDNA_lib']),
        args = "-a -xsmall -no_is -s -e rmblast -dir ./ -gff"
    log:
        "logs/01_satDNA_{GE}"
    threads: 
        config['threads_repeatmasker']
    conda:
        "../envs/repeatmasker.yaml"
    shell:
        """
        mkdir -p {params.outdir} && cd {params.outdir} && 
        RepeatMasker -pa {threads} {params.args} -lib {params.lib} ../{input} > ../{log} 2>&1
        """
