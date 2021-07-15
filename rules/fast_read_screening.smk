rule kraken:
    input:
        read1 = "samples/raw/{sample}_R1.fq",
        read2 = "samples/raw/{sample}_R2.fq",
    output:
        main= "results/{sample}_kraken_output.txt",
        report = "results/{sample}_kraken_report.txt",
    params:
        kraken_db = "/home/groups/CEDAR/callahro/reference_data/",
    threads: 8
    resources:
        mem_mb = 70000 #lambda wildcards, attempt: attempt *3000 + 70000
    conda:
        "../envs/kraken2.yaml"
    shell:
        """
        kraken2 --paired --threads {threads} --db {params.kraken_db} {input.read1} {input.read2} --use-names --report {output.report} --report-zero-counts --output {output.main}
        """

rule mash:
    input:
        read1 = "samples/raw/{sample}_R1.fq",
        read2 = "samples/raw/{sample}_R2.fq",
    params:
        mash_screen_db = "/home/groups/CEDAR/callahro/reference_data/refseq.genomes.k21s1000.msh"
    output:
        "results/{sample}_mash_output.txt"
    threads: 4
    resources:
        mem_mb = 12000
    conda:
        "../envs/trim.yaml"
    shell:
        """
        export PATH=$PATH:/home/groups/CEDAR/callahro/mash-Linux64-v2.2/
        mash screen -p {threads} -w {params.mash_screen_db} {input.read1} >> {output}
        mash screen -p {threads} -w {params.mash_screen_db} {input.read2} >> {output}
        """

