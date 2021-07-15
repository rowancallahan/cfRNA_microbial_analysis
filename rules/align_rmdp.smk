#TODO make sure that this works for files that are zipped

rule trimming:
    input:
        fwd = "samples/raw/{sample}_R1.fq",
        rev = "samples/raw/{sample}_R2.fq"
    output:
        fwd = temp("samples/trimmed/{sample}_R1_t.fq"),
        rev = temp("samples/trimmed/{sample}_R2_t.fq"),
        single = temp("samples/trimmed/{sample}_R1_singletons.fq"),
    threads: 2 
    resources:
        mem_mb = 12000
    run:
        sickle = config["sickle_tool"]

        shell("{sickle} pe -f {input.fwd} -r {input.rev}  -l 40 -q 20 -t sanger  -o {output.fwd} -p {output.rev} -s {output.single} &> {input.fwd}.log")

#TODO ensure that the forward and reverse fastq trimmed files are temporary
rule fastqc:
    input:
        fwd = "samples/trimmed/{sample}_R1_t.fq",
        rev = "samples/trimmed/{sample}_R2_t.fq"
    output:
        fwd = "samples/fastqc/{sample}/{sample}_R1_t_fastqc.zip",
        rev = "samples/fastqc/{sample}/{sample}_R2_t_fastqc.zip"
    threads: 2 
    resources:
        mem_mb = 12000
    conda:
        "../envs/fastqc.yaml"
    shell:
        """fastqc --outdir samples/fastqc/{wildcards.sample} --extract  -f fastq {input.fwd} {input.rev}"""

#TODO does this need to be a part of the pipeline for the microbes?
rule fastqscreen:
    input:
        fwd = "samples/trimmed/{sample}_R1_t.fq",
        rev = "samples/trimmed/{sample}_R2_t.fq"
    output:
        "samples/fastqscreen/{sample}/{sample}_R1_t_screen.html",
        "samples/fastqscreen/{sample}/{sample}_R1_t_screen.png",
        "samples/fastqscreen/{sample}/{sample}_R1_t_screen.txt",
        "samples/fastqscreen/{sample}/{sample}_R2_t_screen.html",
        "samples/fastqscreen/{sample}/{sample}_R2_t_screen.png",
        "samples/fastqscreen/{sample}/{sample}_R2_t_screen.txt"
    params:
        conf = config["conf"]
    threads: 2 
    resources:
        mem_mb = 12000
    conda:
        "../envs/fastqscreen.yaml"
    shell:
        """fastq_screen --aligner bowtie2 --conf {params.conf} --outdir samples/fastqscreen/{wildcards.sample} {input.fwd} {input.rev}"""

#Make sure that this works for zipped reads
rule star_with_multimap:
    input:
        fwd = "samples/trimmed/{sample}_R1_t.fq",
        rev = "samples/trimmed/{sample}_R2_t.fq"
    output:
        "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam",
        "samples/star/{sample}_bam/ReadsPerGene.out.tab",
        "samples/star/{sample}_bam/Log.final.out"
    threads: 12
    resources:
        mem_mb = 29000 
    params:
        gtf=config["gtf_file"]
    run:
        STAR=config["star_tool"],
        pathToGenomeIndex = config["star_index"]


#         shell("""
#                {STAR} --runThreadN {threads} --runMode alignReads --genomeDir {pathToGenomeIndex} \
#                --readFilesIn {input.fwd} {input.rev} \
#                --outFileNamePrefix samples/star/{wildcards.sample}_bam/ \
#                --sjdbGTFfile {params.gtf} --quantMode GeneCounts \
#                --sjdbGTFtagExonParentGene gene_name \
#                --outSAMtype BAM SortedByCoordinate \
#                --readFilesCommand zcat \
#                --twopassMode Basic
#                """)
#
        shell("""
            {STAR} --runThreadN 12 \
             --genomeDir {pathToGenomeIndex} \
             --sjdbGTFfile {params.gtf} \
             --sjdbOverhang 100 \ 
             --readFilesIn {input.fwd} {input.rev} \
             --outFileNamePrefix samples/align_read/{wildcards.sample}/ \
             --outSAMtype BAM Unsorted \
             --winAnchorMultimapNmax 200 \
             --outFilterMultimapNmax 100 \
             --readFilesCommand zcat \
             """)



rule index:
    input:
        "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam"
    output:
        "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam.bai"
    conda:
        "../envs/omic_qc_wf.yaml"
    threads: 2 
    resources:
        mem_mb = 12000
    shell:
        """samtools index {input} {output}"""

rule star_statistics:
    input:
        expand("samples/star/{sample}_bam/Log.final.out",sample=SAMPLES)
    output:
        "results/tables/{project_id}_STAR_mapping_statistics.txt".format(project_id = config["project_id"])
    threads: 2 
    resources:
        mem_mb = 12000
    script:
        "../scripts/compile_star_log.py"


#rule picard:
#    input:
#        "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam"
#    output:
#        temp("samples/genecounts_rmdp/{sample}_bam/{sample}.rmd.bam")
#    params:
#        name="rmd_{sample}",
#    threads: 2 
#    resources:
#        mem_mb = 5300 
#    run:
#      picard=config["picard_tool"]
#
#      shell("java -Xmx3g -jar {picard} \
#      INPUT={input} \
#      OUTPUT={output} \
#      METRICS_FILE=samples/genecounts_rmdp/{wildcards.sample}_bam/{wildcards.sample}.rmd.metrics.text \
#      REMOVE_DUPLICATES=true")
#
#
#rule sort:
#    input:
#      "samples/genecounts_rmdp/{sample}_bam/{sample}.rmd.bam"
#    output:
#      temp("samples/genecounts_rmdp/{sample}_bam/{sample}_sort.rmd.bam"),
#    params:
#      name = "sort_{sample}",
#    threads: 2 
#    resources:
#        mem_mb = 6400
#    conda:
#      "../envs/omic_qc_wf.yaml"
#    shell:
#      """samtools sort -O bam -n {input} -o {output}"""
#
#
#rule samtools_stats:
#    input:
#        "samples/genecounts_rmdp/{sample}_bam/{sample}_sort.rmd.bam"
#    output:
#        "samples/samtools_stats/{sample}.txt"
#    log:
#        "logs/samtools_stats/{sample}_samtools_stats.log"
#    conda:
#        "../envs/omic_qc_wf.yaml"
#    wrapper:
#        "0.17.0/bio/samtools/stats"
#
#
