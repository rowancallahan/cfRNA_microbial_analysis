rule trimming:
	input:
		forward_strand = "samples/raw/{sample}_R1.fastq.gz",
		reverse_strand = "samples/raw/{sample}_R2.fastq.gz"
	output:
		forward_paired = "samples/trim/{sample}_R1_paired.fastq.gz",
		reverse_paired = "samples/trim/{sample}_R2_paired.fastq.gz",
		forward_unpaired = "samples/trim/{sample}_R1_unpaired.fastq.gz",
		reverse_unpaired = "samples/trim/{sample}_R2_unpaired.fastq.gz"
	params:
		adapter = config["ADAPTERS"]
	conda:
		"../envs/trim.yaml"
	threads: 8
	log: "logs/trimming/{sample}.log"
	message:
		"""--- Trimming {wildcards.sample} ---"""
	shell:
		"trimmomatic PE -threads {threads} {input.forward_strand} {input.reverse_strand} \
			{output.forward_paired} {output.forward_unpaired} \
			{output.reverse_paired} {output.reverse_unpaired} \
			ILLUMINACLIP:{params.adapter}:2:30:10 \
			LEADING:3 \
			TRAILING:3 \
			SLIDINGWINDOW:4:15 \
			MINLEN:30"

rule align:
	input:
		forward_paired = rules.trimming.output.forward_paired,
		reverse_paired = rules.trimming.output.reverse_paired
	output:
		"samples/align/bam/{sample}.bam"
	params:
		reference = config["REFERENCE"]
	conda:
		"../envs/bwa.yaml"
	log: "logs/align/{sample}.log"
	message:
		"""--- Aligning {wildcards.sample} to reference."""
	threads: 16
	shell:
		"bwa mem -t {threads} {params.reference} {input.forward_paired} {input.reverse_paired} | samtools view -bS - > {output}"

rule sortbam:
	input:
		rules.align.output
	output:
		"samples/align/sorted/{sample}_sorted.bam",
		"samples/align/sorted/{sample}_sorted.bam.bai"
	conda:
		"../envs/bwa.yaml"
	threads: 4
	message:
		"""--- sorting {wildcards.sample} ---"""
	shell:
		"samtools sort -@ {threads} -m '2G' {input} > {output[0]}; samtools index -b {output[0]}"
