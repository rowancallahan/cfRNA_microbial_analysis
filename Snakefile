#microbial Sequences pipeline
#Rowan Callahan
#Ngo Lab 02/22/2021

#lets locate and find our configuration file
configfile: "config.yaml"

ADAPTERS = config["ADAPTERS"]
REFERENCE = config["REFERENCE"] #The reference to align to

##now lets find the samples and create a list of them
SAMPLES, = glob_wildcards("samples/raw/{sample}_R1.fastq.gz" )
print(SAMPLES)


rule all:
    input:
        # The first rule should define the default target files
        # Subsequent target rules can be specified below. They should start with all_*.
        expand("samples/align/sorted/{sample}_sorted.bam", sample=SAMPLES),
		

include: "rules/quality_and_align.smk"
