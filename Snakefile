#microbial Sequences pipeline
#Rowan Callahan
#Ngo Lab 02/22/2021

import sys
print(sys.version)

#lets locate and find our configuration file
configfile: "config.yaml"

ADAPTERS = config["ADAPTERS"]

##now lets find the samples and create a list of them
SAMPLES, = glob_wildcards("samples/raw/{sample}_R1.fq.gz" )
print(SAMPLES)


rule all:
    input:
        # The first rule should define the default target files
        # Subsequent target rules can be specified below. They should start with all_*.
        #expand("samples/align/sorted/{sample}_sorted.bam", sample=SAMPLES),
        #expand("results/{sample}_LBBC.txt", sample=SAMPLES),
        expand("results/{sample}_kraken_output.txt", sample=SAMPLES),
        expand("results/{sample}_mash_output.txt", sample=SAMPLES),
        expand("results/tables/{project_id}_STAR_mapping_statistics.txt", project_id = config["project_id"]),


#include: "rules/quality_and_align.smk"
#include: "rules/low_biomass_background_correction.smk"
include: "rules/fast_read_screening.smk"
include: "rules/align_rmdp.smk"

