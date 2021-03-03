rule align_bacterial_reference:
    input:
        rules.sortbam.output.bam
        #"samples/align/sorted/{sample}_sorted.bam"
    output:
        "results/{sample}_bacterial_hits.tab"
    shell:
        """ echo test"""

rule max_likelihood:
    input:
        rules.align_bacterial_reference.output
    output:
        "results/{sample}_max_likelihood.txt"
    shell:
        """ echo test"""

rule LBBC:
    input:
        rules.max_likelihood.output
    output:
        "results/{sample}_LBBC.txt"
    shell:
        """echo test"""


