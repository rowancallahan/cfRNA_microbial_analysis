rule kraken:
    input:
        rules.sortbam.output.bam
    output:
        "results/{sample}_kraken_output.txt"
    shell:
        """echo test"""

rule mash:
    input:
        rules.sortbam.output.bam
    output:
        "results/{sample}_mash_output.txt"
    shell:
        """echo test """

