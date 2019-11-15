INDEX = ['1','2','3','4','rev.1','rev.2']

rule all:
    input:
        expand('genome.{index}.bt2', index = INDEX)
                          
rule indexing:
    input:
        '/media/anthony/POULOP/data/{ref}.fa'
    output:
        '{ref}.1.bt2',
        '{ref}.2.bt2',
        '{ref}.3.bt2',
        '{ref}.4.bt2',
        '{ref}.rev.1.bt2',
        '{ref}.rev.2.bt2'
    shell:
        'bowtie2-build {input} {wildcards.ref}'
