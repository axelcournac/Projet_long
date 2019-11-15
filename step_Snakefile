rule indexing:
    input:
        'data/{ref}.fa'
    output:
        'data/{ref}.fa.1.bt2',
        'data/{ref}.fa.2.bt2',
        'data/{ref}.fa.3.bt2',
        'data/{ref}.fa.4.bt2',
        'data/{ref}.fa.rev.1.bt2',
        'data/{ref}.fa.rev.2.bt2'
    shell:
        'bowtie2-build {input} {input.ref}'

rule alignment:
    input:
        left_read='data/{sample}_1.fastq',
        right_read='data/{sample}_2.fastq',
        genome='data/{ref}.fa'
    output:
        left_align='alignment/{sample}_1.sam',
        right_align='alignment/{sample}_2.sam'
    shell:
        """
        bowtie2 -x {input.genome} -p18 --sam-no-sd --sam-no-sq --local --very-sensitive-local -S {output.left_align} {input.left_read}
        bowtie2 -x {input.genome} -p18 --sam-no-sd --sam-no-sq --local --very-sensitive-local -S {output.right_align} {input.right_read}
        """
        
rule sam_extraction:
    input:
        left_old='alignment/{sample}_1.sam',
        right_old='alignment/{sample}_2.sam'
    output:
        left_new='alignment/{sample}_1.sam.0',
        right_new='alignment/{sample}_2.sam.0'
    shell:
        """
        awk '{print $1,$3,$4,$2,$5;}' {input.left_old} > {output.left_new}
        awk '{print $1,$3,$4,$2,$5;}' {input.right_old} > {output.right_new}
        """

rule sam_sorting:
    input:
        left_unsorted='alignment/{sample}_1.sam.0'
        right_unsorted='alignment/{sample}_2.sam.0'
    output:
        left_sorted='alignment/{sample}_1.sam.0.sorted',
        right_sorted='alignment/{sample}_2.sam.0.sorted'
    shell:
        """
        samtools sort -n {input.left_unsorted} > {output.left_sorted}
        samtools sort -n {input.right_unsorted} > {output.right_sorted}
        """
