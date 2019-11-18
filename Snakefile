INDEX = ['1','2']

rule all:
    input:
        '/media/anthony/POULOP/data/alignment/SRR9900851_merged_quality_filtered.dat'
                          
rule indexing:
    input:
        '/media/anthony/POULOP/data/{ref}.fa'
    output:
        '/media/anthony/POULOP/data/{ref}.1.bt2',
        '/media/anthony/POULOP/data/{ref}.2.bt2',
        '/media/anthony/POULOP/data/{ref}.3.bt2',
        '/media/anthony/POULOP/data/{ref}.4.bt2',
        '/media/anthony/POULOP/data/{ref}.rev.1.bt2',
        '/media/anthony/POULOP/data/{ref}.rev.2.bt2'
    shell:
        """
        bowtie2-build {input} {wildcards.ref}
        mv {wildcards.ref}* /media/anthony/POULOP/data
        """

rule alignment:
    input:
        left_read='/media/anthony/POULOP/data/{sample}_1.fastq',
        right_read='/media/anthony/POULOP/data/{sample}_2.fastq',
    output:
        left_align='/media/anthony/POULOP/data/alignment/{sample}_1.sam',
        right_align='/media/anthony/POULOP/data/alignment/{sample}_2.sam'
    shell:
        """
        bowtie2 -x /media/anthony/POULOP/data/genome -p18 --sam-no-hd --sam-no-sq --local --very-sensitive-local -S {output.left_align} {input.left_read}
        bowtie2 -x /media/anthony/POULOP/data/genome -p18 --sam-no-hd --sam-no-sq --local --very-sensitive-local -S {output.right_align} {input.right_read}
        """

rule sam_extraction:
    input:
        left_old ='/media/anthony/POULOP/data/alignment/{sample}_1.sam',
        right_old='/media/anthony/POULOP/data/alignment/{sample}_2.sam'
    output:
        left_new='/media/anthony/POULOP/data/alignment/{sample}_1.sam.0',
        right_new='/media/anthony/POULOP/data/alignment/{sample}_2.sam.0'
    shell:
        """
        awk '{{print $1,$3,$4,$2,$5;}}' {input.left_old} > {output.left_new}
        awk '{{print $1,$3,$4,$2,$5;}}' {input.right_old} > {output.right_new}
        """

rule sam_sorting:
    input:
        left_unsorted='/media/anthony/POULOP/data/alignment/{sample}_1.sam.0',
        right_unsorted='/media/anthony/POULOP/data/alignment/{sample}_2.sam.0'
    output:
        left_sorted='/media/anthony/POULOP/data/alignment/{sample}_1.sam.0.sorted',
        right_sorted='/media/anthony/POULOP/data/alignment/{sample}_2.sam.0.sorted'
    shell:
        """
        sort -V -k1 {input.left_unsorted} > {output.left_sorted}
        sort -V -k1 {input.right_unsorted} > {output.right_sorted}
        """
        
rule pairing:
    input:
        left='/media/anthony/POULOP/data/alignment/{sample}_1.sam.0.sorted',
        right='/media/anthony/POULOP/data/alignment/{sample}_2.sam.0.sorted'
    output:
        '/media/anthony/POULOP/data/alignment/{sample}.1_{sample}.2_merged'
    shell:
        """
         paste {input.left} {input.right} > {output}
         rm {input.left} {input.right}
         """

rule quality_filtering:
    input:
        '/media/anthony/POULOP/data/alignment/{sample}.1_{sample}.2_merged'
    output:
        '/media/anthony/POULOP/data/alignment/{sample}_merged_quality_filtered.dat'
    shell:
        """
        awk '{{if($1==$6 && $5>= 30 && $10 >= 30) print $2,$3,$4,$7,$8,$9}}'  {input}  > {output}
        rm {input}
        """
        
rule fragment_attribution:
    input:
        genome='/media/anthony/POULOP/data/',
        align='/media/anthony/POULOP/data/alignment/{sample}_merged_quality_filtered.dat'
    script:
        "examples_codes/fragment_attribution.py {input.genome} {input.align}"
