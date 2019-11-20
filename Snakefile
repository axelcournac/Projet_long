INDEX = ['1','2','3','4','rev.1','rev.2']

rule all:
    input:
        '/media/anthony/POULOP/HIC/results/SRR9900851.bed2'

rule alignment:
    input:
        left_read='/media/anthony/POULOP/HIC/data/fastq/{sample}_1.fastq',
        right_read='/media/anthony/POULOP/HIC/data/fastq/{sample}_2.fastq',
    output:
        left_align='/media/anthony/POULOP/HIC/alignment/{sample}_1.sam',
        right_align='/media/anthony/POULOP/HIC/alignment/{sample}_2.sam'
    shell:
        """
        cat /media/anthony/POULOP/HIC/data/ref/NC*.fasta > /media/anthony/POULOP/HIC/data/ref/genome.fa
        bowtie2-build /media/anthony/POULOP/HIC/data/ref/genome.fa genome
        mv genome* /media/anthony/POULOP/HIC/data/ref/
        bowtie2 -x /media/anthony/POULOP/HIC/data/ref/genome -p18 --sam-no-hd --sam-no-sq --local --very-sensitive-local -S {output.left_align} {input.left_read}
        bowtie2 -x /media/anthony/POULOP/HIC/data/ref/genome -p18 --sam-no-hd --sam-no-sq --local --very-sensitive-local -S {output.right_align} {input.right_read}
        mv /media/anthony/POULOP/HIC/data/ref/genome.fa /media/anthony/POULOP/HIC/data/obsolete/
        """

rule sam_extraction:
    input:
        left_old ='/media/anthony/POULOP/HIC/alignment/{sample}_1.sam',
        right_old='/media/anthony/POULOP/HIC/alignment/{sample}_2.sam'
    output:
        left_new='/media/anthony/POULOP/HIC/alignment/{sample}_1.sam.0',
        right_new='/media/anthony/POULOP/HIC/alignment/{sample}_2.sam.0'
    shell:
        """
        awk '{{print $1,$3,$4,$2,$5;}}' {input.left_old} > {output.left_new}
        awk '{{print $1,$3,$4,$2,$5;}}' {input.right_old} > {output.right_new}
        rm {input.left_old} {input.right_old}
        """

rule sam_sorting:
    input:
        left_unsorted='/media/anthony/POULOP/HIC/alignment/{sample}_1.sam.0',
        right_unsorted='/media/anthony/POULOP/HIC/alignment/{sample}_2.sam.0'
    output:
        left_sorted='/media/anthony/POULOP/HIC/alignment/{sample}_1.sam.0.sorted',
        right_sorted='/media/anthony/POULOP/HIC/alignment/{sample}_2.sam.0.sorted'
    shell:
        """
        sort -V -k1 {input.left_unsorted} > {output.left_sorted}
        sort -V -k1 {input.right_unsorted} > {output.right_sorted}
        rm {input.left_unsorted} {input.right_unsorted}
        """

rule pairing:
    input:
        left='/media/anthony/POULOP/HIC/alignment/{sample}_1.sam.0.sorted',
        right='/media/anthony/POULOP/HIC/alignment/{sample}_2.sam.0.sorted'
    output:
        '/media/anthony/POULOP/HIC/alignment/{sample}.merged.dat'
    shell:
        """
         paste {input.left} {input.right} > {output}
         rm {input.left} {input.right}
         """

rule quality_filtering:
    input:
        '/media/anthony/POULOP/HIC/alignment/{sample}.merged.dat'
    output:
        '/media/anthony/POULOP/HIC/alignment/{sample}.merged_qualfilt.dat'
    shell:
        """
        awk '{{if($1==$6 && $5>= 1 && $10 >= 1) print $2,$3,$4,$7,$8,$9}}'  {input}  > {output}
        rm {input}
        """

rule fragment_restriction:
    input:
        genome='/media/anthony/POULOP/HIC/data/ref/',
        align='/media/anthony/POULOP/HIC/alignment/{sample}.merged_qualfilt.dat'
    output:
        '/media/anthony/POULOP/HIC/alignment/{sample}.merged_qualfilt.dat.indices'
    script:
        "examples_codes/fragment_attribution.py"
        
rule event_filtering:
    input:
        file='/media/anthony/POULOP/HIC/alignment/{sample}.merged_qualfilt.dat.indices'
    params:
        uncut_threshold='4',
        loop_threshold='2',
        srr='{sample}'
    output:
        '/media/anthony/POULOP/HIC/alignment/{sample}.merged_qualfilt.dat.indices.filtered'
    script:
        "examples_codes/library_events_ARG.py"
        
rule creating:
    input:
        '/media/anthony/POULOP/HIC/alignment/{sample}.merged_qualfilt.dat.indices.filtered'
    output:
        "/media/anthony/POULOP/HIC/alignment/{sample}.merged_qualfilt.dat.indices.filtered.bed2"
    script:
        "examples_codes/convert_pairs_bed2d.py"

rule post_processing:
    input:
        bad_name='/media/anthony/POULOP/HIC/alignment/{sample}.merged_qualfilt.dat.indices.filtered.bed2'
    output:
        good_name='/media/anthony/POULOP/HIC/results/{sample}.bed2'
    shell:
        """
        mv {input.bad_name} {output.good_name}
        mv *.npz /media/anthony/POULOP/HIC/results
        mv *.png /media/anthony/POULOP/HIC/results
        rm /media/anthony/POULOP/HIC/alignment/*.dat
        rm /media/anthony/POULOP/HIC/alignment/*.indices
        """

