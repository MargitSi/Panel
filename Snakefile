SAMPLES = ["IAD143293_241_Designed"]

rule all:
    input:
        annot=expand("final_files/annotated_panel_{sample}.bed", sample=SAMPLES),
        disease=expand("final_files/panel_disease_{sample}.bed", sample=SAMPLES),
        miss=expand("final_files/non_target_{sample}.bed", sample=SAMPLES)


rule intersect:
    input:
        bed=expand("initial_files/{sample}.bed", sample=SAMPLES),
        genes="initial_files/ensembl_genes.gz",
        exons="initial_files/ensembl_exons.gz"
    output:
        indexed=expand("intermediate_files/indexed_{sample}.bed", sample=SAMPLES),
        genes=expand("intermediate_files/intersect_genes_{sample}.bed", sample=SAMPLES),
        exons=expand("intermediate_files/intersect_exons_{sample}.bed", sample=SAMPLES)
    shell:
        "awk '{{print $0\"\t\" NR}}' {input.bed} > {output.indexed} && "
        "bedtools intersect -a {output.indexed} -b {input.genes} -wa -wb > {output.genes} && "
        "bedtools intersect -a {output.indexed} -b {input.exons} -wa -wb > {output.exons}"

rule nontarget:
    input:
        bed=expand("initial_files/{sample}.bed", sample=SAMPLES),
        fa="initial_files/hg19.fa"
    output:
        fasta=expand("intermediate_files/panel_{sample}.fa", sample=SAMPLES),
        blast=expand("intermediate_files/blast_results_{sample}.txt", sample=SAMPLES)
    shell:
        "bedtools getfasta -fi {input.fa} -bed {input.bed} > {output.fasta} && "
        "makeblastdb -in {input.fa} -dbtype nucl -out intermediate_files/DB && "  
        "blastn -query {output.fasta} -db intermediate_files/DB -outfmt 6 > {output.blast}"

rule process_annotate:
    input:
        genes=expand("intermediate_files/intersect_genes_{sample}.bed", sample=SAMPLES),
        exons=expand("intermediate_files/intersect_exons_{sample}.bed", sample=SAMPLES),
        indexed=expand("intermediate_files/indexed_{sample}.bed", sample=SAMPLES),
        clingen="initial_files/Clingen-Gene-Disease-Summary-2024-06-16.csv",        
        blast=expand("intermediate_files/blast_results_{sample}.txt", sample=SAMPLES) 
    output:
        annot=expand("final_files/annotated_panel_{sample}.bed", sample=SAMPLES),
        disease=expand("final_files/panel_disease_{sample}.bed", sample=SAMPLES),
        miss=expand("final_files/non_target_{sample}.bed", sample=SAMPLES)
    script:
        "process_annotate.py"
