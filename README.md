Description and answers:

Pipeline:

-Snakemake pipeline with Snakefile was used by command snakemake --cores 1

Main bash tools in Snakefile are bedtools intersect, bedtools getfasta, blastn

Snakefile also uses process_annotate.py script for subsequent data processing

-environment.yml provides description of environment


1) initial_files:

- IAD143293_241_Designed.bed is initial panel file with hg19 as reference
- ensembl_genes.gz and ensembl_exons.gz were exported from ucsc table browser Ensembl Genes track for hg19 assembly
- hg19.fa.gz was exported from ucsc with command wget -P /tank/projects/margarita/Panel https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz and gunzipped
  hg19.fa.gz is not provided in repository as it exceeds Github file size limit
- Clingen-Gene-Disease-Summary-2024-06-16.csv was downloaded from Clingen website


2) final_files provide answers for the task1 and task2:

- annotated_panel_IAD143293_241_Designed.bed is an initial bed file for which additional information about intersected genes and exons was added (gene_name, gene_id, exon columns)
- panel_disease_IAD143293_241_Designed.bed contains genes from panel bed file with annotation from Clingen
  This file shows that this panel is most likely to be used for testing Monogenic Diabetes also providing some insides in Pulmonary Hypertension and diseases assosiated with abnormal metabolism of fatty acids and amino acids
- non_target_IAD143293_241_Designed.bed demonstrates regions which 100% mapped to nontarget regions, showing full homology. 
  This file contains column showing quantity of nontarget regions as well as coordinates for both target and nontarget regions 




Overall conclusions:
- This panel is most likely to be used for testing Monogenic Diabetes also providing some insides in Pulmonary Hypertension and diseases assosiated with abnormal metabolism of fatty acids and amino acids
- The regions which demonstrate 100% homology with nontarget regions might be excluded from the panel
