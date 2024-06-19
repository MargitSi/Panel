#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt
import numpy as np

from gseapy import Biomart
bm = Biomart()

# Task1 - Annotation

ensembl_annot = bm.query(dataset='hsapiens_gene_ensembl', attributes = ["ensembl_gene_id", "ensembl_transcript_id", "external_gene_name", "go_id" ])

intersect_genes = pd.read_csv(snakemake.input[0], header=None, sep="\t")
intersect_exons = pd.read_csv(snakemake.input[1], header=None, sep="\t")
indexed = pd.read_csv(snakemake.input[2], sep="\t")


id = "6"
target = "10"



def colnaming(dataset):
    dataset.columns=dataset.columns.astype("str")
    return dataset



intersect_genes = colnaming(intersect_genes)
intersect_exons = colnaming(intersect_exons)



def annot_mapping(feature1, feature2, dataset, colname, new_col):
    map_dict = dict(ensembl_annot[[feature1,feature2]].values)
    dataset[new_col] = dataset[colname].map(map_dict)
    return dataset



intersect_genes = annot_mapping("ensembl_transcript_id","external_gene_name", intersect_genes, target, "gene_name")
intersect_genes = annot_mapping("ensembl_transcript_id","ensembl_gene_id", intersect_genes, target, "gene_id")


def annot(final_data, mapping_data, colname):
    map_dict = mapping_data[-mapping_data[colname].isna()].groupby(id)[colname].unique()
    final_data[colname]=final_data['1'].map(map_dict)
    final_data[colname] = final_data.fillna("-")[colname].apply(lambda x : ', '.join(x))
    return final_data


indexed = annot(indexed, intersect_genes, "gene_name")
indexed = annot(indexed, intersect_genes, "gene_id")
indexed = annot(indexed, intersect_exons, target)



indexed = indexed.drop("1", axis=1).rename(columns={"10":"exon"})
indexed.to_csv(snakemake.output[0], sep="\t")




# Task1 - Diagnostics


clingen = pd.read_csv(snakemake.input[3])
clingen.columns = clingen.iloc[3]
clingen = clingen.iloc[5:].reset_index(drop=True)


unique_genes = np.unique(indexed.gene_name[indexed.gene_name != "-"].values)

panel_disease = clingen[clingen["GENE SYMBOL"].isin(unique_genes)][["GENE SYMBOL", "GCEP"]]

panel_disease.to_csv(snakemake.output[1], sep="\t", index=False)


# Task2 - Nontarget mapping regions


blast_result = pd.read_csv(snakemake.input[4], sep='\t', header=None)
blast_result.columns = pd.Index(["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"])
blast_result["mapped"]  = blast_result["sseqid"] + ":" + pd.Series(blast_result["sstart"]-1).astype(str) + "-" + blast_result["send"].astype(str) 


non_target = blast_result[(blast_result["pident"]==100) & (blast_result["qseqid"]!=blast_result["mapped"])]

count_map = non_target["qseqid"].value_counts()
mapped_map = non_target.groupby("qseqid")["mapped"].unique()

non_target["nontarget_counts"] = non_target["qseqid"].map(count_map)
non_target = non_target[["qseqid","nontarget_counts"]].drop_duplicates()
non_target["nontarget_mapped"] = non_target["qseqid"].map(mapped_map)
non_target["nontarget_mapped"] = non_target["nontarget_mapped"].apply(lambda x : ', '.join(x))


non_target = non_target.rename(columns={"qseqid":"target_mapped"})
non_target.to_csv(snakemake.output[2], sep="\t", index=False)

