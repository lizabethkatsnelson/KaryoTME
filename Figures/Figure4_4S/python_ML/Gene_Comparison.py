# %%
import os
import pandas as pd
from collections import defaultdict


def count_genes_in_csv_files(csv_files, gene_col_name):
    # create a list to store gene symbol lists from each file
    gene_lists = []

    # create a dictionary to store the file names for each gene symbol
    gene_files = defaultdict(list)

    # iterate over each CSV file
    for csv_file in csv_files:
        # read in the gene symbol data from the specified column
        data = pd.read_csv(csv_file["path"] + "/" + csv_file["name"])
        gene_list = data[gene_col_name].tolist()

        # store the gene symbols and the file name in the dictionary
        for gene in set(gene_list):
            gene_files[gene].append(csv_file["name"])

        gene_lists.append(gene_list)

    # combine the lists into a single list
    all_genes = [gene for gene_list in gene_lists for gene in gene_list]

    # count the occurrences of each gene symbol
    gene_counts = defaultdict(int)
    for gene in all_genes:
        gene_counts[gene] += 1

    # sort the gene counts by frequency of appearance
    sorted_counts = sorted(gene_counts.items(), key=lambda x: x[1], reverse=True)

    # calculate the percentage of each gene symbol and print the gene counts and file names
    total_genes = len(all_genes)
    for gene, count in sorted_counts:
        percent = (count / total_genes) * 100
        file_names = ", ".join(gene_files[gene])
        print(f'{gene}: {count} ({percent:.2f}%) - {file_names}')


if __name__ == "__main__":
    csv_files  = [
        {"path": "E:/Immune_Evasion_Project/Controls", "name": "pos_controls_parameters.csv"},
        {"path": "E:/Immune_Evasion_Project/Controls", "name": "neg_controls_parameters.csv"},
        {"path": "E:/Immune_Evasion_Project/Mario/Cancer_Types/HNSC", "name": "HNSC_pos_controls_parameters.csv"},
        {"path": "E:/Immune_Evasion_Project/Mario/Cancer_Types/HNSC", "name": "HNSC_neg_controls_parameters.csv"}
    ] #Use this code if you have a list of specific csv files that you want to use. You need to specify the path and the name of each csv file. Adjust accordingly.

    # call the function to count genes in the specified files
    count_genes_in_csv_files(csv_files, "Gene-Symbol")
# %%
