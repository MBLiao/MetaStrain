import pandas as pd


# Read file 1, extract gene names and editing methods
file_1_data = pd.read_csv(r"D:\PycharmProjects\SDanalysis\models\dimensionality_reduction\eciML1515_Tryptophan_gluc-2_2_ecFSEOF_result.csv")
gene_names_and_methods = file_1_data[['gene', 'action']]

# Read file 2, extract gene list
file_2_data = pd.read_csv(r"D:\PycharmProjects\SDanalysis\models\dimensionality_reduction\all_drawpro_genes_eciML1515.csv")
gene_list = file_2_data['gene_name'].tolist()

# Filter entries where gene names and editing methods are in the gene list
filtered_genes = gene_names_and_methods[gene_names_and_methods['gene'].isin(gene_list)]

# Save results in the format of file 3
output_data = pd.DataFrame()
output_data['gene_name'] = filtered_genes['gene']
# output_data['operation'] = filtered_genes['actions']


output_data.to_csv(r"D:\PycharmProjects\SDanalysis\models\dimensionality_reduction\ecFSEOF_MAPPING_eciML1515.csv", index=False)

# import pandas as pd
#
# # Read the first dimensionality reduction result file
# file_1_data = pd.read_csv(r"C:\Users\Wenb1n\PycharmProjects\SDanalysis\models\dimensionality_reduction\fseof_result_AA.csv")
# genes_1 = file_1_data['gene']
#
# # Read the second dimensionality reduction result file (assuming there's a similar file, you need to provide the specific path)
# file_1_data_2 = pd.read_csv(r"C:\Users\Wenb1n\Desktop\乙醇做碳源的FSEOF靶点.csv")
# genes_2 = file_1_data_2['gene ID']
#
# # Get the union of genes from both files
# combined_genes = pd.concat([genes_1, genes_2]).drop_duplicates()
# print(combined_genes.shape)
#
# # Read target file, extract gene list
# target_data = pd.read_csv(r"C:\Users\Wenb1n\PycharmProjects\SDanalysis\models\dimensionality_reduction\all_drawpro_genes.csv")
# target_genes = target_data['gene_name'].tolist()
#
# # Find genes in the union that are also in the target gene list
# common_genes = combined_genes[combined_genes.isin(target_genes)]
#
# # Create output dataframe
# output_data = pd.DataFrame()
# output_data['gene_name'] = common_genes
#
# # Save results
# output_data.to_csv(r"C:\Users\Wenb1n\PycharmProjects\SDanalysis\models\dimensionality_reduction\ecFSEOF_MAPPING_AA_tailored.csv", index=False)
#