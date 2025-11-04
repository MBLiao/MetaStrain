import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re
import warnings

def convert_individual(individual_str):
    # Remove brackets and extra spaces from string, and convert to integer array
    individual_str = individual_str.strip('[]').split()
    individual_array = np.array([int(x) for x in individual_str])
    return individual_array


warnings.filterwarnings("ignore")
pd.set_option('display.width', 1000)
pd.set_option('display.max_columns', 1000)
np.set_printoptions(linewidth=400, threshold=200)  # np.inf represents positive infinity

# file_path = r"C:\Users\Wenb1n\Desktop\expdata\ResultAnalysis\binGA_best_individual.csv"
#
# # Read the last row of each file and extract encoding
# encodings = []
#
# df = pd.read_csv(file_path)
# # Assuming the optimal individual encoding is in the last row of the first column
# encoding_str = df['binGA_f2']
# print(encoding_str)
# converted_individuals = encoding_str.apply(convert_individual)
# print(converted_individuals)
#
# for ind in converted_individuals:
#     print(np.count_nonzero(ind))


# # Convert all individuals to a 2D NumPy array
# individuals_array = np.stack(converted_individuals.values)
#
# # Count the number of times each gene is selected
# gene_selection_counts = np.sum(individuals_array, axis=0)
#
# # Create a new DataFrame where each column is a gene and each row corresponds to an experiment result
# genes_df = pd.DataFrame(individuals_array, columns=[f"Gene{i+1}" for i in range(individuals_array.shape[1])])
#
# # Add statistical results as the last row
# stats_row = pd.DataFrame([gene_selection_counts], columns=genes_df.columns)
# genes_df = pd.concat([genes_df, stats_row], ignore_index=True)
#
# # Add label to the last row
# genes_df.index = list(range(1, len(genes_df))) + ["Count"]


# # Save results to a new CSV file
# output_file_path = r"C:\Users\Wenb1n\Desktop\expdata\gene_selection_statistics.csv"
# genes_df.to_csv(output_file_path, index=False)

file_path = r"C:\Users\Wenb1n\Desktop\expdata\ResultAnalysis\trans_JADE_trp_best_individual20_2-3_without_b1260-1264.csv"
data = pd.read_csv(file_path)

# Initialize a dictionary to store the counts for each gene
stats = {}

# Iterate through each column (gene)
for gene in data.columns:
    # Count the occurrences of each value (0, 1, 2, 3)
    counts = data[gene].value_counts().reindex([0, 1, 2, 3], fill_value=0)
    stats[gene] = counts

# Convert the dictionary to a DataFrame
stats_df = pd.DataFrame(stats)

# Transpose for a more intuitive format (optional)
stats_df = stats_df.transpose()
stats_df.columns = ["Count_0", "Count_1", "Count_2", "Count_3"]

# Save the results to a new CSV file
output_file_path = r"C:\Users\Wenb1n\Desktop\expdata\ResultAnalysis\trans_JADE_trp_count20_2-3_without_b1260-1264.csv"
stats_df.to_csv(output_file_path)