import pandas as pd
import numpy as np

def analyze_gene_flux_simple():
    # 1. Read files
    print("Reading files...")
    # Read gene list (assuming only one column with column name 'genename')
    genes_df = pd.read_csv(r"D:\PycharmProjects\SDanalysis\models\dimensionality_reduction\ecFSEOF_MAPPING_eciML1515_without_b1260-1264.csv")
    GENElist = genes_df['gene_name'].tolist()
    print(GENElist)

    # Read gene-reaction mapping
    df = pd.read_excel(".\\models\\eciML1515_batch.xls")
    gene_protein_map = {gene: enz for gene, enz in zip(df.iloc[4825:6084, 3], df.iloc[4825:6084, 0])}
    # df = pd.read_excel("./models/ecModel_batch.xls")
    # gene_protein_map = {gene: enz for gene, enz in zip(df.iloc[7175:8143, 3], df.iloc[7175:8143, 0])}

    # Read reference fluxes (assuming first column is reaction name, second column is flux value)
    fluxes_df = pd.read_csv(r"C:\Users\Wenb1n\Desktop\expdata\eciML1515_data\without_b1260-1264\all_ref_1.csv")

    # 2. Process data
    # Create a dictionary to store results
    results = {
        'zero_flux': {'reactions': [], 'genes': []},
        'nonzero_flux': {'reactions': [], 'genes': []}
    }

    # Create reverse mapping from enzyme to gene for easy lookup
    enzyme_to_gene = {enz: gene for gene, enz in gene_protein_map.items()}

    # 3. Filter reactions with zero and non-zero flux
    print("Filtering reactions with zero and non-zero flux...")
    count = 0
    for index, row in fluxes_df.iterrows():
        rxn_id = str(row.iloc[0])  # First column is reaction ID
        flux = row.iloc[1]  # Second column is flux value
        count += 1
        # Only consider enzyme-related reactions (assuming these reactions start with 'draw_prot_')
        if rxn_id.startswith('draw_prot'):
            gene = enzyme_to_gene[rxn_id]
            # Check if this enzyme is in our mapping
            if gene in GENElist:
                if flux == 0:
                    results['zero_flux']['reactions'].append(rxn_id)
                    results['zero_flux']['genes'].append(gene)
                else:
                    results['nonzero_flux']['reactions'].append(rxn_id)
                    results['nonzero_flux']['genes'].append(gene)
    print(count)

    # 5. Output results
    print(f"\nNumber of reactions with zero flux: {len(results['zero_flux']['reactions'])}")
    print(f"Number of related genes: {len(results['zero_flux']['genes'])}")

    print(f"\nNumber of reactions with non-zero flux: {len(results['nonzero_flux']['reactions'])}")
    print(f"Number of related genes: {len(results['nonzero_flux']['genes'])}")

    # 6. Save results to CSV
    # Reactions with zero flux and related genes
    zero_flux_data = []
    for rxn in results['zero_flux']['reactions']:
        gene = enzyme_to_gene.get(rxn, "Unknown")
        zero_flux_data.append({"Reaction": rxn, "Gene": gene, "Flux": "Zero"})

    # Reactions with non-zero flux and related genes
    nonzero_flux_data = []
    for rxn in results['nonzero_flux']['reactions']:
        gene = enzyme_to_gene.get(rxn, "Unknown")
        nonzero_flux_data.append({"Reaction": rxn, "Gene": gene, "Flux": "Non-Zero"})

    # Combine and save
    all_data = zero_flux_data + nonzero_flux_data
    results_df = pd.DataFrame(all_data)
    results_df.to_csv("filtered_gene_ref_trp_without_b1260-1264.csv", index=False)

    print("\nResults saved to 'filtered_gene_ref_trp_without_b1260-1264.csv'")


# Read the first CSV file - containing reaction, gene and flux information
def process_files(flux_file_path, gene_file_path, output_file_path):
    # Read the first file - containing Reaction, Gene, Flux information
    df_flux = pd.read_csv(flux_file_path)

    # Find genes with zero flux
    zero_flux_genes = df_flux[df_flux['Flux'] == 'Zero']['Gene'].unique().tolist()
    print(f"Number of genes with zero flux found: {len(zero_flux_genes)}")
    print(f"First 10 genes with zero flux: {zero_flux_genes[:10] if len(zero_flux_genes) > 10 else zero_flux_genes}")

    # Read the second file - last row contains gene names
    df_gene = pd.read_csv(gene_file_path)

    # Assuming last row contains gene names
    genes_row = df_gene.iloc[-1].tolist()

    # Create new row, mark genes with zero flux as 0, others as 1
    new_row = [0 if gene in zero_flux_genes else 1 for gene in genes_row]

    # Add new row to DataFrame
    df_gene.loc[len(df_gene)] = new_row

    # Save modified file
    df_gene.to_csv(output_file_path, index=False)
    print(f"Results saved to {output_file_path}")



def process_gene_file(input_file_path, output_file_path):
    # Read CSV file without specifying header row
    df = pd.read_csv(input_file_path, header=None)

    # Get total number of rows and columns
    total_rows, total_cols = df.shape

    # Get the last row (label row) and second-to-last row (gene names row)
    labels_row = df.iloc[total_rows - 1].tolist()
    gene_names_row = df.iloc[total_rows - 2].tolist()

    print(f"Total rows: {total_rows}, Total columns: {total_cols}")
    print(f"Label row sample: {labels_row[:5]} ... {labels_row[-5:]}")
    print(f"Gene names row sample: {gene_names_row[:5]} ... {gene_names_row[-5:]}")

    # Find column indices with label 0 (among first 62 columns)
    zero_label_columns = []
    for i in range(total_cols - 1):  # Exclude last column (fitness)
        if labels_row[i] == str(0):
            zero_label_columns.append(i)

    print(f"Number of columns with label 0 found: {len(zero_label_columns)}")

    # For columns with label 0, convert values of 2 to 3
    changed_count = 0
    for col_idx in zero_label_columns:
        mask = (df.iloc[:total_rows - 2, col_idx] == str(2))
        changed_count += mask.sum()
        df.iloc[:total_rows - 2, col_idx] = df.iloc[:total_rows - 2, col_idx].replace(str(2), 3)

    print(f"Total converted {changed_count} values from 2 to 3")

    # Save processed file
    df.to_csv(output_file_path, index=False, header=False)
    print(f"Results saved to {output_file_path}")


# Main function
def main():
    analyze_gene_flux_simple()

    flux_file_path = r"D:\PycharmProjects\SDanalysis\filtered_gene_ref_trp_without_b1260-1264.csv"
    gene_file_path = r"C:\Users\Wenb1n\Desktop\expdata\eciML1515_data\without_b1260-1264\trp20_summary_ind_without_b1260-1264.csv"
    output_file_path = 'label_0_gene_trp_without_b1260-1264.csv'

    process_files(flux_file_path, gene_file_path, output_file_path)

    input_file_path = r"D:\PycharmProjects\SDanalysis\individual_trans\filter_ind_summary_trp20_without_b1260-1264.csv"
    output_file_path = r'individual_trans\trans_ind_summary_trp20_2-3_without_b1260-1264.csv'

    try:
        process_gene_file(input_file_path, output_file_path)
        # You can choose which function to use based on the actual format of the second file
    except Exception as e:
        print(f"Error processing file with standard method: {e}")
        print("Trying alternative method...")

if __name__ == "__main__":
    main()
