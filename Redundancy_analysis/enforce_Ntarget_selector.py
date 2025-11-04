import os
import numpy as np
import matplotlib.pyplot as plt
import time
from itertools import combinations
import pandas as pd

from ModelOperation_eciML1515 import *
from JADE_operator import *

np.set_printoptions(precision=15, suppress=True)


def fitness_sym_eciML1515_readref(individual, model0, target_table, GeneProteinMap, tartgetRxn, ref, metrxn):
    with model0 as model:
        for i, gene in enumerate(individual):
            if gene != 0:
                target_gene = target_table.iloc[i]["gene_name"]
                target_enzyme = GeneProteinMap[target_gene]
                # Adjust model according to operation
                if gene == 1:
                    enzymeUsage = ref[target_enzyme]
                    if enzymeUsage <= 0.00000001:
                        model.reactions.get_by_id(target_enzyme).lower_bound = 0.00000004
                    else:
                        model.reactions.get_by_id(target_enzyme).lower_bound = enzymeUsage * 4
                    # print('E' * 12)
                    # print('E' * 12)
                    # print(target_gene, target_enzyme)
                    # print(model.reactions.get_by_id(target_enzyme).bounds)
                elif gene == 2:
                    enzymeUsage = ref[target_enzyme]
                    model.reactions.get_by_id(target_enzyme).upper_bound = enzymeUsage * 0.5
                    # print('D' * 12)
                    # print('D' * 12)
                    # print(target_gene, target_enzyme)
                    # print(model.reactions.get_by_id(target_enzyme).bounds)

                elif gene == 3:
                    model.reactions.get_by_id(target_enzyme).upper_bound = 0
                    # print('O' * 12)
                    # print('O' * 12)
                    # print(target_gene, target_enzyme)
                    # print(model.reactions.get_by_id(target_enzyme).bounds)
        try:
            ft1 = time.time()
            solution_m = MPMA.moma(model=model, reference_fluxes=metrxn, linear=False)
            ft2 = time.time()

            status = solution_m.status
            if status == 'optimal':
                mutantYield_m = solution_m.fluxes[tartgetRxn] / solution_m.fluxes['EX_glc__D_e_REV']
            else:
                mutantYield_m = 0.0001
        except:
            mutantYield_m = 0.0001

    return mutantYield_m


def evaluate_individual(individual, model0, target_table, GeneProteinMap, tartgetRxn, ref, metrxn):
    fitness = fitness_sym_eciML1515_readref(individual=individual, model0=model0, target_table=target_table,
                                            GeneProteinMap=GeneProteinMap, tartgetRxn=tartgetRxn, ref=ref,
                                            metrxn=metrxn)
    return -fitness


def create_combination_strategy(original_strategy, selected_positions):
    combination_strategy = [0] * len(original_strategy)
    for pos in selected_positions:
        combination_strategy[pos] = original_strategy[pos]
    return combination_strategy


def analyze_single_experiment_combinations(experiment_id, N, model0, target_table, GeneProteinMap, tartgetRxn,
                                           ref, metrxn, original_strategy, max_combinations=1000):
    # Find positions that need editing
    edit_positions = [i for i, gene in enumerate(original_strategy) if gene != 0]
    total_edits = len(edit_positions)

    if total_edits < N:
        print(f"exp {experiment_id}: location ({total_edits}) less than ({N}), skipped")
        return []

    print(f"exp {experiment_id}: contains {total_edits} editing positions")

    # Calculate original fitness
    original_fitness = -evaluate_individual(original_strategy, model0, target_table, GeneProteinMap,
                                            tartgetRxn, ref, metrxn)

    # Get all strategies with N targets
    all_combinations = list(combinations(edit_positions, N))
    total_combinations = len(all_combinations)

    # Limiting the number of targets
    if total_combinations > max_combinations:
        import random
        all_combinations = random.sample(all_combinations, max_combinations)

    # Save the results
    combination_results = []

    # Analyze each combination
    for idx, combination in enumerate(all_combinations):
        if (idx + 1) % 10 == 0:
            print(f"  - processing: {idx + 1}/{len(all_combinations)}")

        # Create combination strategy
        combination_strategy = create_combination_strategy(original_strategy, combination)

        # Evaluate combination strategy
        combination_fitness = -evaluate_individual(combination_strategy, model0, target_table, GeneProteinMap,
                                                   tartgetRxn, ref, metrxn)

        # Calculate relative fitness percentage
        fitness_percentage = (combination_fitness / original_fitness) * 100 if original_fitness != 0 else 0

        # Get gene information
        gene_info = []
        for pos in combination:
            gene_name = target_table.iloc[pos]["gene_name"]
            operation_names = {1: "OE", 2: "KD", 3: "KO"}
            op_type = operation_names.get(original_strategy[pos], f"Unknown({original_strategy[pos]})")
            gene_info.append(f"{gene_name}({op_type})")

        combination_results.append({
            'Experiment_ID': experiment_id,
            'Original_Edits': total_edits,
            'Original_Fitness': original_fitness,
            'Combination_Positions': str(combination),
            'Combination_Genes': " + ".join(gene_info),
            'Combination_Fitness': combination_fitness,
            'Fitness_Percentage': fitness_percentage
        })

    return combination_results


def analyze_all_experiments_combinations(N, model0, target_table, GeneProteinMap, tartgetRxn,
                                         all_strategies_list, max_combinations_per_exp=1000):
    print(f"\n========== Starting analysis of {N}-target combinations for all experiments ==========\n")
    all_results = []
    # Analyze each experiment
    for exp_id in range(1, len(all_strategies_list) + 1):
        print(f"\nProcessing experiment {exp_id}/{len(all_strategies_list)}...")
        try:
            read_fluxes = pd.read_csv(f"C:/Users/Wenb1n/Desktop/expdata/eciML1515_data/ref/ref_{exp_id}.csv",
                                      header=None, index_col=0)
            read_all_fluxes = pd.read_csv(f"C:/Users/Wenb1n/Desktop/expdata/eciML1515_data/ref/all_ref_{exp_id}.csv",
                                          header=None, index_col=0)
            ref_fluxes = pd.Series(read_fluxes[1])
            all_ref_fluxes = pd.Series(read_all_fluxes[1])
        except:
            print(f"  - Unable to read reference data for experiment {exp_id}, skipped")
            continue

        # Get strategy for this experiment
        original_strategy = all_strategies_list[exp_id - 1]

        # Analyze combinations for this experiment
        exp_results = analyze_single_experiment_combinations(
            experiment_id=exp_id,
            N=N,
            model0=model0,
            target_table=target_table,
            GeneProteinMap=GeneProteinMap,
            tartgetRxn=tartgetRxn,
            ref=all_ref_fluxes,
            metrxn=ref_fluxes,
            original_strategy=original_strategy,
            max_combinations=max_combinations_per_exp
        )
        all_results.extend(exp_results)

    return all_results


def create_summary_analysis(results_df, N):

    print(f"\n========== {N}-target combination summary analysis ==========\n")

    # 1. Overall statistics
    total_combinations = len(results_df)
    total_experiments = results_df['Experiment_ID'].nunique()

    print(f"Analyzed {total_combinations} combinations from {total_experiments} experiments in total")

    # 2. Find best combination for each experiment
    best_per_experiment = results_df.loc[results_df.groupby('Experiment_ID')['Fitness_Percentage'].idxmax()]
    best_per_experiment = best_per_experiment.sort_values('Fitness_Percentage', ascending=False)

    print(f"\nBest {N}-target combination for each experiment:")
    print("-" * 100)
    print(f"{'Strategy ID':<8}{'origin_edit_count':<12}{'origin_fitness':<12}{'combined_fitness':<12}{'percentage':<10}{'Gene_combination'}")
    print("-" * 100)

    for _, row in best_per_experiment.iterrows():
        print(f"{row['Experiment_ID']:<8}{row['Original_Edits']:<12}{row['Original_Fitness']:<12.4f}"
              f"{row['Combination_Fitness']:<12.4f}{row['Fitness_Percentage']:<10.1f}%  {row['Combination_Genes'][:50]}...")

    # 3. Find globally best combinations
    print(f"\nTop {N} best {N}-target combinations globally:")
    print("-" * 100)
    top_global = results_df.nlargest(N, 'Fitness_Percentage')

    for idx, (_, row) in enumerate(top_global.iterrows(), 1):
        print(f"{idx}. Experiment {row['Experiment_ID']} - {row['Fitness_Percentage']:.1f}% - {row['Combination_Genes']}")

    # 4. Create visualization
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # 4.1 Best combination performance for each experiment
    ax1 = axes[0, 0]
    best_per_experiment_sorted = best_per_experiment.sort_values('Experiment_ID')
    ax1.bar(best_per_experiment_sorted['Experiment_ID'], best_per_experiment_sorted['Fitness_Percentage'])
    ax1.set_xlabel('Strategy ID')
    ax1.set_ylabel('Best combination fitness percentage (%)')
    ax1.set_title(f'Fitness for {N} targets for each strategy')
    ax1.axhline(y=90, color='red', linestyle='--', alpha=0.7, label='90%')
    ax1.axhline(y=80, color='orange', linestyle='--', alpha=0.7, label='80%')
    ax1.legend()

    # 4.2 Fitness percentage distribution
    ax2 = axes[0, 1]
    ax2.hist(results_df['Fitness_Percentage'], bins=30, alpha=0.7, edgecolor='black')
    ax2.set_xlabel('Fitness_Percentage (%)')
    ax2.set_ylabel('Number of combinations')
    ax2.set_title(f'Fitness distribution for all {N}-target combinations')
    ax2.axvline(x=results_df['Fitness_Percentage'].mean(), color='red', linestyle='--',
                label=f'average: {results_df["Fitness_Percentage"].mean():.1f}%')
    ax2.legend()

    # 4.3 Number of combinations in different performance ranges
    ax3 = axes[1, 0]
    performance_bins = [0, 50, 60, 70, 80, 90, 100]
    performance_labels = ['<50%', '50-60%', '60-70%', '70-80%', '80-90%', '≥90%']
    bin_counts = pd.cut(results_df['Fitness_Percentage'], bins=performance_bins).value_counts().sort_index()
    ax3.bar(range(len(bin_counts)), bin_counts.values)
    ax3.set_xticks(range(len(bin_counts)))
    ax3.set_xticklabels(performance_labels, rotation=45)
    ax3.set_ylabel('Number of combinations')
    ax3.set_title('Distribution of combinations in different performance ranges')

    # 4.4 Original edits vs best combination performance
    ax4 = axes[1, 1]
    ax4.scatter(best_per_experiment['Original_Edits'], best_per_experiment['Fitness_Percentage'], s=100, alpha=0.7)
    ax4.set_xlabel('Original number of edits')
    ax4.set_ylabel('Best combination fitness percentage (%)')
    ax4.set_title('Relationship between original edits and best combination performance')

    # Add trend line
    z = np.polyfit(best_per_experiment['Original_Edits'], best_per_experiment['Fitness_Percentage'], 1)
    p = np.poly1d(z)
    ax4.plot(best_per_experiment['Original_Edits'], p(best_per_experiment['Original_Edits']),
             "r--", alpha=0.8, label=f'Trend line: y={z[0]:.2f}x+{z[1]:.2f}')
    ax4.legend()

    plt.tight_layout()
    plt.savefig(f'all_experiments_{N}targets_summary.png', dpi=300, bbox_inches='tight')
    plt.close()

    # 5. Save summary statistics
    summary_stats = {
        'Total_Experiments': total_experiments,
        'Total_Combinations': total_combinations,
        'Mean_Fitness_Percentage': results_df['Fitness_Percentage'].mean(),
        'Std_Fitness_Percentage': results_df['Fitness_Percentage'].std(),
        'Max_Fitness_Percentage': results_df['Fitness_Percentage'].max(),
        'Min_Fitness_Percentage': results_df['Fitness_Percentage'].min(),
        'Combinations_Above_90': len(results_df[results_df['Fitness_Percentage'] >= 90]),
        'Combinations_Above_80': len(results_df[results_df['Fitness_Percentage'] >= 80]),
        'Combinations_Above_70': len(results_df[results_df['Fitness_Percentage'] >= 70])
    }

    summary_df = pd.DataFrame([summary_stats])
    summary_df.to_csv(f'all_experiments_{N}targets_summary_stats.csv', index=False)

    print(f"\nSummary statistics:")
    for key, value in summary_stats.items():
        print(f"{key}: {value:.2f}" if isinstance(value, float) else f"{key}: {value}")


if __name__ == '__main__':
    warnings.filterwarnings("ignore")
    pd.set_option('display.width', 1000)
    pd.set_option('display.max_columns', 1000)
    np.set_printoptions(linewidth=400, threshold=200)

    plt.rc("font", family='YouYuan')
    # Initialize model
    ec_iML1515 = load_ec_iML1515()

    GeneProteinMap = {}
    df = pd.read_excel(".\\models\\eciML1515_batch.xls")
    GeneProteinMap = {gene: enz for gene, enz in zip(df.iloc[4825:6084, 3], df.iloc[4825:6084, 0])}
    ecFSEOF_tabel = pd.read_csv(".\\models\\dimensionality_reduction\\ecFSEOF_MAPPING_eciML1515.csv")

    rxnID = []
    for i, x in enumerate(ec_iML1515.reactions):
        rxnID.append(x.id)

    for re in rxnID:
        if ec_iML1515.reactions.get_by_id(re).upper_bound == np.inf:
            ec_iML1515.reactions.get_by_id(re).upper_bound = 1000.0

    ec_iML1515.reactions.get_by_id('BIOMASS_Ec_iML1515_core_75p37M').lower_bound = 0.1
    ec_iML1515.reactions.get_by_id('EX_trp__L_e').lower_bound = 0.01

    print(ec_iML1515.reactions.get_by_id('BIOMASS_Ec_iML1515_core_75p37M').bounds)
    print(ec_iML1515.reactions.get_by_id('EX_glc__D_e_REV').bounds)
    print(ec_iML1515.reactions.get_by_id('EX_trp__L_e').bounds)

    # Read all strategies
    ind_readall = pd.read_csv(
        r"C:\Users\Wenb1n\Desktop\expdata\ResultAnalysis\trans_JADE_trp_best_individual20_2-3.csv",
        header=0, index_col=False)
    gene_columns = [col for col in ind_readall.columns if col.startswith('Gene_')]
    strategies_df = ind_readall[gene_columns]
    strategies_array = strategies_df.values.astype(int)

    all_strategies_list = []
    for iter_time in range(20):
        stra_list = strategies_array[iter_time].tolist()
        all_strategies_list.append(stra_list)

    # Set parameters
    N = int(input("Please enter the number of targets in the combination (recommend 3 or 4): "))
    max_combinations_per_exp = int(input("Please enter the maximum number of combinations per experiment (recommend 100-1000): "))

    # Analyze all experiments
    all_results = analyze_all_experiments_combinations(
        N=N,
        model0=ec_iML1515,
        target_table=ecFSEOF_tabel,
        GeneProteinMap=GeneProteinMap,
        tartgetRxn='EX_trp__L_e',
        all_strategies_list=all_strategies_list,
        max_combinations_per_exp=max_combinations_per_exp
    )

    # Convert to DataFrame
    if all_results:
        results_df = pd.DataFrame(all_results)

        # Save detailed results
        output_filename = f'all_experiments_{N}targets_combinations.csv'
        results_df.to_csv(output_filename, index=False)
        print(f"\nDetailed results saved to: {output_filename}")

        # Create summary analysis
        create_summary_analysis(results_df, N)

        print("\nAnalysis completed!")
    else:
        print("\nNo valid results generated!")




