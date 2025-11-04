import os
import numpy as np
import matplotlib.pyplot as plt
import time

from ModelOperation_eciML1515 import *
from JADE_operator import *

# Set numeric display format
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
                # print("---------------------------------------------MOMA----------------------------------------------")
                # print("TIME: ", ft2 - ft1)
                # growth = solution_m.fluxes['BIOMASS_Ec_iML1515_core_75p37M']
                # product = solution_m.fluxes[tartgetRxn]
                # glucose = solution_m.fluxes['EX_glc__D_e_REV']
                # print(growth, '--', product, '--', glucose, '--', mutantYield_m)
                # print('-----------------------------------------------------------------------------------------------')
            else:
                mutantYield_m = 0.0001
        except:
            mutantYield_m = 0.0001

    return mutantYield_m


def evaluate_individual(individual, model0, target_table, GeneProteinMap, tartgetRxn, ref, metrxn):
    fitness = fitness_sym_eciML1515_readref(individual=individual, model0=model0, target_table=target_table,
                                            GeneProteinMap=GeneProteinMap, tartgetRxn=tartgetRxn, ref=ref, metrxn=metrxn)
    return -fitness


def analyze_redundancy(strategy, model0, target_table, GeneProteinMap, tartgetRxn, ref, metrxn, threshold=0.95):
    print("\nStarting redundancy analysis...")

    original_fitness = -evaluate_individual(strategy, model0, target_table, GeneProteinMap,
                                            tartgetRxn, ref, metrxn)
    print(f"Verified original strategy fitness: {original_fitness:.6f}")

    edit_positions = [i for i, gene in enumerate(strategy) if gene != 0]
    edit_count = len(edit_positions)

    print(f"Original number of edits: {edit_count}")

    # Record redundant and essential positions
    redundant_positions = []
    essential_positions = []

    # Test the necessity of each edit position
    for pos in edit_positions:
        # Create a copy with this edit removed
        test_strategy = strategy.copy()
        test_strategy[pos] = 0  # Set to no operation

        # Evaluate fitness
        test_fitness = -evaluate_individual(test_strategy, model0, target_table, GeneProteinMap,
                                            tartgetRxn, ref, metrxn)

        # Determine if redundant
        is_redundant = test_fitness >= original_fitness * threshold

        gene_name = target_table.iloc[pos]["gene_name"]
        target_enzyme = GeneProteinMap[gene_name]

        operation_names = {1: "OE", 2: "KD", 3: "KO"}
        op_type = operation_names.get(strategy[pos], f"Unknown({strategy[pos]})")

        if is_redundant:
            redundant_positions.append((pos, gene_name, target_enzyme, op_type, test_fitness))
            print(
                f"Position {pos} ({gene_name}, {op_type}) is redundant - fitness after removal: {test_fitness:.6f} ({test_fitness / original_fitness:.2%})")
        else:
            essential_positions.append((pos, gene_name, target_enzyme, op_type))
            print(
                f"Position {pos} ({gene_name}, {op_type}) is essential - fitness after removal: {test_fitness:.6f} ({test_fitness / original_fitness:.2%})")

    # Create minimal strategy (keep only essential edits)
    minimal_strategy = strategy.copy()
    for pos, _, _, _, _ in redundant_positions:
        minimal_strategy[pos] = 0

    # Evaluate minimal strategy
    minimal_fitness = -evaluate_individual(minimal_strategy, model0, target_table, GeneProteinMap,
                                           tartgetRxn, ref, metrxn)

    remaining_edits = sum(1 for gene in minimal_strategy if gene != 0)

    print(f"\nRedundancy analysis results:")
    print(f"Original number of edits: {edit_count}")
    print(f"Number of redundant edits: {len(redundant_positions)}")
    print(f"Number of essential edits: {len(essential_positions)}")
    print(f"Minimal strategy number of edits: {remaining_edits}")
    print(f"Minimal strategy fitness: {minimal_fitness:.6f} ({minimal_fitness / original_fitness:.2%} relative to original)")

    return minimal_strategy, essential_positions, redundant_positions, original_fitness, minimal_fitness


if __name__ == '__main__':
    warnings.filterwarnings("ignore")
    pd.set_option('display.width', 1000)
    pd.set_option('display.max_columns', 1000)
    np.set_printoptions(linewidth=400, threshold=200)  # np.inf means positive infinity

    ec_iML1515 = load_ec_iML1515()

    GeneProteinMap = {}
    df = pd.read_excel(".\\models\\eciML1515_batch.xls")
    GeneProteinMap = {gene: enz for gene, enz in zip(df.iloc[4825:6084, 3], df.iloc[4825:6084, 0])}
    ecFSEOF_tabel = pd.read_csv(".\\models\\dimensionality_reduction\\ecFSEOF_MAPPING_eciML1515_without_b1260-1264.csv")

    rxnID = []
    geneAll = []
    for i, x in enumerate(ec_iML1515.reactions):
        rxnID.append(x.id)
    for gene in ec_iML1515.genes:
        geneAll.append(gene.id)
    rxnID_protein_draw = [x for i, x in enumerate(rxnID) if 'draw_prot_' in x]

    for re in rxnID:
        if ec_iML1515.reactions.get_by_id(re).upper_bound == np.inf:
            ec_iML1515.reactions.get_by_id(re).upper_bound = 1000.0

    ec_iML1515.reactions.get_by_id('BIOMASS_Ec_iML1515_core_75p37M').lower_bound = 0.1
    ec_iML1515.reactions.get_by_id('EX_trp__L_e').lower_bound = 0.01

    print(ec_iML1515.reactions.get_by_id('BIOMASS_Ec_iML1515_core_75p37M').bounds)
    print(ec_iML1515.reactions.get_by_id('EX_glc__D_e_REV').bounds)
    print(ec_iML1515.reactions.get_by_id('EX_trp__L_e').bounds)

    ind_readall = pd.read_csv(
        r"C:\Users\Wenb1n\Desktop\expdata\ResultAnalysis\trans_JADE_trp_best_individual20_2-3_without_b1260-1264.csv",
        header=0, index_col=False)
    gene_columns = [col for col in ind_readall.columns if col.startswith('Gene_')]
    strategies_df = ind_readall[gene_columns]
    strategies_array = strategies_df.values.astype(int)

    all_strategies_list = []
    all_check_fitness = []
    all_op_counts = []
    for iter_time in range(20):
        read_fluxes = pd.read_csv(f"C:/Users/Wenb1n/Desktop/expdata/eciML1515_data/without_b1260-1264/ref_{iter_time+1}.csv",
                                     header=None, index_col=0)
        read_all_fluxes = pd.read_csv(f"C:/Users/Wenb1n/Desktop/expdata/eciML1515_data/without_b1260-1264/all_ref_{iter_time+1}.csv",
                                      header=None, index_col=0)
        read_fluxes = pd.Series(read_fluxes[1])
        read_all_fluxes = pd.Series(read_all_fluxes[1])

        ind_data = pd.read_csv(f"C:/Users/Wenb1n/Desktop/expdata/eciML1515_data/without_b1260-1264/JADE_ind_trp_{iter_time+1}.csv", index_col=False)
        ind_ref_fitness = ind_data['Best Fitness'].values[0]

        op_counts = {
            "0": np.sum(strategies_array[iter_time] == 0),
            "OE": np.sum(strategies_array[iter_time] == 1),
            "KD": np.sum(strategies_array[iter_time] == 2),
            "KO": np.sum(strategies_array[iter_time] == 3)
        }
        all_op_counts.append(op_counts)
        print(f"Experiment {iter_time+1}, strategy:{strategies_array[iter_time]}:")
        print(f"  Operation statistics: {op_counts}")
        print(f"  Fitness value: {ind_ref_fitness}")

        stra_list = strategies_array[iter_time].tolist()
        all_strategies_list.append(stra_list)
        # check_fitness = evaluate_individual(stra_list, ec_iML1515, ecFSEOF_tabel, GeneProteinMap,
        #                                     'EX_trp__L_e', read_all_fluxes, read_fluxes)
        # all_check_fitness.append(check_fitness)

    all_minimal_strategies = []
    all_essential_positions = []
    all_redundant_positions = []
    all_minimal_fitness = []

    print("\n========== Starting redundancy analysis ==========")
    for iter_time in range(len(all_strategies_list)):
        print(f"\nAnalyzing redundancy for experiment {iter_time + 1}/20...")

        # Get data for this experiment
        strategy = all_strategies_list[iter_time]
        read_fluxes = pd.read_csv(f"C:/Users/Wenb1n/Desktop/expdata/eciML1515_data/without_b1260-1264/ref_{iter_time + 1}.csv",
                                  header=None, index_col=0)
        read_all_fluxes = pd.read_csv(f"C:/Users/Wenb1n/Desktop/expdata/eciML1515_data/without_b1260-1264/all_ref_{iter_time + 1}.csv",
                                      header=None, index_col=0)
        ref_fluxes = pd.Series(read_fluxes[1])
        all_ref_fluxes = pd.Series(read_all_fluxes[1])

        # Get number of edits in original strategy
        edit_count = sum(1 for gene in strategy if gene != 0)
        if edit_count <= 1:
            print(f"Strategy has only {edit_count} edit(s), skipping redundancy analysis")
            all_minimal_strategies.append(strategy)
            all_essential_positions.append([])
            all_redundant_positions.append([])
            continue

        # Perform redundancy analysis
        minimal_strategy, essential, redundant, origin_f, minimal_f = analyze_redundancy(
            strategy, ec_iML1515, ecFSEOF_tabel, GeneProteinMap,
            'EX_trp__L_e', all_ref_fluxes, ref_fluxes, threshold=0.95
        )

        all_minimal_strategies.append(minimal_strategy)
        all_essential_positions.append(essential)
        all_redundant_positions.append(redundant)
        all_check_fitness.append(origin_f)
        all_minimal_fitness.append(minimal_f)

    # Save results of all minimal strategies
    print("\n========== Redundancy analysis results summary ==========")
    results = []
    for i in range(len(all_minimal_strategies)):
        original_edits = sum(1 for gene in all_strategies_list[i] if gene != 0)
        minimal_edits = sum(1 for gene in all_minimal_strategies[i] if gene != 0)
        original_fitness = all_check_fitness[i]
        minimal_fitness = all_minimal_fitness[i]

        results.append({
            'Strategy_ID': i + 1,
            'Original_Edits': original_edits,
            'Minimal_Edits': minimal_edits,
            'Reduced_By': original_edits - minimal_edits,
            'Original_Fitness': original_fitness,
            'Minimal_Fitness': minimal_fitness,
            'Fitness_Ratio': minimal_fitness / original_fitness if original_fitness != 0 else 0,
            'Essential_Count': len(all_essential_positions[i]),
            'Redundant_Count': len(all_redundant_positions[i])
        })
    # Convert to DataFrame and sort
    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values('Minimal_Fitness', ascending=False).reset_index(drop=True)

    # Print results
    print("\nRedundancy analysis results for all strategies (sorted by minimal fitness):")
    print(results_df.to_string())

    # Save results
    results_df.to_csv("redundancy_analysis_results_without_b1260-1264.csv", index=False)

    # Combine all minimal strategies into one file
    print("\nCreating summary file...")
    all_minimal_combined = []

    for i in range(len(all_minimal_strategies)):
        strategy_id = i + 1
        minimal_strategy = all_minimal_strategies[i]
        minimal_fitness = results_df[results_df['Strategy_ID'] == strategy_id]['Minimal_Fitness'].values[0]

        # Save detailed information of essential edits
        essential = all_essential_positions[i]
        if essential:
            essential_df = pd.DataFrame(essential, columns=['Position', 'Gene_Name', 'Enzyme', 'Operation'])
            essential_df.to_csv(f"essential_separate/essential_edits_exp{strategy_id}_without_b1260-1264.csv", index=False)

        # Create a dictionary containing strategy ID and fitness
        strategy_dict = {
            'Strategy_ID': strategy_id,
            'Fitness': minimal_fitness
        }

        # Add all gene edits
        for j, gene in enumerate(minimal_strategy):
            gene_col = f"Gene_{j + 1}"
            strategy_dict[gene_col] = gene

        all_minimal_combined.append(strategy_dict)

    # Convert to DataFrame and sort
    combined_df = pd.DataFrame(all_minimal_combined)
    combined_df = combined_df.sort_values('Fitness', ascending=False).reset_index(drop=True)

    # Move ID to front, fitness to end
    cols = combined_df.columns.tolist()
    cols.remove('Strategy_ID')
    cols.remove('Fitness')
    cols = ['Strategy_ID'] + cols + ['Fitness']
    combined_df = combined_df[cols]

    # Save to CSV
    combined_df.to_csv("all_minimal_strategies_without_b1260-1264.csv", index=False)
    print(f"Saved all {len(all_minimal_strategies)} minimal strategies to all_minimal_strategies_without_b1260-1264.csv")

    # Find the best minimal strategy
    best_idx = int(results_df.iloc[0]['Strategy_ID'] - 1)

    best_minimal_strategy = all_minimal_strategies[best_idx]
    best_essential = all_essential_positions[best_idx]

    print(f"\nBest minimal strategy is from experiment {best_idx + 1}")
    print(f"Original number of edits: {results_df.iloc[0]['Original_Edits']}")
    print(f"Minimal number of edits: {results_df.iloc[0]['Minimal_Edits']}")
    print(
        f"Reduced by {results_df.iloc[0]['Reduced_By']} edit(s) ({results_df.iloc[0]['Reduced_By'] / results_df.iloc[0]['Original_Edits']:.2%})")
    print(
        f"Minimal fitness: {results_df.iloc[0]['Minimal_Fitness']:.6f} ({results_df.iloc[0]['Fitness_Ratio']:.2%} relative to original)")

    print("\nDetailed display of best minimal strategy:")
    print("-" * 80)
    print(f"{'Position':<8}{'Gene Name':<15}{'Enzyme':<20}{'Operation Type'}")
    print("-" * 80)

    for pos, gene_name, enzyme, op_type in best_essential:
        print(f"{pos:<8}{gene_name:<15}{enzyme:<20}{op_type}")

    print("\nRedundancy analysis completed!")


