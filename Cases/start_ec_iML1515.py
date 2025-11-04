import os

import numpy as np
import ray
from ModelOperation_eciML1515 import *
import time
import matplotlib.pyplot as plt
from JADE_operator import *

# Set numerical display format
np.set_printoptions(precision=15, suppress=True)

# Ray initialization
ray.init()


def fitness_sym_eciML1515(individual, model0, target_table, GeneProteinMap, tartgetRxn, ref, metrxn):
    with model0 as model:
        for i, gene in enumerate(individual):
            if gene != 0:
                target_gene = target_table.iloc[i]["gene_name"]
                target_enzyme = GeneProteinMap[target_gene]
                # Adjust model according to operation
                if gene == 1:
                    enzymeUsage = ref.fluxes[target_enzyme]
                    if enzymeUsage <= 0.00000001:
                        model.reactions.get_by_id(target_enzyme).lower_bound = 0.00000004
                    else:
                        model.reactions.get_by_id(target_enzyme).lower_bound = enzymeUsage * 4
                    # print('E' * 12)
                    # print('E' * 12)
                    # print(model.reactions.get_by_id(target_enzyme).bounds)
                elif gene == 2:
                    enzymeUsage = ref.fluxes[target_enzyme]
                    model.reactions.get_by_id(target_enzyme).upper_bound = enzymeUsage * 0.5
                    # print('D' * 12)
                    # print('D' * 12)
                    # print(model.reactions.get_by_id(target_enzyme).bounds)

                elif gene == 3:
                    model.reactions.get_by_id(target_enzyme).upper_bound = 0
                    # print('O' * 12)
                    # print('O' * 12)
                    # print(model.reactions.get_by_id(target_enzyme).bounds)
        try:
            ft1 = time.time()
            solution_m = MPMA.moma(model=model, reference_fluxes=metrxn, linear=False)
            ft2 = time.time()

            status = solution_m.status
            if status == 'optimal':
                mutantYield_m = solution_m.fluxes[tartgetRxn] / solution_m.fluxes['EX_glc__D_e_REV']
                print("---------------------------------------------MOMA----------------------------------------------")
                print("TIME: ", ft2 - ft1)
                growth = solution_m.fluxes['BIOMASS_Ec_iML1515_core_75p37M']
                product = solution_m.fluxes[tartgetRxn]
                glucose = solution_m.fluxes['EX_glc__D_e_REV']
                print(growth, '--', product, '--', glucose, '--', mutantYield_m)
                print('-----------------------------------------------------------------------------------------------')
            else:
                mutantYield_m = 0.0001
        except:
            mutantYield_m = 0.0001

    return mutantYield_m

@ray.remote
def evaluate_individual(individual, model0, target_table, GeneProteinMap, tartgetRxn, ref, metrxn):
    fitness = fitness_sym_eciML1515(individual=individual, model0=model0, target_table=target_table,
                                    GeneProteinMap=GeneProteinMap, tartgetRxn=tartgetRxn, ref=ref, metrxn=metrxn)
    return -fitness


def map_individual(individual):
    """
    Map continuous values [0, 4] to discrete values [0, 1, 2, 3]:
    0: No adjustment
    1: Overexpression
    2: Knockdown
    3: Knockout
    """
    # bins=[1, 2, 3]
    discrete_individual = np.digitize(individual, bins=[3.4, 3.6, 3.8], right=True)
    return discrete_individual


def evaluate(pop, model0, target_table, GeneProteinMap, tartgetRxn, ref, metrxn):
    """
    Evaluate the fitness values of the population.
    For each individual, first map to discrete values, then call fitness function.
    """
    pop_separete = []
    for individual in pop:
        discrete_individual = map_individual(individual)
        pop_separete.append(discrete_individual)
        # print(discrete_individual)
    # print(pop_separete)

    futures = [evaluate_individual.remote(individual, model0, target_table, GeneProteinMap, tartgetRxn, ref, metrxn)
               for individual in pop_separete]
    results = ray.get(futures)

    return results


# Run main function
if __name__ == '__main__':
    warnings.filterwarnings("ignore")
    pd.set_option('display.width', 1000)
    pd.set_option('display.max_columns', 1000)
    np.set_printoptions(linewidth=400, threshold=200)  # np.inf represents positive infinity

    ec_iML1515 = load_ec_iML1515()

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

    max_growth(ec_iML1515)
    max_product(ec_iML1515)
    WildSolution = Wild_Growth_Trp(ec_iML1515)

    # step = np.arange(1, 4, 0.1)
    # for oe in step:
    #     print(oe)
    #     with ec_iML1515 as m:
    #         m.reactions.get_by_id('draw_prot_P0A877').lower_bound = oe*WildSolution.fluxes['draw_prot_P0A877']
    #         m.reactions.get_by_id('draw_prot_P0A879').lower_bound = oe*WildSolution.fluxes['draw_prot_P0A879']
    #         m.reactions.get_by_id('draw_prot_P00909').lower_bound = oe*WildSolution.fluxes['draw_prot_P00909']
    #         m.reactions.get_by_id('draw_prot_P00904').lower_bound = oe*WildSolution.fluxes['draw_prot_P00904']
    #         m.reactions.get_by_id('draw_prot_P00895').lower_bound = oe*WildSolution.fluxes['draw_prot_P00895']
    #         Wild_Growth_Trp(m)

    ec_iML1515.reactions.get_by_id('BIOMASS_Ec_iML1515_core_75p37M').lower_bound = 0.1
    ec_iML1515.reactions.get_by_id('EX_trp__L_e').lower_bound = 0.01

    print(ec_iML1515.reactions.get_by_id('BIOMASS_Ec_iML1515_core_75p37M').bounds)
    print(ec_iML1515.reactions.get_by_id('EX_glc__D_e_REV').bounds)
    print(ec_iML1515.reactions.get_by_id('EX_trp__L_e').bounds)

    wildYield = WildSolution.fluxes['EX_trp__L_e'] / WildSolution.fluxes['EX_glc__D_e_REV']

    metabolic_solution = pd.Series()
    enzyme_solution = pd.Series()
    arm_solution = pd.Series()
    protein_pool = pd.Series()
    all_solution = pd.Series()
    for rid in WildSolution.fluxes.index:
        all_solution[rid] = WildSolution.fluxes[rid]
        if rid.startswith('draw_'):
            enzyme_solution[rid] = WildSolution.fluxes[rid]
        else:
            metabolic_solution[rid] = WildSolution.fluxes[rid]
        if rid.startswith('arm_'):
            arm_solution[rid] = WildSolution.fluxes[rid]
    protein_pool['prot_pool_exchange'] = WildSolution.fluxes['prot_pool_exchange']

    gene_with_Kcat = [x.replace('draw_prot_', '') for x in rxnID_protein_draw]

    # Save the reference fluxes
    metabolic_solution.to_csv("ref_check.csv", header=False, index=True)
    all_solution.to_csv("all_ref_check——1.csv", header=False, index=True)

    GeneProteinMap = {}
    df = pd.read_excel(".\\models\\eciML1515_batch.xls")
    GeneProteinMap = {gene: enz for gene, enz in zip(df.iloc[4825:6084, 3], df.iloc[4825:6084, 0])}

    ecFSEOF_tabel = pd.read_csv(".\\models\\dimensionality_reduction\\ecFSEOF_MAPPING_eciML1515_nob1260-1264.csv")

    print('#####################--MODEL INFORMATION--#####################')
    print('model tolerance:       ', ec_iML1515.tolerance)
    print('rxnID:                 ', len(rxnID))
    print('geneAll:               ', len(geneAll))
    print('rxnID_protein_draw:    ', len(rxnID_protein_draw))
    print('gene_with_Kcat:        ', len(gene_with_Kcat))
    print('metabolic reactions:   ', len(metabolic_solution))
    print('enzyme reactions:      ', len(enzyme_solution))
    print('arm reactions:         ', len(arm_solution))
    print('prot_pool_exchange:    ', len(protein_pool))
    print('exchanges:             ', len(ec_iML1515.exchanges))
    print('demands:               ', len(ec_iML1515.demands))
    print('sinks:                 ', len(ec_iML1515.sinks))
    print('boundary:              ', len(ec_iML1515.boundary))
    print('medium:                ', ec_iML1515.medium)
    print('wild yield:            ', wildYield)
    print('wild yield regulated:  ', 0)
    print('###############################################################')

    n = 57
    popsize = 100
    totalTime = 1

    start_time = time.time()

    print("JADE Algorithm")

    # Parameter and problem setup
    lu = np.array([[0] * n, [4] * n])

    all_result = []
    all_best_individual = []

    time_iter = 1

    while time_iter <= totalTime:
        # seed = random.randint(0, len(ecYeast.exchanges))
        # np.random.seed(int(time.time()))

        popold = lu[0] + np.random.rand(popsize, n) * (lu[1] - lu[0])

        valParents = evaluate(popold, ec_iML1515, ecFSEOF_tabel, GeneProteinMap,
                              'EX_trp__L_e', WildSolution, metabolic_solution)
        # print(valParents)

        c = 0.1
        p = 0.05
        CRm = 0.5
        Fm = 0.5
        Afactor = 1
        fes_v = np.array([0.01] + list(0.1 * np.arange(1, 11)))
        fj = 1

        FES = 0
        best_fitness = []
        goodCR, goodF = [], []

        pNP = max(int(round(p * popsize)), 2)

        # Initialize archive
        archive = {
            'NP': int(Afactor * popsize),
            'pop': np.zeros((0, n)),
            'funvalues': np.zeros(0)
        }

        FES += popsize

        # Sort and get best individual
        indBest = np.argsort(valParents)
        valBest = np.sort(valParents)

        min_val = np.min(valParents)
        best_fitness.append(-min_val)

        # Output information for this generation
        print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
        print("Time iter: ", time_iter)
        print("FES: ", FES)
        print("CRm: ", CRm)
        print("Fm: ", Fm)
        print(len(popold))
        print(len(archive['pop']), len(archive['funvalues']))
        print(pNP)
        print(f"Generation {0}: "
              f"Max Fitness = {-min_val}\n"
              f"Solvable Num of individual: {sum(i != -0.0001 for i in valParents)}")
        print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')

        while FES < 500:
            ti1 = time.time()
            pop = popold.copy()
            if FES > 1 and len(goodCR) > 0 and np.sum(goodF) > 0:
                CRm = (1 - c) * CRm + c * np.mean(goodCR)
                Fm = (1 - c) * Fm + c * np.sum(goodF ** 2) / np.sum(goodF)

            # Generate F and CR
            F, CR = randFCR(popsize, CRm, 0.1, Fm, 0.1)

            # Generate r1 and r2 indices
            popAll = np.vstack([pop, archive['pop']])
            r1, r2 = gnR1R2(popsize, popAll.shape[0])

            # Select p-best solutions
            randindex = np.clip(np.ceil(np.random.rand(popsize) * pNP).astype(int), 1, pNP) - 1
            pbest = pop[indBest[randindex], :]

            # Mutation operation
            vi = pop + F[:, None] * (pbest - pop + pop[r1, :] - popAll[r2, :])
            vi = boundConstraint(vi, pop, lu)

            # Crossover operation
            mask = np.random.rand(popsize, n) > CR[:, None]
            rows = np.arange(popsize)
            cols = np.floor(np.random.rand(popsize) * n).astype(int)
            jrand = (rows, cols)
            mask[jrand] = False
            ui = vi.copy()
            ui[mask] = pop[mask]

            # Calculate offspring fitness
            valOffspring = evaluate(ui, ec_iML1515, ecFSEOF_tabel, GeneProteinMap,
                                    'EX_trp__L_e', WildSolution, metabolic_solution)
            FES += popsize

            # Selection operation
            combined = np.vstack([valParents, valOffspring]).T
            valParents, I = np.min(combined, axis=1), np.argmin(combined, axis=1)
            archive = updateArchive(archive, popold[I == 1, :], valParents[I == 1])
            popold[I == 1, :] = ui[I == 1, :]

            goodCR = CR[I == 1]
            goodF = F[I == 1]

            # Record the best fitness value of current iteration
            min_val = np.min(valParents)
            best_fitness.append(-min_val)

            valBest = np.sort(valParents)
            indBest = np.argsort(valParents)
            ti2 = time.time()

            # Output information for this generation
            print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
            print("Time iter: ", time_iter)
            print("FES: ", FES)
            print("CRm: ", CRm)
            print("Fm: ", Fm)
            print("CR: ", CR)
            print("F: ", F)
            print(len(pop))
            print(len(popold))
            print(len(popAll))
            print(len(archive['pop']), len(archive['funvalues']))
            print(pNP)
            print(f"Generation {int((FES / popsize) - 1)}: "
                  f"Max Fitness = {-min_val}\n"

                  f"Time: {ti2 - ti1}\n"
                  f"Solvable Num of individual: {sum(i != -0.0001 for i in valParents)}")
            print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')

        all_result.append(best_fitness)
        all_best_individual.append(popold[indBest[0]])
        time_iter += 1

    end_time = time.time()
    print(f"Total Time: {end_time - start_time:.2f} seconds")

    mapped_best_individuals = [map_individual(individual) for individual in all_best_individual]

    # Convert best individuals and fitness values to DataFrame
    df_ind = pd.DataFrame(mapped_best_individuals, columns=[f'Gene_{i + 1}' for i in range(n)])
    df_val = pd.DataFrame(all_result, columns=[f'Iteration_{i}' for i in range(len(all_result[0]))])

    # Add fitness values to individual DataFrame
    df_ind['Best Fitness'] = [max(fitnesses) for fitnesses in all_result]

    # Save DataFrame to CSV file
    df_ind.to_csv("JADE_ind_trp.csv", index=False)
    df_val.to_csv("JADE_val_trp.csv", index=False)

    ray.shutdown()

    # Plot convergence curve
    average_fitness = np.mean(all_result, axis=0)
    plt.plot(average_fitness)
    plt.xlabel('Iteration')
    plt.ylabel('Fitness')
    plt.title('Convergence Curve')
    plt.show()

    print("JADE algorithm completed.")
