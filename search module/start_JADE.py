import ray
from GA_Operator import *
from ModelOperation import *
import time
import matplotlib.pyplot as plt
from JADE_operator import *

# Set numerical display format
np.set_printoptions(precision=15, suppress=True)

# Ray initialization
ray.init()


@ray.remote
def evaluate_individual(individual, model0, target_table, GeneProteinMap, tartgetRxn, ref, metrxn):
    fitness = fitness_sym(individual=individual, model0=model0, target_table=target_table,
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

    ecYeast = loadModel()

    rxnID = []
    for i, x in enumerate(ecYeast.reactions):
        rxnID.append(x.id)
    rxnID_protein_draw = [x for i, x in enumerate(rxnID) if 'draw_prot_' in x]
    geneAll = []
    for gene in ecYeast.genes:
        geneAll.append(gene.id)

    '''
    Set the max upper bound to 1000 instead of inf
    '''
    for re in rxnID:
        if ecYeast.reactions.get_by_id(re).upper_bound == np.inf:
            ecYeast.reactions.get_by_id(re).upper_bound = 1000.0

    wildSolution_moma = WildModel_Growth(ecYeast)

    ecYeast.reactions.get_by_id('r_2111').lower_bound = 0.1
    ecYeast.reactions.get_by_id('r_1589').lower_bound = 0.1
    ecYeast.reactions.get_by_id('r_1714_REV').lower_bound = 0.1

    print(ecYeast.reactions.get_by_id('r_2111').bounds)
    print(ecYeast.reactions.get_by_id('r_1589').bounds)
    print(ecYeast.reactions.get_by_id('r_1714_REV').bounds)

    wildYield = wildSolution_moma.fluxes['r_1589'] / wildSolution_moma.fluxes['r_1714_REV']

    metabolic_solution = pd.Series()
    enzyme_solution = pd.Series()
    arm_solution = pd.Series()
    protein_pool = pd.Series()
    all_solution = pd.Series()
    for rid in wildSolution_moma.fluxes.index:
        all_solution[rid] = wildSolution_moma.fluxes[rid]
        if rid.startswith('r_'):
            metabolic_solution[rid] = wildSolution_moma.fluxes[rid]
        if rid.startswith('draw_'):
            enzyme_solution[rid] = wildSolution_moma.fluxes[rid]
        if rid.startswith('arm_'):
            arm_solution[rid] = wildSolution_moma.fluxes[rid]
    protein_pool['prot_pool_exchange'] = wildSolution_moma.fluxes['prot_pool_exchange']

    # get the gene list for the related proteins with kcat number
    gene_with_Kcat = [x.replace('draw_prot_', '') for x in rxnID_protein_draw]

    # Save the reference fluxes
    metabolic_solution.to_csv(r"C:\Users\Wenb1n\Desktop\ref_check.csv", header=False, index=True)
    all_solution.to_csv(r"C:\Users\Wenb1n\Desktop\all_ref_check.csv", header=False, index=True)

    GeneProteinMap = {}
    df = pd.read_excel(r"C:\Users\Wenb1n\StrainDesign\models\ecModel_batch.xls")
    GeneProteinMap = {gene: enz for gene, enz in zip(df.iloc[7175:8143, 3], df.iloc[7175:8143, 0])}

    ecFSEOF_tabel = pd.read_csv(r"C:\Users\Wenb1n\StrainDesign\models\dimensionality_reduction\ecFSEOF_MAPPING.csv")

    '''
    Show the basic information of the ecYeast model
    '''
    print('#####################--MODEL INFORMATION--#####################')
    print('model tolerance:       ', ecYeast.tolerance)
    print('rxnID:                 ', len(rxnID))
    print('geneAll:               ', len(geneAll))
    print('rxnID_protein_draw:    ', len(rxnID_protein_draw))
    print('gene_with_Kcat:        ', len(gene_with_Kcat))
    print('metabolic reactions:   ', len(metabolic_solution))
    print('enzyme reactions:      ', len(enzyme_solution))
    print('arm reactions:         ', len(arm_solution))
    print('prot_pool_exchange:    ', len(protein_pool))
    print('exchanges:             ', len(ecYeast.exchanges))
    print('demands:               ', len(ecYeast.demands))
    print('sinks:                 ', len(ecYeast.sinks))
    print('boundary:              ', len(ecYeast.boundary))
    print('wild yield:            ', wildYield)
    print('wild yield regulated:  ', 0)
    print('###############################################################')

    n = 52
    popsize = 100
    totalTime = 2

    start_time = time.time()

    print("JADE Algorithm")

    # Parameter and problem setup
    lu = np.array([[0] * n, [4] * n])

    all_result = []
    all_best_individual = []

    # edit_ratio = 0.5  # Control random initialization ratio, the rest are all-zero individuals

    # pr_ind = individual_filling()
    # num_existing = len(pr_ind)
    # num_random = int(popsize * edit_ratio)

    time_iter = 1

    while time_iter <= totalTime:
        np.random.seed(int(time.time()))

        # # Initialize population
        # popold = np.zeros((popsize, n))  # Initialize all individuals to all zeros
        #
        # if num_random > 0:
        #     random_indices = np.random.choice(popsize, num_random, replace=False)  # Randomly select some individuals
        #     popold[random_indices] = lu[0] + np.random.rand(num_random, n) * (lu[1] - lu[0])  # Replace some individuals with random values
        # popold[:num_existing] = pr_ind  # Fill the first part with existing individuals

        popold = lu[0] + np.random.rand(popsize, n) * (lu[1] - lu[0])

        valParents = evaluate(popold, ecYeast, ecFSEOF_tabel, GeneProteinMap,
                              'r_1589', wildSolution_moma, metabolic_solution)
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

        while FES < 200:
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
            valOffspring = evaluate(ui, ecYeast, ecFSEOF_tabel, GeneProteinMap,
                                    'r_1589', wildSolution_moma, metabolic_solution)
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
    df_ind.to_csv(r"C:\Users\Wenb1n\Desktop\expdata\JADE_ind.csv", index=False)
    df_val.to_csv(r"C:\Users\Wenb1n\Desktop\expdata\JADE_val.csv", index=False)

    # Plot convergence curve
    average_fitness = np.mean(all_result, axis=0)
    plt.plot(average_fitness)
    plt.xlabel('Iteration')
    plt.ylabel('Fitness')
    plt.title('Convergence Curve')
    plt.show()

    print("JADE algorithm completed.")
    