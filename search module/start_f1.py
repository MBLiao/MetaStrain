import pandas as pd
import ray
from GA_Operator import *
from ModelOperation import *


'''
Global variables
'''
POPULATION_SIZE = 50  # population size
GENE_LENGTH = 52      # gene length
ACTIVE_GENES = 5       # active genes
MAX_GENERATIONS = 500  # max generation
CROSSOVER_RATE = 0.7   # cross rate
MUTATION_RATE = 0.003   # mutation rate
# Convergence detection parameters
CONVERGENCE_THRESHOLD = 0.0001  # Define the threshold for fitness change
CONVERGENCE_STEPS = 1000        # Consider converged when fitness change is below threshold for this many consecutive generations

# Ray initialization
ray.init()


@ray.remote
def evaluate_individual(individual, model0, target_table, GeneProteinMap, tartgetRxn, ref, metrxn):
    fitness = fitness_f1(individual=individual, model0=model0, target_table=target_table,
                          GeneProteinMap=GeneProteinMap, tartgetRxn=tartgetRxn, ref=ref, metrxn=metrxn)
    return fitness


def parallel_fitness_evaluation(pop, model0, target_table, GeneProteinMap, tartgetRxn, ref, metrxn):
    futures = [evaluate_individual.remote(individual, model0, target_table, GeneProteinMap, tartgetRxn, ref, metrxn)
               for individual in pop]
    results = ray.get(futures)
    return results


if __name__ == '__main__':
    warnings.filterwarnings("ignore")
    pd.set_option('display.width', 1000)
    pd.set_option('display.max_columns', 1000)
    np.set_printoptions(linewidth=400, threshold=200)  # np.inf represents positive infinity

    start_time = time.time()  # Global start time

    ecYeast = loadModel()
    # s = S_Matrix(ecYeast)  # Get the S matrix


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
        # print(ecYeast.reactions.get_by_id(re).upper_bound)

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
    metabolic_solution.to_csv(r"C:\Users\Wenb1n\Desktop\expdata\reffluxes_GA01_50_500_0.7_0.003_m1.csv", header=False, index=True)
    all_solution.to_csv(r"C:\Users\Wenb1n\Desktop\expdata\all_reffluxes_GA01_50_500_0.7_0.003_m1.csv", header=False, index=True)

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
    print('wild yield regulated:  ', wildYield + 52 / 10000)
    print('###############################################################')

    '''
    Start GA searching
    '''
    time_all_start = time.time()

    best_fitness_over_generations = []
    best_individual_over_generations = []
    mean_fitness_over_generations = []
    # Initialize convergence detection related variables
    last_fitness = 0
    convergence_counter = 0
    num_elites = int(POPULATION_SIZE * 0.05)

    ts_init = time.time()
    population = initialize_pop_loose(POPULATION_SIZE, GENE_LENGTH)
    fitness_scores = parallel_fitness_evaluation(population, ecYeast, ecFSEOF_tabel, GeneProteinMap,
                                                 'r_1589', wildSolution_moma, metabolic_solution)

    # Bind individuals with their fitness values and sort by fitness in descending order
    pop_with_fitness = sorted(zip(population, fitness_scores), key=lambda x: x[1], reverse=True)

    # Get best fitness and best individual
    best_individual, max_fitness = pop_with_fitness[0]
    mean_fitness = np.mean(fitness_scores)

    # Save best fitness and best individual
    best_fitness_over_generations.append(max_fitness)
    best_individual_over_generations.append(best_individual)
    mean_fitness_over_generations.append(mean_fitness)
    te_init = time.time()

    # Output information for this generation
    print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
    print(num_elites)
    print(len(population))
    print(len(fitness_scores))
    print(f"Generation {0}: "
          f"Max Fitness = {max_fitness}\n"
          f"Mean Fitness = {mean_fitness}\n"
          f"Time: {te_init - ts_init}\n"
          f"Solvable Num of individual: {sum(i != 0.0001 for i in fitness_scores)}")
    print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')

    for generation in range(MAX_GENERATIONS):
        ti1 = time.time()

        new_population = []
        new_fitness_scores = []

        # Keep elite individuals
        elite_indices = np.argsort(fitness_scores)[-num_elites:]
        elites = [population[i] for i in elite_indices]
        elite_fitness_scores = [fitness_scores[i] for i in elite_indices]

        population = [population[i] for i in range(len(population)) if i not in elite_indices]
        fitness_scores = [fitness_scores[i] for i in range(len(fitness_scores)) if i not in elite_indices]

        # Selection operation
        remaining_size = POPULATION_SIZE - num_elites
        selected_indices = list(roulette_wheel_selection_POPSIZE(population, fitness_scores, remaining_size))
        # print(remaining_size)
        # print(len(selected_indices))
        if int(len(selected_indices) / 2) != len(selected_indices) / 2:
            selected_indices.append(selected_indices[0])
        # print(selected_indices)

        c_pop = []
        c_fitness_scores = []
        c_pop_temp = []
        c_fit_temp = []
        s_population = []
        for i in selected_indices:
            s_population.append(population[i])

        for i in range(0, len(selected_indices), 2):
            parent1, parent2 = s_population[i], s_population[i + 1]
            child1, child2, fg = crossover(parent1, parent2, CROSSOVER_RATE, GENE_LENGTH)
            if fg == 1:
                c_pop_temp.extend([child1, child2])
            else:
                c_pop.append(population[selected_indices[i]])
                c_pop.append(population[selected_indices[i + 1]])
                c_fitness_scores.append(fitness_scores[selected_indices[i]])
                c_fitness_scores.append(fitness_scores[selected_indices[i + 1]])

        c_fit_temp = parallel_fitness_evaluation(c_pop_temp, ecYeast, ecFSEOF_tabel, GeneProteinMap,
                                                 'r_1589', wildSolution_moma, metabolic_solution)
        c_pop.extend(c_pop_temp)
        c_fitness_scores.extend(c_fit_temp)

        m_pop_temp = []
        m_fit_temp = []

        for i in range(len(c_pop)):
            mutated_ind, mflag = mutation(c_pop[i], MUTATION_RATE, GENE_LENGTH)
            if mflag == 1:
                m_pop_temp.append(mutated_ind)
            else:
                new_population.append(c_pop[i])
                new_fitness_scores.append(c_fitness_scores[i])

        m_fit_temp = parallel_fitness_evaluation(m_pop_temp, ecYeast, ecFSEOF_tabel, GeneProteinMap,
                                                 'r_1589', wildSolution_moma, metabolic_solution)
        new_population.extend(m_pop_temp)
        new_fitness_scores.extend(m_fit_temp)

        new_population.extend(elites)
        new_fitness_scores.extend(elite_fitness_scores)

        # Update population!!!!!!!!!!
        population, fitness_scores = select_fixed_size_population(new_population, new_fitness_scores, POPULATION_SIZE)

        print('-' * 120)
        print('-' * 120)

        # Bind individuals with their fitness values and sort by fitness in descending order
        pop_with_fitness = sorted(zip(population, fitness_scores), key=lambda x: x[1], reverse=True)

        # Get best fitness and best individual
        best_individual, max_fitness = pop_with_fitness[0]
        mean_fitness = np.mean(fitness_scores)

        # Save best fitness and best individual
        best_fitness_over_generations.append(max_fitness)
        best_individual_over_generations.append(best_individual)
        mean_fitness_over_generations.append(mean_fitness)

        ti2 = time.time()

        # Output information for this generation
        print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
        print(num_elites)
        print(len(c_pop_temp))
        print(len(c_pop))
        print(len(m_pop_temp))
        print(len(new_population))
        print(len(new_fitness_scores))
        print(len(population))
        print(len(fitness_scores))
        print(f"Generation {generation + 1}: "
              f"Max Fitness = {max_fitness}\n"
              f"Mean Fitness = {mean_fitness}\n"
              f"Time: {ti2 - ti1}\n"
              f"Solvable Num of individual: {sum(i != 0.0001 for i in fitness_scores)}")
        print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')

        # Convergence detection
        # If the best fitness does not change significantly for several consecutive generations, the algorithm can be terminated early
        fitness_change = abs(max_fitness - last_fitness)
        if fitness_change < CONVERGENCE_THRESHOLD:
            convergence_counter += 1
            if convergence_counter >= CONVERGENCE_STEPS:
                print(f"Convergence reached at generation {generation}.")
                break
        else:
            convergence_counter = 0
        last_fitness = max_fitness

    # After the loop ends, output or save the best results of the entire process
    print(f"Best Fitness over all generations: {max(best_fitness_over_generations)}")
    best_overall = \
        best_individual_over_generations[best_fitness_over_generations.index(max(best_fitness_over_generations))]
    print(f"Best Individual over all generations: {best_overall}")

    data = {'individual': best_individual_over_generations, 'best fitness': best_fitness_over_generations,
            'mean fitness': mean_fitness_over_generations}
    df = pd.DataFrame(data)

    # Save DataFrame to CSV file
    df.to_csv(r"C:\Users\Wenb1n\Desktop\expdata\GA01_50_500_0.7_0.003_m1.csv", index=False)

    # Plot the change in maximum fitness
    plt.plot(best_fitness_over_generations, color='red')

    plt.title('Best Fitness Over Generations')
    plt.xlabel('Generation')
    plt.ylabel('Fitness')
    plt.show()

    plt.plot(mean_fitness_over_generations, color='blue')
    plt.title('Mean Fitness Over Generations')
    plt.xlabel('Generation')
    plt.ylabel('Fitness')
    plt.show()

    time_all_end = time.time()
    print('Searching Time:', time_all_end - time_all_start)
    end_time = time.time()  # Global end time
    print('Total Time:', end_time - start_time)
    ray.shutdown()
