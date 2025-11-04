import pandas as pd
import ray
from GA_Operator import *
from ModelOperation import *


'''
Global variables
'''
POPULATION_SIZE = 20  # population size
GENE_LENGTH = 52      # gene length
ACTIVE_GENES = 5       # active genes
MAX_GENERATIONS = 5  # max generation
CROSSOVER_RATE = 0.7   # cross rate
MUTATION_RATE = 0.003   # mutation rate
# 收敛检测参数
CONVERGENCE_THRESHOLD = 0.001  # 定义适应度变化的阈值
CONVERGENCE_STEPS = 1000        # 连续多少代适应度变化小于阈值时认为已收敛

# Ray初始化
ray.init()


@ray.remote
def evaluate_individual(individual, model0, target_table, GeneProteinMap, tartgetRxn, ref, metrxn):
    fitness = fitness_fun(individual=individual, model0=model0, target_table=target_table,
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
    np.set_printoptions(linewidth=400, threshold=200)  # np.inf表示正无穷

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

    # ecYeast.reactions.get_by_id('r_1654_REV').upper_bound = 1000
    # ecYeast.reactions.get_by_id('r_1714_REV').upper_bound = 1
    # ecYeast.reactions.get_by_id('r_1832_REV').upper_bound = 1000
    # ecYeast.reactions.get_by_id('r_1861_REV').upper_bound = 1000
    # ecYeast.reactions.get_by_id('r_1992_REV').upper_bound = 1000
    # ecYeast.reactions.get_by_id('r_2005_REV').upper_bound = 1000
    # ecYeast.reactions.get_by_id('r_2020_REV').upper_bound = 1000
    # ecYeast.reactions.get_by_id('r_2049_REV').upper_bound = 1000
    # ecYeast.reactions.get_by_id('r_2060_REV').upper_bound = 1000
    # ecYeast.reactions.get_by_id('r_2100_REV').upper_bound = 1000
    # ecYeast.reactions.get_by_id('r_4593_REV').upper_bound = 1000
    # ecYeast.reactions.get_by_id('r_4594_REV').upper_bound = 1000
    # ecYeast.reactions.get_by_id('r_4595_REV').upper_bound = 1000
    # ecYeast.reactions.get_by_id('r_4596_REV').upper_bound = 1000
    # ecYeast.reactions.get_by_id('r_4597_REV').upper_bound = 1000
    # ecYeast.reactions.get_by_id('r_4600_REV').upper_bound = 1000

    # wildSolution_moma = WildModel_MOMA(ecYeast)
    wildSolution_moma = WildModel_Growth(ecYeast)
    # wildSolution_fba = WildModel_FBA(ecYeast)

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
    print(rxnID_protein_draw[: 5])
    print(gene_with_Kcat[:5])

    # Save the reference fluxes
    metabolic_solution.to_csv(r"C:\Users\Wenb1n\Desktop\expdata\reffluxes.csv", header=False, index=True)
    all_solution.to_csv(r"C:\Users\Wenb1n\Desktop\expdata\all_reffluxes.csv", header=False, index=True)

    '''
    Mapping the genes and enzyme
    '''
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
    # 初始化收敛检测相关变量
    last_fitness = 0
    convergence_counter = 0
    num_elites = int(POPULATION_SIZE * 0.05)

    ts_init = time.time()
    population = initialize_pop_loose(POPULATION_SIZE, GENE_LENGTH)
    fitness_scores = parallel_fitness_evaluation(population, ecYeast, ecFSEOF_tabel, GeneProteinMap,
                                                 'r_1589', wildSolution_moma, metabolic_solution)
    # 第0代的最优适应度和最优个体
    max_fitness = max(fitness_scores)
    best_individual = population[fitness_scores.index(max_fitness)]
    mean_fitness = np.mean(fitness_scores)

    # 保存最优适应度和最优个体
    best_fitness_over_generations.append(max_fitness)
    best_individual_over_generations.append(best_individual)
    mean_fitness_over_generations.append(mean_fitness)
    te_init = time.time()

    # 输出这一代的信息
    print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
    print(num_elites)
    print(len(population))
    print(len(fitness_scores))
    print(f"Generation {0}: "
          f"Max Fitness = {max_fitness}\n"
          f"Mean Fitness = {mean_fitness}\n"
          f"Time: {te_init - ts_init}\n"
          f"Solvable Num of individual: {sum(i != 0.01 for i in fitness_scores)}")
    print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')

    for generation in range(MAX_GENERATIONS):
        ti1 = time.time()

        new_population = []
        new_fitness_scores = []

        # 保留精英个体
        elite_indices = np.argsort(fitness_scores)[-num_elites:]
        elites = [population[i] for i in elite_indices]
        elite_fitness_scores = [fitness_scores[i] for i in elite_indices]
        # print(len(elites))
        # print(len(elite_fitness_scores))
        # print(elites)
        # print(elite_fitness_scores)

        population = [population[i] for i in range(len(population)) if i not in elite_indices]
        fitness_scores = [fitness_scores[i] for i in range(len(fitness_scores)) if i not in elite_indices]
        # print(len(population))
        # print(len(fitness_scores))
        # print(population)
        # print(fitness_scores)

        # 选择操作
        remaining_size = POPULATION_SIZE - num_elites
        selected_indices = list(roulette_wheel_selection_POPSIZE(population, fitness_scores, remaining_size))
        # print(remaining_size)
        # print(len(selected_indices))
        if int(len(selected_indices) / 2) != len(selected_indices) / 2:
            selected_indices.append(selected_indices[0])
        # print(selected_indices)

        c_pop = []
        c_fitness_scores = []
        s_population = []
        c_count = 0
        for i in selected_indices:
            s_population.append(population[i])

        for i in range(0, len(selected_indices), 2):
            parent1, parent2 = s_population[i], s_population[i + 1]
            # print(parent1)
            # print(parent2)
            # print(population[selected_indices[i]])
            # print(population[selected_indices[i+1]])
            child1, child2, fg = crossover(parent1, parent2, CROSSOVER_RATE, GENE_LENGTH)
            if fg == 1:
                c_pop.extend([child1, child2])
                c_count += 2
                cs1 = fitness_fun(child1, ecYeast, ecFSEOF_tabel,
                                  GeneProteinMap, 'r_1589', wildSolution_moma, metabolic_solution)
                cs2 = fitness_fun(child2, ecYeast, ecFSEOF_tabel,
                                  GeneProteinMap, 'r_1589', wildSolution_moma, metabolic_solution)
                c_fitness_scores.append(cs1)
                c_fitness_scores.append(cs2)

            else:
                # location_i = selected_indices[i]
                c_pop.append(population[selected_indices[i]])
                c_pop.append(population[selected_indices[i + 1]])
                c_fitness_scores.append(fitness_scores[selected_indices[i]])
                c_fitness_scores.append(fitness_scores[selected_indices[i + 1]])

        m_pop = []
        m_fitness_scores = []
        m_count = 0

        for i in range(len(c_pop)):
            mutated_ind, mflag = mutation(c_pop[i], MUTATION_RATE, GENE_LENGTH)
            if mflag == 1:
                m_count += 1
                m_pop.append(mutated_ind)
                ms = fitness_fun(mutated_ind, ecYeast, ecFSEOF_tabel,
                                 GeneProteinMap, 'r_1589', wildSolution_moma, metabolic_solution)
                m_fitness_scores.append(ms)
            else:
                m_pop.append(c_pop[i])
                m_fitness_scores.append(c_fitness_scores[i])

        new_population.extend(elites)
        # new_population.extend(c_pop)
        new_population.extend(m_pop)

        new_fitness_scores.extend(elite_fitness_scores)
        # new_fitness_scores.extend(c_fitness_scores)
        new_fitness_scores.extend(m_fitness_scores)

        # 更新种群!!!!!!!!!!
        population, fitness_scores = select_fixed_size_population(new_population, new_fitness_scores, POPULATION_SIZE)

        print('-' * 120)
        print('-' * 120)

        # 找出这一代的最优适应度和最优个体
        max_fitness = max(fitness_scores)
        best_individual = population[fitness_scores.index(max_fitness)]
        mean_fitness = np.mean(fitness_scores)

        # 保存最优适应度和最优个体
        best_fitness_over_generations.append(max_fitness)
        best_individual_over_generations.append(best_individual)
        mean_fitness_over_generations.append(mean_fitness)

        ti2 = time.time()

        # 输出这一代的信息
        print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
        print(num_elites)
        print(c_count)
        print(m_count)
        print(len(c_pop))
        print(len(m_pop))
        print(len(new_population))
        print(len(new_fitness_scores))
        print(len(population))
        print(len(fitness_scores))
        print(f"Generation {generation + 1}: "
              f"Max Fitness = {max_fitness}\n"
              f"Mean Fitness = {mean_fitness}\n"
              f"Time: {ti2 - ti1}\n"
              f"Solvable Num of individual: {sum(i != 0.01 for i in fitness_scores)}")
        print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')

        # 收敛检测
        # 如果连续几代的最优适应度没有显著变化，可以提前终止算法
        fitness_change = abs(max_fitness - last_fitness)
        if fitness_change < CONVERGENCE_THRESHOLD:
            convergence_counter += 1
            if convergence_counter >= CONVERGENCE_STEPS:
                print(f"Convergence reached at generation {generation}.")
                break
        else:
            convergence_counter = 0
        last_fitness = max_fitness

    # 在循环结束后，可以输出或保存整个进程的最优结果
    print(f"Best Fitness over all generations: {max(best_fitness_over_generations)}")
    best_overall = \
        best_individual_over_generations[best_fitness_over_generations.index(max(best_fitness_over_generations))]
    print(f"Best Individual over all generations: {best_overall}")

    data = {'individual': best_individual_over_generations, 'best fitness': best_fitness_over_generations,
            'mean fitness': mean_fitness_over_generations}
    df = pd.DataFrame(data)

    # 将DataFrame保存到CSV文件中
    df.to_csv(r"C:\Users\Wenb1n\Desktop\expdata\standardGA.csv", index=False)

    # 绘制最大适应度的变化
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

    print('1')

    # fitness_s = []
    # for individual in population:
    #     fitness = fitness_fun(individual=individual, model0=ecYeast,
    #                           target_table=ecFSEOF_tabel, GeneProteinMap=GeneProteinMap,
    #                           tartgetRxn='r_1589', ref=wildSolution_moma, metrxn=metabolic_solution)
    #     fitness_s.append(fitness)
    #
    # readpd = pd.read_csv(r"C:\Users\Wenb1n\Desktop\expdata\reffluxes.csv", header=None, index_col=0)
    # readse = pd.Series(readpd[1])
    #
    # fitn = parallel_fitness_evaluation(population, ecYeast, ecFSEOF_tabel, GeneProteinMap,
    #                                              'r_1589', wildSolution_moma, readse)
    #
    # print(fitness_scores)
    # print(fitn)







