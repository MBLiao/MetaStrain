import pandas as pd
from GA_Operator import *
from ModelOperation import *


def f_f(individual, model0, target_table, GeneProteinMap, tartgetRxn, ref, metrxn):
    with model0 as model:
        # 遍历基因编码
        count_zero = np.count_nonzero(individual == 0)/10000
        print(count_zero)
        mutantyield_m = 0

        for i, gene in enumerate(individual):
            if gene == 1:

                # 获取目标基因/酶的名称和调整方式
                target_gene = target_table.iloc[i]["gene_name"]
                operation = target_table.iloc[i]["operation"]
                target_enzyme = GeneProteinMap[target_gene]
                # print(target_gene, target_enzyme)
                # 根据操作调整模型
                if operation == "OE":
                    # print('E' * 12)
                    # print('E' * 12)
                    enzymeUsage = ref[target_enzyme]
                    if enzymeUsage <= 0.00000001:
                        model.reactions.get_by_id(target_enzyme).lower_bound = 0.00000004
                    else:
                        model.reactions.get_by_id(target_enzyme).lower_bound = enzymeUsage * 4
                elif operation == "KD":
                    # print('D' * 12)
                    # print('D' * 12)
                    enzymeUsage = ref[target_enzyme]
                    model.reactions.get_by_id(target_enzyme).upper_bound = enzymeUsage * 0.5
                    # print(model.reactions.get_by_id(target_enzyme).bounds)

                elif operation == "KO":
                    # print('O' * 12)
                    # print('O' * 12)
                    model.reactions.get_by_id(target_enzyme).upper_bound = 0
                    # print(model.reactions.get_by_id(target_enzyme).bounds)
        try:
            ft1 = time.time()
            solution_m = MPMA.moma(model=model, reference_fluxes=metrxn, linear=False)
            ft2 = time.time()

            status = solution_m.status
            if status == 'optimal':
                mutantyield_origin = solution_m.fluxes[tartgetRxn] / solution_m.fluxes['r_1714_REV']
                mutantyield_m = mutantyield_origin + count_zero
                print("---------------------------------------------MOMA----------------------------------------------")
                print("TIME: ", ft2 - ft1)
                growth = solution_m.fluxes['r_2111']
                product = solution_m.fluxes[tartgetRxn]
                glucose = solution_m.fluxes['r_1714_REV']
                print(growth, '|--|', product, '|--|', glucose, '|--|', mutantyield_origin, '|--|', mutantyield_m)
                print('-----------------------------------------------------------------------------------------------')

                # calculate the product and biomass fold change in mutant strain compared with wild strain

            else:
                mutantyield_m = 0.0001
        except:
            mutantyield_m = 0.0001

    return mutantyield_m


def convert_individual(individual_str):
    # 去掉字符串中的方括号和多余空格，并将其转换为整数数组
    individual_str = individual_str.strip('[]').split()
    individual_array = np.array([int(x) for x in individual_str])
    return individual_array


if __name__ == '__main__':
    warnings.filterwarnings("ignore")
    pd.set_option('display.width', 1000)
    pd.set_option('display.max_columns', 1000)
    np.set_printoptions(linewidth=400, threshold=200)  # np.inf表示正无穷

    start_time = time.time()  # Global start time
    ecYeast = loadModel()

    GeneProteinMap = {}
    df = pd.read_excel(r'C:\Users\Wenb1n\StrainDesign\models\ecModel_batch.xls')
    GeneProteinMap = {gene: enz for gene, enz in zip(df.iloc[7175:8143, 3], df.iloc[7175:8143, 0])}

    ecFSEOF_tabel = pd.read_csv(r"C:\Users\Wenb1n\StrainDesign\models\dimensionality_reduction\ecFSEOF_MAPPING.csv")

    rxnID = []
    for i, x in enumerate(ecYeast.reactions):
        rxnID.append(x.id)
    rxnID_protein_draw = [x for i, x in enumerate(rxnID) if 'draw_prot_' in x]
    geneAll = []
    for gene in ecYeast.genes:
        # print(gene.id)
        geneAll.append(gene.id)

    '''
    Set the max upper bound to 1000 instead of inf
    '''
    for re in rxnID:
        if ecYeast.reactions.get_by_id(re).upper_bound == np.inf:
            ecYeast.reactions.get_by_id(re).upper_bound = 1000.0

    ecYeast.reactions.get_by_id('r_2111').lower_bound = 0.1
    ecYeast.reactions.get_by_id('r_1589').lower_bound = 0.1
    ecYeast.reactions.get_by_id('r_1714_REV').lower_bound = 0.1

    read_fluxes_pd = pd.read_csv(r"C:\Users\Wenb1n\Desktop\expdata\binGASOD\remote_result\remote_ref\reffluxes_ne3.csv",
                                 header=None, index_col=0)
    read_all_fluxes = pd.read_csv(r"C:\Users\Wenb1n\Desktop\expdata\binGASOD\remote_result\remote_ref\all_reffluxes_ne3.csv",
                                  header=None, index_col=0)
    read_fluxes = pd.Series(read_fluxes_pd[1])
    read_all_fluxes = pd.Series(read_all_fluxes[1])

    result_data = pd.read_csv(r"C:\Users\Wenb1n\Desktop\expdata\binGASOD\remote_result\expdata\GASD01_ne3.csv")
    read_individuals = result_data['individual']
    unique_individuals = read_individuals.drop_duplicates()
    print(len(unique_individuals))
    converted_individuals = unique_individuals.apply(convert_individual)

    # 检查转换后的个体
    print(converted_individuals)

    max_fit = result_data['best fitness']
    mean_fit = result_data['mean fitness']

    fitness_cheking = []
    for individual in converted_individuals:
        fc = f_f(individual, ecYeast, ecFSEOF_tabel, GeneProteinMap,
                 'r_1589', read_all_fluxes, read_fluxes)
        fitness_cheking.append(fc)

    print(max_fit)
    print(mean_fit)
    print(fitness_cheking)


