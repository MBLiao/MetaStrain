import numpy as np
import pandas as pd
import ray
from GA_Operator import *
from ModelOperation import *

GENE_LENGTH = 86        # Gene length


if __name__ == '__main__':
    pd.set_option('display.width', 1000)
    pd.set_option('display.max_columns', 1000)
    np.set_printoptions(linewidth=400, threshold=200)  # np.inf represents positive infinity
    warnings.filterwarnings("ignore")

    ecYeast = loadModel()

    wildSolution_moma = Wild_Growth_sp(ecYeast)

    ecYeast.reactions.get_by_id('r_2111').lower_bound = 0.1
    ecYeast.reactions.get_by_id('r_2051').lower_bound = 0.1
    ecYeast.reactions.get_by_id('r_1714_REV').lower_bound = 0.1

    print(ecYeast.reactions.get_by_id('r_2111').bounds)
    print(ecYeast.reactions.get_by_id('r_2051').bounds)
    print(ecYeast.reactions.get_by_id('r_1714_REV').bounds)

    wildYield = wildSolution_moma.fluxes['r_2051'] / wildSolution_moma.fluxes['r_1714_REV']

    rxnID = []
    for i, x in enumerate(ecYeast.reactions):
        # print(x)
        rxnID.append(x.id)
    rxnID_protein_draw = [x for i, x in enumerate(rxnID) if 'draw_prot_' in x]
    geneAll = []
    for gene in ecYeast.genes:
        # print(gene.id)
        geneAll.append(gene.id)

    metabolic_solution = pd.Series()
    for rid in wildSolution_moma.fluxes.index:
        if rid.startswith('r_'):
            metabolic_solution[rid] = wildSolution_moma.fluxes[rid]

    enzyme_solution = pd.Series()
    for rid in wildSolution_moma.fluxes.index:
        if rid.startswith('draw_'):
            enzyme_solution[rid] = wildSolution_moma.fluxes[rid]

    arm_solution = pd.Series()
    for rid in wildSolution_moma.fluxes.index:
        if rid.startswith('arm_'):
            arm_solution[rid] = wildSolution_moma.fluxes[rid]

    all_reff = pd.Series()
    for rid in wildSolution_moma.fluxes.index:
        all_reff[rid] = wildSolution_moma.fluxes[rid]


    # get the gene list for the related proteins with kcat number
    gene_with_Kcat = [x.replace('draw_prot_', '') for x in rxnID_protein_draw]
    print('----------------------------MODEL INFORMATION----------------------------')
    print('rxnID: ', len(rxnID))
    print('geneAll: ', len(geneAll))
    print('rxnID_protein_draw: ', len(rxnID_protein_draw))
    print('gene_with_Kcat: ', len(gene_with_Kcat))
    print('metabolic reactions: ', len(metabolic_solution))
    print('enzyme reactions: ', len(enzyme_solution))
    print('arm reactions: ', len(arm_solution))
    print('WILD YEILD: ', wildYield)
    print(rxnID_protein_draw[:5])
    print(gene_with_Kcat[:5])
    print('-------------------------------------------------------------------------')

    GeneProteinMap = {}
    df = pd.read_excel(r"C:\Users\Wenb1n\StrainDesign\models\ecModel_batch.xls")
    GeneProteinMap = {gene: enz for gene, enz in zip(df.iloc[7175:8143, 3], df.iloc[7175:8143, 0])}

    ecFSEOF_tabel = pd.read_csv(r"C:\Users\Wenb1n\StrainDesign\models\dimensionality_reduction\ecFSEOF_MAPPING_sp.csv")

    time_all_start = time.time()

    population = ini_pop(GENE_LENGTH)
    # print(len(population))
    # print(len(population[0]))c
    # for i in population:
    #     print(i)
    xlength = []
    ind_fit_all = []

    for x in range(1, 87):
        xlength.append(x)

    for ind in population:
        fit = fitness_f1(ind, ecYeast, ecFSEOF_tabel,
                        GeneProteinMap, 'r_2051', wildSolution_moma, metabolic_solution)
        ind_fit_all.append(fit)

        print('-' * 120)


    data = {'individual': population, 'fitness': ind_fit_all}
    df = pd.DataFrame(data)

    # Save DataFrame to CSV file
    df.to_csv(r"C:\Users\Wenb1n\Desktop\moma_single_sp.csv", index=False)

    # Plot the change in maximum fitness
    plt.plot(xlength, ind_fit_all, color='red')

    plt.title('Best Fitness Over Generations')
    plt.xlabel('Generation')
    plt.ylabel('Fitness')
    plt.show()

    time_all_end = time.time()
    print('Total Time:', time_all_end - time_all_start)

    