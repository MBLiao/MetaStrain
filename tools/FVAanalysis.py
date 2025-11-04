import pandas as pd
import ray
from GA_Operator import *
from ModelOperation import *


if __name__ == '__main__':
    warnings.filterwarnings("ignore")
    pd.set_option('display.width', 1000)
    pd.set_option('display.max_columns', 1000)
    np.set_printoptions(linewidth=400, threshold=200)

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

    # wildSolution_moma = WildModel_Growth(ecYeast)

    # bound1
    # ecYeast.reactions.get_by_id('r_2111').lower_bound = 0.1
    # ecYeast.reactions.get_by_id('r_1589').lower_bound = 0.1
    # ecYeast.reactions.get_by_id('r_1714_REV').lower_bound = 0.1
    #
    # print(ecYeast.reactions.get_by_id('r_2111').bounds)
    # print(ecYeast.reactions.get_by_id('r_1589').bounds)
    # print(ecYeast.reactions.get_by_id('r_1714_REV').bounds)

    # bound2
    # ecYeast.reactions.get_by_id('r_1589').lower_bound = 2

    # bound3
    # ecYeast.reactions.get_by_id('r_1589').lower_bound = 2
    # ecYeast.reactions.get_by_id('r_1714_REV').upper_bound = 20

    # bound4
    # ecYeast.reactions.get_by_id('r_1714_REV').upper_bound = 1

    # bound5
    # ecYeast.reactions.get_by_id('r_1714_REV').upper_bound = 5

    # bound6
    # ecYeast.reactions.get_by_id('r_1714_REV').upper_bound = 10


    print('#####################--MODEL INFORMATION--#####################')
    print('model tolerance:       ', ecYeast.tolerance)
    print('rxnID:                 ', len(rxnID))
    print('geneAll:               ', len(geneAll))
    print('rxnID_protein_draw:    ', len(rxnID_protein_draw))
    print('exchanges:             ', len(ecYeast.exchanges))
    print('demands:               ', len(ecYeast.demands))
    print('sinks:                 ', len(ecYeast.sinks))
    print('boundary:              ', len(ecYeast.boundary))
    print('###############################################################')

    with ecYeast as model0:
        fva_result = cobra.flux_analysis.flux_variability_analysis(model0, fraction_of_optimum=1)
        print(fva_result)
        fva_result.to_csv(r"C:\Users\Wenb1n\Desktop\SparseFluxOpt\fva_100000%_1E-9.csv", index=True)








