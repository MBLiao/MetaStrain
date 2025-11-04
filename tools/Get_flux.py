from GA_Operator import *
from ModelOperation import *
import time
import matplotlib.pyplot as plt
from JADE_operator import *

np.set_printoptions(precision=15, suppress=True)


if __name__ == '__main__':
    warnings.filterwarnings("ignore")
    pd.set_option('display.width', 1000)
    pd.set_option('display.max_columns', 1000)
    np.set_printoptions(linewidth=400, threshold=200)

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

    # STEP 1
    # ecYeast.reactions.get_by_id('r_2051').bounds = (0.1, 1)  # set the initial rate of the production
    ecYeast.objective = 'r_2111'
    solution1 = ecYeast.optimize()

    all_solution = pd.Series()
    for rid in solution1.fluxes.index:
        all_solution[rid] = solution1.fluxes[rid]

    all_solution.to_csv(r"C:\Users\Wenb1n\Desktop\Get_flux_1E-9_3.csv", header=False, index=True)


