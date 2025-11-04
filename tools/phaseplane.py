import cobra
import straindesign as sd
import pandas as pd
import numpy as np
import warnings
import time
import os
from GA_Operator import *
from ModelOperation import *
from ModelOperation_eciML1515 import *


if __name__ == '__main__':
    warnings.filterwarnings("ignore")
    pd.set_option('display.width', 1000)
    pd.set_option('display.max_columns', 1000)
    np.set_printoptions(linewidth=400, threshold=200)

    start_time = time.time()  # Global start time

    # ecModel = load_ec_iML1515()
    ecModel = loadModel()
    # s = S_Matrix(ecYeast)  # Get the S matrix

    rxnID = []
    for i, x in enumerate(ecModel.reactions):
        rxnID.append(x.id)
    rxnID_protein_draw = [x for i, x in enumerate(rxnID) if 'draw_prot_' in x]
    geneAll = []
    for gene in ecModel.genes:
        geneAll.append(gene.id)

    '''
    Set the max upper bound to 1000 instead of inf
    '''
    for re in rxnID:
        if ecModel.reactions.get_by_id(re).upper_bound == np.inf:
            ecModel.reactions.get_by_id(re).upper_bound = 1000.0
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

    print('#####################--MODEL INFORMATION--#####################')
    print('model tolerance:       ', ecModel.tolerance)
    print('rxnID:                 ', len(rxnID))
    print('geneAll:               ', len(geneAll))
    print('rxnID_protein_draw:    ', len(rxnID_protein_draw))
    print('exchanges:             ', len(ecModel.exchanges))
    print('demands:               ', len(ecModel.demands))
    print('sinks:                 ', len(ecModel.sinks))
    print('boundary:              ', len(ecModel.boundary))
    print('###############################################################')
    cons = ['r_1714_REV<=20']  # 'r_1654_REV<=5' 'EX_trp__L_e' EX_glc__D_e_REV
    # cons = []
    datapoints, triangulation, plot_obj = sd.plot_flux_space(ecModel, ('r_1589', 'r_2111', 'r_1714_REV'),
                       constraints=cons, plt_backend='TKAgg', transparent=False)
    print('#####################--PHASE PLANE--#####################')












