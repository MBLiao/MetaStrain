import pandas as pd
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import escher
import json
import geckopy

from os.path import join, pardir

import cobra
from cobra.io import load_matlab_model, read_sbml_model, load_json_model
from cobra.flux_analysis import moma
from cobra.flux_analysis import pfba
from cobra.util import create_stoichiometric_matrix
import gurobipy as gp

import warnings
import time
import os
import sys

import MPMA


def frange(start, stop, step):
    """
    This function is like range, step can be float value, like 0.1
    :param start:
    :param stop:
    :param step:
    :return:
    """
    i = start
    while i < stop:
        yield i
        i += step


def load_ec_iML1515():
    cobra.Configuration().tolerance = 1e-9

    data_dir = ".\\models\\fixed_eciML1515_batch.mat"
    # data_dir = ".\\models\\fixed_eciML1515_batch_MGHT_10_1.xml"
    ec_model = load_matlab_model(data_dir)

    # ec_model = read_sbml_model(data_dir)
    # data_dir = r".\models\eciML1515.xml.gz"
    # data_dir = "D:\Downloads\eciML1515.xml"
    # data_dir = r".\models\eciML1515_batch.xml"
    # print(cobra.io.sbml.validate_sbml_model(data_dir))

    # ec_model = load_json_model(data_dir)

    ec_model.solver = 'gurobi'
    ec_model.solver.configuration.verbosity = 0
    print(ec_model.objective)

    # Medium
    medium = ec_model.medium
    medium["EX_tyr__L_e_REV"] = 1000.
    medium["EX_phe__L_e_REV"] = 1000.
    medium["EX_cit_e_REV"] = 1000.
    medium["EX_fe3_e_REV"] = 0
    medium["EX_so3_e_REV"] = 0
    # medium["EX_ni2_e_REV"] = 0
    # medium["EX_mobd_e_REV"] = 0
    medium["EX_so2_e_REV"] = 0
    # EX_fe3_e_REV
    # EX_so3_e_REV
    # EX_ni2_e_REV
    # EX_mobd_e_REV
    # EX_so2_e_REV
    ec_model.medium = medium

    ec_model.reactions.get_by_id('draw_prot_P07023').upper_bound = 0  # b2600
    ec_model.reactions.get_by_id('draw_prot_P0A853').upper_bound = 0  # b3708
    ec_model.reactions.get_by_id('draw_prot_P0A9J8').upper_bound = 0  # b2599

    return ec_model


def saveExcel(infile, outfile):
    writer = pd.ExcelWriter(outfile)
    infile.to_excel(writer, 'Sheet1')
    writer.save()


def S_Matrix(Model):
    with Model as model:
        s = create_stoichiometric_matrix(model)
        plt.figure(dpi=300)
        plt.spy(s, markersize=.5)
        plt.show()
    return s


def Wild_Growth_Trp(Model):
    with Model as model:
        model.reactions.get_by_id('EX_trp__L_e').bounds = (0.02, 0.021)  # set the initial rate of the production
        model.reactions.get_by_id('EX_glc__D_e_REV').bounds = (1, 1.8)  # Glucose

        model.objective = 'BIOMASS_Ec_iML1515_core_75p37M'
        solution1 = model.optimize()
        growth = float(solution1.objective_value)

        print("-------------REFERENCE-------------")
        print('Biomass: ', solution1.fluxes['BIOMASS_Ec_iML1515_core_75p37M'])
        print('Tryptophan: ', solution1.fluxes['EX_trp__L_e'])
        print('glucose uptake: ', solution1.fluxes['EX_glc__D_e_REV'])
        print('enzymeUsage: ', solution1.fluxes['prot_pool_exchange'])
        print('growth: ', growth)
        print('-----------------------------------')

    return solution1


def max_growth(Model):
    with Model as model:
        # model.reactions.get_by_id('EX_trp__L_e').bounds = (0, 1000)  # set the initial rate of the production
        # model.reactions.get_by_id('EX_glc__D_e_REV').bounds = (0, 0.8)  # Glucose

        model.objective = 'BIOMASS_Ec_iML1515_core_75p37M'
        solution1 = model.optimize()
        growth = float(solution1.objective_value)

        print("-------------Max Growth-------------")
        print('Biomass: ', solution1.fluxes['BIOMASS_Ec_iML1515_core_75p37M'])
        print('Tryptophan: ', solution1.fluxes['EX_trp__L_e'])
        print('glucose uptake: ', solution1.fluxes['EX_glc__D_e_REV'])
        print('enzymeUsage: ', solution1.fluxes['prot_pool_exchange'])
        print('growth: ', growth)
        print('------------------------------------')

    return None

def max_product(Model):
    with Model as model:
        # model.reactions.get_by_id('EX_trp__L_e').bounds = (0, 1000)  # set the initial rate of the production
        # model.reactions.get_by_id('EX_glc__D_e_REV').bounds = (0, 0.8)  # Glucose

        model.objective = 'EX_trp__L_e'
        solution2 = model.optimize()
        product = float(solution2.objective_value)

        print("-------------Max Product-------------")
        print('Biomass: ', solution2.fluxes['BIOMASS_Ec_iML1515_core_75p37M'])
        print('Tryptophan: ', solution2.fluxes['EX_trp__L_e'])
        print('glucose uptake: ', solution2.fluxes['EX_glc__D_e_REV'])
        # print('ethanol uptake: ', solution1.fluxes['r_1761_REV'])
        print('enzymeUsage: ', solution2.fluxes['prot_pool_exchange'])
        print('product: ', product)
        print('-------------------------------------')
    return None



def WildModel_Growth(Model):
    with Model as model:
        # model.reactions.get_by_id("r_1714_REV").bounds = (0, 5)  # open the glucose uptake rate

        # STEP 1
        model.reactions.get_by_id('r_1589').bounds = (2, 2.5)  # set the initial rate of the production
        model.objective = 'r_2111'
        t_g = time.time()
        solution1 = model.optimize()
        t_ge = time.time()
        growth0 = float(solution1.objective_value)

        print("-------------STEP 1/1-------------")
        print('TIME: ', t_ge - t_g)
        print('fluxes[r_2111]: ', solution1.fluxes['r_2111'])
        print('fluxes[r_1589]: ', solution1.fluxes['r_1589'])
        print('glucose uptake: ', solution1.fluxes['r_1714_REV'])
        print('enzymeUsage: ', solution1.fluxes['prot_pool_exchange'])
        print('growth0: ', growth0)
        print('----------------------------------')

    return solution1
