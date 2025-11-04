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


def loadModel():
    cobra.Configuration().tolerance = 1e-9

    data_dir = r".\models\ecYeastGEM_batch.mat"
    model = load_matlab_model(data_dir)
    model.solver = 'gurobi'
    model.solver.configuration.verbosity = 0

    return model


def load_ec_iML1515():
    cobra.Configuration().tolerance = 1e-9

    data_dir = "D:\\PycharmProjects\\SDanalysis\models\\fixed_eciML1515_batch.mat"
    ec_model = load_matlab_model(data_dir)

    # data_dir = r".\models\eciML1515.xml.gz"
    # data_dir = "D:\下载\eciML1515.xml"
    # data_dir = r".\models\eciML1515_batch.xml"
    # print(cobra.io.sbml.validate_sbml_model(data_dir))
    # ec_model = read_sbml_model(data_dir)
    # ec_model = load_json_model(data_dir)

    ec_model.solver = 'gurobi'
    ec_model.solver.configuration.verbosity = 0

    return ec_model


def upgradeModel():
    cobra.Configuration().tolerance = 1e-9

    data_dir = "./models/ecYeastGEM_batch.mat"
    model = load_matlab_model(data_dir)
    model.solver = 'gurobi'
    model.solver.configuration.verbosity = 0

    # First reaction
    reaction = cobra.Reaction('r_Amorpha_diene_sys')
    reaction.name = 'Amorpha_diene systhesis'
    reaction.subsystem = 'Cell cytoplasm Biosynthesis'
    reaction.lower_bound = 0.
    reaction.upper_bound = 1000.

    Amorpha_diene = cobra.Metabolite(
        'Amorpha_diene',
        formula='C15H24',
        name='Amorpha_diene',
        compartment='c')

    # Enzyme metabolite
    '''prot_ADS = Metabolite(
        'prot_ADS',
        formula='',
        name='ADS enzyme',
        compartment='c')'''

    reaction.add_metabolites({
        model.metabolites.get_by_id('s_0190'): -1.0,  # FPP
        Amorpha_diene: 1.0,  # Amorpha-4,11-diene
        model.metabolites.get_by_id('s_0633'): 1.0,  # Pyrophosphate
        # prot_ADS: -1.9305/3600                    # Here we temporarily make it consistent with the following kcat
    })
    model.add_reactions([reaction])

    # Second reaction
    reaction = cobra.Reaction('r_artemisinic_alcohol_sys')
    reaction.name = 'artemisinic_alcohol systhesis'
    reaction.subsystem = 'Cell cytoplasm Biosynthesis'
    reaction.lower_bound = 0.
    reaction.upper_bound = 1000.

    artemisinic_alcohol = cobra.Metabolite(
        'artemisinic_alcohol',
        formula='C15H24O',
        name='artemisinic_alcohol',
        compartment='c')

    # Enzyme metabolite
    '''prot_ALDH1 = Metabolite(
        'prot_ALDH1',
        formula='',
        name='prot_ALDH1',
        compartment='c')
    '''
    reaction.add_metabolites({
        Amorpha_diene: -1.0,  # Amorpha-4,11-diene
        model.metabolites.get_by_id('s_1275'): -1.0,  # O2
        model.metabolites.get_by_id('s_1212'): -1.0,  # NADPH
        model.metabolites.get_by_id('s_1207'): 1.0,  # NADP+
        model.metabolites.get_by_id('s_0794'): -1.0,  # H+
        model.metabolites.get_by_id('s_0803'): 1.0,  # H2O
        artemisinic_alcohol: 1.0,  # Artemisinic alcohol
        # prot_ALDH1: -0.65359/3600/4                        # 1/kcat for step 2 enzyme
    })
    model.add_reactions([reaction])

    # Third reaction
    reaction = cobra.Reaction('r_artemisinic_aldehyde_sys')
    reaction.name = 'artemisinic_aldehyde systhesis'
    reaction.subsystem = 'Cell cytoplasm Biosynthesis'
    reaction.lower_bound = 0.
    reaction.upper_bound = 1000.

    artemisinic_aldehyde = cobra.Metabolite(
        'artemisinic_aldehyde',
        formula='C15H22O',
        name='artemisinic_aldehyde',
        compartment='c')

    # Enzyme metabolite
    '''prot_ALDH1 = Metabolite(
        'prot_ADS',
        formula='',
        name='ADS enzyme',
        compartment='c')'''

    reaction.add_metabolites({
        artemisinic_alcohol: -1.0,  # Artemisinic alcohol
        model.metabolites.get_by_id('s_1275'): -1.0,  # O2
        model.metabolites.get_by_id('s_1212'): -1.0,  # NADPH
        model.metabolites.get_by_id('s_1207'): 1.0,  # NADP+
        model.metabolites.get_by_id('s_0794'): -1.0,  # H+
        model.metabolites.get_by_id('s_0803'): 2.0,  # H2O
        artemisinic_aldehyde: 1.0,  # Artemisinic aldehyde
        # prot_ALDH1: -0.65359/3600/4                                # 1/kcat for step 3 enzyme
    })
    model.add_reactions([reaction])

    # Fourth reaction
    reaction = cobra.Reaction('r_artemisinic_acid_sys')
    reaction.name = 'artemisinic_acid systhesis'
    reaction.subsystem = 'Cell cytoplasm Biosynthesis'
    reaction.lower_bound = 0.
    reaction.upper_bound = 1000.

    artemisinic_acid = cobra.Metabolite(
        'artemisinic_acid',
        formula='C15H21O2',
        name='artemisinic_acid',
        compartment='c')

    # Enzyme metabolite
    '''prot_ALDH1 = Metabolite(
        'prot_ADS',
        formula='',
        name='ADS enzyme',
        compartment='c')'''

    reaction.add_metabolites({
        artemisinic_aldehyde: -1.0,  # Artemisinic aldehyde
        model.metabolites.get_by_id('s_1275'): -1.0,  # O2
        model.metabolites.get_by_id('s_1212'): -1.0,  # NADPH
        model.metabolites.get_by_id('s_1207'): 1.0,  # NADP+
        model.metabolites.get_by_id('s_0794'): -1.0,  # H+
        model.metabolites.get_by_id('s_0803'): 1.0,  # H2O
        artemisinic_acid: 1.0,  # Artemisinic acid
        # prot_ALDH1: -0.65359/3600/4                        # 1/kcat for step 4 enzyme
    })
    model.add_reactions([reaction])

    model.add_boundary(model.metabolites.get_by_id("artemisinic_acid"), type="sink")

    medium = model.medium
    medium["r_1714_REV"] = 0
    medium["r_1761_REV"] = 1000.
    # model.reactions.get_by_id('draw_prot_Q05521').upper_bound = 0

    model.medium = medium

    return model


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


def wild_glucose_constraint(Model):
    with Model as model:
        t_g = time.time()
        model.objective = 'r_2111'
        solution1 = model.optimize()
        growth0 = float(solution1.objective_value)
        model.reactions.get_by_id('r_2111').bounds = (0.90 * growth0, growth0)
        model.objective = 'r_1589'
        solution2 = model.optimize()
        t_ge = time.time()
        growth1 = float(solution2.objective_value)

        print("-------------STEP 1/1-------------")
        print('TIME: ', t_ge - t_g)
        print('fluxes[r_2111]: ', solution2.fluxes['r_2111'])
        print('fluxes[r_1589]: ', solution2.fluxes['r_1589'])
        print('glucose uptake: ', solution2.fluxes['r_1714_REV'])
        print('enzymeUsage: ', solution2.fluxes['prot_pool_exchange'])
        print('growth1: ', growth1)
        print('----------------------------------')

    return solution2


def Wild_Growth_sp(Model):
    with Model as model:
        # model.reactions.get_by_id("r_1714_REV").bounds = (0, 5)  # open the glucose uptake rate

        # STEP 1
        model.reactions.get_by_id('r_2051').bounds = (0.1, 1)  # set the initial rate of the production
        model.objective = 'r_2111'
        t_g = time.time()
        solution1 = model.optimize()
        t_ge = time.time()
        growth0 = float(solution1.objective_value)

        print("-------------STEP 1/1-------------")
        print('TIME: ', t_ge - t_g)
        print('fluxes[r_2111]: ', solution1.fluxes['r_2111'])
        print('fluxes[r_2051]: ', solution1.fluxes['r_2051'])
        print('glucose uptake: ', solution1.fluxes['r_1714_REV'])
        print('enzymeUsage: ', solution1.fluxes['prot_pool_exchange'])
        print('growth0: ', growth0)
        print('----------------------------------')

    return solution1


def Wild_Growth_AA(Model):
    with Model as model:
        # STEP 1
        model.reactions.get_by_id('r_artemisinic_acid_sys').bounds = (0.05, 1000)  # set the initial rate of the production
        model.reactions.get_by_id('r_1761_REV').bounds = (0, 0.8)
        model.objective = 'r_2111'
        t_g = time.time()
        solution1 = model.optimize()
        t_ge = time.time()
        growth0 = float(solution1.objective_value)

        print("-------------STEP 1/1-------------")
        print('TIME: ', t_ge - t_g)
        print('fluxes[r_2111]: ', solution1.fluxes['r_2111'])
        print('fluxes[r_artemisinic_acid_sys]: ', solution1.fluxes['r_artemisinic_acid_sys'])
        print('glucose uptake: ', solution1.fluxes['r_1714_REV'])
        print('ethanol uptake: ', solution1.fluxes['r_1761_REV'])
        print('enzymeUsage: ', solution1.fluxes['prot_pool_exchange'])
        print('growth0: ', growth0)
        print('----------------------------------')

    return solution1


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


def WildModel_MOMA(Model):
    with Model as model:
        print('--------MODEL INFORMATION--------')
        print('Reactions:', len(model.reactions))
        print('Metabolites:', len(model.metabolites))
        print('Genes', len(model.genes))
        print('Glucose uptake bound: ', model.reactions.get_by_id("r_1714_REV").bounds)
        print('-' * 40)

        # Basic model solving
        print("--------BASIC RESULT1--------")
        model.objective = 'r_2111'
        solution1 = model.optimize()
        print(solution1)
        print('fluxes[r_2111]: ', solution1.fluxes['r_2111'])
        print('fluxes[r_1589]: ', solution1.fluxes['r_1589'])
        print('glucose uptake: ', solution1.fluxes['r_1714_REV'])
        print('enzymeUsage: ', solution1.fluxes['prot_pool_exchange'])
        print('enzymeUsage(P11154): ', solution1.fluxes['draw_prot_P11154'])
        print('---------------------------------')

        print("--------BASIC RESULT2--------")
        model.objective = 'r_1589'
        solution2 = model.optimize()
        print(solution2)
        print('fluxes[r_2111]: ', solution2.fluxes['r_2111'])
        print('fluxes[r_1589]: ', solution2.fluxes['r_1589'])
        print('glucose uptake: ', solution2.fluxes['r_1714_REV'])
        print('enzymeUsage: ', solution2.fluxes['prot_pool_exchange'])
        print('enzymeUsage(P11154): ', solution2.fluxes['draw_prot_P11154'])
        print('---------------------------------')

        print("--------BASIC RESULT3--------")
        model.objective = {model.reactions.prot_pool_exchange: -1}
        solution3 = model.optimize()
        print(solution3)
        print('fluxes[r_2111]: ', solution3.fluxes['r_2111'])
        print('fluxes[r_1589]: ', solution3.fluxes['r_1589'])
        print('glucose uptake: ', solution3.fluxes['r_1714_REV'])
        print('enzymeUsage: ', solution3.fluxes['prot_pool_exchange'])
        print('enzymeUsage(P11154): ', solution3.fluxes['draw_prot_P11154'])
        print('---------------------------------')

        # model.reactions.get_by_id("r_1714_REV").bounds = (0, 5)  # open the glucose uptake rate

        # STEP 1
        model.reactions.get_by_id('r_1589').bounds = (2, 2.5)  # set the initial rate of the production
        model.objective = 'r_2111'
        solution4 = model.optimize()
        growth0 = float(solution4.objective_value)
        print("--------STEP 1--------")
        print('fluxes[r_2111]: ', solution4.fluxes['r_2111'])
        print('fluxes[r_1589]: ', solution4.fluxes['r_1589'])
        print('glucose uptake: ', solution4.fluxes['r_1714_REV'])
        print('enzymeUsage: ', solution4.fluxes['prot_pool_exchange'])

        print('growth0: ', growth0)
        print('---------------------------------')

        # STEP 2
        # fix growth0 and minimization the glucose uptake
        model.reactions.get_by_id('r_2111').bounds = (0.90 * growth0, growth0)
        model.objective = {model.reactions.get_by_id('r_1714_REV'): -1}
        solution5 = model.optimize()
        glucose0 = -float(solution5.objective_value)
        print("--------STEP 2--------")
        print('fluxes[r_2111]: ', solution5.fluxes['r_2111'])
        print('fluxes[r_1589]: ', solution5.fluxes['r_1589'])
        print('glucose uptake: ', solution5.fluxes['r_1714_REV'])
        print('enzymeUsage: ', solution5.fluxes['prot_pool_exchange'])

        print('glucose0: ', glucose0)
        print('---------------------------------')

        # STEP 3
        # fix glucose uptake and minimization the protein pool
        model.reactions.get_by_id('r_1714_REV').bounds = (glucose0, 1.1 * glucose0)
        model.objective = {model.reactions.prot_pool_exchange: -1}
        solution6 = model.optimize()
        print("--------STEP 3--------")
        print('fluxes[r_2111]: ', solution6.fluxes['r_2111'])
        print('fluxes[r_1589]: ', solution6.fluxes['r_1589'])
        print('glucose uptake: ', solution6.fluxes['r_1714_REV'])
        print('enzymeUsage: ', solution6.fluxes['prot_pool_exchange'])
        print('---------------------------------')

        print('*' * 40)
        print('WILD: ', solution6.fluxes['r_1589'], solution6.fluxes['r_1714_REV'])
        wildYield = solution6.fluxes['r_1589'] / solution6.fluxes['r_1714_REV']
        BPCY = wildYield * solution6.fluxes['r_2111']
        print('wildYield: ', wildYield)
        print('BPCY_W: ', BPCY)
        print('*' * 40)

        print(model.reactions.get_by_id('r_1714_REV').bounds)
        print(model.reactions.get_by_id('r_2111').bounds)

    return solution6


def WildModel_FBA(Model):
    model = Model.copy()
    with Model as model:
        # model.reactions.get_by_id("r_1714_REV").bounds = (0, 5)  # open the glucose uptake rate

        # STEP 1
        model.objective = 'r_2111'
        solution4 = model.optimize()
        growth0 = float(solution4.objective_value)
        print("--------STEP 1--------")
        print('fluxes[r_2111]: ', solution4.fluxes['r_2111'])
        print('fluxes[r_1589]: ', solution4.fluxes['r_1589'])
        print('glucose uptake: ', solution4.fluxes['r_1714_REV'])
        print('enzymeUsage: ', solution4.fluxes['prot_pool_exchange'])

        print('growth0: ', growth0)
        print('---------------------------------')

        # STEP 2
        # fix growth0 and minimization the glucose uptake
        model.reactions.get_by_id('r_2111').bounds = (0.80 * growth0, growth0)
        model.objective = {model.reactions.get_by_id('r_1714_REV'): -1}
        solution5 = model.optimize()
        glucose0 = -float(solution5.objective_value)
        print("--------STEP 2--------")
        print('fluxes[r_2111]: ', solution5.fluxes['r_2111'])
        print('fluxes[r_1589]: ', solution5.fluxes['r_1589'])
        print('glucose uptake: ', solution5.fluxes['r_1714_REV'])
        print('enzymeUsage: ', solution5.fluxes['prot_pool_exchange'])

        print('glucose0: ', glucose0)
        print('---------------------------------')

        # STEP 3
        # fix glucose uptake and minimization the protein pool
        model.reactions.get_by_id('r_1714_REV').bounds = (glucose0, 1.2 * glucose0)
        model.objective = {model.reactions.prot_pool_exchange: -1}
        solution6 = model.optimize()
        protein_pool0 = -float(solution6.objective_value)
        print("--------STEP 3--------")
        print('fluxes[r_2111]: ', solution6.fluxes['r_2111'])
        print('fluxes[r_1589]: ', solution6.fluxes['r_1589'])
        print('glucose uptake: ', solution6.fluxes['r_1714_REV'])
        print('enzymeUsage: ', solution6.fluxes['prot_pool_exchange'])
        print('---------------------------------')

        # STEP 4
        # fix the protein pool and max product
        model.reactions.get_by_id('prot_pool_exchange').bounds = (protein_pool0, 1.2 * protein_pool0)
        model.objective = {model.reactions.get_by_id('r_1589'): 1}
        solution7 = model.optimize()
        print("--------STEP 4--------")
        print('fluxes[r_2111]: ', solution7.fluxes['r_2111'])
        print('fluxes[r_1589]: ', solution7.fluxes['r_1589'])
        print('glucose uptake: ', solution7.fluxes['r_1714_REV'])
        print('enzymeUsage: ', solution7.fluxes['prot_pool_exchange'])
        print('---------------------------------')


        print('*' * 40)
        print('WILD: ', solution7.fluxes['r_1589'], solution7.fluxes['r_1714_REV'])
        wildYield = solution7.fluxes['r_1589'] / solution7.fluxes['r_1714_REV']
        BPCY = wildYield * solution7.fluxes['r_2111']
        print('wildYield: ', wildYield)
        print('BPCY_W: ', BPCY)
        print('*' * 40)

        print(model.reactions.get_by_id('r_1714_REV').bounds)
        print(model.reactions.get_by_id('r_2111').bounds)

    return solution7
