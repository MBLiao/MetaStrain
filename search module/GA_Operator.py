import numpy as np
import random
import matplotlib.pyplot as plt
import pandas as pd
import MPMA
import time


def double_hook_function(x, a, b):
    return a * x + (b / x)


def initialize_pop_loose(population_size, gene_length):
    """
    Initialization
    """
    population = []
    for _ in range(population_size):
        individual = np.zeros(gene_length, dtype=int)
        for i in range(gene_length):
            if random.random() < 0.09:
                individual[i] = 1
        population.append(individual)
    return population


def fitness_fun(individual, model0, target_table, GeneProteinMap, tartgetRxn, ref, metrxn):
    with model0 as model:
        # Iterate through gene encoding
        count_zero = np.count_nonzero(individual == 0)/10000
        print(count_zero)
        mutantyield_m = 0

        for i, gene in enumerate(individual):
            if gene == 1:

                # Get target gene/enzyme name and adjustment method
                target_gene = target_table.iloc[i]["gene_name"]
                operation = target_table.iloc[i]["operation"]
                target_enzyme = GeneProteinMap[target_gene]
                # Adjust model according to operation
                if operation == "OE":
                    enzymeUsage = ref.fluxes[target_enzyme]
                    if enzymeUsage <= 0.00000001:
                        model.reactions.get_by_id(target_enzyme).lower_bound = 0.00000004
                    else:
                        model.reactions.get_by_id(target_enzyme).lower_bound = enzymeUsage * 4
                elif operation == "KD":
                    enzymeUsage = ref.fluxes[target_enzyme]
                    model.reactions.get_by_id(target_enzyme).upper_bound = enzymeUsage * 0.5
                elif operation == "KO":
                    model.reactions.get_by_id(target_enzyme).upper_bound = 0
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
            else:
                mutantyield_m = 0.0001
        except:
            mutantyield_m = 0.0001

    return mutantyield_m


def fitness_f1(individual, model0, target_table, GeneProteinMap, tartgetRxn, ref, metrxn):
    with model0 as model:
        # Iterate through gene encoding
        mutantyield_m = 0
        for i, gene in enumerate(individual):
            if gene == 1:
                # Get target gene/enzyme name and adjustment method
                target_gene = target_table.iloc[i]["gene_name"]
                operation = target_table.iloc[i]["operation"]
                target_enzyme = GeneProteinMap[target_gene]
                # Adjust model according to operation
                if operation == "OE":
                    enzymeUsage = ref.fluxes[target_enzyme]
                    if enzymeUsage <= 0.00000001:
                        model.reactions.get_by_id(target_enzyme).lower_bound = 0.00000004
                    else:
                        model.reactions.get_by_id(target_enzyme).lower_bound = enzymeUsage * 4
                elif operation == "KD":
                    enzymeUsage = ref.fluxes[target_enzyme]
                    model.reactions.get_by_id(target_enzyme).upper_bound = enzymeUsage * 0.5
                elif operation == "KO":
                    model.reactions.get_by_id(target_enzyme).upper_bound = 0
        try:
            ft1 = time.time()
            solution_m = MPMA.moma(model=model, reference_fluxes=metrxn, linear=False)
            ft2 = time.time()

            status = solution_m.status
            if status == 'optimal':
                mutantyield_m = solution_m.fluxes[tartgetRxn] / solution_m.fluxes['r_1714_REV']
                print("---------------------------------------------MOMA----------------------------------------------")
                print("TIME: ", ft2 - ft1)
                growth = solution_m.fluxes['r_2111']
                product = solution_m.fluxes[tartgetRxn]
                glucose = solution_m.fluxes['r_1714_REV']
                print(growth, '|--|', product, '|--|', glucose, '|--|', mutantyield_m)
                print('-----------------------------------------------------------------------------------------------')
            else:
                mutantyield_m = 0.0001
        except:
            mutantyield_m = 0.0001

    return mutantyield_m


def fitness_double_hook(individual, model0, target_table, GeneProteinMap, tartgetRxn, ref, metrxn):
    with model0 as model:
        count_non_zero = np.count_nonzero(individual)
        # print(count_non_zero)
        mutantyield_m = 0

        if count_non_zero == 0:
            return 0.0001

        reg = double_hook_function(count_non_zero, 1, 25)/10000
        # print(reg)

        for i, gene in enumerate(individual):
            if gene == 1:
                # Get target gene/enzyme name and adjustment method
                target_gene = target_table.iloc[i]["gene_name"]
                operation = target_table.iloc[i]["operation"]
                target_enzyme = GeneProteinMap[target_gene]
                # Adjust model according to operation
                if operation == "OE":
                    enzymeUsage = ref.fluxes[target_enzyme]
                    if enzymeUsage <= 0.00000001:
                        model.reactions.get_by_id(target_enzyme).lower_bound = 0.00000004
                    else:
                        model.reactions.get_by_id(target_enzyme).lower_bound = enzymeUsage * 4
                elif operation == "KD":
                    enzymeUsage = ref.fluxes[target_enzyme]
                    model.reactions.get_by_id(target_enzyme).upper_bound = enzymeUsage * 0.5

                elif operation == "KO":
                    model.reactions.get_by_id(target_enzyme).upper_bound = 0
        try:
            ft1 = time.time()
            solution_m = MPMA.moma(model=model, reference_fluxes=metrxn, linear=False)
            ft2 = time.time()

            status = solution_m.status
            if status == 'optimal':
                mutantyield_origin = solution_m.fluxes[tartgetRxn] / solution_m.fluxes['r_1714_REV']
                mutantyield_m = mutantyield_origin - reg
                mutantyield_m = max(mutantyield_m, 0.0001)
                print("---------------------------------------------MOMA----------------------------------------------")
                print("TIME: ", ft2 - ft1)
                growth = solution_m.fluxes['r_2111']
                product = solution_m.fluxes[tartgetRxn]
                glucose = solution_m.fluxes['r_1714_REV']
                print(growth, '|--|', product, '|--|', glucose, '|--|', mutantyield_origin, '|--|', mutantyield_m)
                print('-----------------------------------------------------------------------------------------------')
            else:
                mutantyield_m = 0.0001
        except:
            mutantyield_m = 0.0001

    return mutantyield_m


def fitness_sym(individual, model0, target_table, GeneProteinMap, tartgetRxn, ref, metrxn):
    with model0 as model:
        for i, gene in enumerate(individual):
            if gene != 0:
                target_gene = target_table.iloc[i]["gene_name"]
                target_enzyme = GeneProteinMap[target_gene]
                # Adjust model according to operation
                if gene == 1:
                    enzymeUsage = ref.fluxes[target_enzyme]
                    if enzymeUsage <= 0.00000001:
                        model.reactions.get_by_id(target_enzyme).lower_bound = 0.00000004
                    else:
                        model.reactions.get_by_id(target_enzyme).lower_bound = enzymeUsage * 4
                    # print('E' * 12)
                    # print('E' * 12)
                    # print(model.reactions.get_by_id(target_enzyme).bounds)
                elif gene == 2:
                    enzymeUsage = ref.fluxes[target_enzyme]
                    model.reactions.get_by_id(target_enzyme).upper_bound = enzymeUsage * 0.5
                    # print('D' * 12)
                    # print('D' * 12)
                    # print(model.reactions.get_by_id(target_enzyme).bounds)

                elif gene == 3:
                    model.reactions.get_by_id(target_enzyme).upper_bound = 0
                    # print('O' * 12)
                    # print('O' * 12)
                    # print(model.reactions.get_by_id(target_enzyme).bounds)
        try:
            ft1 = time.time()
            solution_m = MPMA.moma(model=model, reference_fluxes=metrxn, linear=False)
            ft2 = time.time()

            status = solution_m.status
            if status == 'optimal':
                mutantYield_m = solution_m.fluxes[tartgetRxn] / solution_m.fluxes['r_1714_REV']
                print("---------------------------------------------MOMA----------------------------------------------")
                print("TIME: ", ft2 - ft1)
                growth = solution_m.fluxes['r_2111']
                product = solution_m.fluxes[tartgetRxn]
                glucose = solution_m.fluxes['r_1714_REV']
                print(growth, '--', product, '--', glucose, '--', mutantYield_m)
                print('-----------------------------------------------------------------------------------------------')
            else:
                mutantYield_m = 0.0001
        except:
            mutantYield_m = 0.0001

    return mutantYield_m


def fitness_sym_f2(individual, model0, target_table, GeneProteinMap, tartgetRxn, ref, metrxn):
    with model0 as model:
        count_zero = np.count_nonzero(individual == 0)/10000
        print(count_zero)
        mutantyield_m = 0

        for i, gene in enumerate(individual):
            if gene != 0:
                target_gene = target_table.iloc[i]["gene_name"]
                target_enzyme = GeneProteinMap[target_gene]
                # Adjust model according to operation
                if gene == 1:
                    enzymeUsage = ref.fluxes[target_enzyme]
                    if enzymeUsage <= 0.00000001:
                        model.reactions.get_by_id(target_enzyme).lower_bound = 0.00000004
                    else:
                        model.reactions.get_by_id(target_enzyme).lower_bound = enzymeUsage * 4
                    # print('E' * 12)
                    # print('E' * 12)
                    # print(model.reactions.get_by_id(target_enzyme).bounds)
                elif gene == 2:
                    enzymeUsage = ref.fluxes[target_enzyme]
                    model.reactions.get_by_id(target_enzyme).upper_bound = enzymeUsage * 0.5
                    # print('D' * 12)
                    # print('D' * 12)
                    # print(model.reactions.get_by_id(target_enzyme).bounds)

                elif gene == 3:
                    model.reactions.get_by_id(target_enzyme).upper_bound = 0
                    # print('O' * 12)
                    # print('O' * 12)
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
            else:
                mutantyield_m = 0.0001
        except:
            mutantyield_m = 0.0001

    return mutantyield_m


def roulette_wheel_selection_POPSIZE(population, fitness_scores, fixed_size):
    # Ensure population and fitness scores list are of the same length
    assert len(population) == len(fitness_scores), "Population and fitness scores must be of the same length."
    assert len(population) == fixed_size, "Population and fixed_size must be of the same length."

    total_fitness = sum(fitness_scores)
    selection_probs = [f / total_fitness for f in fitness_scores]
    selected_indices = np.random.choice(range(len(population)), size=fixed_size, replace=True, p=selection_probs)

    return selected_indices


# Crossover operation
def crossover(parent1, parent2, crossover_rate, gene_length):
    fg = 0
    if random.random() < crossover_rate:
        crossover_point = random.randint(1, gene_length - 1)
        child1 = np.concatenate([parent1[:crossover_point], parent2[crossover_point:]])
        child2 = np.concatenate([parent2[:crossover_point], parent1[crossover_point:]])
        fg = 1
        return child1, child2, fg
    return parent1, parent2, fg


# Mutation operation
def mutation(individual, mutation_rate, gene_length, per_gene_variation=True):
    mflag = 0
    if per_gene_variation:
        for i in range(len(individual)):
            if random.random() < mutation_rate:
                individual[i] = 1 if individual[i] == 0 else 0
                mflag = 1
    else:
        if random.random() < mutation_rate:
            mp = random.randint(0, gene_length - 1)
            individual[mp] = 1 if individual[mp] == 0 else 0
            mflag = 1

    return individual, mflag


def select_fixed_size_population(population, fitness_scores, fixed_size):
    # Ensure population and fitness scores list are of the same length
    assert len(population) == len(fitness_scores), "Population and fitness scores must be of the same length."
    population_with_scores = list(zip(population, fitness_scores))

    sorted_population_with_scores = sorted(population_with_scores, key=lambda x: x[1], reverse=True)
    selected_population_with_scores = sorted_population_with_scores[:fixed_size]
    selected_population, selected_fitness_scores = zip(*selected_population_with_scores)

    return list(selected_population), list(selected_fitness_scores)


def ini_rand_operation(population_size, gene_length):
    population = []
    for _ in range(population_size):
        individual = np.zeros(gene_length, dtype=int)
        for i in range(gene_length):
            o = random.random()
            if o < 0.03:
                individual[i] = 1
            elif o < 0.06:
                individual[i] = 2
            elif o < 0.09:
                individual[i] = 3
        population.append(individual)
    return population


def mutation_rand(individual, mutation_rate, gene_length, per_gene_variation=False):
    mflag = 0
    if random.random() < mutation_rate:
        for i in range(len(individual)):
            om = random.random()
            if om < 0.01:
                individual[i] = 1
            elif om < 0.02:
                individual[i] = 2
            elif om < 0.03:
                individual[i] = 3
            elif om < 0.05:
                individual[i] = 0
        mflag = 1

    return individual, mflag


def ini_pop(GENE_LENGTH):
    population = []
    for i in range(GENE_LENGTH):
        individual = np.zeros(GENE_LENGTH, dtype=int)
        individual[i] = 1
        population.append(individual)

    return population
