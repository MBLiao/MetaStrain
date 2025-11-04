import pandas as pd
from Benchmark import *


def randFCR(popsize, CRm, CRsigma, Fm, Fsigma):
    CR = np.clip(CRm + CRsigma * np.random.randn(popsize), 0, 1)
    F = np.minimum(1, np.abs(Fm + Fsigma * np.random.standard_cauchy(size=popsize)))
    while np.any(F <= 0):
        F[F <= 0] = np.minimum(1, np.abs(Fm + Fsigma * np.random.standard_cauchy(size=np.sum(F <= 0))))

    return F, CR


def gnR1R2(NP1, NP2):
    r1 = np.random.randint(0, NP1, NP1)
    r2 = np.random.randint(0, NP2, NP1)
    for i in range(NP1):
        while r1[i] == i:
            r1[i] = np.random.randint(0, NP1)
        while r2[i] == r1[i] or r2[i] == i:
            r2[i] = np.random.randint(0, NP2)
    return r1, r2


def boundConstraint(vi, pop, lu):
    xl = np.tile(lu[0], (pop.shape[0], 1))
    xu = np.tile(lu[1], (pop.shape[0], 1))
    vi = np.where(vi < xl, (pop + xl) / 2, vi)
    vi = np.where(vi > xu, (pop + xu) / 2, vi)
    return vi


# def evaluate(pop):
#     eva_scores = []
#     for individual in pop:
#         add = ackley(individual)
#         eva_scores.append(add)
#     return eva_scores


def updateArchive(archive, pop, funvalues):
    popAll = np.vstack([archive['pop'], pop])
    funvaluesAll = np.hstack([archive['funvalues'], funvalues])
    unique_pop, indices = np.unique(popAll, axis=0, return_index=True)
    archive['pop'] = unique_pop
    archive['funvalues'] = funvaluesAll[indices]
    if len(archive['pop']) > archive['NP']:
        indices = np.random.choice(len(archive['pop']), archive['NP'], replace=False)
        archive['pop'] = archive['pop'][indices]
        archive['funvalues'] = archive['funvalues'][indices]
    return archive


def convert_individual(individual_str):
    # Remove brackets and extra spaces from string, and convert to integer array
    individual_str = individual_str.strip('[]').split()
    individual_array = np.array([int(x) for x in individual_str])
    return individual_array


def reverse_map_individual(individual):
    """
    Map discrete values [0, 1, 2, 3] back to random values in continuous interval [0, 4].
    """
    reverse_bins = [0, 1, 2, 3, 4]
    continuous_values = np.array([np.random.uniform(reverse_bins[val], reverse_bins[val + 1])
                                  for val in individual])
    return continuous_values
