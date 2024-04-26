import sys
import time
from absl import app
import numpy as np
import pandas as pd
import evaluate_module
import os

path = "/Users/leonivandijk/Desktop/thesis/pyfiles"

# load in gene/wgcna options
#TOM_HD = pd.read_csv("./data/tom_hd.csv", delimiter=';')
start_module = evaluate_module.start_module

# initialise
np.random.seed(42)
f = evaluate_module.f

# parameters settings
n_variabes = len(start_module)  # might change
pop_size = 100
n_generations = 100
mutation_rate = 1 / n_variabes
mutation_rate_init = 5 * mutation_rate
crossover_probability = 0.6


def mutation(p, rate):
    """
    Random mutation. (De-)Selects every gene with an equal probability equal to the mutation rate.
    :param p:
    :param rate:
    :return:
    """
    for i in range(len(p)):
        if np.random.uniform(0, 1) < rate:
            # bit flipping
            p[i] = 1 - p[i]


def guided_mutation(p):
    # TODO guided_mutation: write method body. This method uses information from the TOM matrix to guide the mutation process.
    p


def init_population(func, guide=False):
    """
    Constructs the initial population for the Genetic Algorithm. The population consists of binary vectors where
    every element corresponds to a gene. An element gets the value ‘1’ if the gene is part of the given module and
    ‘0’ if not. We want the initial population to explore the neighbourhood of the AD module in the HD brain. We can
    do this by randomly adding and removing genes in the module.
    :return:
    """
    parents = [start_module]  # start with only the AD module
    parents_f = [func(start_module)]

    while len(parents) != pop_size:
        parent = start_module.copy()
        if not guide:
            mutation(parent, mutation_rate_init)
        else:
            guided_mutation(parent)
        if evaluate_module.is_valid(parent):
            parents.append(parent)
            parents_f.append(func(parent))

    return parents, parents_f


def roulette_wheel_selection(parent, parent_f):
    # Plusing 0.001 to avoid dividing 0
    f_min = min(parent_f)
    f_sum = sum(parent_f) - (f_min - 0.001) * len(parent_f)

    rw = [(parent_f[0] - f_min + 0.001) / f_sum]
    for i in range(1, len(parent_f)):
        rw.append(rw[i - 1] + (parent_f[i] - f_min + 0.001) / f_sum)

    select_parent = []
    for i in range(len(parent)):
        r = np.random.uniform(0, 1)
        index = 0
        while r > rw[index]:
            index = index + 1

        select_parent.append(parent[index].copy())
    return select_parent


def onepoint_crossover(p1, p2):
    if np.random.uniform(0, 1) < crossover_probability:
        idx = np.random.randint(n_variabes)
        t = p1[idx:]
        p1[idx:] = p2[idx:]
        p2[idx:] = t


def genetic_algorithm(func, generations_left=None):
    if generations_left is None:
        generations_left = n_generations

    f_opt = sys.float_info.min
    x_opt = None

    # construct the initial population
    parents, parents_f = init_population(func)

    while generations_left > 0:
        offspring = roulette_wheel_selection(parents, parents_f)

        for i in range(0, pop_size - (pop_size % 2), 2):
            onepoint_crossover(offspring[i], offspring[i + 1])

        for i in range(pop_size):
            mutation(offspring[i], rate=mutation_rate)

        parents = offspring.copy()
        for i in range(pop_size):
            parents_f[i] = func(parents[i])
            if parents_f[i] > f_opt:
                f_opt = parents_f[i]
                x_opt = parents[i].copy()
        generations_left = generations_left - 1

    print(f_opt, x_opt)


    return f_opt, x_opt


def main(argv):
    del argv  # the function doesn't use command-line arguments
    f_opt_result = 0
    # We run the algorithm 100 independent times.
    n_runs = 100
    for _ in range(n_runs):
        print("start run ", _)
        f_opt, _ = genetic_algorithm(f)
        f_opt_result += f_opt
        f.reset()
    print("The average f_opt is %s" % (f_opt_result / n_runs))


if __name__ == '__main__':
    start = time.time()
    app.run(main)
    end = time.time()
    print("The program takes %s seconds" % (end - start))
