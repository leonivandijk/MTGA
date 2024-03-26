import sys
import time
from absl import app
import numpy as np
import pandas as pd
import evaluate_module

# load in gene/wgcna options
TOM_HD = pd.read_csv("./data/tom_hd.csv", delimiter=';')
AD_module = np.loadtxt("./data/lightcyan1.csv", delimiter=';', usecols=1, skiprows=1)

# initialise
np.random.seed(42)
f = evaluate_module.f

# parameters settings
pop_size = 1000
n_generations = 1000
tournament_k = 2
mutation_rate = 1/len(AD_module)
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
    #TODO guided_mutation: write method body. This method uses information from the TOM matrix to guide the mutation process.


def init_population(func, guide = False):
    """
    Constructs the initial population for the Genetic Algorithm. The population consists of binary vectors where
    every element corresponds to a gene. An element gets the value ‘1’ if the gene is part of the given module and
    ‘0’ if not. We want the initial population to explore the neighbourhood of the AD module in the HD brain. We can
    do this by randomly adding and removing genes in the module.
    :return:
    """
    parents = [AD_module] #start with only the AD module
    parents_f = [func(AD_module)]

    while len(parents) != pop_size:
        parent = AD_module.copy()
        if not guide:
            mutation(parent, mutation_rate_init)
        else:
            guided_mutation(parent)
        if evaluate_module.is_valid(parent):
            parents.append(parent)
            parents_f.append(func(parent))

    return parents, parents_f


def roulette_wheel_selection(parent1, parent2):
    #TODO roulette_wheel_selection: write method body
    return selected_parent


def onepoint_crossover(param, param1):
    #TODO onepoint_crossover(): write complete method body
    pass


def genetic_algorithm(func, generations_left = None):
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
            mutation(offspring[i])

        parents = offspring.copy()
        for i in range(pop_size):
            parents_f[i] = func(parents[i])
            generations_left = generations_left - 1
            if parents_f[i] > f_opt:
                f_opt = parents_f[i]
                x_opt = parents[i].copy()

    print(f_opt, x_opt)

    return f_opt, x_opt



def main(argv):
    del argv  # the function doesn't use command-line arguments
    f_opt_result = 0
    # We run the algorithm 100 independent times.
    n_runs = 100
    for _ in range(n_runs):
        f_opt, _ = genetic_algorithm(f)
        f_opt_result += f_opt
        f.reset()
    print("The average f_opt is %s" % (f_opt_result / n_runs))

if __name__ == '__main__':
    start = time.time()
    app.run(main)
    end = time.time()
    print("The program takes %s seconds" % (end - start))