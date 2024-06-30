import argparse
import random
import sys
import time

import numpy as np
from absl import app

import evaluate_module


def mutation(x, rate):
    """
    In-place random mutation. (De-)Selects every gene with an equal probability equal to the mutation rate.
    :param x: A binary array representing a module
    :param rate: The mutation rate
    """
    for i in range(n_variables):
        if np.random.uniform(0, 1) < rate:
            # bit flipping
            x[i] = 1 - x[i]


def guided_mutation(x, rate, tom):
    """
    In-place guided mutation. (De-)Selects genes with a probability proportional to its connectivity to one of the
    genes already in the module.
    :param x: A binary array representing a module
    :param rate: The mutation rate
    :param tom: The topological overlap matrix of the disease network
    """

    n_mutate = int(np.round(rate * n_variables))

    for i in range(n_mutate):
        module_genes = np.flatnonzero(x)
        module_size = sum(x)
        # select a node in the module
        selected_index = int(np.round(float(np.random.uniform(0, 1) * (module_size - 1))))
        selected_gene = int(module_genes[selected_index])

        # choice between add/remove depends on the size of the current module
        # select gene to be added; higher tom-similarity to node get higher selection probability
        if np.random.uniform(0, 1) > (2 * module_size / n_variables):
            row = tom[selected_gene, :].copy()
            row[selected_gene] = 0
            prob_add = row / np.sum(row)
            add = np.random.choice(n_variables, p=prob_add)
            x[add] = 1

        # or select a node to remove; lower tom-similarity to node gets higher deselection probability
        else:
            to_module = tom[x.astype(bool).squeeze()][:, x.astype(bool).squeeze()]
            row = to_module[selected_index, :].copy()
            row[selected_index] = 1
            row = 1 - row
            prob_remove = row / np.sum(row)
            remove = np.random.choice(module_size, p=prob_remove)
            remove_gene = int(module_genes[remove])
            x[remove_gene] = 0


def init_population(func, start_module, tom, guide=False):
    """
    Constructs the initial population for the Genetic Algorithm. The population consists of binary vectors where
    every element corresponds to a gene. An element gets the value ‘1’ if the gene is part of the given module and
    ‘0’ if not. We want the initial population to explore the neighbourhood of the AD module in the HD brain. We can
    do this by randomly adding and removing genes in the module.
    :param func: The fitness function as defined in evaluate_module.py
    :param start_module: The module found by WGCNA which we want to optimize
    :param tom: The topological overlap matrix of the disease network
    :param guide: Method of mutation; "False" is random mutation, "True" is guided mutation.
    :return: The initial parent population with its function values
    """
    parents = [start_module]  # start with only the AD module
    parents_f = [func(start_module)]  # compute fitness value

    while len(parents) != pop_size:
        parent = start_module.copy()
        if not guide:
            mutation(parent, mutation_rate_init)
        else:
            guided_mutation(parent, rate=mutation_rate_init, tom=tom)
        if evaluate_module.is_valid(parent):
            parents.append(parent)
            parents_f.append(func(parent))

    return parents, parents_f


def roulette_wheel_selection(parent, parent_f):
    """
    Selection operator. The probability of selection is proportional to each individual's fitness, allowing higher
    fitness individuals a greater chance of being selected.
    :param parent: The population
    :param parent_f: The fitness values of the population
    :return: The population after the selection operator
    """
    f_min = min(parent_f)
    f_sum = sum(parent_f) - (f_min - 0.001) * pop_size  # Plusing 0.001 to avoid dividing 0

    # compute probabilities for every individual
    rw = [(parent_f[0] - f_min + 0.001) / f_sum]
    for i in range(1, pop_size):
        rw.append(rw[i - 1] + (parent_f[i] - f_min + 0.001) / f_sum)

    select_parent = []
    for i in range(pop_size):
        r = np.random.uniform(0, 1)
        index = 0
        while r > rw[index]:
            index = index + 1

        select_parent.append(parent[index].copy())
    return select_parent


def tournament_selection(parent, parent_f):
    """
    Selection operator. Selects individuals from a population by running several "tournaments" among a few individuals
    chosen at random. The winner of each tournament, typically the individual with the highest fitness,
    is selected for the next generation.
    :param parent: The population
    :param parent_f: The fitness values of the population
    :return: The population after the selection operator
    """
    select_parent = []
    for i in range(pop_size):
        # select five individuals at random
        pre_select = np.random.choice(pop_size, 5, replace=False)
        max_f = parent_f[pre_select[0]]
        index = pre_select[0]
        for p in pre_select:
            if parent_f[p] > max_f:
                index = p
                max_f = parent_f[p]
        select_parent.append(parent[index].copy())
    return select_parent


def npoint_crossover(n, x1, x2, crossover_probability):
    """
    In-place crossover operator. Combines genetic material from two parent individuals by selecting N crossover points,
    swapping segments between the parents to produce offspring.
    :param n: The number of crossover points
    :param x1: A binary array representing module 1
    :param x2: A binary array representing module 2
    """
    if np.random.uniform(0, 1) < crossover_probability:
        idxs = sorted(random.sample(range(n_variables), n))
        for idx in idxs:
            t = x1[idx:]
            x1[idx:] = x2[idx:]
            x2[idx:] = t


def uniform_crossover(x1, x2, crossover_probability):
    """
    In-place crossover operator. Combines genetic material from two parent individuals by randomly selecting genes
    from each parent with a fixed probability.
    :param x1: A binary array representing module 1
    :param x2: A binary array representing module 2
    """
    if np.random.uniform(0, 1) < crossover_probability:
        for i in range(n_variables):
            if np.random.uniform(0, 1) < 0.5:
                t = x1[i]
                x1[i] = x2[i]
                x2[i] = t


def genetic_algorithm(func, start_module, tom, generations_left=None):
    """
    Executes the genetic algorithm.
    :param func: The fitness function as defined in evaluate_module.py
    :param start_module: The module found by WGCNA which we want to optimize
    :param tom: The topological overlap matrix of the disease network
    :param generations_left: The duration of the algorithm
    :return: The best found solution
    """

    if generations_left is None:
        generations_left = n_generations

    f_opt = sys.float_info.min
    x_opt = None

    # construct the initial population
    parents, parents_f = init_population(func, start_module, tom)

    while generations_left > 0:
        # 1. selection
        offspring = tournament_selection(parents, parents_f)
        # 2. crossover: higher probability when population gets more diverse
        crossover_probability = crossover_probability_start * ((n_generations - generations_left) / n_generations * 0.6)
        for i in range(0, pop_size - (pop_size % 2), 2):
            uniform_crossover(offspring[i], offspring[i + 1], crossover_probability)
        # 3. mutation
        for i in range(pop_size):
            guided_mutation(offspring[i], rate=mutation_rate, tom=tom)

        parents = offspring.copy()
        for i in range(pop_size):
            # calculate the fitness scores of the new population
            parents_f[i] = func(parents[i])
            if parents_f[i] > f_opt:
                f_opt = parents_f[i]
                x_opt = parents[i].copy()
        generations_left = generations_left - 1

    print(f_opt, x_opt)

    return f_opt, x_opt


def main():
    # save variables as local variables
    start_module = evaluate_module.start_module
    tom = evaluate_module.tom

    f = evaluate_module.f
    f_opt_result = 0

    # We run the algorithm n_runs independent times.
    for _ in range(n_runs):
        print("start run", _)
        start_subrun = time.time()
        f_opt, _ = genetic_algorithm(func=f, start_module=start_module, tom=tom)
        f_opt_result += f_opt
        f.reset()
        end_subrun = time.time()
        print((end_subrun - start_subrun), "seconds")
    print("The average f_opt is %s" % (f_opt_result / n_runs))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--disease',
                        help="Choose between AD or HD",
                        required=True)
    parser.add_argument('--seed',
                        default=42,
                        help="Choose a random seed. Default: 42")
    args, unknown = parser.parse_known_args()

    # initialise data and parameters
    evaluate_module.load_data(args.disease)
    n_variables = len(evaluate_module.search_space)
    np.random.seed(args.seed)
    pop_size = 150
    n_generations = 250
    crossover_probability_start = 1 / pop_size
    mutation_rate = 1 / n_variables
    mutation_rate_init = 2 * mutation_rate
    n_runs = 5

    print("Starting", n_runs, "runs of the algorithm with n-pop:", pop_size,
          "n-gen:", n_generations,
          "p-cross-start:", crossover_probability_start,
          "p-mut:", mutation_rate,
          "p-mut-init:", mutation_rate_init)

    start = time.time()
    app.run(main())
    end = time.time()
    print("Done after %s seconds" % (end - start))
