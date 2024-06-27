import sys
import time
from absl import app
import numpy as np
import evaluate_module
import random
from numba import cuda, float32, int32
import argparse
import os
# Set this at the beginning of your script
os.environ["CUDA_VISIBLE_DEVICES"] = "0"

# initialise
np.random.seed(42)

pop_size = 250
n_generations = 400
crossover_probability = 0.6

# we use 3 cuda kernels according to paper by M. Abbasi.
## kernel 1: constructs initial population and evaluates its fitness. 1 individual per thread.
@cuda.jit
def init_population_kernel(func, parents, parents_f, start_module, mutation_rate_init, n_variables, tom, guide=True):
    idx = cuda.grid(1) # 1dimensional grid as our populations are arrays
    if idx < n_variables-1:
        individual = start_module.copy()

    # Perform mutation or guided mutation
    if not guide:
        mutation(individual, mutation_rate_init, n_variables)
    else:
        guided_mutation(individual, mutation_rate_init, n_variables, tom)

    if evaluate_module.is_valid(individual):
        # Store the valid individual and its fitness score
        parents[idx] = individual
        parents_f[idx] = func(individual)


def init_population(func, start_module, mutation_rate_init, n_variables, guide=False):
    """
    Constructs the initial population for the Genetic Algorithm. The population consists of binary vectors where
    every element corresponds to a gene. An element gets the value ‘1’ if the gene is part of the given module and
    ‘0’ if not. We want the initial population to explore the neighbourhood of the AD module in the HD brain. We can
    do this by randomly adding and removing genes in the module.
    :return:
    """

    # Initialize the population and fitness scores arrays
    population = np.zeros((pop_size, n_variables), dtype=np.int32)
    fitness_scores = np.zeros(pop_size, dtype=np.float32)

    # Allocate device memory and copy data from host to device
    d_population = cuda.to_device(population)
    d_fitness_scores = cuda.to_device(fitness_scores)
    d_start_module = cuda.to_device(np.array(start_module, dtype=np.int32))
    d_tom = cuda.to_device(np.array(evaluate_module.tom, dtype=np.int32))

    # Define the number of threads per block and the number of blocks per grid
    threads_per_block = 512
    blocks_per_grid = (pop_size + threads_per_block - 1) // threads_per_block

    # Launch the kernel
    init_population_kernel[blocks_per_grid, threads_per_block](func, d_population, d_fitness_scores, d_start_module,
                                                               mutation_rate_init, n_variables, d_tom, guide)

    # Copy the results back to the host
    d_population.copy_to_host(population)
    d_fitness_scores.copy_to_host(fitness_scores)

    return population, fitness_scores

# kernel 2: forms new children by applying multiple operators and evaluates their fitness
@cuda.jit
def fitness_kernel( ):
    return []
# kernel 3: selects the new population



def mutation(p, rate, n_variables):
    """
    Random mutation. (De-)Selects every gene with an equal probability equal to the mutation rate.
    :param p:
    :param rate:
    :return:
    """
    for i in range(n_variables):
        if np.random.uniform(0, 1) < rate:
            # bit flipping
            p[i] = 1 - p[i]


def guided_mutation(p, rate, n_variables, tom):
    # randomly select a number of nodes according to the rate

    n_nodes = int(np.round(rate * n_variables))
    for i in range(n_nodes):
        module_genes = np.flatnonzero(p)
        module_size = sum(p)
        # select a node in the module
        selected_node = int(np.round(float(np.random.uniform(0, 1) * (module_size - 1))))
        selected_gene = int(module_genes[selected_node])

        # mutate in the neighborhood of the selected gene: interchange low connected node with highly connected node
        # first select gene to be added; high connection higher selection probability
        row = tom[selected_gene,:].copy()
        row[selected_gene] = 0
        prob_add = row/np.sum(row)
        add = np.random.choice(n_variables, p=prob_add)

        # then select the node to remove
        # low connection higher selection probability
        if np.random.uniform(0, 1) > 0.3:
            to_module = tom[p.astype(bool).squeeze()][:, p.astype(bool).squeeze()]
            row = to_module[selected_node,:].copy()
            row[selected_node] = 1
            row = 1 - row
            prob_remove = row/np.sum(row)
            remove = np.random.choice(module_size, p=prob_remove)
            remove_gene = int(module_genes[remove])
            p[remove_gene] = 0

        p[add] = 1

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

def tournament_selection(parent, parent_f):
    select_parent = []
    for i in range(len(parent)) :
        pre_select = np.random.choice(len(parent_f), 5, replace = False)
        max_f = parent_f[pre_select[0]]
        index = pre_select[0]
        for p in pre_select:
            if parent_f[p] > max_f:
                index = p
                max_f = parent_f[p]
        select_parent.append(parent[index].copy())
    return select_parent


def onepoint_crossover(p1, p2, n_variables):
    if np.random.uniform(0, 1) < crossover_probability:
        idx = np.random.randint(n_variables)
        t = p1[idx:]
        p1[idx:] = p2[idx:]
        p2[idx:] = t


def npoint_crossover(n, p1, p2, n_variables):
    if (np.random.uniform(0, 1) < crossover_probability):
        idxs = sorted(random.sample(range(n_variables), n))
        for idx in idxs:
            t = p1[idx:]
            p1[idx:] = p2[idx:]
            p2[idx:] = t

def uniform_crossover(p1, p2, n_variables):
    if (np.random.uniform(0, 1) < crossover_probability):
        for i in range(n_variables):
            if np.random.uniform(0, 1) < 0.5:
                t = p1[i]
                p1[i] = p2[i]
                p2[i] = t


def uniform_crossover(p1, p2, n_variables):
    if (np.random.uniform(0, 1) < crossover_probability):
        for i in range(n_variables):
            if np.random.uniform(0, 1) < 0.5:
                t = p1[i]
                p1[i] = p2[i]
                p2[i] = t



def genetic_algorithm(func, start_module, generations_left=None):
    # parameters settings
    n_variables = len(start_module)
    mutation_rate = 2 / n_variables
    mutation_rate_init = mutation_rate

    if generations_left is None:
        generations_left = n_generations

    f_opt = sys.float_info.min
    x_opt = None

    # construct the initial population
    parents, parents_f = init_population(func, start_module, mutation_rate_init, n_variables, guide=True)
    cuda.synchronize()

    while generations_left > 0:

        offspring = tournament_selection(parents, parents_f)

        for i in range(0, pop_size - (pop_size % 2), 2):
            npoint_crossover(10, offspring[i], offspring[i + 1], n_variables)


        for i in range(pop_size):
            #mutation(offspring[i], rate=mutation_rate, n_variables=n_variables)
            guided_mutation(offspring[i], rate=mutation_rate, n_variables=n_variables, tom=evaluate_module.tom)

        parents = offspring.copy()
        for i in range(pop_size):
            parents_f[i] = func(parents[i])
            if parents_f[i] > f_opt:
                f_opt = parents_f[i]
                x_opt = parents[i].copy()
        generations_left = generations_left - 1

    print(f_opt, x_opt)

    return f_opt, x_opt


def main(disease):
    start_module = evaluate_module.start_module
    f = evaluate_module.f

    f_opt_result = 0
    # We run the algorithm 100 independent times.
    n_runs = 25
    for _ in range(n_runs):
        print("start run", _)
        start_subrun = time.time()
        f_opt, _ = genetic_algorithm(func=f, start_module=start_module)
        f_opt_result += f_opt
        f.reset()
        end_subrun = time.time()
        print((end_subrun - start_subrun), "seconds")
    print("The average f_opt is %s" % (f_opt_result / n_runs))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--disease',
                        help="Choose between AD or HD")
    args, unknown = parser.parse_known_args()
    evaluate_module.load_data(args.disease)

    start = time.time()
    app.run(main(args.disease))
    end = time.time()
    print("The program takes %s seconds" % (end - start))
