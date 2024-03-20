import sys
import time
from absl import app
import numpy as np
import pandas as pd
import evaluate_module

# load in gene/wgcna options
TOM_HD = pd.read_csv("./data/tom_hd.csv", delimiter=';')

# initialise
np.random.seed(42)
f = evaluate_module.f

# parameters settings
pop_size = 1000
n_generations = 1000
tournament_k = 2
mutation_rate = 0.0001
crossover_probability = 0.6


def genetic_algorithm(func, generations = None):
    if generations is None:
        generations = n_generations

    f_opt = sys.float_info.min
    x_opt = None

    # construct the initial population
    parent = []
    parent_f = []




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