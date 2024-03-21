import numpy as np
import pandas as pd
import os
from ioh import problem, OptimizationType, get_problem, logger

# EXPRMAT_AD = pd.read_csv("./data/expr_ad.csv", delimiter=';')
# EXPRMAT_AD_geneMap = np.array(EXPRMAT_AD.columns)
# EXPRMAT_AD = np.array(EXPRMAT_AD)

PHENO_HD = np.loadtxt("./data/pheno_hd.csv", delimiter=';', usecols=1, skiprows=1)
EXPRMAT_HD = pd.read_csv("./data/expr_hd.csv", delimiter=';')
EXPRMAT_HD_geneMap = np.array(EXPRMAT_HD.columns)
EXPRMAT_HD = np.array(EXPRMAT_HD)

# module parameters
min_size = 30
min_gene_overlap = .5
min_func_overlap = .75

def is_valid(x):
    """
    method that checks a module for module requirements:
    1. the module should be formatted as a binary vector with one element for every gene in the co-expression matrix
    2. the minimal size of the module is 30
    3. the minimal 1-1 gene overlap with the reference module is 0.5
    4. the minimal functional overlap with the reference module is 0.75
    :param x: a module, i.e. a group of genes, in the HD network
    :return: "True" if the module meets all requirements, "False" otherwise.
    """
    if len(np.unique(x)) > 2:
        return False
    if sum(x) < min_size:
        return False
    # TODO: add method body for requirements 3 and 4


def module_eigengene(x):
    """
    computes the module eigengene for a group of genes in the HD network.
    :param x: a module, i.e. a group of genes, in the HD network
    :return: the first principal component of the expression matrix of the corresponding module
    """
    mod_expr = EXPRMAT_HD[:, x.astype(bool).squeeze()]
    scaled_mod_expr = ((mod_expr - np.mean(mod_expr, axis=0))/np.std(mod_expr, axis=0)).T

    u, d, v = np.linalg.svd(scaled_mod_expr)
    v = v[0:1, :]
    pc = v[0].tolist()
    return pc


def fitness(x):
    """
    method that computes the "quality" of a given module in the HD network
    :param x: a module, i.e. a group of genes, in the HD network
    :return: the fitness score of the module
    """
    if not is_valid(x):
        return 0
    me = module_eigengene(x)
    result = np.corrcoef(me, PHENO_HD)[0, 1]
    # TODO: add more statistics into the fitness formula
    return result


problem.wrap_integer_problem(fitness, 'module_fitness', optimization_type=OptimizationType.MAX)

# Call get_problem to instantiate a version of this problem
f = get_problem('module_fitness', instance=0, dimension=26, problem_type="Integer")
logger = logger.Analyzer(
    root=os.getcwd(),
    # Store data in the current working directory
    folder_name="results",
    algorithm_name="GA-transfer-learning",
    # meta-data for the algorithm used to generate these results.
    algorithm_info="type-current-test-settings-here",
    store_positions=True  # store x-variables in the logged files
)
f.attach_logger(logger)
