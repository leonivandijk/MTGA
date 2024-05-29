import numpy as np
import pandas as pd
import gc
import os
from ioh import problem, OptimizationType, get_problem, logger, ProblemClass

path = "/Users/leonivandijk/Desktop/thesis/pyfiles"

# hd data
# PHENO_HD   = np.loadtxt(path + "/data/PHENO_HD.csv", delimiter=';', usecols=1, skiprows=1)
# EXPRMAT_HD = pd.read_csv(path + "/data/expr_hd_turquoise.csv", delimiter=';')
# EXPRMAT_HD = np.array(EXPRMAT_HD)

# ad data
EXPRMAT_AD = pd.read_csv(path + "/data/EXPRMAT_AD.csv", delimiter=';')
AD_DEGS = pd.read_csv(path + "/data/results_deseq2_15385.csv", delimiter='\t', index_col=0)
TOM_AD = pd.read_csv(path + "/data/TOM_AD.csv", delimiter=';')
PHENO_AD = np.loadtxt(path + "/data/PHENO_AD.csv", delimiter=';', usecols=1, skiprows=1, dtype=int)

AD_MODULE = np.array(pd.read_table(path + "/data/saddlebrown.txt", dtype=str))
AD_GENEMEMBERSHIP = pd.read_csv(path + "/data/geneModuleMembership.csv", delimiter=';', usecols=range(3), index_col=0)

# algorithm parameters
min_size = 30
min_gene_overlap = .5
# min_func_overlap = .75
deg_threshold = .05
member_threshold = .5


def create_searchspace(geneMemberShip, degs):
    degs = degs[degs["padj"] < deg_threshold].index
    members = geneMemberShip[abs(geneMemberShip["Correlation"]) > member_threshold].index
    searchspace = degs.union(members).to_numpy(dtype=str)
    return searchspace


# problem specific settings
search_space = create_searchspace(geneMemberShip=AD_GENEMEMBERSHIP, degs=AD_DEGS)
disease_data = [EXPRMAT_AD, TOM_AD, PHENO_AD]

# make data-matrices smaller for efficient memory use
expr_mat = np.array(disease_data[0].loc[:, search_space])
tom = np.array(disease_data[1].loc[search_space, search_space])
pheno = disease_data[2]

del disease_data
gc.collect()

start_module = np.zeros(shape=search_space.shape, dtype=int)
for i in AD_MODULE:
    index = np.where(search_space == i)
    start_module[index] = 1


def is_valid(x):
    """
    method that checks a module for module requirements:
    1. the module should be formatted as a binary vector with one element for every gene in the co-expression matrix
    2. the minimal size of the module is 30
    3. the minimal 1-1 gene overlap ratio with respect to the reference module is 0.5
    4. the minimal functional overlap with the reference module is 0.75
    :param x: a module, i.e. a group of genes, in the HD network
    :return: "True" if the module meets all requirements, "False" otherwise.
    """
    if len(np.unique(x)) > 2:
        return False
    if sum(x) < min_size:
        return False
    if 1 - (sum(start_module) - sum(x[start_module.astype(bool)])) / sum(start_module) < min_gene_overlap:
        return False
    # TODO is_valid(x): add method body for requirement 4 (dev-stage 2)
    return True


def module_eigengene(x):
    """
    computes the module eigengene for a group of genes in the HD network.
    :param x: a module, i.e. a group of genes, in the HD network
    :return: the first principal component of the expression matrix of the corresponding module
    """
    mod_expr = expr_mat[:, x.astype(bool).squeeze()]
    scaled_mod_expr = ((mod_expr - np.mean(mod_expr, axis=0)) / np.std(mod_expr, axis=0)).T

    u, d, v = np.linalg.svd(scaled_mod_expr)
    v = v[0:1, :]
    pc = v[0].tolist()
    return pc


def fitness(x):
    """
    method that computes the "quality" of a given module in the HD network
    module quality is based on:
    1. the (differential) expression of the genes in the module with respect to the phenotype
    2. the density of the interactions between the genes in the module
    :param x: a module, i.e. a group of genes, in the network
    :return: the fitness score of the module
    """
    if not is_valid(x):
        return 0

    # quality metric 1
    # might be less strong than metric 2, as this metric will probably never reach >.8
    me = module_eigengene(x)
    signal = np.corrcoef(me, pheno)[0, 1]
    signal = abs(signal)

    # quality metric 2
    # might favour smaller modules due to sparseness of our graph.
    # if so, think of an extra tradeoff metric that favours larger modules
    edges = np.sum(np.triu(tom[start_module.astype(bool).squeeze()][:, start_module.astype(bool).squeeze()], k=1)) #exclude diagonal
    # scale max number of possible edges by upper bound on strength of a connection
    max_tom = np.max(np.triu(tom, k=1))
    possible_edges = (sum(x) * (sum(x) - 1) / 2) * max_tom

    return signal + (edges / possible_edges)


problem.wrap_integer_problem(fitness,
                             name='module_fitness',
                             optimization_type=OptimizationType.MAX,
                             lb=0)

# Call get_problem to instantiate a version of this problem
f = get_problem('module_fitness',
                problem_class=ProblemClass.INTEGER,
                instance=0,
                dimension=len(search_space))

logger = logger.Analyzer(
    root=os.getcwd(),
    # Store data in the current working directory
    folder_name="results",
    algorithm_name="GA-transfer-learning",
    # meta-data for the algorithm used to generate these results.
    algorithm_info="AD-tests-saddlebrown",
    store_positions=True  # store x-variables in the logged files
)
f.attach_logger(logger)
