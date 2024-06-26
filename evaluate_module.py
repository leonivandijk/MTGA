import numpy as np
import pandas as pd
import os
from ioh import problem, OptimizationType, get_problem, logger, ProblemClass

# global variables
search_space = None
expr_mat = None
tom = None
max_tom = None
pheno = None
AD_MODULE = None
start_module = None
f = None

path = "/data/s3035158/data/"

# algorithm parameters
min_size = 30
min_gene_overlap = .5
# min_func_overlap = .75
deg_threshold = .05
member_threshold = .5

# initial solution
AD_MODULE = np.array(pd.read_table(path + "saddlebrown.txt", dtype=str))


def computeGeneModuleMembership():
    module_expression = np.array(expr_mat.loc[:, AD_MODULE[:, 0]])
    me = pd.DataFrame(module_eigengene(x=None, mod_expr=module_expression), index=expr_mat.index)
    expr_mat['me'] = me
    return expr_mat.corr()['me'].iloc[:-1]


def create_searchspace(degs):
    degs = degs[degs["pvalue"] < deg_threshold].index
    gene_module_membership = computeGeneModuleMembership()
    members = gene_module_membership[abs(gene_module_membership) > member_threshold].index
    searchspace = degs.union(members).to_numpy(dtype=str)
    return searchspace

def module_eigengene(x=None, mod_expr=None):
    """
    computes the module eigengene for a group of genes in the HD network.
    :param mod_expr: if mod_expr is already computed elsewhere
    :param x: a module, i.e. a group of genes, in the HD network
    :return: the first principal component of the expression matrix of the corresponding module
    """
    if x is not None:
        mod_expr = expr_mat[:, x.astype(bool).squeeze()]
    else:
        mod_expr = mod_expr
    scaled_mod_expr = ((mod_expr - np.mean(mod_expr, axis=0)) / np.std(mod_expr, axis=0)).T

    u, d, v = np.linalg.svd(scaled_mod_expr)
    v = v[0:1, :]
    pc = v[0].tolist()
    return pc

def load_data(disease):
    print("loading data of", disease, "network")
    global search_space, expr_mat, tom, max_tom, pheno, AD_MODULE, start_module, f

    if disease == "AD":
        # ad data
        expr_mat = pd.read_csv(path + "EXPRMAT_AD.csv", delimiter=';')
        tom = pd.read_csv(path + "TOM_AD.csv", delimiter=';')
        pheno = np.loadtxt(path + "PHENO_AD.csv", delimiter=';', usecols=1, skiprows=1, dtype=int)
        # extra results
        degs = pd.read_csv(path + "results_deseq2_15385.csv", delimiter='\t', index_col=0)

    elif disease == "HD":
        # hd data
        expr_mat = pd.read_csv(path + "EXPRMAT_HD.csv", delimiter=';')
        tom = pd.read_csv(path + "TOM_HD.csv", delimiter=';')
        pheno = np.loadtxt(path + "PHENO_HD.csv", delimiter='\t', usecols=1, skiprows=1, dtype=int)
        # extra results
        degs = pd.read_csv(path + "result_deseq2_15984.csv", delimiter='\t', index_col=0)

    else:
        raise ValueError("Unsupported disease code: " + disease)

    search_space = create_searchspace(degs=degs)
    # make data-matrices smaller for efficient memory use
    expr_mat = np.array(expr_mat.loc[:, search_space])
    tom = np.array(tom.loc[search_space, search_space])
    max_tom = np.max(np.triu(tom, k=1))

    start_module = np.zeros(shape=search_space.shape, dtype=int)
    for i in AD_MODULE:
        index = np.where(search_space == i)
        start_module[index] = 1
    print("loading done")

    # Call get_problem to instantiate a version of this problem
    f = get_problem('module_fitness',
                    problem_class=ProblemClass.INTEGER,
                    instance=0,
                    dimension=len(search_space))
    f.attach_logger(logger)


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
    me = module_eigengene(x=x, mod_expr=None)
    signal = np.corrcoef(me, pheno)[0, 1]
    signal = abs(signal)

    # quality metric 2
    # might favour smaller modules due to sparseness of our graph.
    # if so, think of an extra tradeoff metric that favours larger modules
    module_size = sum(x)
    edges = np.sum(np.triu(tom[x.astype(bool).squeeze()][:, x.astype(bool).squeeze()], k=1)) #exclude diagonal
    # scale max number of possible edges by upper bound on strength of a connection
    possible_edges = (module_size * (module_size - 1) / 2) * max_tom

    #print("correlation:", signal)
    #print("connectivity:", (edges / possible_edges))

    return signal * (edges / possible_edges)


problem.wrap_integer_problem(fitness,
                             name='module_fitness',
                             optimization_type=OptimizationType.MAX,
                             lb=0)


logger = logger.Analyzer(
    root=os.getcwd(),
    # Store data in the current working directory
    folder_name="results",
    algorithm_name="GA-HD-saddlebrown",
    # meta-data for the algorithm used to generate these results.
    algorithm_info="HD-tests-saddlebrown",
    store_positions=True  # store x-variables in the logged files
)
