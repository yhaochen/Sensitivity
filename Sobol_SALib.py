# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
%reset -f

import numpy as np
import math
import matplotlib.pyplot as plt

from SALib.sample import latin
from SALib.sample import sobol_sequence
from SALib.sample import saltelli
from SALib.analyze import sobol


#step 1: define models (unnecessary if directly have outputs)
def Reliability(x1, x2):
    return float(math.floor(x1*x2/0.75))

#step 2: define problem input variables and boundaries
problem = {
    'num_vars': 2,
    'names': ['x2', 'x1'],
    'bounds': [[0, 1],[0, 1]]
}

# sample
N=2**17
param_values = latin.sample(problem,N*6)
param_values = sobol_sequence.sample(N,2)
param_values = saltelli.sample(problem,N)

# evaluate
y = np.array([Reliability(*params) for params in param_values])

# analyze: "problem" contains input parameters, "Y" contains output
# if (by default) second-order indices are calculated, then sample size must be N(2k+2)
# if only need to calculate first order indices, sample size must be N(k+2)
sobol_indices = sobol.analyze(problem, y, print_to_console=True,calc_second_order=False)

np.savetxt("param_values.txt", param_values)
np.savetxt("y.txt", y)
param_values = np.loadtxt("param_values.txt", float)
y = np.loadtxt("y.txt", float)
