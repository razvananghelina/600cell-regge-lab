"""
EXP-153: Log-Periodic Corrections to Gravity from 600-Cell Greens Function
============================================================================
CLASSIFICATION: To be determined by results.
"""

import numpy as np
from itertools import permutations, product as iprod
import scipy.linalg as la
from scipy.optimize import curve_fit

PHI = (1 + np.sqrt(5)) / 2
print('test')
