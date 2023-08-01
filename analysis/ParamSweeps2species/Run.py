import os
import sys

home_direc = os.path.dirname(os.path.realpath(__file__))
sys.path.append(home_direc)
sys.path.append(home_direc + '/..')
save_direc = home_direc + '/../../../ModelData/advection_two_species/'

from Model2species import Model
from paramspace2d import ParamSpace2D
import numpy as np
import copy

print(sys.argv[1])


def calc_asi(m):
    ant = np.mean(m[:30])
    post = np.mean(m[-30:])
    asi = abs((ant - post) / (2 * (ant + post)))
    return asi


"""

"""

# Base parameter set
BaseM = Model(Dm=0.1, kon1=1, kon2=1, koff1=1, koff2=1, psi=0.174, tot=1.56, xsteps=100)

# Parameter range boundaries
p1_boundaries = (0.001, 0.25, 0.5, 0.75, 0.999)  # fraction in fast component
p2_boundaries = (-4, -3, -2, -1, 0)  # log koff_slow

# Split parameter ranges
param_range_groups = []
for i in range(len(p1_boundaries) - 1):
    for j in range(len(p2_boundaries) - 1):
        p1_range = [p1_boundaries[i], p1_boundaries[i + 1]]
        p2_range = [p2_boundaries[j], p2_boundaries[j + 1]]
        param_range_groups.append([p1_range, p2_range])
len_param_range_groups = len(param_range_groups)
# print(len_param_range_groups)


"""
ASI after flow 500s

"""

if int(sys.argv[1]) in range(0, len_param_range_groups):
    prange = param_range_groups[int(sys.argv[1])]
    p1_range = prange[0]
    p2_range = prange[1]


    def func(fraction, log_koffS):
        m = copy.deepcopy(BaseM)

        # Parameters
        m.koff1 = 10 ** log_koffS  # slow component off rate
        m.kon1 = 1 - fraction  # slow component on rate
        m.kon2 = fraction  # fast component on rate

        # Run
        soln, _, _, _ = m.run(Tmax=500)

        # ASI
        asi = calc_asi(soln[0] + soln[1])
        return asi

"""
Ratio at membrane (whole embryo)

"""

if int(sys.argv[1]) in range(len_param_range_groups, 2 * len_param_range_groups):
    prange = param_range_groups[int(sys.argv[1]) - len_param_range_groups]
    p1_range = prange[0]
    p2_range = prange[1]


    def func(fraction, log_koffS):
        m = copy.deepcopy(BaseM)

        # Parameters
        m.koff1 = 10 ** log_koffS  # slow component off rate
        m.kon1 = 1 - fraction  # slow component on rate
        m.kon2 = fraction  # fast component on rate

        # Run
        soln, _, _, _ = m.run(Tmax=500)

        # Ratio at membrane
        mem_fast_fraction = np.mean(soln[1]) / (np.mean(soln[0] + soln[1]))
        return mem_fast_fraction

"""
Ratio at membrane (anterior)

"""

if int(sys.argv[1]) in range(2 * len_param_range_groups, 3 * len_param_range_groups):
    prange = param_range_groups[int(sys.argv[1]) - 2 * len_param_range_groups]
    p1_range = prange[0]
    p2_range = prange[1]


    def func(fraction, log_koffS):
        m = copy.deepcopy(BaseM)

        # Parameters
        m.koff1 = 10 ** log_koffS  # slow component off rate
        m.kon1 = 1 - fraction  # slow component on rate
        m.kon2 = fraction  # fast component on rate

        # Run
        soln, _, _, _ = m.run(Tmax=500)

        # Ratio at membrane
        mem_fast_fraction = np.mean(soln[1][:30]) / (np.mean(soln[0][:30] + soln[1][:30]))
        return mem_fast_fraction

###############################################################################################


ParamSpace2D(func, p1_range=p1_range, p2_range=p2_range, cores=32, resolution0=11, direc=save_direc + sys.argv[1],
             parallel=True, replace=True).run()
