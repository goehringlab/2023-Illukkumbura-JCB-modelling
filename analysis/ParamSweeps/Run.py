import os
import sys

home_direc = os.path.dirname(os.path.realpath(__file__))
sys.path.append(home_direc)
sys.path.append(home_direc + '/..')
sys.path.append(home_direc + '/../Modules/polaritymodel')
save_direc = home_direc + '/../../../ModelData/advection_v2/'

from Model import Model
from polaritypde import ParamSpace2D
import numpy as np
import copy

print(sys.argv[1])


def calc_asi(m):
    ant = np.mean(m[:30])
    post = np.mean(m[-30:])
    asi = abs((ant - post) / (2 * (ant + post)))
    return asi


"""
Wider D range for square plots

"""

# Base parameter set
BaseM = Model(Dm=0, kon=0, koff=0, psi=0.174, tot=1.56, xsteps=100)

# Parameter range boundaries
p1_boundaries = (-3, -2, -1, 0, 1)  # log D
p2_boundaries = (-4, -3, -2, -1, 0)  # log koff and kon

# Split parameter ranges
param_range_groups = []
for i in range(len(p1_boundaries) - 1):
    for j in range(len(p2_boundaries) - 1):
        p1_range = [p1_boundaries[i], p1_boundaries[i + 1]]
        p2_range = [p2_boundaries[j], p2_boundaries[j + 1]]
        param_range_groups.append([p1_range, p2_range])
len_param_range_groups = len(param_range_groups)
# print(len_param_range_groups) = 12

"""
ASI after flow 500s

"""

if int(sys.argv[1]) in range(0, len_param_range_groups):
    prange = param_range_groups[int(sys.argv[1])]
    p1_range = prange[0]
    p2_range = prange[1]


    def func(log_D, log_koff):
        m = copy.deepcopy(BaseM)

        # Parameters
        D = 10 ** log_D
        koff = 10 ** log_koff
        m.Dm = D
        m.koff = koff
        m.kon = koff

        # Run
        soln, _, _, _ = m.run(Tmax=500)

        # ASI
        asi = calc_asi(soln[0])
        return asi

"""
ASI after flow 1000s

"""

if int(sys.argv[1]) in range(len_param_range_groups, 2 * len_param_range_groups):
    prange = param_range_groups[int(sys.argv[1]) - len_param_range_groups]
    p1_range = prange[0]
    p2_range = prange[1]


    def func(log_D, log_koff):
        m = copy.deepcopy(BaseM)

        # Parameters
        D = 10 ** log_D
        koff = 10 ** log_koff
        m.Dm = D
        m.koff = koff
        m.kon = koff

        # Run
        soln, _, _, _ = m.run(Tmax=1000)

        # ASI
        asi = calc_asi(soln[0])
        return asi

"""
Timescale to decrease ASI by half after 500s flow

"""

if int(sys.argv[1]) in range(2 * len_param_range_groups, 3 * len_param_range_groups):
    prange = param_range_groups[int(sys.argv[1]) - 2 * len_param_range_groups]
    p1_range = prange[0]
    p2_range = prange[1]


    def func(log_D, log_koff):
        m = copy.deepcopy(BaseM)

        # Parameters
        D = 10 ** log_D
        koff = 10 ** log_koff
        m.Dm = D
        m.koff = koff
        m.kon = koff
        m.Tmax = 1000

        # Run with flow
        soln, time, solns, times = m.run(Tmax=500)

        # ASI
        asi_orig = calc_asi(soln[0])

        # ASI threshold
        def killfunc(X, asi_thresh=asi_orig * 0.5):
            current_asi = calc_asi(X[0])
            if current_asi < asi_thresh:
                return True
            return False

        # Run without flow
        soln, time, _, _ = m.run(Tmax=10000, start=soln[0], flow=False, killfunc=killfunc, maxstep=0.01)
        return time

"""
Relative depletion of posterior cortex after 500s

"""

if int(sys.argv[1]) in range(3 * len_param_range_groups, 4 * len_param_range_groups):
    prange = param_range_groups[int(sys.argv[1]) - 3 * len_param_range_groups]
    p1_range = prange[0]
    p2_range = prange[1]


    def func(log_D, log_koff):
        m = copy.deepcopy(BaseM)

        # Parameters
        D = 10 ** log_D
        koff = 10 ** log_koff
        m.Dm = D
        m.koff = koff
        m.kon = koff

        # Run
        soln, _, _, _ = m.run(Tmax=500)

        # Relative depletion
        # starting_p = (m.kon * m.tot) / (m.koff + m.psi * m.kon)
        starting_p = np.mean(soln[0])
        current_p = np.mean(soln[0][-30:])
        relative_depletion = (starting_p - current_p) / starting_p
        return relative_depletion

"""
Relative depletion of posterior cortex after 1000s

"""

if int(sys.argv[1]) in range(4 * len_param_range_groups, 5 * len_param_range_groups):
    prange = param_range_groups[int(sys.argv[1]) - 4 * len_param_range_groups]
    p1_range = prange[0]
    p2_range = prange[1]


    def func(log_D, log_koff):
        m = copy.deepcopy(BaseM)

        # Parameters
        D = 10 ** log_D
        koff = 10 ** log_koff
        m.Dm = D
        m.koff = koff
        m.kon = koff

        # Run
        soln, _, _, _ = m.run(Tmax=1000)

        # Relative depletion
        # starting_p = (m.kon * m.tot) / (m.koff + m.psi * m.kon)
        starting_p = np.mean(soln[0])
        current_p = np.mean(soln[0][-30:])
        relative_depletion = (starting_p - current_p) / starting_p
        return relative_depletion

"""
Timescale to decrease relative depletion of posterior by half after 500s with flow

"""

if int(sys.argv[1]) in range(5 * len_param_range_groups, 6 * len_param_range_groups):
    prange = param_range_groups[int(sys.argv[1]) - 5 * len_param_range_groups]
    p1_range = prange[0]
    p2_range = prange[1]


    def func(log_D, log_koff):
        m = copy.deepcopy(BaseM)

        # Parameters
        D = 10 ** log_D
        koff = 10 ** log_koff
        m.Dm = D
        m.koff = koff
        m.kon = koff

        # Run
        soln, _, _, _ = m.run(Tmax=500)

        # Relative depletion
        # starting_p = (m.kon * m.tot) / (m.koff + m.psi * m.kon)
        starting_p = np.mean(soln[0])
        current_p = np.mean(soln[0][-30:])
        relative_depletion_orig = (starting_p - current_p) / starting_p

        # Relative depletion threshold
        def killfunc(X, rd_thresh=relative_depletion_orig * 0.5):
            current_p = np.mean(X[0][-30:])
            current_rd = (starting_p - current_p) / starting_p
            if current_rd < rd_thresh:
                return True
            return False

        # Run without flow
        soln, time, _, _ = m.run(Tmax=10000, start=soln[0], flow=False, killfunc=killfunc, maxstep=0.01)
        return time

###############################################################################################


ParamSpace2D(func, p1_range=p1_range, p2_range=p2_range, cores=32, resolution0=11, direc=save_direc + sys.argv[1],
             parallel=True, replace=True).run()
