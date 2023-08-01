import numpy as np
import os
import sys

home_direc = os.path.dirname(os.path.realpath(__file__))
sys.path.append(home_direc)
from pde_rk import pde_rk


class Model:
    def __init__(self, Dm, kon, koff, deltat=0.01, xsteps=100, psi=0.174, tot=1.56, noise=0):

        # Diffusion
        self.Dm = Dm

        # Flow profile
        def generate_flow_profile(A=74, B=391, C=1000, D=100):
            X = np.linspace(0, 60, xsteps)
            part1 = ((X / A) * np.exp(-(X ** 2) / B))
            part2 = (((X - 60) / C) * np.exp(-((X - 60) ** 2) / D))
            Y = part1 + part2
            Y -= Y[0] + ((Y[-1] - Y[0]) / 60) * X  # adding this term to force the ends to zero
            return Y[::-1]

        self.flow_profile = generate_flow_profile()
        self.noise = noise

        # Membrane exchange
        self.kon = kon
        self.koff = koff

        # Misc
        self.tot = tot
        self.psi = psi
        self.xsteps = int(xsteps)
        self.deltat = deltat
        self.deltax = 60 / self.xsteps

    def diffusion(self, concs, dx):
        concs_ = np.r_[concs[0], concs, concs[-1]]
        d = concs_[:-2] - 2 * concs_[1:-1] + concs_[2:]
        return d / (dx ** 2)

    def flow(self, concs, dx):
        # Calculate gradient in both directions, take average
        return (np.r_[0, np.diff(concs * self.flow_profile)] +
                np.r_[np.diff(concs * self.flow_profile), 0]) / (2 * dx)

    def flow_noisy(self, concs, dx):
        noise = np.random.normal(0, self.noise, size=self.xsteps)
        return (np.r_[0, np.diff(concs * (self.flow_profile + noise))] +
                np.r_[np.diff(concs * (self.flow_profile + noise)), 0]) / (2 * dx)

    def dxdt(self, X):
        m = X[0]
        c = self.tot - self.psi * np.mean(m)

        if self.noise == 0:
            flow = self.flow(m, self.deltax)
        else:
            flow = self.flow_noisy(m, self.deltax)
            flow[0], flow[-1] = [0, 0]

        dm = (self.kon * c) - (self.koff * m) + (self.Dm * self.diffusion(m, self.deltax)) + flow
        return [dm, ]

    def dxdt_no_flow(self, X):
        m = X[0]
        c = self.tot - self.psi * np.mean(m)

        dm = (self.kon * c) - (self.koff * m) + (self.Dm * self.diffusion(m, self.deltax))
        return [dm, ]

    def run(self, Tmax, t_eval=None, start=None, flow=True, killfunc=None, maxstep=None, rk=True):

        # Evaluation times
        if t_eval is None:
            t_eval = np.arange(0, Tmax + 0.0001, Tmax)

        # Starting conditions
        if start is None:
            start = (self.kon * self.tot) / (self.koff + self.psi * self.kon)

        # Flow regime
        if flow:
            func = self.dxdt
        else:
            func = self.dxdt_no_flow

        soln, time, solns, times = pde_rk(dxdt=func, X0=[np.ones([self.xsteps]) * start, ],
                                         Tmax=Tmax, deltat=self.deltat, t_eval=t_eval, killfunc=killfunc,
                                         maxstep=maxstep, rk=rk)
        return soln, time, solns, times


def calc_asi(m):
    ant = np.mean(m[:50])
    post = np.mean(m[50:])
    asi = abs((ant - post) / (2 * (ant + post)))
    return asi
