from pdetools import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

"""
Parameters

"""

D = 0.1
koff = 0.005

"""
Model functions

"""


def generate_flow_profile(A=74, B=391, C=1000, D=100):
    X = np.linspace(0, 60, 100)
    v = (((60 - X) / A) * np.exp(-((60 - X) ** 2) / B)) - ((X / C) * np.exp(-(X ** 2) / D))
    return v


class Model:
    def __init__(self, D, kon, koff, psi=0.174, tot=1.56):

        # Diffusion
        self.D = D  # diffusion coefficient

        # Flow profile
        self.flow_profile = generate_flow_profile()

        # Membrane exchange
        self.kon = kon  # membrane binding rate
        self.koff = koff  # membrane unbinding rate

        # Misc
        self.tot = tot  # total amount of protein
        self.psi = psi  # surface area to volume ratio
        self.xsteps = 100  # number of positions to split the spatial dimension into
        self.L = 60  # system length
        self.deltax = self.L / self.xsteps  # spatial step

    def diffusion(self, concs, dx):
        concs_ = np.r_[concs[0], concs, concs[-1]]  # Dirichlet boundary conditions
        d = concs_[:-2] - 2 * concs_[1:-1] + concs_[2:]
        return d / (dx ** 2)

    def flow(self, concs, dx):
        # Calculates the gradient in both directions and takes the average
        return (np.r_[0, np.diff(concs * self.flow_profile)] +
                np.r_[np.diff(concs * self.flow_profile), 0]) / (2 * dx)

    def dxdt(self, X):
        m = X[0]
        c = self.tot - self.psi * np.mean(m)  # calculate uniform cytoplasmic concentration
        flow = self.flow(m, self.deltax)
        dm = (self.kon * c) - (self.koff * m) + (self.D * self.diffusion(m, self.deltax)) + flow
        return [dm, ]

    def dxdt_no_flow(self, X):
        m = X[0]
        c = self.tot - self.psi * np.mean(m)  # calculate uniform cytoplasmic concentration
        dm = (self.kon * c) - (self.koff * m) + (self.D * self.diffusion(m, self.deltax))
        return [dm, ]

    def run(self, Tmax, t_eval=None, start=None, flow=True, killfunc=None, maxstep=None, rk=True):

        # Specify evaluation times
        if t_eval is None:
            t_eval = np.arange(0, Tmax + 0.0001, Tmax)

        # Specify starting conditions
        if start is None:
            # Start from uniform equilibrium membrane concentration, calculated analytically
            start = (self.kon * self.tot) / (self.koff + self.psi * self.kon)
            start *= np.ones([self.xsteps])

        # Specify flow regime
        if flow:
            func = self.dxdt
        else:
            func = self.dxdt_no_flow

        # Run simulation
        soln, time, solns, times = pdeRK(dxdt=func, X0=[start, ], Tmax=Tmax, deltat=0.01, t_eval=t_eval,
                                         killfunc=killfunc, maxstep=maxstep)
        # Return results
        return soln, time, solns, times


"""
Run simulation

"""

if __name__ == '__main__':
    # Build class
    m = Model(D=D, koff=koff, kon=koff)

    # Run with flow
    soln_seg, time_seg, solns_seg, times_seg = m.run(Tmax=500, flow=True, t_eval=np.arange(0, 501, 10))

    # Run without flow
    soln_rel, time_rel, solns_rel, times_rel = m.run(start=soln_seg[0], flow=False, Tmax=1000,
                                                     t_eval=np.arange(0, 1001, 10))

    # Append results
    solns = [np.r_[i, j] for i, j in zip(solns_seg, solns_rel)]
    times = np.r_[times_seg, times_rel + time_seg]
    flows = np.r_[np.tile(generate_flow_profile(), (solns_seg[0].shape[0], 1)),
                  np.tile(np.zeros(100), (solns_rel[0].shape[0], 1))]

    # Create animation
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    fig.set_size_inches(4, 3)
    fig.tight_layout()
    fig.subplots_adjust(right=0.85)


    def update(Time=0):
        ax1.clear()
        ax2.clear()
        tpoint = np.argmin(abs(times - Time))
        ax1.plot(np.linspace(0, 60, 100), solns[0][tpoint])
        ax2.plot(np.linspace(0, 60, 100), flows[tpoint], c='0.8', linestyle='--')
        ax1.set_ylim(0, 3)
        ax2.set_ylim(-0.005, 0.125)
        ax1.set_xlabel('Position (x)')
        ax1.set_ylabel('Membrane concentration')
        ax2.set_ylabel('Flow velocity (µm/s)', fontsize=8, color='0.6')
        ax1.tick_params(axis='both', labelsize=8)
        ax2.tick_params(axis='both', labelsize=8)
        ax2.spines['right'].set_color('0.6')
        ax2.tick_params(axis='y', colors='0.6')


    frames = np.r_[0, np.arange(0, max(times), 10)]
    anim = animation.FuncAnimation(fig, update, frames=iter(frames), save_count=len(frames))
    writer = animation.writers['ffmpeg']
    writer = writer(fps=24, bitrate=2000)
    anim.save('animation.gif', writer=writer, dpi=300)
