import numpy as np
import os
import multiprocessing
import shutil


class ParamSpace2D:
    def __init__(self, func, p1_range, p2_range, resolution, direc, parallel=False, cores=None, args=[], replace=False):
        """

        Class to perform parameter space sweeps

        Run by calling run function
        - performs analysis, saves results to a .csv file, and final array as a .txt
        - progress is saved, if interrupted can be resumed without loss by calling the run function again

        Computation parameters
        :param func: Function to evaluate. Must take two parameter values as inputs, and return a single number
        representing a measure of model behaviour. This can be a float, in the case of a quantitative measure
        (e.g. asymmetry index, aPAR domain size), or an integer in the case of a quantitative measure
        (e.g. 0/1 = polarised/not polarised).
        :param p1_range: range for parameter 1 (lower, upper)
        :param p2_range: range for parameter 2 (lower, upper)
        :param resolution: grid resolution
        :param parallel: if True, will run in parallel using number of cores specified
        :param cores: number of cores on machine to use in parallel
        :param args: optional additional arguments for func
        :param replace: If True, any results already stored in direc will be deleted and replaced. If False, any results
        already stored in direc will be imported and not replaced.

        # Saving parameters
        :param direc: directory to save results. Directory must already exist

        """

        # Computation
        self.func = func
        self.p1_range = p1_range
        self.p2_range = p2_range
        self.resolution = resolution
        self.parallel = parallel
        self.args = args

        if cores is None:
            self.cores = multiprocessing.cpu_count()
        else:
            self.cores = cores

        # Saving
        self.direc = direc
        self.replace = replace

        # Results
        self.res = None
        self.n_sims = None

    def single_eval(self, p1val_p2val):
        """
        Single function call for given p1 and p2 values, save result

        """

        # Run function
        state = self.func(*[float(i) for i in p1val_p2val.split(',')] + self.args)

        # Save state
        with open(self.direc + '/' + 'Sims.csv', 'a') as f:
            f.write(p1val_p2val + ',' + str(state) + '\n')

    def batch_eval(self, pcombs):
        """
        Evaluate parameter sets in bulk
        pcombs is list of strings 'p1val,p2val'

        """
        if self.parallel:
            pool = multiprocessing.Pool(self.cores)
            pool.map(self.single_eval, pcombs)
        else:
            for k in iter(pcombs):
                self.single_eval(k)

    def import_res(self):
        """
        Import all results from current iteration, load into self.res

        """

        with open(self.direc + '/' + 'Sims.csv') as g:
            for line in g:
                p1, p2, val = line[:-1].split(',')
                xind = ((float(p1) - self.p1_range[0]) * (self.n_sims - 1)) / (self.p1_range[1] - self.p1_range[0])
                yind = ((float(p2) - self.p2_range[0]) * (self.n_sims - 1)) / (self.p2_range[1] - self.p2_range[0])
                if '.' in val:
                    self.res[round(xind), round(yind)] = float(val)
                else:
                    self.res[round(xind), round(yind)] = int(val)

    def run(self):
        """
        Run algorithm, save figure

        """

        # Clear results directory if replace is True
        if self.replace:
            if os.path.isdir(self.direc):
                shutil.rmtree(self.direc)

        # Make results directory, if it doesn't already exist
        if not os.path.isdir(self.direc):
            os.mkdir(self.direc)

        # Grid
        self.n_sims = self.resolution
        run_bool = np.ones([self.n_sims, self.n_sims])

        # Parameter combinations
        sims_array_ind = np.nonzero(run_bool)
        p1vals = self.p1_range[0] + sims_array_ind[0] * (self.p1_range[1] - self.p1_range[0]) / (self.n_sims - 1)
        p2vals = self.p2_range[0] + sims_array_ind[1] * (self.p2_range[1] - self.p2_range[0]) / (self.n_sims - 1)
        pcombs = ["{:.12f}".format(p1vals[i]) + ',' + "{:.12f}".format(p2vals[i]) for i in range(len(p1vals))]

        # Remove parameters already tested (if algorithm run before)
        if os.path.isfile(self.direc + '/' + 'Sims.csv'):
            with open(self.direc + '/' + 'Sims.csv') as f:
                for line in f:
                    p = line.split(',')[0] + ',' + line.split(',')[1]
                    if p in pcombs:
                        pcombs.remove(p)

        # Run
        self.batch_eval(pcombs)

        # Import results
        self.res = np.nan * np.zeros([self.n_sims, self.n_sims])
        self.import_res()

        # Specify int/float format
        f = (self.res % 1) == 0

        # Save row by row
        with open(self.direc + '/Res.txt', 'w') as fh:
            for i, row in enumerate(self.res):
                formats = f[i]
                line = ' '.join(
                    "{:.0f}".format(value) if formats[j] else "{:.12f}".format(value) for j, value in enumerate(row))
                fh.write(line + '\n')
