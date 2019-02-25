#!env python
import sys
import numpy as np
import scipy as sp
import dit
import scipy.io
from dit.algorithms.distribution_optimizers import BROJABivariateOptimizer

def pid_broja_dist(dist, sources, target, rv_mode=None):
    broja = BROJABivariateOptimizer(dist, sources, target, rv_mode)
    broja.optimize(niter=10)
    opt_dist = broja.construct_dist()
    return opt_dist

fname = sys.argv[1]
dat = sp.io.loadmat(fname)
P = np.squeeze(dat['P'])
s = P.shape
Nvar = len(P.shape)
if Nvar != 3:
    raise ValueError, "not supported"

d = dit.Distribution(*zip(*np.ndenumerate(P)))

bd = pid_broja_dist(d, [[0],[1]], [2])
bd.make_dense()

Pb = bd.pmf.reshape(s,order='C')
dat['Pb'] = Pb
sp.io.savemat(fname,dat)
