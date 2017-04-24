#!env python
import sys
import numpy as np
import scipy as sp
import dit
import scipy.io
from dit.algorithms.scipy_optimizers import pid_broja_dist


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
