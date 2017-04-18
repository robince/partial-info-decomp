#!env python
import numpy as np
import scipy as sp
import dit
import scipy.io

dat = sp.io.loadmat('pyP.mat')
P = np.squeeze(dat['P'])
s = P.shape
Nvar = len(P.shape)
if Nvar != 4:
    raise ValueError, "not supported"

d = dit.Distribution(*zip(*np.ndenumerate(P)))

from dit.algorithms.scipy_optimizers import maxent_dist
# pairwise predictor-target and full 3 way predictor-predictor
me = maxent_dist(d, [ [0,3], [1,3], [2,3], [0,1,2] ])
me.make_dense()

Pme = me.pmf.reshape(s,order='C')
dat['Pme'] = Pme
sp.io.savemat('pyP.mat',dat)

