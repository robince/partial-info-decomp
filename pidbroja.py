#!env python
import numpy as np
import scipy as sp
import dit
import scipy.io
from dit.algorithms.scipy_optimizers import pid_broja

dat = sp.io.loadmat('pyP.mat')
P = np.squeeze(dat['P'])
s = P.shape
Nvar = len(P.shape)
if Nvar != 3:
    raise ValueError, "not supported"

d = dit.Distribution(*zip(*np.ndenumerate(P)))

x = pid_broja(d, [[0],[1]], [2])

pid = np.zeros(4)
pid[0] = x.R
pid[1] = x.U0
pid[2] = x.U1
pid[3] = x.S

dat['pid'] = pid
sp.io.savemat('pyP.mat',dat)
