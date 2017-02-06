#!env python
import numpy as np
import scipy as sp
import dit
import scipy.io

dat = sp.io.loadmat('pyP.mat')
P = np.squeeze(dat['P'])
s = P.shape
Nvar = len(P.shape)
if Nvar != 3:
    raise ValueError, "not supported"

d = dit.Distribution(*zip(*np.ndenumerate(P)))

ui, rdn, syn, mi_orig, mi_opt = dit.algorithms.pid_broja.unique_informations(d,[[0], [1]], [2], rv_mode='indices', tol=1e-4, verbose=True)

pid = np.zeros(4)
pid[0] = rdn
pid[1] = ui[0]
pid[2] = ui[1]
pid[3] = syn

dat['pid'] = pid
sp.io.savemat('pyP.mat',dat)
