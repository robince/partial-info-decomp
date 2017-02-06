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

#outcomes = [];
#probs = [];
#for xi in range(P.shape[0]):
    #for yi in range(P.shape[1]):
        #for zi in range(P.shape[2]):
            #outcomes.append('{:X}{:X}{:X}'.format(xi,yi,zi))
            #probs.append(P[xi,yi,zi])

#d = dit.Distribution(outcomes, probs)
d = dit.Distribution(*zip(*np.ndenumerate(P)))
#d.set_rv_names('XYZ')
me = dit.algorithms.marginal_maxent_dists(d,k_max=2,maxiters=2000,tol=1e-4,verbose=True);
me1 = me[1]
me2 = me[2]
me1.make_dense()
me2.make_dense()

P1 = me1.pmf.reshape(s,order='C')
P2 = me2.pmf.reshape(s,order='C')
dat['P1'] = P1
dat['P2'] = P2
sp.io.savemat('pyP.mat',dat)

