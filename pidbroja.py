#!env python
import sys
import numpy as np
import scipy as sp
import dit
import scipy.io
from dit.pid import PID_BROJA


fname = sys.argv[1]
dat = sp.io.loadmat(fname)
P = np.squeeze(dat['P'])
s = P.shape
Nvar = len(P.shape)
if Nvar != 3:
    raise ValueError, "not supported"

d = dit.Distribution(*zip(*np.ndenumerate(P)))

x = PID_BROJA(d, [[0],[1]], [2])

pid = np.zeros(4)
pid[0] = x.get_partial(((0,), (1,)))
pid[1] = x.get_partial(((0,),))
pid[2] = x.get_partial(((1,),))
pid[3] = x.get_partial(((0,1),))

dat['pid'] = pid
sp.io.savemat(fname,dat)
