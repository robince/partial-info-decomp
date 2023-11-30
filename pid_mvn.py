import numpy as np

class Lattice2D:
    def __init__(self):
        # Build 2D redundancy lattice
        self.Nx = 2
        self.Nnodes = 4
        self.Nelements = 3
        self.elements = [1, 2, [1, 2]]
        
        self.A = [[1, 2], [1], [2], [[1, 2]]]
        self.level = [1, 2, 2, 3]
        self.Nlevels = max(self.level)
        self.labels = ['{1}{2}', '{1}', '{2}', '{12}']

        self.nodes = {}
        for ni, label in enumerate(self.labels):
            self.nodes[label] = ni + 1
        
        parents = [[]] * self.Nnodes
        parents[self.nodes['{1}{2}'] - 1] = [self.nodes['{1}'], self.nodes['{2}']]
        parents[self.nodes['{1}'] - 1] = [self.nodes['{12}']]
        parents[self.nodes['{2}'] - 1] = [self.nodes['{12}']]
        self.parents = parents

        # Build children for bidirectional link
        children = [[]] * self.Nnodes
        for ni in range(1, self.Nnodes + 1):
            for parent in parents[ni - 1]:
                if children[parent - 1] == []:
                    children[parent - 1] = []
                children[parent - 1].append(ni)
        self.children = children

        self.top = self.nodes['{12}']
        self.bottom = self.nodes['{1}{2}']
        self.Icap = [float('nan')] * self.Nnodes
        self.PI = [float('nan')] * self.Nnodes

def logmvnpdf(x, mu, Sigma):
    N, D = x.shape
    const = -0.5 * D * np.log(2 * np.pi)

    if len(x) == len(mu):
        xc = np.array([x[i] - mu[i] for i in range(N)])
    else:
        xc = x - mu

    term1 = -0.5 * np.sum((xc @ np.linalg.pinv(np.atleast_2d(Sigma))) * xc, axis=1)
    term2 = const - 0.5 * logdet(Sigma)

    logp = term1 + term2
    return logp


    # Copyright (c) 2011, Benjamin Dichter
    # All rights reserved.
    #
    # Redistribution and use in source and binary forms, with or without
    # modification, are permitted provided that the following conditions are
    # met:
    #
    # * Redistributions of source code must retain the above copyright
        # notice, this list of conditions and the following disclaimer.
    # * Redistributions in binary form must reproduce the above copyright
        # notice, this list of conditions and the following disclaimer in
        # the documentation and/or other materials provided with the distribution
    #
    # THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    # AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    # IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    # ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
    # LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
    # CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
    # SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
    # INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
    # CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    # ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    # POSSIBILITY OF SUCH DAMAGE.

def logdet(A):
    U = np.linalg.cholesky(A)
    y = 2 * np.sum(np.log(np.diag(U)))
    return y

def Iccs_mvn(A, Cfull, varsizes):
    # Calculate redundancy between a set of Gaussian sources
    # from pointwise common change in surprise
    # using conditional independent joint-source distributions
    
    if np.sum(varsizes) != Cfull.shape[0]:
        raise ValueError('Wrong number of variables specified')
    if len(varsizes) != 3:
        raise ValueError('Only 2 variables supported')

    NA = len(A)
    NVs = varsizes[-1]
    Nx = len(varsizes) - 1
    NVx = varsizes[:-1]
    varstart = np.cumsum(varsizes)
    varstart = np.insert(varstart[:-1], 0, 0)
    
    sidx = np.arange(varstart[-1], varstart[-1] + NVs)
    Cs = Cfull[sidx, sidx]

    AC = []
    for ai in range(NA):
        thsA = A[ai]
        if isinstance(thsA, int):
            thsA = [thsA]
        aidxfull = [None] * len(thsA)
        aidx = [None] * len(thsA)
        thsvstart = 0
        for vi in range(len(thsA)):
            aidxfull[vi] = np.arange(varstart[thsA[vi]-1] , varstart[thsA[vi] - 1] + NVx[thsA[vi]-1])
            thsL = len(aidxfull[vi])
            aidx[vi] = np.arange(thsvstart, thsvstart + thsL)
            thsvstart = thsvstart + thsL
        thsNv = sum(len(a) for a in aidx)
        Ca = np.zeros((thsNv, thsNv))
        
        # fill in blocks
        # diagonal
        for vi in range(len(thsA)):
            Ca[aidx[vi], aidx[vi]] = Cfull[aidxfull[vi], aidxfull[vi]]
        # off diagonal
        for vi in range(len(thsA)):
            for vj in range(len(thsA)):
                if vi == vj:
                    continue
                Ca[aidx[vi], aidx[vj]] = Cfull[aidxfull[vi], aidxfull[vj]]

        Cas = np.zeros((thsNv + NVs, thsNv + NVs))
        Cas[:thsNv, :thsNv] = Ca
        # joint with S
        # diagonal
        thssidx = np.arange(thsNv, thsNv + NVs)
        Cas[thssidx, thssidx] = Cs
        # off diagonal
        for vi in range(len(thsA)):
            Cas[aidx[vi], thssidx] = Cfull[aidxfull[vi], sidx] # 0,2
            Cas[thssidx, aidx[vi]] = Cfull[sidx, aidxfull[vi]] # 2,0

        Casoff = Cas[:thsNv, thssidx]
        CXYY1 = Casoff @ np.linalg.pinv(np.atleast_2d(Cs))
        Cacs = Ca - CXYY1 @ Cas[thssidx, :thsNv]
        MacsF = CXYY1

        AC.append({
            'Ca': Ca,
            'Cas': Cas,
            'Cacs': Cacs,
            'Casoff': Casoff,
            'MacsF': CXYY1,
            'Nv': thsNv
        })

    if NA == 1:
        # use closed form expression
        chA = np.linalg.cholesky(np.atleast_2d(AC[0]['Ca']))
        chS = np.linalg.cholesky(np.atleast_2d(Cs))
        chAS = np.linalg.cholesky(np.atleast_2d(AC[0]['Cas']))
        # normalizations cancel for information
        HA = np.sum(np.log(np.diag(chA)))
        HS = np.sum(np.log(np.diag(chS)))
        HAS = np.sum(np.log(np.diag(chAS)))
        Iccs = (HA + HS - HAS) / np.log(2)

    if NA == 2:
        # Covariance for Pind(A1,A2)
        thsNv = AC[0]['Nv'] + AC[1]['Nv'] + NVs
        ANv = AC[0]['Nv'] + AC[1]['Nv']
        a1idx = slice(0,AC[0]['Nv'])
        a2idx = slice(AC[0]['Nv'], AC[0]['Nv'] + AC[1]['Nv'])
        a12idx = slice(0, AC[0]['Nv'] + AC[1]['Nv'])

        # P2 == full gaussian covariance
        C = Cfull[a12idx,a12idx]

        Caasoff = Cfull[a12idx, sidx]
        CXYY1 = Caasoff @ np.linalg.pinv(np.atleast_2d(Cs))
        Caacs = C - CXYY1 @ Cfull[sidx, a12idx]
        MaacsF = CXYY1

        thssidx = slice(AC[0]['Nv']+AC[1]['Nv'], AC[0]['Nv']+AC[1]['Nv']+NVs)

        # Monte Carlo Integration for Iccs
        # 100,000 works well for 5d space.
        # 10,000 might be ok but some variance
        Nmc = 100000

        mcx = np.random.multivariate_normal(mean=np.zeros(ANv + NVs), cov=Cfull, size=Nmc)
        
        px = logmvnpdf(mcx[:, a1idx], np.zeros(AC[0]['Nv']), AC[0]['Ca'])
        pxcs = logmvnpdf(mcx[:, a1idx], (AC[0]['MacsF'] @ mcx[:, thssidx].T).T.squeeze(), AC[0]['Cacs'])
        
        py = logmvnpdf(mcx[:, a2idx], np.zeros(AC[1]['Nv']), AC[1]['Ca'])
        pycs = logmvnpdf(mcx[:, a2idx], (AC[1]['MacsF'] @ mcx[:, thssidx].T).T.squeeze(), AC[1]['Cacs'])
        
        pxy = logmvnpdf(mcx[:, a12idx], np.zeros(ANv), C)
        pxycs = logmvnpdf(mcx[:,a12idx], (MaacsF@mcx[:,thssidx].T).T, Caacs)
        
        dhx = pxcs - px
        dhy = pycs - py
        dhxy = pxycs - pxy

        lnii = dhx + dhy - dhxy
        keep = np.logical_and(np.sign(dhx) == np.sign(lnii),
                              np.logical_and(np.sign(dhx) == np.sign(dhy),
                                             np.sign(dhx) == np.sign(dhxy)))
        lnii[~keep] = 0
        
        Iccs = np.nanmean(lnii)
        # convert to bits
        Iccs = Iccs / np.log(2)

    return Iccs


def calc_pi_mvn(lat, Cfull, varsizes, Icap, forcenn=False, normlevels=False):
    # Calculate PI on a redundancy lattice using Williams and Beer summation

    # if only lat provided calculate PI using existing Icap
    # otherwise, recalculate Icap
    if len(Cfull) != sum(varsizes) or len(varsizes) != 3:
        raise ValueError('Wrong number of variables specified')

    for ni in range(0, lat.Nnodes):
        lat.Icap[ni] = Icap(lat.A[ni], Cfull, varsizes)

    if lat.Nx > 3:
        raise ValueError('calc_pi: too many variables')

    # use equation (7) from Williams and Beer to calculate PI at each node
    lat.PI = np.full_like(lat.Icap, np.nan)
    # raw PI before non-disjoint normalization
    lat.PIraw = np.full_like(lat.Icap, np.nan)

    # ascend through levels of the lattice
    Nlevels = np.max(lat.level)
    for li in range(1, Nlevels):
        nodes = np.where(np.array(lat.level) == li)[0]
        for ni in nodes:
            lat = calc_pi_node(lat, ni, forcenn, normlevels)

    # don't enforce non-negativity for top node
    lat = calc_pi_node(lat, lat.top - 1, False, normlevels)

    return lat

def calc_pi_node(lat, ni, nonneg=False, normlevels=False):

    children = lat.children[ni]

    if not children:
        # no children
        thsPI = lat.Icap[ni]

        if nonneg:
            thsPI = max(thsPI, 0)

        lat.PI[ni] = thsPI
        lat.PIraw[ni] = thsPI
        return lat

    all_children = recurse_children(lat, ni, [])

    if normlevels:
        PIchildren = normalise_levels(lat, all_children)
    else:
        PIchildren = lat.PI[np.array(all_children)-1]

    thsPI = lat.Icap[ni] - np.sum(PIchildren)

    if nonneg:
        thsPI = max(thsPI, 0)

    lat.PI[ni] = thsPI
    lat.PIraw[ni] = thsPI

    if ni == lat.top - 1:
        lat.PI[np.array(all_children)-1] = PIchildren

    return lat


def normalise_levels(lat, children):
    # normalize to correct for non-additivity of non-disjoint nodes

    # values for this set of children
    PIraw = lat.PIraw[children]
    levels = lat.level[children]
    labels = lat.labels[children]
    A = lat.A[children]
    normPI = PIraw.copy()

    for li in range(1, lat.Nlevels + 1):
        nodes = np.where(levels == li)[0]
        levelPI = PIraw[nodes]
        posPInodes = nodes[np.abs(levelPI) > 1e-15]
        posPIvars = np.concatenate(A[posPInodes])
        posPIvars = np.concatenate([np.array(var) for var in posPIvars])

        if len(posPIvars) != len(np.unique(posPIvars)):
            # have non-disjoint positive PI contributions at this level

            # using structure of 3rd order lattice (might need more logic to
            # determine pairwise disjoint-ness for higher order lattices)
            if li == 4:
                # special case level 4 for 3 variable lattice
                # one node contains all variables
                fullnode = np.where(labels == '{12}{13}{23}')[0]
                if not fullnode or PIraw[fullnode] < 1e-15:
                    # all sources at this level are disjoint so no
                    # normalization required
                    continue
                elif len(posPInodes) == 1:
                    # only {12}{13}{23} is non-zero so no normalization
                    # required
                    continue

                # only normalize by 2 here even if more posPInodes because
                # there are only 2 disjoint copies at this level
                normPI[posPInodes] = PIraw[posPInodes] / 2
            else:
                normPI[posPInodes] = PIraw[posPInodes] / len(posPInodes)

    return normPI

def recurse_children(lat, ni, children):
    children.extend(lat.children[ni])

    for ci in lat.children[ni]:
        children = recurse_children(lat, ci-1, children)

    return list(set(children))
