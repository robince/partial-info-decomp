# local_interaction

This folder implements some pointwise approaches to decomposing multivariate entropy and mutual information into redundant, unique and synergistic terms. 

These functions make use of the Python [discrete information theory (dit) toolbox](https://github.com/dit/dit). 
To use these you will need Python on your path with `dit` installed. 

### Partial Entropy Decomposition based on common surprisal (new)

Implementation to accompany the manuscript:

RAA Ince, **The Partial Entropy Decomposition: Decomposing multivariate entropy and mutual information via pointwise common surprisal**  
[arXiv:1702.01591](http://arxiv.org/abs/1702.01591) [cs.IT]

Functions: `calc_pe`, `Hcs`, `compare_ped`
`pid_from_ped`, `monopid_from_ped`, `redpid_from_ped`

Functions with a `mech_` prefix calculate the respective information decomposition from PED but separating source and mechansitic redundancy.

Examples from paper: `examples_ped.m`

### Partial Information Decomposition based on common change in surprisal

Implementation to accompany the manuscript:

RAA Ince, **Measuring multivariate redundant information with pointwise common change in surprisal**  
[arXiv:1602.05063](http://arxiv.org/abs/1602.05063) [cs.IT]

Functions: `calc_pi`, `Iccs`, `Iccs_Pind`, `Iccs_fulljoint`,  `compare`, `calc_pi_mvn`, `Iccs_mvn`

Examples from paper: `examples_2d.m`, `examples_3d.m`, `examples_2dmvn.m`, `discrete_pred_pred.m`

Please note v2 of the arXiv paper made some changes to the method which are now reflected in the code. 
Different marginal constraints are used for the maximum entropy solution, and zero thresholding and normalisation are no longer applied to the lattice.
Thresholding and normalisation are still available as options to the `calc_pi` function, but do not occur by default.

### Other measures


For comparison purposes the original `Imin` of [Williams and Beer (2010)](http://arxiv.org/abs/1004.2515) is also implemented.

The measure Ibroja [Bertschinger et al. (2014)](http://www.mdpi.com/1099-4300/16/4/2161) is calculated via the [discrete information theory (dit) toolbox](https://github.com/dit/dit). 
To use these you will need Python on your path with `dit` installed.
This is also used for calculating maximum entropy distributions over discrete spaces subject to pairwise marginal equality constraints (used in `Hcs.m`).

Functions: `Imin`, `pid_broja`, `pid_broja_dist`, `marg_maxent2`, `marg_maxent_3pred`

Common functions: `lattice2d`, `lattice3d`


### Usage

The lattice is implemented as a Matlab structure:

`lat = lattice2d();` or `lat = lattice3d()` for 2 or 3 variable decomposition respectively [(Williams and Beer, 2010; Figure 2)](http://arxiv.org/abs/1004.2515).

The joint distribution for the system of interest should be structured in an array with an axis for each variable with the privileged (target) variable as the last axis.

Then the redundancy for any element can be calculated by providing a cell array of sources, where each source is an array specifying the constituent variables:

```matlab
% 3 way redundancy (bottom node) for 3 variable case
R = Icap({1 2 3}, Pxxxy);
% full system information (top node) for 3 variable case
R = Icap({[1 2 3]}, Pxxxy);
% other nodes
R = Icap({1 [2 3]}, Pxxxy);
```

where `Icap` is one of the redundancy functions (information or entropy).

The PID can be calulated by passing a lattice, full joint distribution, and function handle to a redundancy function to one of the `calc_pi` functions. 
These implement different summation strategies over the lattice:

```matlab
lat = lattice2d();
latmin = calc_pi(lat, Pxxy, @Imin);
latccs = calc_pi(lat, Pxxy, @Iccs);
lat.labels % labels for each element in the PID
latmin.PI % Williams and Beer PID
latccs.PI % CCS based PID
```

A range of examples for the PED with Hcs are implemented in the script [`examples_ped.m`](examples_ped.m).
A range of examples for the PID with Iccs are implemented in the scripts [`examples_2d.m`](examples_2d.m) and [`examples_3d.m`](examples_3d.m). 
The output of these scripts is in [`examples2d_output.txt`](examples2d_output.txt) and [`examples3d_output.txt`](examples3d_output.txt) respectively.

Similarly, for the Partial Entropy Decomposition of the same system, and the derived PIDs from the PED:

```matlab
lat = lattice3d();
lat = calc_pe(lat,Pxxy,@Hcs);
lat.labels % each node in the PED
lat.PI % Partial Entropy Values
pifull = pid_from_ped(lat); % full PID from PED
pimono = monopid_from_ped(lat); % monosemous PID from PED
pimonomech = mech_monopid_from_ped(lat); % monosemous PID from PED with separate mechanistic and source redundancy
pibroja = pid_broja(Pxxy);
```

### Main functions

- [`calc_pi.m`](calc_pi_wb.m) : PID decomposition from [Williams and Beer (2010)](http://arxiv.org/abs/1004.2515)
- [`Iccs.m`](Iccs.m) : Information redundancy measure based on pointwise common change in surprisal with game theoretic operational maxent constraints.
- [`Imin.m`](Imin.m) : Imin function from [Williams and Beer (2010)](http://arxiv.org/abs/1004.2515)

- [`Hcs.m`](Hcs.m) : Entropy redundancy measure based on pointwise common surprisal
- [`calc_pe.m`](calc_pe.m) : Partial Entropy Decomposition (no normalisation or thresholding of values on the lattice
- [`pid_from_ped.m`](pid_from_ped.m) : Full PID from a PED
- [`monopid_from_ped.m`](monopid_from_ped.m) : Monosemous PID from a PED, includes only non-ambiguous partial entropy terms



### Guassian / Multivariate Normal (mvn) functions

The functions below implement the PID for multivariate normal / Gaussian variables (currently only 2 predictor variables). The input format here is the full covariance matrix of the system, following by a vector `varsizes`, which lists the multivariate structure of each considered variable. I.e. if `size(Cfull) = [6 6]` and `varsizes = [1 3 2]`, this means variable `X1` is a univariate Gaussian in the first position in `Cfull`, `X2` is a trivariate Gaussian and `S` is a bivariate Gaussian.

- [`Immi_mvn.m`](Immi_mvn.m) : Redundancy measure based on minumum mutual information (equivalent to Imin for Gaussian sources).
- [`Iccs_mvn.m`](Iccs_mvn.m) : Redundancy measure based on pointwise common change in surprisal.
- [`calc_pi_mvn.m`](calc_pi_mvn.m) : Alternative PID decomposition on the lattice which normalises non-zero PI values across non-disjoint nodes of the same height in the lattice, to avoid overcounting.
- [`examples_2dmvn.m`](examples_2dmvn.m) : Run Gaussian examples.
- [`pid_mvn.py`](pid_mvn.m) : Python package to implement PID using the Iccs redundancy measure. See the [PyPID_mvn](https://github.com/AlessandroCorsini/PyPID_mvn) repo for more info on usage.


### Retired experimental functions

- `Iccs_fulljoint.m` : pointwise common change in surprisal with true joint element distributions.
- `Inlii_indjoint.m` : redundancy measure based on negative local/pointwise interaction information with independent joint element distributions.
- `Hcs_fulljoint.m` : common surprisal entropy redundancy with true joint element distribution
- `Inlii_fulljoint.m` : negative local interaction information with true joint element distributions.
- `calc_pi_backprop_normalisation.m`, `calc_pi_normalise_levels.m` : alternative approaches to sum over the lattice to obtain a consistent non-negative decomposition.

