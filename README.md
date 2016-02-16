# local_interaction

This folder implements a redundancy measure based on pointwise common change of surprisal and a partial information decomposition based on this measure.

For comparison the original `Imin` of [Williams and Beer (2010)](http://arxiv.org/abs/1004.2515) is also implemented.

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

where `Icap` is one of the redundancy functions.

The PID can be calulated by passing a lattice, full joint distribution, and function handle to a redundancy function to one of the `calc_pi` functions. 
These implement different summation strategies over the lattice:

```matlab
lat = lattice2d();
latmin = calc_pi_wb(lat, Pxxy, @Imin);
latccs = calc_pi(lat, Pxxy, @Iccs);
lat.labels % labels for each element in the PID
latmin.PI % Williams and Beer PID
latccs.PI % CCS based PID
```

A range of examples are implemented in the scripts [`examples_2d.m`](examples_2d.m) and [`examples_3d.m`](examples_3d.m). 
The output of these scripts is in [`examples2d_output.txt`](examples2d_output.txt) and [`examples3d_output.txt`](examples3d_output.txt) respectively.


### Main functions

- [`Imin.m`](Imin.m) : Imin function from [Williams and Beer (2010)](http://arxiv.org/abs/1004.2515)
- [`calc_pi_wb.m`](calc_pi_wb.m) : PID decomposition from [Williams and Beer (2010)](http://arxiv.org/abs/1004.2515)

- [`Iccs.m`](Iccs.m) : Redundancy measure based on pointwise common change in surprisal.
- [`calc_pi.m`](calc_pi.m) : Alternative PID decomposition on the lattice which normalises non-zero PI values across non-disjoint nodes of the same height in the lattice, to avoid overcounting.

### Guassian / Multivariate Normal (mvn) functions

The functions below implement the PID for multivariate normal / Gaussian variables (currently only 2 predictor variables). The input format here is the full covariance matrix of the system, following by a vector `varsizes`, which lists the multivariate structure of each considered variable. I.e. if `size(Cfull) = [6 6]` and `varsizes = [1 3 2]`, this means variable `X1` is a univariate Gaussian in the first position in `Cfull`, `X2` is a trivariate Gaussian and `S` is a bivariate Gaussian.

- [`Immi_mvn.m`](Immi_mvn.m) : Redundancy measure based on minumum mutual information (equivalent to Imin for Gaussian sources).
- [`Iccs_mvn.m`](Iccs_mvn.m) : Redundancy measure based on pointwise common change in surprisal.
- [`calc_pi_mvn.m`](calc_pi_mvn.m) : Alternative PID decomposition on the lattice which normalises non-zero PI values across non-disjoint nodes of the same height in the lattice, to avoid overcounting.
- [`examples_2dmvn.m`](examples_2dmvn.m) : Run Gaussian examples.


### Retired experimental functions

- `Iccs_fulljoint.m` : pointwise common change in surprisal with true joint element distributions.
- `Inlii_indjoint.m` : redundancy measure based on negative local/pointwise interaction information with independent joint element distributions.
- `Inlii_fulljoint.m` : negative local interaction information with true joint element distributions.
- `calc_pi_backprop_normalisation.m`, `calc_pi_normalise_levels.m` : alternative approaches to sum over the lattice to obtain a consistent non-negative decomposition.

