# local_interaction

This folder implements a redundancy measure based on pointwise common change of surprisal.
It also implements a partial information decomposition based  n this measure.

For comparison the original `Imin` of [Williams and Beer (2010)](http://arxiv.org/abs/1004.2515) is also implemented.

### Usage

The lattice is implemented as a Matlab structure:

`lat = lattice2d();` or `lat = lattice3d()` for 2 or 3 variable decomposition respectively [(Williams and Beer, 2010; Figure 2)](http://arxiv.org/abs/1004.2515).

The joint distribution for the system of interest should be structured in an array with an axis for each variable with the privileged (independent) variable as the last axis.

Then the redundancy for any source can be calculated by specifying the source as a cell array of elements, where each element is an array of the constituent variables:

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
These implement different summation strategies over the lattice. 

A range of examples are implemented in the scripts [`examples_2d.m`](examples_2d.m) and [`examples_3d.m`](examples_3d.m). 
The output of these scripts is in [`examples2d_output.txt`](examples2d_output.txt) and [`examples3d_output.txt`](examples3d_output.txt) respectively.


### Main functions

- [`Imin.m`](Imin.m) : Imin function from [Williams and Beer (2010)](http://arxiv.org/abs/1004.2515)
- [`calc_pi_wb.m`](calc_pi_wb.m) : PID decomposition from [Williams and Beer (2010)](http://arxiv.org/abs/1004.2515)

- [`Iccs.m`](Iccs.m) : Redundancy measure based on pointwise common change in surprisal.
- [`calc_pi.m`](calc_pi.m) : Alternative PID decomposition on the lattice which normalises non-zero PI values across non-disjoint nodes of the same height in the lattice, to avoid overcounting.


### Retired experimental functions

- `Iccs_fulljoint.m` : pointwise common change in surprisal with true joint element distributions.
- `Inlii_indjoint.m` : redundancy measure based on negative local/pointwise interaction information with independent joint element distributions.
- `Inlii_fulljoint.m` : negative local interaction information with true joint element distributions.
- `calc_pi_backprop_normalisation.m`, `calc_pi_normalise_levels.m` : alternative approaches to sum over the lattice to obtain a consistent non-negative decomposition.

