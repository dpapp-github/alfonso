# alfonso: ALgorithm FOr Non-Symmetric Optimization

## A fast, open-source Matlab solver for conic optimization

`alfonso` is an open-source Matlab package for solving convex optimization problems in conic form, created by Dávid Papp and Sercan Yıldız. It enables optimization over any convex cone as long as a suitable barrier is available for either the cone or its dual. This includes many nonsymmetric cones, for example, hyperbolicity cones and their duals (such as sum-of-squares cones), semidefinite and second-order cone representable cones, power cones, and the exponential cone. It also offers performance advantages for problems whose symmetric cone programming representation requires a large number of auxiliary variables or has a special structure that can be exploited in the barrier computation. An example of this is *sums-of-squares optimization*.

The interfaces of the original version of the software are described in the accompanying software paper:
> D. Papp and S. Yıldız. alfonso: Matlab package for nonsymmetric conic optimization. arXiv:2101.04274 [https://arxiv.org/abs/2101.04274](https://arxiv.org/abs/2101.04274)

The "oracle interface" of the current version is slightly different; see any of the built-in barriers (e.g., gH_LP.m) for details. (*TBD: updating the preprint above.*)

The package also includes an implementation of the sum-of-squares optimization algorithm based on non-symmetric conic optimization and polynomial interpolants presented in:

> D. Papp and S. Yıldız. Sum-of-squares optimization without semidefinite programming. *SIAM Journal on Optimization* 29(1), 2019, pp. 822-851. [https://doi.org/10.1137/17M1160124](https://doi.org/10.1137/17M1160124)

The code is distributed under the [2-Clause BSD License](LICENSE).

## Citing alfonso

To cite alfonso, please mention the [research article](https://doi.org/10.1137/17M1160124) for which the code was originally developed and the [software paper](https://doi.org/10.1287/ijoc.2021.1058). Here's a [BibTeX file](alfonso.bib) with these references.

## What's new?
July 2024: The `alfonso_simple` interface has a built-in barrier function for *very* efficient optimization over sets defined by the types of linear matrix inequalities that are characteristic of *sums-of-squares* optimization and *design of experiments*.

July 2024: The custom barrier functions for `alfonso` are implemented using a slightly different interface: rather than returning the Hessian matrix H and its Cholesky factor L, they return a function handles representing the actions of the inverse Hessian and inverse Cholesky factor. See one of the built-in barrier functions, e.g., `gH_LP.m` for details.

Big thanks to [Giovanni Fantuzzi](https://dcn.nat.fau.eu/giovanni-fantuzzi/) for his helpful discussions and the occasional beta (alpha?) testing prior to this update. 

## Installation

`alfonso` is entirely written in Matlab m-code. To install, unzip the downloaded files in any directory and add the `src` subdirectory to the Matlab path.

The polynomial optimization (sum-of-squares) examples in `examples/poly_opt/` require the following additional software:

1. `Padua2DM` (a Matlab/Octave package by M. Caliari, S. De Marchi, A. Sommariva, and M. Vianello for interpolation and
cubature at the Padua points): [download](http://profs.sci.univr.it/~caliari/software.htm) to any directory and add its directory to the Matlab path.

2. `Chebfun` (an open-source package for numerical computing with functions): [download](http://www.chebfun.org/download/) and install following the installation instructions.
Our code has been developed and tested with Chebfun version 5.5.0.

3. `partitions` (a function by John D'Errico to compute all partitions of an integer): [download](https://www.mathworks.com/matlabcentral/fileexchange/12009-partitions-of-an-integer) to any directory and add its directory to the Matlab path.

## Contents

* `src`: source files for alfonso. Use `alfonso.m` for the oracle interface, and `alfonso_simple` for the simple interface, which lets you specify complicated cones from simple, built-in building blocks.
* `examples`: contains a number of optimization problems solved with `alfonso` or `alfonso_simple`:
  * `examples/random_lp`: solves a random linear programming problem using either interface.
  * `examples/exp_design`: solves an optimal design of experiments problem.
  * `examples/poly_opt`: solves various polynomial optimization problems using sum-of-squares optimization.
  * `examples/portfolio`: solves a mean-risk portfolio optimization problem with a factor risk model and market impact.

## Additional information

This package is based on a infeasible-start primal-dual interior-point algorithm that originally appeared in:

> A. Skajaa and Y. Ye, A homogeneous interior-point algorithm for nonsymmetric convex conic optimization, *Mathematical Programming Ser. A* 150 (2015), pp. 391-422. [https://doi.org/10.1007/s10107-014-0773-1](https://doi.org/10.1007/s10107-014-0773-1)

The alfonso implementation is derived from the corrected analysis of this algorithm published in 2017:

> D. Papp and S. Yıldız. On "A homogeneous interior-point algorithm for nonsymmetric convex conic optimization". [https://arxiv.org/abs/1712.00492](https://arxiv.org/abs/1712.00492)

The Matlab software was first published here:
> D. Papp and S. Yıldız. alfonso: Matlab package for nonsymmetric conic optimization,
*INFORMS Journal on Computing* 34(1) (2022), pp. 11-19. [https://doi.org/10.1287/ijoc.2021.1058](https://doi.org/10.1287/ijoc.2021.1058)

The first version of alfonso was created at North Carolina State University (NCSU) and at the The Statistical and Applied Mathematical Sciences Institute (SAMSI), with partial support from the National Science Foundation. Recent (2023-) improvements of alfonso are supported in part by the Air Force Office of Scientific Research. Grant numbers:
* NSF:
  * DMS-1638521
  * DMS-1719828
  * DMS-1847865
* AFOSR:
  * FA9550-23-1-0370
