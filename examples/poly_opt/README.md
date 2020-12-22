# examples/poly_opt

* `polyOpt.m`
This code formulates and solves the *polynomial optimization* problem
using a sum-of-squares optimization approach and polynomial interpolants
as described in:

> Papp, D; Yildiz, S: Sum-of-squares optimization without semidefinite
> programming. *SIAM Journal on Optimization* 29(1), 2019, pp. 822-851. 
> URL: https://doi.org/10.1137/17M1160124


* `demo_polyOpt.m`
This script demonstrates how to use the provided methods to solve the
polynomial optimization problems described in the above paper.


* `polyEnv.m`
This code formulates and solves the *polynomial envelope* problem
using a sum-of-squares optimization approach and polynomial interpolants
as described in:

> Papp, D; Yildiz, S: Sum-of-squares optimization without semidefinite
> programming. *SIAM Journal on Optimization* 29(1), 2019, pp. 822-851. 
> URL: https://doi.org/10.1137/17M1160124


* `demo_polyEnv.m`
This script demonstrates how to use the provided methods to solve the
polynomial envelope problem also described in the above paper.


* `ChebInterval.m`
This code generates parameters for the interpolant basis representation
of univariate sum-of-squares polynomials.


* `PaduaSquare.m`
This code generates parameters for the interpolant basis representation
of bivariate sum-of-squares polynomials.


* `FeketeCube.m`
This code generates parameters for the interpolant basis representation
of sum-of-squares polynomials with three or more variables. It follows
the approach described in:

> A. Sommariva and M. Vianello, Computing approximate Fekete points by
> QR factorizations of Vandermonde matrices, *Computers & Mathematics*
> *with Applications*, 57 (2009), pp. 1324-1336. URL:
> https://doi.org/10.1016/j.camwa.2008.11.011.


EXTERNAL FUNCTIONS REQUIRED BY THESE EXAMPLES:

* `chebpts`, `chebpolyval` from Chebfun. Chebfun is an open-source package for
numerical computing with functions: http://www.chebfun.org/.
Our code has been developed and tested with Chebfun version 5.5.0.  
The latest version of the Chebfun package can be downloaded from
http://www.chebfun.org/download/.

* `pdpts`, `pdwtsMM` from Padua2DM. Padua2DM is a Matlab package from M. Caliari,
S. De Marchi, A. Sommariva, and M. Vianello for interpolation and
cubature at the Padua points. It can be downloaded from
http://profs.sci.univr.it/~caliari/software.htm.

* `partitions`. The partitions function computes all partitions of an integer.
We use the implementation of John D'Errico from a MATLAB Central File 
Exchange post which is available at
https://www.mathworks.com/matlabcentral/fileexchange/12009-partitions-of-an-integer. 