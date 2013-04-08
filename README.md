Embedded Conic Solver (ECOS)
====

ECOS is a numerical software for solving convex second-order cone programs (SOCPs) of type
```
min  c'*x 
s.t. A*x = b
     G*x <=_K h
```
where the last inequality is generalized, i.e. `h - G*x` belongs to the cone `K`. 
ECOS supports the positive orthant `R_+` and second-order cones `Q_n` defined as
```
Q_n = { (t,x) | t >= || x ||_2 } 
```
In the definition above, t is a scalar and `x` is in `R_{n-1}`. The cone `K` is therefore
a direct product of the positive orthant and second-order cones:
```
K = R_+ x Q_n1 x ... x Q_nN
```

Features of ECOS
----
 
+ *ECOS runs on embedded platforms*. Written in ANSI C, 
  it can be compiled for any platform for which a C compiler is available. Excluding the problem setup
  part, no memory manager is needed for solving problem instances of same structure. 
+ *ECOS is efficient*. Using sparse linear algebra routines, it computes only what is really necessary.
  The interior point algorithm it implements is one of the fastest converging methods that are currently
  in use for solving convex conic problems.
+ *ECOS has a tiny footprint*. The ECOS solver consists of 750 lines of C code (excluding the problem setup 
   code).
+ *ECOS is numerically robust*. Using regularization and iterative refinement coupled with a carefully chosen 
  sparse representation of scaling matrices, millions of problem instances are solved reliably.
+ *ECOS comes with a MATLAB and CVX 2.0 interface*. This way, you can prototype, simulate and verify the performance
  of ECOS before you implement the very same code on your embedded hardware.
+ *ECOS is library-free*. No need to link any external library to ECOS, apart from `AMD` and `sparseLDL`, both 
  from Timothy A. Davis, which are included in this project.


Using ECOS in MATLAB
====


Compiling ECOS for MATLAB
----

ECOS comes with a makefile which resides in the `matlab` subdirectory of the code. To build ECOS for MATLAB:
```matlab
cd <ecos-directory>/matlab
makemex
```
You should now have a binary file `ecos.[ending]`, with a platform-specific ending. This is the solver binary. 
Add the directory `<ecos-directory>/matlab` to your path to be able to call ECOS from any place. The command
```matlab 
makemex clean
```
deletes unecessary files that were produced during compilation.

Calling ECOS from MATLAB
----

You can directly call ECOS from Matlab using its native interface, that basically takes the problem data `c,G,h,A,b`
and some dimension information that is given in the struct `dims`. This structure has the following fields:
```
dims.l - scalar, dimension of positive orthant (LP-cone) R_+
dims.q - vector with dimensions of second order cones
```
The length of `dims.q` determines the number of second order cones. If you do not have a cone in your problem, use
the empty matrix `[ ]` instead.

### Example
