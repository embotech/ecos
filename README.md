Embedded Conic Solver (ECOS)
====

[![Build Status](https://travis-ci.org/embotech/ecos.svg?branch=master)](https://travis-ci.org/embotech/ecos)

ECOS is a numerical software for solving convex second-order cone programs (SOCPs) of type
```
min  c'*x
s.t. A*x = b
     G*x <=_K h
```
where the last inequality is generalized, i.e. `h - G*x` belongs to the cone `K`.
ECOS supports the positive orthant `R_+` and second-order cones `Q_n` defined as
```
Q_n = { (t,x) | t >= || x ||_2 } 
```
In the definition above, t is a scalar and `x` is in `R_{n-1}`. The cone `K` is therefore
a direct product of the positive orthant and second-order cones:
```
K = R_+ x Q_n1 x ... x Q_nN
```

License
----

ECOS is distributed under the [GNU General Public License v3.0](http://www.gnu.org/copyleft/gpl.html). Other licenses may be available upon request from [embotech](http://www.embotech.com).


Features of ECOS
----

+ *ECOS runs on embedded platforms*. Written in ANSI C (except for the timing code),
  it can be compiled for any platform for which a C compiler is available. Excluding the problem setup
  part, no memory manager is needed for solving problem instances of same structure.
+ *ECOS is efficient*. Using sparse linear algebra routines, it computes only what is really necessary.
  The interior point algorithm it implements is one of the fastest converging methods that are currently
  in use for solving convex conic problems.
+ *ECOS has a tiny footprint*. The ECOS solver consists of 750 lines of C code (excluding the problem setup
   code).
+ *ECOS is numerically robust*. Using regularization and iterative refinement coupled with a carefully chosen
  sparse representation of scaling matrices, millions of problem instances are solved reliably.
+ *ECOS comes with a MATLAB and CVX 2.1 interface*, and *it is supported by YALMIP*. With [CVX](http://cvxr.com) and
  [YALMIP](http://users.isy.liu.se/johanl/yalmip/) you can prototype, simulate and verify the performance of ECOS before you implement the very same code on
your embedded hardware.  
+ *ECOS comes with a Python interface*. This interface is built on top of
  [NUMPY](http://numpy.org) and [SCIPY](http://scipy.org/) and uses its sparse data structures.
+ *There is a Julia interface for ECOS*. [Julia](http://julialang.org) is a high-level, high-performance language for technical and scientific computing. You can pull the [Julia interface here](https://github.com/JuliaOpt/ECOS.jl).
+ *ECOS is library-free*. No need to link any external library to ECOS, apart from `AMD` and `sparseLDL`, both
  from Timothy A. Davis, which are included in this project.


Credits
----

The solver is essentially based on Lieven Vandenberghe's [CVXOPT](http://cvxopt.org) [ConeLP](http://www.ee.ucla.edu/~vandenbe/publications/coneprog.pdf) solver, although it differs in the particular way the linear systems are treated.

The following people have been, and are, involved in the development and maintenance of ECOS:

+ Alexander Domahidi (principal developer)
+ Eric Chu (Python interface, unit tests)
+ Stephen Boyd (methods and maths)
+ Michael Grant (CVX interface)
+ Johan Löfberg (YALMIP interface)
+ João Felipe Santos, Iain Dunning (Julia inteface)
+ Han Wang (ECOS branch and bound wrapper)

The main technical idea behind ECOS is described in a short [paper](http://www.stanford.edu/~boyd/papers/ecos.html). More details are given in Alexander Domahidi's [PhD Thesis](http://e-collection.library.ethz.ch/view/eth:7611?q=domahidi) in Chapter 9.

If you find ECOS useful, you can cite it using the following BibTex entry:

```
@INPROCEEDINGS{bib:Domahidi2013ecos,
author={Domahidi, A. and Chu, E. and Boyd, S.},
booktitle={European Control Conference (ECC)},
title={{ECOS}: {A}n {SOCP} solver for embedded systems},
year={2013},
pages={3071-3076}
}
```

Using ECOS with CVXPY
====

[CVXPY](http://cvxpy.org) is a powerful Python modeling framework for convex optimization, similar to the MATLAB counterpart CVX. ECOS is one of the default solvers in CVXPY, so there is nothing special you have to do in order to use ECOS with CVXPY, besides specifying it as a solver. Here is a small [example](http://www.cvxpy.org/en/latest/tutorial/advanced/index.html#solve-method-options) from the CVXPY tutorial:

```
# Solving a problem with different solvers.
x = Variable(2)
obj = Minimize(norm(x, 2) + norm(x, 1))
constraints = [x >= 2]
prob = Problem(obj, constraints)

# Solve with ECOS.
prob.solve(solver=ECOS)
print "optimal value with ECOS:", prob.value
```


Using ECOS in MATLAB
====

Please refer to the repository [ecos-matlab](https://github.com/embotech/ecos-matlab) hosting the ECOS Matlab interface repository.

Using ECOS in Python
====

Please refer to the repository [ecos-python](https://github.com/embotech/ecos-python) hosting the ECOS Python interface repository.

Using ECOS in C
====
ECOS exports 3 functions, see ecos.h. You need to call these in the following sequence:

1. Setup
----

Setup allocates memory for ECOS, computes the fill-in reducing ordering and provides other initialization necessary before solve can start. Use the following function to initialize ECOS:
```
pwork* ECOS_setup(idxint n, idxint m, idxint p, idxint l, idxint ncones, idxint* q,
                   pfloat* Gpr, idxint* Gjc, idxint* Gir,
                   pfloat* Apr, idxint* Ajc, idxint* Air,
                   pfloat* c, pfloat* h, pfloat* b);
```
where you have to pass the following arguments:
* `n` is the number of variables,
* `m` is the number of inequality constraints (dimension 1 of the matrix `G` and the length of the vector `h`),
* `p` is the number of equality constraints (can be 0)
* `l` is the dimension of the positive orthant, i.e. in `Gx+s=h, s in K`, the first `l` elements of `s` are `>=0`
* `ncones` is the number of second-order cones present in `K`
* `q` is an array of integers of length `ncones`, where each element defines the dimension of the cone
* `Gpr`, `Gjc`, `Gir` are the the data, the column index, and the row index arrays, respectively, for the matrix `G` represented in column compressed storage (CCS) format (Google it if you need more information on this format, it is one of the standard sparse matrix representations)
* `Apr`, `Ajc`, `Air` is the CCS representation of the matrix `A` (can be all `NULL` if no equalities are present)
* `c` is an array of type `pfloat` of size `n`
* `h` is an array of type `pfloat` of size `m`
* `b` is an array of type `pfloat` of size `p` (can be `NULL` if no equalities are present)
The setup function returns a struct of type ```pwork```, which you need to define first.

2. Solve
----
After the initialization is done, you can call the solve function, which contains the actual interior point method, by
```idxint ECOS_solve(pwork* w);```
The return value is an integer, see below.

3. Cleanup
----
Call
```void ECOS_cleanup(pwork* w, idxint keepvars);```
to free all allocated memory.

Exitcodes
====
ECOS defines a number of exitcodes that indicate the quality of the returned solution. In general, positive values indicate that the solver has converged within the given tolerance. More specifically,
+ 0: optimal
+ 1: primal infeasible
+ 2: dual infeasible

An exact definition of when these flags are returned can be found in Alexander Domahidi's [PhD Thesis](http://e-collection.library.ethz.ch/view/eth:7611?q=domahidi) in Chapter 9.4.2. (pp. 163).

Negative numbers indicate that the problem could not be solved to the required accuracy (the returned iterate might still be satisfactory - please do check the duality gap etc.)
+ -1: maximum number of iterations reached
+ -2: numerical problems (unreliable search direction)
+ -3: numerical problems (slacks or multipliers became exterior)
+ -7: unknown problem in solver

It is in general good practice to check the exitcode, in particular when solving optimization problems in an unsupervised, automated fashion (in a batch job, for example). Please report optimization problems for which ECOS struggles to converge to one of the authors.


EXTENSIONS
====

Mixed Integer SOCP

ECOS_BB is a mixed integer extension of ECOS with support for boolean variables only (full integer support is planned). ECOS_BB solves convex second-order cone programs (SOCPs) of type

```
min  c'*x
s.t. A*x = b
     G*x <=_K h
some x_i in {0,1}
```
where the last inequality is generalized, i.e. `h - G*x` belongs to the cone `K`.
ECOS_BB supports the positive orthant `R_+` and second-order cones `Q_n` defined as
```
Q_n = { (t,x) | t >= || x ||_2 } 
```
In the definition above, t is a scalar and `x` is in `R_{n-1}`. The cone `K` is therefore
a direct product of the positive orthant and second-order cones:
```
K = R_+ x Q_n1 x ... x Q_nN
```
As with ECOS, ECOS_BB is designed for embedded systems so it has hard bounds on runtime and memory footprint.

#### Runtime Constraints

ECOS_BB will call ecos_solve() on a relaxed subproblem at most `(MI_MAXITER * MAXIT)` times where `MI_MAXITER` is the maximum number of iterations allowed by the branch and bound wrapper and `MAXIT` is the maximum number of ECOS iterations allowed per sub problem. `MI_MAXITER` is set within the ecos_bb.h header and `MAXIT` is set within the ecos.h header.

#### Memory Constraints

ECOS_BB stores each active constraint per boolean var as a *char*. Thus, ECOS_BB requires `num_bool_vars` *char*s to store the full set of constraints for one branch and bound node. Since branch and bound generates two nodes per iteration and one node can reuse the parent node's memory, ECOS_BB preallocates enough memory for `MI_MAXITER` nodes. Total memory requirements for ECOS_BB is `MI_MAXITER * num_bool_vars` *bytes* + memory necessary for the relaxed ECOS subproblem.

Features of ECOS_BB
----

+ *ECOS_BB runs on embedded platforms*. Written in ANSI C (except for the timing code),
  it can be compiled for any platform for which a C compiler is available. Excluding the problem setup
  part, no memory manager is needed for solving problem instances of same structure.
+ *ECOS_BB is efficient*. 
    The branch and bound algorithm is a direct translation of Stephan Boyd's [lecture slides](http://stanford.edu/class/ee364b/lectures/bb_slides.pdf) from EE364. And has proven excellent performance in small to medium problems.
+ *ECOS_BB has a tiny footprint*. The ECOS_BB solver consists of 200 lines of C code on top of ECOS's 750 (excluding the problem setup
   code).
+ *ECOS is numerically robust*. Using regularization and iterative refinement coupled with a carefully chosen
  sparse representation of scaling matrices, millions of problem instances are solved reliably.
+ *ECOS is library-free*. No need to link any external library to ECOS, apart from `AMD` and `sparseLDL`, both
  from Timothy A. Davis, which are included in this project.


Using ECOS_BB in C
====

ECOS_BB exports 3 functions, see ecos_bb.h. You need to call these in the following sequence:

1. Setup
----

Setup allocates memory for ECOS_BB, computes the fill-in reducing ordering and provides other initialization necessary before solve can start. Use the following function to initialize ECOS_BB:
```
pwork* ecos_bb_setup(idxint n, idxint m, idxint p, idxint l, idxint ncones, idxint* q,
                   pfloat* Gpr, idxint* Gjc, idxint* Gir,
                   pfloat* Apr, idxint* Ajc, idxint* Air,
                   pfloat* c, pfloat* h, pfloat* b, idxint num_bool_vars);
```
where you have to pass the following arguments:

* `n` is the number of variables,
* `m` is the number of inequality constraints (dimension 1 of the matrix `G` and the length of the vector `h`),
* `p` is the number of equality constraints (can be 0)
* `l` is the dimension of the positive orthant, i.e. in `Gx+s=h, s in K`, the first `l` elements of `s` are `>=0`
* `ncones` is the number of second-order cones present in `K`
* `q` is an array of integers of length `ncones`, where each element defines the dimension of the cone
* `Gpr`, `Gjc`, `Gir` are the the data, the column index, and the row index arrays, respectively, for the matrix `G` represented in column compressed storage (CCS) format (Google it if you need more information on this format, it is one of the standard sparse matrix representations)
* `Apr`, `Ajc`, `Air` is the CCS representation of the matrix `A` (can be all `NULL` if no equalities are present)
* `c` is an array of type `pfloat` of size `n`
* `h` is an array of type `pfloat` of size `m`
* `b` is an array of type `pfloat` of size `p` (can be `NULL` if no equalities are present)
* `num_bool_vars` is the number of boolean variables in this problem. ECOS_bb will assume that the first num_bool_vars variables are boolean (i.e. x[0] to x[num_bool_vars-1] ).

The setup function returns a struct of type ```ecos_bb_pwork```, which you need to define first.

2. Solve
----
After the initialization is done, you can call the solve function, which contains the branch and bound wrapper around the internal point method, by
```idxint ecos_bb_solve(ecos_bb_pwork* w);```
The return value is an integer, see below.

3. Cleanup
----
Call
```void ecos_bb_cleanup(ecos_bb_pwork* w);```
to free all allocated memory.

Exitcodes
====
ECOS_BB defines a number of exitcodes that indicate the quality of the returned solution. In general, positive values indicate that the solver has converged within the given tolerance. More specifically,

+ 0: optimal
+ 1: max iterations reached no solution has been found
+ 2: problem is proven infeasible

Negative numbers indicate that the problem could not be solved to the required accuracy (the returned iterate might still be satisfactory - please do check the duality gap etc.)

+ -1: maximum number of iterations reached, but a feasible solution has been discovered

It is in general good practice to check the exitcode, in particular when solving optimization problems in an unsupervised, automated fashion (in a batch job, for example). Please report optimization problems for which ECOS_BB struggles to converge to one of the authors.                   



