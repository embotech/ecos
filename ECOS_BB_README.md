Using ECOS_BB in C
====
ECOS_BB is a mixed integer extension of ECOS with support for boolean variables only (full integer support is planned).

As with ECOS, ECOS_BB is designed for embedded systems so it has hard bounds on runtime and memory footprint.

ECOS_BB will call ecos_solve() on a relaxed subproblem at most (MI_MAXITER * MAXIT) times where MI_MAXITER is the maximum number of iterations allowed by the branch and bound wrapper and MAXIT is the maximum number of ECOS iterations allowed per sub problem. MI_MAXITER is set within the ecos_bb.h header and MAXIT is set within the ecos.h header.

Also, the memory required for ecos_bb is O(MI_MAXITER * num_bool_vars).

The branch and bound algorithm is a direct translation of Stephan Boyd's [lecture slides](http://stanford.edu/class/ee364b/lectures/bb_slides.pdf) from EE364.

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
