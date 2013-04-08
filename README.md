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

### Example: L1 minimization (Linear Programming)

In the following, we show how to solve a L1 minimization problem, which arises for example in sparse signal
reconstruction problems (compressed sensing):
```
    minimize  ||x||_1         (L1)
  subject to  Ax = b
```
where `x` is in `R^n`, `A` in `R^{m x n}` with `m <= n`. We use the epigraph reformulation to express the L1-norm of `x`, 
```
    x <= u
   -x <= v
```
where `u,v` are in `R^n`, and we minimize `sum(u) + sum(v)`. Hence the optimization variables are stacked as follows:
```
   z = [x; u; v]
```
With this reformulation, (L1) can be written as linear program (LP),
``` 
  minimize   c'*z
  subject to Atilde*z = b;    (LP)
             Gx <= h
```
where the inequality is w.r.t. the positive orthant. The following MATLAB code generates a random instance of this problem
and calls ECOS to solve the problem:
```
% set dimensions and sparsity of A
n = 1000; 
m = 10;
density = 0.01;

% linear term
c = [zeros(n,1); ones(2*n,1)];

% equality constraints
A = sprandn(m,n,density);
Atilde = [A, zeros(m,2*n)];
b = randn(m,1);

% linear inequality constraints
I = speye(n); O = zeros(n);
G = [  I -I  O;
      -I  O -I];
h = zeros(2*n,1);

% cone dimensions (LP cone only)
dims.l = 2*n;
dims.q = [];

% call solver
fprintf('Calling solver...');
z = ecos(c,G,h,dims,Atilde,b);
x = z(1:n);
u = z(n+1:2*n);
v = z(2*n+1:3*n);
nnzx = sum(abs(x) > 1e-8);

% print sparsity info
fprintf('Optimal x has %d/%d (%4.2f%%) non-zero (>1e-8 in abs. value) entries.\n', nnzx , n,  nnzx/n*100);
```

### Example: Quadratic Programming

In this example, we consider problems of form
```
  minimize 0.5*x'*H*x + f'*x
subject to A*x <= b                (QP)
           Aeq*x = beq
           lb <= x <= ub
```
where we assume that `H` is positive definite. This is the standard formulation that also MATLAB's built-in solver
`quadprog` uses. To deal with the quadratic objective, you have to reformulate it into a second-order cone constraint to
directly call ECOS. We do provide a MATLAB interface called `ecosqp` that automatically does this transformation for you,
and has the exact same interface as `quadprog`. Hence you can just use
```
[x,fval,exitflag,output,lambda,t] = ecosqp(H,f,A,b,Aeq,beq,lb,ub)
```
to solve (QP). See `help ecosqp` for more details. The last output argument, `t`, gives the solution time.
