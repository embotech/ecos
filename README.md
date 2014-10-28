Java/Scala drivers to Embedded Conic Solver (ECOS) and comparisons with Proximal and Interior Point Methods
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

Details of ECOS internals are described in the original repository https://github.com/ifa-ethz/ecos

ECOS Credits
----

Core ECOS solver is essentially based on Lieven Vandenberghe's [CVXOPT](http://cvxopt.org) [ConeLP](http://www.ee.ucla.edu/~vandenbe/publications/coneprog.pdf) solver, although it differs in the particular way the linear systems are treated.

The following people have been, and are, involved in the development and maintenance of ECOS:

+ Alexander Domahidi (principal developer)
+ Eric Chu (Python interface, unit tests)
+ Stephen Boyd (methods and maths)
+ Michael Grant (CVX interface)
+ Johan Löfberg (YALMIP interface)
+ Karanveer Mohan (Julia interface)
+ Debasish Das (Java/Scala interfaces)

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

You can directly call ECOS from Matlab using its native interface:
```
[x,y,info,s,z] = ecos(c,G,h,dims,A,b)
```
It takes the problem data `c,G,h,A,b` and some dimension information that is given in the struct `dims`. Note that
`A` and `G` have to be given in sparse format. The equality constraints defined by `A` and `b` are optional and can be
omitted. The `dims` structure has the following fields:
```
dims.l - scalar, dimension of positive orthant (LP-cone) R_+
dims.q - vector with dimensions of second order cones
```
The length of `dims.q` determines the number of second order cones. If you do not have a cone in your problem, use
the empty matrix `[ ]` instead, for example `dims.q = [ ]` if you do not have second-order cones. After a solve,
ECOS returns the following variables
```
  x: primal variables
  y: dual variables for equality constraints
  s: slacks for Gx + s <= h, s \in K
  z: dual variables for inequality constraints s \in K
```
In addition, the struct `info` is returned which contains the following fields:
```
    exitflag: 0=OPTIMAL, 1=PRIMAL INFEASIBLE, etc. (see exitcodes section in this readme)
  infostring: gives information about the status of solution
       pcost: value of primal objective
       dcost: value of dual objective
        pres: primal residual on inequalities and equalities
        dres: dual residual
        pinf: primal infeasibility measure
        dinf: dual infeasibility measure
     pinfres: NaN
     dinfres: 3.9666e+15
         gap: duality gap
      relgap: relative duality gap
          r0: ???
      numerr: numerical error?
        iter: number of iterations
      timing: struct with timing information
```

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
   -x <= u
```
where `u` is in `R^n`, and we minimize `sum(u)`. Hence the optimization variables are stacked as follows:
```
   z = [x; u]
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
c = [zeros(n,1); ones(n,1)];

% equality constraints
A = sprandn(m,n,density);
Atilde = [A, zeros(m,n)];
b = randn(m,1);

% linear inequality constraints
I = speye(n);
G = [  I -I;
      -I -I];
h = zeros(2*n,1);

% cone dimensions (LP cone only)
dims.l = 2*n;
dims.q = [];

% call solver
fprintf('Calling solver...');
z = ecos(c,G,h,dims,Atilde,b);
x = z(1:n);
u = z(n+1:2*n);
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

Using ECOS in Java/Scala
====

1. MacOSX and Linux JNI libraries are included and published in the jar
2. amd and ldl JNI libraries are licensed under LGPL
3. ecos JNI library is licensed under GPL
4. Java packages are licensed under Apache
+ Java driver for ECOS SocpSolver is com.github.ecos.RunECOS
+ Scala driver for Quadratic Programming Solver is com.github.ecos.QpSolver

Matlab comparisons of Quadratic Minimization solvers:
====

1. ECOS based Quadratic Programming
2. MOSEK based Quadratic Programming
3. ADMM and Accelerated ADMM based Quadratic Minimization
4. PDCO based Quadratic Programming



