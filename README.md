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
+ *ECOS comes with a MATLAB and CVX 2.0 interface*, and *it is supported by YALMIP*. With [CVX](http://cvxr.com) and
  [YALMIP](http://users.isy.liu.se/johanl/yalmip/) you can prototype, simulate and verify the performance of ECOS before you implement the very same code on 
your embedded hardware.  
+ *ECOS comes with a Python interface*. This interface is built on top of 
  [NUMPY](http://numpy.org) and [SCIPY](http://scipy.org/) and uses its sparse data structures.
+ *ECOS is library-free*. No need to link any external library to ECOS, apart from `AMD` and `sparseLDL`, both 
  from Timothy A. Davis, which are included in this project.


Credits
----

The solver is essentially based on Lieven Vandenberghe's [CVXOPT](http://cvxopt.org) [ConeLP](http://www.ee.ucla.edu/~vandenbe/publications/coneprog.pdf) solver, although it differs in the particular way the linear systems are treated.

The following people have been, and are, involved in the development and maintenance of ECOS:

+ Alexander Domahidi (principal developer)
+ Eric Chu (Python interface, testing, testing, testing)
+ Stephen Boyd (methods and maths)
+ Michael Grant (CVX interface)
+ Johan Löfberg (YALMIP interface)

The main technical idea behind ECOS is described in a short [paper](http://www.stanford.edu/~boyd/papers/ecos.html). More details are given in Alexander Domahidi's [PhD Thesis](http://e-collection.library.ethz.ch/view/eth:7611?q=domahidi).

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


Using ECOS with CVX
====

The simplest way to use ECOS is to install a CVX 2.0 shim. For this to work, you must have the latest version of 
[CVX](http://cvxr.com) installed in MATLAB. Once CVX is installed, open MATLAB and run,

     cd <ecos-directory>/matlab
     cvx_install_ecos

This will automatically build ECOS and install the CVX shim. Please report any error messages to us.

Once the ECOS shim is installed, the CVX solver can be switched using the `cvx_solver` command. For instance,

     cvx_begin
          cvx_solver ecos     % without this line, CVX will use its default solver
          variable x(n)
          
          minimize sum_square(A*x - b)
          subject to
               x >= 0
     cvx_end
     
*IMPORTANT*: Not all of CVX's atoms are SOCP-representable. Some of the atoms implemented in CVX require the use of SDP cones. Additionally, other atoms that could be implemented with a second-order cone are instead implemented as SDPs. See 
[Issue #8](https://github.com/ifa-ethz/ecos/issues/8) for more information.

Using ECOS with YALMIP
====
As of release R20130628, [YALMIP](http://users.isy.liu.se/johanl/yalmip/) supports ECOS as a solver - simply use the command
```
sdpsettings('solver','ecos');
```
to select ECOS as the solver for your problem. Below is a toy example:
```
% Solve 1000 SOCPs
x = sdpvar(3,1);
Ufo= [norm(x) <= 2, norm(x+1) <= 2];
plot(Ufo,x,'y',1000,sdpsettings('solver','ecos'))
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
    exitflag: 0=OPTIMAL, 1=PRIMAL INFEASIBLE, 2=DUAL INFEASIBLE, -1=MAXIT REACHED
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

Using ECOS in Python
====

Compiling ECOS for Python
----
To create the Python interface, the following lines of code should work:
```
cd <ecos-directory>/python
python setup.py install
```
You may need `sudo` privileges for a global installation. 

**Important note**: Even if you do not have CVXOPT installed, this code will likely compile without errors. However, you will run into dynamic link issues when attempting to use the module.

Calling ECOS from Python
----
After installing the ECOS interface, you must import the module with
```
import ecos
```
This module provides a single function `ecos` with one of the following calling sequences:
```
solution = ecos.solve(c,G,h,dims)
solution = ecos.solve(c,G,h,dims,A,b)
```
The arguments `c`, `h`, and `b` are NUMPY arrays (i.e., matrices with a single column). 
The arguments `G` and `A` are SCIPY *sparse* matrices in CSR format; if they are
not of the proper format, ECOS will attempt to convert them. 
The argument `dims` is a dictionary with two fields, `dims['l']` and `dims['q']`. These are the same fields as in the Matlab case. If the fields are omitted or empty, they default to 0. The arguments `A` and `b` are optional.

The returned object is a dictionary containing the fields `solution['x']`, `solution['y']`, `solution['s']`, `solution['z']`, and `solution['info']`. 
The first four are NUMPY arrays containing the relevant solution. The last field contains a dictionary with the same fields as the `info` struct in the MATLAB interface.
