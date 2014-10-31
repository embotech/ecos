Java/Scala drivers to Second Order Cone Programming (ECOS) Solver and Matlab comparisons of 
ECOS with Proximal Algorithms (ADMM) and Interior Point Methods (MOSEK)

====
LICENSE

All contributions by Debasish Das copyright © 2014 Verizon, and are licensed per GNU GPL v3.  
See COPYING file for details.

====
What's available in the project

1. MacOSX and Linux JNI libraries for ECOS are licensed under GPL based on original ECOS license.
2. amd and ldl JNI libraries are licensed under LGPL honoring original Tim Davis's license
3. Java packages are licensed under Apache
+ Java driver for ECOS SocpSolver is com.github.ecos.RunECOS
+ Scala driver for Quadratic Programming Solver is com.github.ecos.QpSolver

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

====

Example problem formulations in ECOS

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

Matlab comparisons of Quadratic Minimization solvers:
====

For the matlab comparisons we are focused on the following quadratic minimization formulation
which covers interest machine learning use-cases. 

Formulation1. Quadratic minimization with bounds
```
minimize 0.5*x'*H*x + f'*x
subject to lb <= x <= ub
```

Formulation2. Quadratic minimization with elastic net regularization
```
minimize 0.5*x'*H*x + f'*x + lambda*beta*||x||_1 + lambda*(1-beta)*||x||_2^{2}
subject to ||x||_1 <= c
```

Formulation3. Quadratic minimization with equality and bounds
```
minimize 0.5*x'*H*x + f'*x
subject to Ax = b, lb <= x <= ub
```

We did Matlab comparisons as MOSEK, ECOS and PDCO can be called from Matlab. ADMM based proximal algorithm 
and accelerated ADMM algorithm using Nesterov's idea are also implemented in Matlab. MOSEK and ECOS calls mex file 
and therefore should be more efficient implementations. The runtimes are as reported by tic/toc.

### Formulation1 comparisons

The driver script is in matlab/admm/qprandom.m. The input to the script are rank, lambda, beta and alpha for ADMM 
over-relaxation. Badly conditioned gram matrices are generated in matlab through rand('state', 10). We run with 
lambda=0 for running Formulation1. The last line indicate the runtime and iterations for ADMM solvers compared to 
mosek and ecos. mosek solution is considered golden.

Rank=100

qprandom(100, 0.0, 0.99, 1.0, true)

variables 100 lambdaL1 0 lambdaL2 0 beta 0.99 alpha 1 rho 53.1261

mosek-ecos norm 0.000197934

mosek-admm norm 9.26015e-06

mosek-admmfast norm 2.0384e-08

mosek 0.043701 ecos 0.0351161 admm 0.07394 iters 322 admmfast 0.22126 accIters 1084

rank=500

qprandom(500, 0.0, 0.99, 1.0, true)

variables 500 lambdaL1 0 lambdaL2 0 beta 0.99 alpha 1 rho 258.788

mosek-ecos norm 0.0016086

mosek-admm norm 3.34802e-05

mosek-admmfast norm 3.34802e-05

mosek 0.398574 ecos 1.29301 admm 0.264346 iters 400 admmfast 0.413664 accIters 635

rank=1000

qprandom(1000, 0.0, 0.99, 1.0, true)

variables 1000 lambdaL1 0 lambdaL2 0 beta 0.99 alpha 1 rho 514.312

mosek-ecos norm 4.57589e-05

mosek-admm norm 3.16294e-05

mosek-admmfast norm 3.16294e-05

mosek 1.89482 ecos 8.31283 admm 1.63412 iters 418 admmfast 2.60943 accIters 679

###Conclusion

ADMM based Proximal Algorithm is at par with MOSEK without over-relaxation. ECOS is slower than both MOSEK and ADMM.

###Formulation2 comparisons

The driver script is in matlab/admm/qprandom.m. The input to the script are rank, lambda, beta and alpha for ADMM over-relaxation. 
Badly conditioned gram matrices are generated in matlab through rand('state', 10). We run with non-zero lambda for running Formulation1. 
L1 regularization is computed as lambda*beta. L2 regularization is computed as lambda*(1-beta). 

rank=100

qprandom(100, 1.0, 0.99, 1.0, false)

variables 100 lambdaL1 0.99 lambdaL2 0.01 beta 0.99 alpha 1 rho 53.1261

mosek-ecos norm 15.151

mosek-admm norm 5.05056e-06

mosek-admmfast norm 8.29321e-09

mosek 0.245771 ecos 0.0895326 admm 0.080307 iters 499 admmfast 0.119292 accIters 771

rank=500

qprandom(500, 1.0, 0.99, 1.0, false)

variables 500 lambdaL1 0.99 lambdaL2 0.01 beta 0.99 alpha 1 rho 258.788

mosek-ecos norm 2.94544

mosek-admm norm 1.40105e-06

mosek-admmfast norm 5.53243e-08

mosek 1.78544 ecos 3.23721 admm 0.305345 iters 460 admmfast 0.527242 accIters 806

rank=1000

qprandom(1000, 1.0, 0.99, 1.0, false)

variables 1000 lambdaL1 0.99 lambdaL2 0.01 beta 0.99 alpha 1 rho 514.312

mosek-ecos norm 1.84793

mosek-admm norm 2.87457e-06

mosek-admmfast norm 2.87457e-06

mosek 10.0573 ecos 23.868 admm 1.0004 iters 390 admmfast 1.44243 accIters 569

###Conclusions
ADMM based proximal algorithm is ~10X faster than MOSEK. ECOS here failed since ecosqp matlab interface 
assume positive definite gram matrix which is not true for L1 formulation as a Quadratic Program 

###Formulation3 comparisons

We picked pdco along with mosek and ecos for the equality and bound formulation comparisons. PDCO was suggested by
Xiangrui Meng at our Spark Summit 2014 talk http://spark-summit.org/2014/talk/quadratic-programing-solver-for-non-negative-matrix-factorization-with-spark. 
The script matlab/pdco4/code/pdcotestQP.m takes m equalities and n as rank. It generates dense gram matrix of size nxn and 
equality constraint Ax = b along with lower and upper bounds randomly with bad conditioning. The script also takes an 
alpha for over-relaxation of ADMM but it is fixed at 1.0

###PLSA and SVM formulations

rank=100, equality=1

pdcotestQP(1, 100, 1.0)

mosek-ecos norm 0.00506569

mosek-pdco norm 0.000303848

mosek-admm norm 3.52694e-05 iters 212

mosek-accadmm norm 5.58596e-08 iters 1677

mosek 0.039641 ecos 0.0341901 pdco 0.08 admm 0.041568 admmacc 0.239477

rank=500, equality=1

pdcotestQP(1, 500, 1.0)

mosek-ecos norm 0.220284

mosek-pdco norm 9.08018

mosek-admm norm 7.23129e-05 iters 987

mosek-accadmm norm 1.31676e-05 iters 5545

mosek 0.323458 ecos 2.89661 pdco 4.44 admm 0.647283 admmacc 3.56523

rank=1000, equality=1

pdcotestQP(1, 1000, 1.0)

mosek-ecos norm 6.28512

mosek-pdco norm 8.9423

mosek-admm norm 0.00140746 iters 2041

mosek-accadmm norm 0.00140746 iters 20000

mosek 1.58106 ecos 18.3085 pdco 28.95 admm 12.7529 admmacc 120.245

###Conclusions
MOSEK is at par with ADMM for low ranks but faster when ranks are high. Both PDCO and ECOS produce bad quality results. 
ADMM is closest to MOSEK in terms of solution quality. PDCO at rank=500+ failed to converge as well. 
It is worth noting that both ADMM and PDCO are running LU factorization.

###Generic Quadratic Program with Equality and bound constraints

rank=1000 equality=50

pdcotestQP(50, 1000, 1.0)

mosek-ecos norm 9.6798

mosek-pdco norm 13.2054

mosek-admm norm 0.00474487 iters 1356

mosek-accadmm norm 0.00474487 iters 5356

mosek 2.08932 ecos 20.3385 pdco 32.13 admm 7.01939 admmacc 27.5484

###Conclusions
ADMM is close to MOSEK both in terms of soluton accuracy and runtime.

###Spark integration of ADMM based Proximal Algorithm

We added ADMM based Proximal Algorithm as a Quadratic Minimization solver for Spark mllib
ALS based recommendation. Each ALS subproblem is solved as an instance of 
Quadratic Minimization solver with positivity, bounds, elastic-net or equality with positivity.
Further details on machine learning use-cases and dataset experiments are presented at:

1. Spark JIRA https://issues.apache.org/jira/browse/SPARK-2426
2. Spark PR https://github.com/apache/spark/pull/2705

