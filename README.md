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

We did Matlab comparisons as MOSEK, ECOS and PDCO can be called from Matlab. ADMM based proximal algorithm and accelerated ADMM algorithm using Nesterov's idea are also implemented in Matlab. MOSEK and ECOS calls mex file and therefore should be more efficient implementations. The runtimes are as reported by tic/toc

### Formulation1 comparisons

The driver script is in matlab/admm/qprandom.m. The input to the script are rank, lambda, beta and alpha for ADMM over-relaxation. Badly conditioned gram matrices are generated in matlab through rand('state', 10). We run with lambda=0 for running Formulation1. The last line indicate the runtime and iterations for ADMM solvers compared to mosek and ecos. mosek solution is considered golden.

Rank=100

qprandom(100, 0.0, 0.99, 1.0, true)

eigenMinApprox 0.000141167 eigenMaxApprox 2822.39 rho 53.1261 eigenMin 0.000141167 eigenMax 2481.31
variables 100 lambdaL1 0 lambdaL2 0 beta 0.99 alpha 1 rho 53.1261
mosek done
ECOSQP: Time for Cholesky: 0.00 seconds
ecos done
mosek-ecos norm 0.000197934
mosek-admm norm 9.26015e-06
mosek-admmfast norm 2.0384e-08
mosek 0.043701 ecos 0.0351161 admm 0.07394 iters 322 admmfast 0.22126 accIters 1084

rank=500

qprandom(500, 0.0, 0.99, 1.0, true)

eigenMinApprox 0.000135883 eigenMaxApprox 66971.1 rho 258.788 eigenMin 0.000135883 eigenMax 62618.1
variables 500 lambdaL1 0 lambdaL2 0 beta 0.99 alpha 1 rho 258.788
mosek done
ECOSQP: Time for Cholesky: 0.01 seconds
ecos done
mosek-ecos norm 0.0016086
mosek-admm norm 3.34802e-05
mosek-admmfast norm 3.34802e-05
mosek 0.398574 ecos 1.29301 admm 0.264346 iters 400 admmfast 0.413664 accIters 635

rank=1000

qprandom(1000, 0.0, 0.99, 1.0, true)

eigenMinApprox 8.77638e-07 eigenMaxApprox 264517 rho 514.312 eigenMin 8.77638e-07 eigenMax 250499
variables 1000 lambdaL1 0 lambdaL2 0 beta 0.99 alpha 1 rho 514.312
mosek done
ECOSQP: Time for Cholesky: 0.03 seconds
ecos done
mosek-ecos norm 4.57589e-05
mosek-admm norm 3.16294e-05
mosek-admmfast norm 3.16294e-05
mosek 1.89482 ecos 8.31283 admm 1.63412 iters 418 admmfast 2.60943 accIters 679

###Conclusion

ADMM based Proximal Algorithm is at par with MOSEK without over-relaxation. ECOS is slower than both MOSEK and ADMM.

###Formulation2 comparisons

The driver script is in matlab/admm/qprandom.m. The input to the script are rank, lambda, beta and alpha for ADMM over-relaxation. Badly conditioned gram matrices are generated in matlab through rand('state', 10). We run with non-zero lambda for running Formulation1. L1 regularization is computed as lambda*beta. L2 regularization is computed as lambda*(1-beta). 

rank=100

qprandom(100, 1.0, 0.99, 1.0, false)
eigenMinApprox 0.000141167 eigenMaxApprox 2822.39 rho 53.1261 eigenMin 0.000141167 eigenMax 2481.31
variables 100 lambdaL1 0.99 lambdaL2 0.01 beta 0.99 alpha 1 rho 53.1261
mosek done
Warning: Hessian not positive definite, using sqrt(H) instead of chol 
  In ecosqp at 135
  In qpdriver at 44
  In qprandom at 19 
ecos done
mosek-ecos norm 15.151
mosek-admm norm 5.05056e-06
mosek-admmfast norm 8.29321e-09
mosek 0.245771 ecos 0.0895326 admm 0.080307 iters 499 admmfast 0.119292 accIters 771

rank=500

qprandom(500, 1.0, 0.99, 1.0, false)

eigenMinApprox 0.000135883 eigenMaxApprox 66971.1 rho 258.788 eigenMin 0.000135883 eigenMax 62618.1
variables 500 lambdaL1 0.99 lambdaL2 0.01 beta 0.99 alpha 1 rho 258.788
mosek done
Warning: Hessian not positive definite, using sqrt(H) instead of chol 
  In ecosqp at 135
  In qpdriver at 44
  In qprandom at 19 
ecos done
mosek-ecos norm 2.94544
mosek-admm norm 1.40105e-06
mosek-admmfast norm 5.53243e-08
mosek 1.78544 ecos 3.23721 admm 0.305345 iters 460 admmfast 0.527242 accIters 806

rank=1000

>> qprandom(1000, 1.0, 0.99, 1.0, false)
eigenMinApprox 8.77638e-07 eigenMaxApprox 264517 rho 514.312 eigenMin 8.77638e-07 eigenMax 250499
variables 1000 lambdaL1 0.99 lambdaL2 0.01 beta 0.99 alpha 1 rho 514.312
mosek done
Warning: Hessian not positive definite, using sqrt(H) instead of chol 
  In ecosqp at 135
  In qpdriver at 44
  In qprandom at 19 
ecos done
mosek-ecos norm 1.84793
mosek-admm norm 2.87457e-06
mosek-admmfast norm 2.87457e-06
mosek 10.0573 ecos 23.868 admm 1.0004 iters 390 admmfast 1.44243 accIters 569

###Conclusions
ADMM based proximal algorithm is ~10X faster than MOSEK. ECOS here failed since ecosqp matlab interface 
assume positive definite gram matrix which is not true for L1 formulation as a Quadratic Program 

###Formulation3 comparisons

We picked pdco along with mosek and ecos for the equality and bound formulation. PDCO was suggested by
Xiangrui Meng at our Spark Summit 2014 talk. The script matlab/pdco4/code/pdcotestQP.m takes m equalities
and n as rank. It generates dense gram matrix of size nxn and equality constraint Ax = b along with lower
and upper bounds randomly with bad conditioning. The script also takes an alpha for over-relaxation of
ADMM but it is fixed at 1.0

###PLSA and SVM formulations

rank=100, equality=1

Generating Dense Hessian and Inequalities
ECOSQP: Time for Cholesky: 0.00 seconds
eigenMinApprox 0.00420724 eigenMaxApprox 2824.68 rho 53.1477 eigenMin 0.00420724 eigenMax 2485.63
rho 53.1477 alpha 1

   --------------------------------------------------------
   pdco.m                            Version of 23 Nov 2013
   Primal-dual barrier method to minimize a convex function
   subject to linear constraints Ax + r = b,  bl <= x <= bu
                                                           
   Michael Saunders       SOL and ICME, Stanford University
   Contributors:     Byunggyoo Kim (SOL), Chris Maes (ICME)
                     Santiago Akle (ICME), Matt Zahr (ICME)
   --------------------------------------------------------

The objective is defined by a function handle: 
   @(x)deal(0.5*(x'*H*x)+c'*x,H*x+c,H)
The matrix A is an explicit dense matrix

m        =        1     n        =      100      nnz(A)  =      100
max |b | =  49.2086     max |x0| =  1.0e-02      xsize   =  5.0e+00
max |y0| =        0     max |z0| =  1.0e-02      zsize   =  8.0e+00

x0min    =        1     featol   =  1.0e-06      d1max   =  1.0e-03
z0min    =        1     opttol   =  1.0e-06      d2max   =  1.0e-03
mu0      =  0.0e+00     steptol  =     0.99     bigcenter=     1000

LSMR/MINRES:
atol1    =  1.0e-10     atol2    =  1.0e-15      btol    =  0.0e+00
conlim   =  1.0e+12     itnlim   =       10      show    =        0

Method   =       21     (1=chol  2=QR  3=LSMR  4=MINRES  21=SQD(LU)  22=SQD(MA57))


Bounds:
  [0,inf]  [-inf,0]  Finite bl  Finite bu  Two bnds   Fixed    Free
        0         0        100        100       100       0       0
  [0, bu]  [bl,  0]  excluding fixed variables
      100         0

Itn   mu stepx stepz  Pinf  Dinf  Cinf   Objective    nf  center      SQD
  0                    1.6   3.2   0.3  3.0949013e+06        2.0
  1 -0.8 0.522 0.522   1.3   2.9   0.1  7.0914626e+05  1    74.1     5151 x2
  2 -0.8 0.349 0.349   1.1   2.7  -0.1  3.0275466e+05  1   149.3
  3 -0.8 0.205 0.205   1.0   2.6  -0.2  1.9377586e+05  1   206.1
  4 -0.8 0.166 0.166   0.9   2.6  -0.2  1.3815738e+05  1   152.7
  5 -0.8 0.089 0.089   0.9   2.5  -0.2  1.1804173e+05  1   170.3
  6 -0.8 0.121 0.121   0.8   2.5  -0.3  9.6983202e+04  1   104.2
  7 -0.8 0.123 0.123   0.8   2.4  -0.3  8.2356043e+04  1   476.7
  8 -0.8 0.116 0.116   0.7   2.4  -0.4  7.2732990e+04  1   156.8
  9 -0.8 0.153 0.153   0.6   2.3  -0.4  6.3754370e+04  1   156.9
 10 -0.8 0.160 0.160   0.6   2.2  -0.5  5.7025526e+04  1   114.4
 11 -0.8 0.186 0.186   0.5   2.1  -0.5  5.1151984e+04  1   200.5
 12 -0.8 0.476 0.476   0.2   1.8  -0.6  4.0158049e+04  1    68.8
 13 -0.8 1.000 1.000 -15.3 -12.6   0.1  2.9924585e+04  1    19.4
 14 -1.0 1.000 1.000 -14.6 -12.6  -0.6  2.9726372e+04  1     2.6
 15 -1.8 1.000 1.000 -14.3 -12.4  -1.5  2.9411142e+04  1     3.0
 16 -2.7 1.000 1.000 -14.5 -12.4  -2.3  2.9360812e+04  1     2.4
 17 -3.6 1.000 1.000 -15.1 -12.5  -3.4  2.9354244e+04  1     1.5
 18 -5.0 1.000 1.000 -15.0 -12.5  -4.9  2.9353349e+04  1     1.1
 19 -7.0 0.999 0.999 -14.7 -12.6  -6.9  2.9353308e+04  1     1.0
   Converged

max |x| =     2.000    max |y| =   152.872    max |z| =   159.117   scaled
max |x| =    10.000    max |y| =  1222.973    max |z| =  1272.939 unscaled
PDitns  =        19    SQDitns =         0    cputime =       0.1

Distribution of vector     x         z
[  1e+03,  1e+04 )         0        21
[    100,  1e+03 )         0        65
[     10,    100 )         0         5
[      1,     10 )         8         1
[    0.1,      1 )         2         0
[   0.01,    0.1 )         0         0
[  0.001,   0.01 )         0         0
[ 0.0001,  0.001 )         0         0
[  1e-05, 0.0001 )         0         0
[      0,  1e-05 )        90         8 
Elapsed time is 0.079090 seconds.
mosek-ecos norm 0.00506569
mosek-pdco norm 0.000303848
mosek-admm norm 3.52694e-05 iters 212
mosek-accadmm norm 5.58596e-08 iters 1677
mosek 0.039641 ecos 0.0341901 pdco 0.08 admm 0.041568 admmacc 0.239477

rank=500, equality=1

Generating Dense Hessian and Inequalities
ECOSQP: Time for Cholesky: 0.00 seconds
eigenMinApprox 0.000226668 eigenMaxApprox 66962.1 rho 258.77 eigenMin 0.000226668 eigenMax 62605.2
rho 258.77 alpha 1

   --------------------------------------------------------
   pdco.m                            Version of 23 Nov 2013
   Primal-dual barrier method to minimize a convex function
   subject to linear constraints Ax + r = b,  bl <= x <= bu
                                                           
   Michael Saunders       SOL and ICME, Stanford University
   Contributors:     Byunggyoo Kim (SOL), Chris Maes (ICME)
                     Santiago Akle (ICME), Matt Zahr (ICME)
   --------------------------------------------------------

The objective is defined by a function handle: 
   @(x)deal(0.5*(x'*H*x)+c'*x,H*x+c,H)
The matrix A is an explicit dense matrix

m        =        1     n        =      500      nnz(A)  =      500
max |b | =  249.573     max |x0| =  2.0e-03      xsize   =  5.0e+00
max |y0| =        0     max |z0| =  2.0e-03      zsize   =  8.0e+00

x0min    =        1     featol   =  1.0e-06      d1max   =  1.0e-03
z0min    =        1     opttol   =  1.0e-06      d2max   =  1.0e-03
mu0      =  0.0e+00     steptol  =     0.99     bigcenter=     1000

LSMR/MINRES:
atol1    =  1.0e-10     atol2    =  1.0e-15      btol    =  0.0e+00
conlim   =  1.0e+12     itnlim   =       10      show    =        0

Method   =       21     (1=chol  2=QR  3=LSMR  4=MINRES  21=SQD(LU)  22=SQD(MA57))


Bounds:
  [0,inf]  [-inf,0]  Finite bl  Finite bu  Two bnds   Fixed    Free
        0         0        500        500       500       0       0
  [0, bu]  [bl,  0]  excluding fixed variables
      500         0

Itn   mu stepx stepz  Pinf  Dinf  Cinf   Objective    nf  center      SQD
  0                    2.3   4.6   0.3  3.9103891e+08        2.0
  1 -0.8 0.388 0.388   2.1   4.4   0.1  1.4646501e+08  1    81.1   125751 x2
  2 -0.8 0.169 0.169   2.0   4.3   0.1  1.0122174e+08  1   343.9
  3 -0.8 0.176 0.176   1.9   4.2  -0.0  6.8755817e+07  1   263.1
  4 -0.8 0.144 0.144   1.9   4.2  -0.1  5.0379486e+07  1   173.9
  5 -0.8 0.124 0.124   1.8   4.1  -0.1  3.8657877e+07  1   220.5
  6 -0.8 0.100 0.100   1.8   4.1  -0.1  3.1367784e+07  1   387.0
  7 -0.8 0.065 0.065   1.7   4.0  -0.2  2.7424928e+07  1   180.9
  8 -0.8 0.061 0.061   1.7   4.0  -0.2  2.4209530e+07  1   139.8
  9 -0.8 0.043 0.043   1.7   4.0  -0.2  2.2180595e+07  1   386.3
 10 -0.8 0.041 0.041   1.7   4.0  -0.2  2.0415808e+07  1   184.8
 11 -0.8 0.035 0.035   1.6   4.0  -0.2  1.9051144e+07  1   190.1
 12 -0.8 0.037 0.037   1.6   3.9  -0.2  1.7735483e+07  1   264.9
 13 -0.8 0.024 0.024   1.6   3.9  -0.2  1.6946439e+07  1   287.1
 14 -0.8 0.021 0.021   1.6   3.9  -0.3  1.6296062e+07  1   283.9
 15 -0.8 0.016 0.016   1.6   3.9  -0.3  1.5819877e+07  1   483.2
 16 -0.8 0.018 0.018   1.6   3.9  -0.3  1.5313816e+07  1   224.0
 17 -0.8 0.019 0.019   1.6   3.9  -0.3  1.4800320e+07  1   258.3
 18 -0.8 0.016 0.016   1.6   3.9  -0.3  1.4410316e+07  1   268.5
 19 -0.8 0.015 0.015   1.6   3.9  -0.3  1.4053937e+07  1   341.7
 20 -0.8 0.020 0.020   1.6   3.9  -0.3  1.3599758e+07  1   131.9
 21 -0.8 0.005 0.005   1.6   3.9  -0.3  1.3487642e+07  1   700.1
 22 -0.8 0.013 0.013   1.6   3.9  -0.3  1.3230133e+07  1   361.9
 23 -0.8 0.014 0.014   1.5   3.9  -0.3  1.2955680e+07  1   133.9
 24 -0.8 0.014 0.014   1.5   3.9  -0.3  1.2705153e+07  1   193.8
 25 -0.8 0.014 0.014   1.5   3.9  -0.3  1.2469233e+07  1   189.6
 26 -0.8 0.012 0.012   1.5   3.9  -0.3  1.2276199e+07  1   623.1
 27 -0.8 0.013 0.013   1.5   3.8  -0.3  1.2079342e+07  1   164.8
 28 -0.8 0.013 0.013   1.5   3.8  -0.3  1.1897226e+07  1  6161.2
 29 -0.8 0.014 0.014   1.5   3.8  -0.3  1.1718930e+07  1   354.9
 30 -0.8 0.013 0.013   1.5   3.8  -0.3  1.1553821e+07  1   155.8
   Too many iterations

max |x| =     2.158    max |y| =   493.734    max |z| =   456.611   scaled
max |x| =    10.788    max |y| =  3949.869    max |z| =  3652.890 unscaled
PDitns  =        30    SQDitns =         0    cputime =       4.4

Distribution of vector     x         z
[  1e+04,  1e+05 )         0       500
[  1e+03,  1e+04 )         0         0
[    100,  1e+03 )         0         0
[     10,    100 )        13         0
[      1,     10 )        51         0
[    0.1,      1 )         8         0
[   0.01,    0.1 )        26         0
[  0.001,   0.01 )       298         0
[ 0.0001,  0.001 )       104         0
[      0, 0.0001 )         0         0 
Elapsed time is 3.358827 seconds.
mosek-ecos norm 0.220284
mosek-pdco norm 9.08018
mosek-admm norm 7.23129e-05 iters 987
mosek-accadmm norm 1.31676e-05 iters 5545
mosek 0.323458 ecos 2.89661 pdco 4.44 admm 0.647283 admmacc 3.56523

rank=1000, equality=1

>> pdcotestQP(1, 1000, 1.0)
Generating Dense Hessian and Inequalities
ECOSQP: Time for Cholesky: 0.03 seconds
eigenMinApprox 1.19736e-06 eigenMaxApprox 264516 rho 514.311 eigenMin 1.19736e-06 eigenMax 250493
rho 514.311 alpha 1

   --------------------------------------------------------
   pdco.m                            Version of 23 Nov 2013
   Primal-dual barrier method to minimize a convex function
   subject to linear constraints Ax + r = b,  bl <= x <= bu
                                                           
   Michael Saunders       SOL and ICME, Stanford University
   Contributors:     Byunggyoo Kim (SOL), Chris Maes (ICME)
                     Santiago Akle (ICME), Matt Zahr (ICME)
   --------------------------------------------------------

The objective is defined by a function handle: 
   @(x)deal(0.5*(x'*H*x)+c'*x,H*x+c,H)
The matrix A is an explicit dense matrix

m        =        1     n        =     1000      nnz(A)  =     1000
max |b | =  490.417     max |x0| =  1.0e-03      xsize   =  5.0e+00
max |y0| =        0     max |z0| =  1.0e-03      zsize   =  8.0e+00

x0min    =        1     featol   =  1.0e-06      d1max   =  1.0e-03
z0min    =        1     opttol   =  1.0e-06      d2max   =  1.0e-03
mu0      =  0.0e+00     steptol  =     0.99     bigcenter=     1000

LSMR/MINRES:
atol1    =  1.0e-10     atol2    =  1.0e-15      btol    =  0.0e+00
conlim   =  1.0e+12     itnlim   =       10      show    =        0

Method   =       21     (1=chol  2=QR  3=LSMR  4=MINRES  21=SQD(LU)  22=SQD(MA57))


Bounds:
  [0,inf]  [-inf,0]  Finite bl  Finite bu  Two bnds   Fixed    Free
        0         0       1000       1000      1000       0       0
  [0, bu]  [bl,  0]  excluding fixed variables
     1000         0

Itn   mu stepx stepz  Pinf  Dinf  Cinf   Objective    nf  center      SQD
  0                    2.6   5.2   0.3  3.1300828e+09        2.0
  1 -0.8 0.358 0.358   2.4   5.0   0.1  1.2885906e+09  1    82.6   501501 x2
  2 -0.8 0.108 0.108   2.4   5.0   0.1  1.0243840e+09  1  1063.7
  3 -0.8 0.144 0.144   2.3   4.9   0.0  7.5111452e+08  1   438.6
  4 -0.8 0.134 0.134   2.2   4.8  -0.0  5.6318206e+08  1   252.8
  5 -0.8 0.106 0.106   2.2   4.8  -0.0  4.4979895e+08  1   232.7
  6 -0.8 0.110 0.110   2.1   4.7  -0.1  3.5629320e+08  1   173.5
  7 -0.8 0.074 0.074   2.1   4.7  -0.1  3.0558299e+08  1   271.6
  8 -0.8 0.092 0.092   2.0   4.7  -0.1  2.5211230e+08  1   576.7
  9 -0.8 0.057 0.057   2.0   4.6  -0.2  2.2406357e+08  1   346.6
 10 -0.8 0.053 0.053   2.0   4.6  -0.2  2.0082254e+08  1   236.0
 11 -0.8 0.037 0.037   2.0   4.6  -0.2  1.8637662e+08  1   312.5
 12 -0.8 0.039 0.039   2.0   4.6  -0.2  1.7234863e+08  1   310.3
 13 -0.8 0.033 0.033   1.9   4.6  -0.2  1.6129356e+08  1   399.5
 14 -0.8 0.028 0.028   1.9   4.6  -0.2  1.5253778e+08  1   283.0
 15 -0.8 0.024 0.024   1.9   4.6  -0.2  1.4528429e+08  1   248.1
 16 -0.8 0.016 0.016   1.9   4.5  -0.2  1.4073026e+08  1   455.3
 17 -0.8 0.022 0.022   1.9   4.5  -0.3  1.3458332e+08  1   246.4
 18 -0.8 0.014 0.014   1.9   4.5  -0.3  1.3088041e+08  1   204.8
 19 -0.8 0.016 0.016   1.9   4.5  -0.3  1.2693311e+08  1   253.5
 20 -0.8 0.014 0.014   1.9   4.5  -0.3  1.2343544e+08  1   158.9
 21 -0.8 0.013 0.013   1.9   4.5  -0.3  1.2027249e+08  1   267.7
 22 -0.8 0.012 0.012   1.9   4.5  -0.3  1.1753837e+08  1   291.6
 23 -0.8 0.007 0.007   1.9   4.5  -0.3  1.1599386e+08  1   777.5
 24 -0.8 0.014 0.014   1.9   4.5  -0.3  1.1304516e+08  1   193.7
 25 -0.8 0.008 0.008   1.9   4.5  -0.3  1.1142308e+08  1   924.5
 26 -0.8 0.007 0.007   1.9   4.5  -0.3  1.1001115e+08  1   258.8
 27 -0.8 0.007 0.007   1.9   4.5  -0.3  1.0870063e+08  1  7597.6
 28 -0.8 0.013 0.013   1.9   4.5  -0.3  1.0623539e+08  1   243.9
 29 -0.8 0.004 0.004   1.9   4.5  -0.3  1.0542355e+08  1   411.4
 30 -0.8 0.008 0.008   1.8   4.5  -0.3  1.0394653e+08  1   450.4
   Too many iterations

max |x| =     2.159    max |y| =   341.786    max |z| =   305.990   scaled
max |x| =    10.797    max |y| =  2734.288    max |z| =  2447.922 unscaled
PDitns  =        30    SQDitns =         0    cputime =      29.0

Distribution of vector     x         z
[  1e+05,  1e+06 )         0      1000
[  1e+04,  1e+05 )         0         0
[  1e+03,  1e+04 )         0         0
[    100,  1e+03 )         0         0
[     10,    100 )         4         0
[      1,     10 )       175         0
[    0.1,      1 )        34         0
[   0.01,    0.1 )        89         0
[  0.001,   0.01 )       668         0
[      0,  0.001 )        30         0 
Elapsed time is 23.448384 seconds.
mosek-ecos norm 6.28512
mosek-pdco norm 8.9423
mosek-admm norm 0.00140746 iters 2041
mosek-accadmm norm 0.00140746 iters 20000
mosek 1.58106 ecos 18.3085 pdco 28.95 admm 12.7529 admmacc 120.245

###Conclusions
MOSEK is at par with ADMM for low ranks but faster when ranks are high. Both PDCO and ECOS produce bad quality results and ADMM is closest to MOSEK in terms of solution quality. PDCO at rank=500+ failed to converge as well. It is worth noting that both ADMM and PDCO are running LU factorization.

###Generic Quadratic Program with Equality and bound constraints

rank=1000 equality=50

>> pdcotestQP(50, 1000, 1.0)
Generating Dense Hessian and Inequalities
ECOSQP: Time for Cholesky: 0.03 seconds
eigenMinApprox 0.0729589 eigenMaxApprox 264565 rho 514.358 eigenMin 0.0641894 eigenMax 250555
rho 514.358 alpha 1

   --------------------------------------------------------
   pdco.m                            Version of 23 Nov 2013
   Primal-dual barrier method to minimize a convex function
   subject to linear constraints Ax + r = b,  bl <= x <= bu
                                                           
   Michael Saunders       SOL and ICME, Stanford University
   Contributors:     Byunggyoo Kim (SOL), Chris Maes (ICME)
                     Santiago Akle (ICME), Matt Zahr (ICME)
   --------------------------------------------------------

The objective is defined by a function handle: 
   @(x)deal(0.5*(x'*H*x)+c'*x,H*x+c,H)
The matrix A is an explicit dense matrix

m        =       50     n        =     1000      nnz(A)  =    50000
max |b | =  520.321     max |x0| =  1.0e-03      xsize   =  5.0e+00
max |y0| =        0     max |z0| =  1.0e-03      zsize   =  8.0e+00

x0min    =        1     featol   =  1.0e-06      d1max   =  1.0e-03
z0min    =        1     opttol   =  1.0e-06      d2max   =  1.0e-03
mu0      =  0.0e+00     steptol  =     0.99     bigcenter=     1000

LSMR/MINRES:
atol1    =  1.0e-10     atol2    =  1.0e-15      btol    =  0.0e+00
conlim   =  1.0e+12     itnlim   =      500      show    =        0

Method   =       21     (1=chol  2=QR  3=LSMR  4=MINRES  21=SQD(LU)  22=SQD(MA57))


Bounds:
  [0,inf]  [-inf,0]  Finite bl  Finite bu  Two bnds   Fixed    Free
        0         0       1000       1000      1000       0       0
  [0, bu]  [bl,  0]  excluding fixed variables
     1000         0

Itn   mu stepx stepz  Pinf  Dinf  Cinf   Objective    nf  center      SQD
  0                    2.6   5.2   0.3  3.1308620e+09        2.0
  1 -0.8 0.060 0.060   2.6   5.2   0.3  2.7648605e+09  1    97.5   551775 x2
  2 -0.8 0.012 0.012   2.6   5.2   0.3  2.7016333e+09  1  4664.1
  3 -0.8 0.017 0.017   2.6   5.2   0.3  2.6095717e+09  1   767.5
  4 -0.8 0.025 0.025   2.6   5.2   0.3  2.4823822e+09  1   702.4
  5 -0.8 0.024 0.024   2.6   5.2   0.2  2.3637270e+09  1   307.0
  6 -0.8 0.023 0.023   2.5   5.1   0.2  2.2552283e+09  1   443.4
  7 -0.8 0.023 0.023   2.5   5.1   0.2  2.1515851e+09  1   388.5
  8 -0.8 0.026 0.026   2.5   5.1   0.2  2.0430638e+09  1   267.6
  9 -0.8 0.015 0.015   2.5   5.1   0.2  1.9837909e+09  1  1730.2
 10 -0.8 0.019 0.019   2.5   5.1   0.2  1.9091248e+09  1   411.8
 11 -0.8 0.018 0.018   2.5   5.1   0.2  1.8421745e+09  1   892.2
 12 -0.8 0.014 0.014   2.5   5.1   0.2  1.7924956e+09  1   352.0
 13 -0.8 0.019 0.019   2.5   5.1   0.2  1.7269175e+09  1   403.4
 14 -0.8 0.013 0.013   2.5   5.1   0.2  1.6848739e+09  1   264.1
 15 -0.8 0.008 0.008   2.5   5.1   0.2  1.6585594e+09  1  9629.8
 16 -0.8 0.014 0.014   2.5   5.1   0.2  1.6144739e+09  1   349.4
 17 -0.8 0.010 0.010   2.5   5.1   0.2  1.5848748e+09  1  4756.5
 18 -0.8 0.010 0.010   2.5   5.1   0.2  1.5540082e+09  1   525.6
 19 -0.8 0.007 0.007   2.5   5.1   0.2  1.5349144e+09  1   546.3
 20 -0.8 0.012 0.012   2.5   5.1   0.2  1.4992293e+09  1   680.3
 21 -0.8 0.007 0.007   2.5   5.1   0.1  1.4797578e+09  1   439.3
 22 -0.8 0.007 0.007   2.5   5.0   0.1  1.4601235e+09  1  1311.0
 23 -0.8 0.005 0.005   2.4   5.0   0.1  1.4465142e+09  1 15952.2
 24 -0.8 0.010 0.010   2.4   5.0   0.1  1.4213101e+09  1   573.3
 25 -0.8 0.007 0.007   2.4   5.0   0.1  1.4038884e+09  1  4270.1
 26 -0.8 0.007 0.007   2.4   5.0   0.1  1.3878290e+09  1   787.5
 27 -0.8 0.007 0.007   2.4   5.0   0.1  1.3715322e+09  1 11524.6
 28 -0.8 0.006 0.006   2.4   5.0   0.1  1.3562288e+09  1   840.8
 29 -0.8 0.004 0.004   2.4   5.0   0.1  1.3475177e+09  1  1051.7
 30 -0.8 0.006 0.006   2.4   5.0   0.1  1.3342844e+09  1  2966.0
   Too many iterations

max |x| =     2.644    max |y| =    96.789    max |z| =   394.699   scaled
max |x| =    13.219    max |y| =   774.314    max |z| =  3157.590 unscaled
PDitns  =        30    SQDitns =         0    cputime =      32.1

Distribution of vector     x         z
[  1e+05,  1e+06 )         0      1000
[  1e+04,  1e+05 )         0         0
[  1e+03,  1e+04 )         0         0
[    100,  1e+03 )         0         0
[     10,    100 )       176         0
[      1,     10 )       212         0
[    0.1,      1 )        59         0
[   0.01,    0.1 )       208         0
[  0.001,   0.01 )       294         0
[      0,  0.001 )        51         0 
Elapsed time is 26.611777 seconds.
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

