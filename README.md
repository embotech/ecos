Embedded Conic Solver (ECOS)
====

[![Build Status](https://travis-ci.org/embotech/ecos.svg?branch=master)](https://travis-ci.org/embotech/ecos)
[![Windows Build](https://ci.appveyor.com/api/projects/status/github/embotech/ecos)](https://ci.appveyor.com/api/projects/status/github/embotech/ecos)
[![Coverage Status](https://coveralls.io/repos/embotech/ecos/badge.svg?branch=develop)](https://coveralls.io/r/embotech/ecos)

ECOS is a numerical software for solving convex second-order cone programs (SOCPs) of type
```
min  c'*x
s.t. A*x = b
     G*x <=_K h
```
where the last inequality is generalized, i.e. `h - G*x` belongs to the cone `K`.
ECOS supports the positive orthant `R_+`, second-order cones `Q_n` defined as
```
Q_n = { (t,x) | t >= || x ||_2 } 
```
with t a scalar and `x` in `R_{n-1}`, 
and the exponential cone `K_e` defined as 
```
K_e = closure{(x,y,z) | exp(x/z) <= y/z, z>0}
```
where `(x,y,z)` is in `R_3`. 
The cone `K` is therefore
a direct product of the positive orthant, second-order, and exponential cones:
```
K = R_+ x Q_n1 x ... x Q_nN x K_e x ... x K_e
```


Mixed-Integer SOCPs (ECOS_BB)
----
Through a recent extension by Han Wang, ECOS now comes with a branch-and-bound procedure (a direct translation of Stephen Boyd's [lecture slides](http://stanford.edu/class/ee364b/lectures/bb_slides.pdf)) called `ECOS_BB` for solving mixed-integer or mixed-boolean programs of the form

```
min  c'*x
s.t. A*x = b
     G*x <=_K h
     some x_i in {0,1}
     some x_j integer
```

Note: the branch-and-bound module has been designed to solve small problems at acceptable speeds and with minimum added code complexity (ca. 200 lines of code on top of ECOS). 

Interfaces
----

ECOS has numerous interfaces, each hosted in a separate git repository. The core ECOS solver (this repository) is included in the interface repositories as a submodule. You should run `git submodule init` and `git submodule update` after cloning the interface repositories.

* [MATLAB interface](https://github.com/embotech/ecos-matlab)
* [Python interface](https://github.com/embotech/ecos-python)
* [Julia interface](https://github.com/jump-dev/ECOS.jl)
* [R interface](https://github.com/bnaras/ECOSolveR). Also on [CRAN](https://cran.r-project.org/package=ECOSolveR).


Please refer to the corresponding repositories or the [wiki](https://github.com/embotech/ecos/wiki) for information on how to install and use ECOS through these interfaces.


License
----

ECOS is distributed under the [GNU General Public License v3.0](http://www.gnu.org/copyleft/gpl.html). Other licenses for the core solver may be available upon request from [embotech](http://www.embotech.com).


Documentation
----
The current home of the documentation is [here](https://github.com/embotech/ecos/wiki). If you find something is missing, feel free to open an issue and describe what you'd like to be documented better.


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
+ Santiago Akle (Exponential cone method and extension)

The main technical idea behind ECOS is described in a short [paper](http://www.stanford.edu/~boyd/papers/ecos.html). More details are given in Alexander Domahidi's [PhD Thesis](http://e-collection.library.ethz.ch/view/eth:7611?q=domahidi) in Chapter 9.

The algorithm for exponential cones can be found in Santiago Akle's 
[PhD Thesis](http://searchworks.stanford.edu/view/10769798) in Chapter 10.

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
