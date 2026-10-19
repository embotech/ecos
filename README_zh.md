嵌入式锥求解器（ECOS）
====

[![Build Status](https://travis-ci.org/embotech/ecos.svg?branch=master)](https://travis-ci.org/embotech/ecos)
[![Windows Build](https://ci.appveyor.com/api/projects/status/github/embotech/ecos)](https://ci.appveyor.com/api/projects/status/github/embotech/ecos)
[![Coverage Status](https://coveralls.io/repos/embotech/ecos/badge.svg?branch=develop)](https://coveralls.io/r/embotech/ecos)

[English](README.md) | 简体中文

ECOS 是一个用于求解凸二阶锥规划（SOCP）的数值软件，问题形式如下：
```
min  c'*x
s.t. A*x = b
     G*x <=_K h
```
其中最后一个不等式是广义不等式，即 `h - G*x` 属于锥 `K`。
ECOS 支持正锥（非负象限）`R_+`、二阶锥 `Q_n`，其定义为
```
Q_n = { (t,x) | t >= || x ||_2 }
```
其中 t 为标量，`x` 属于 `R_{n-1}`，
以及指数锥 `K_e`，其定义为
```
K_e = closure{(x,y,z) | exp(x/z) <= y/z, z>0}
```
其中 `(x,y,z)` 属于 `R_3`。
因此，锥 `K` 是正锥、若干二阶锥和若干指数锥的直积：
```
K = R_+ x Q_n1 x ... x Q_nN x K_e x ... x K_e
```


混合整数 SOCP（ECOS_BB）
----

得益于 Han Wang 的扩展，ECOS 现在包含一个分支定界算法（Stephen Boyd 的[课程讲义](http://stanford.edu/class/ee364b/lectures/bb_slides.pdf)的直接实现），称为 `ECOS_BB`，用于求解如下形式的混合整数或混合布尔规划：

```
min  c'*x
s.t. A*x = b
     G*x <=_K h
     some x_i in {0,1}
     some x_j integer
```

注意：分支定界模块的设计目标是以可接受的速度求解小规模问题，并尽量少地增加代码复杂度（在 ECOS 之上仅增加约 200 行代码）。

接口
----

ECOS 拥有众多接口，每个接口托管在独立的 git 仓库中。核心 ECOS 求解器（即本仓库）以子模块的形式包含在各接口仓库中。克隆接口仓库后，请运行 `git submodule init` 和 `git submodule update`。

* [MATLAB 接口](https://github.com/embotech/ecos-matlab)
* [Python 接口](https://github.com/embotech/ecos-python)
* [Julia 接口](https://github.com/jump-dev/ECOS.jl)
* [R 接口](https://github.com/bnaras/ECOSolveR)，也可在 [CRAN](https://cran.r-project.org/package=ECOSolveR) 上获取。


有关如何通过这些接口安装和使用 ECOS 的信息，请参阅相应的仓库或 [wiki](https://github.com/embotech/ecos/wiki)。


许可证
----

ECOS 基于 [GNU General Public License v3.0](http://www.gnu.org/copyleft/gpl.html) 发布。核心求解器的其他许可证可向 [embotech](http://www.embotech.com) 申请获取。


文档
----

文档目前位于[这里](https://github.com/embotech/ecos/wiki)。如果你发现缺少某些内容，欢迎提交 issue 并说明你希望补充哪些文档。


致谢
----

该求解器本质上基于 Lieven Vandenberghe 的 [CVXOPT](http://cvxopt.org) [ConeLP](http://www.ee.ucla.edu/~vandenbe/publications/coneprog.pdf) 求解器，但在线性方程组的处理方式上有所不同。

以下人员曾经并仍在参与 ECOS 的开发与维护：

+ Alexander Domahidi（主要开发者）
+ Eric Chu（Python 接口、单元测试）
+ Stephen Boyd（方法与数学）
+ Michael Grant（CVX 接口）
+ Johan Löfberg（YALMIP 接口）
+ João Felipe Santos、Iain Dunning（Julia 接口）
+ Han Wang（ECOS 分支定界封装）
+ Santiago Akle（指数锥方法与扩展）

ECOS 背后的主要技术思想在一篇简短[论文](http://www.stanford.edu/~boyd/papers/ecos.html)中有描述。更多细节见 Alexander Domahidi 的[博士论文](http://e-collection.library.ethz.ch/view/eth:7611?q=domahidi)第 9 章。

指数锥的算法见 Santiago Akle 的[博士论文](http://searchworks.stanford.edu/view/10769798)第 10 章。

如果你觉得 ECOS 有用，可以使用以下 BibTeX 条目引用它：

```
@INPROCEEDINGS{bib:Domahidi2013ecos,
author={Domahidi, A. and Chu, E. and Boyd, S.},
booktitle={European Control Conference (ECC)},
title={{ECOS}: {A}n {SOCP} solver for embedded systems},
year={2013},
pages={3071-3076}
}
```
