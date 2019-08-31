# What is this?
This project contains scripts to reproduce experiments from the paper 
[Algorithmic Analysis and Statistical Estimation of SLOPE via Approximate Message Passing](https://arxiv.org/abs/1907.07502)
by Zhiqi Bu, Jason Klusowski, Cynthia Rush, Weijie Su.

# The Problem of Interest
SLOPE is a relatively new convex optimization procedure for high-dimensional linear regression, which includes LASSO as a special case. SLOPE penalizes the solution via the sorted $\ell_1$ penalty: the larger the rank of the fitted coefficient, the larger the penalty. In this paper, we develop an iterative algorithm, known as approximate message passing (AMP), for SLOPE problem which provably converges to the true minimizer and numerical simulations show that the convergence is surprisingly fast. In addition, AMP allows us to conduct inference on SLOPE minimizer in the asymptotic manner.

# Description of Files

You need to install R-package ['SLOPE'](https://cran.r-project.org/web/packages/SLOPE/index.html) to run the following codes.

## [lambda_to_alpha.py](lambda_to_alpha.R) 

Creates numpy archives (.npz) and matlab (.mat) files with (y,x,A) for the sparse linear problem y=Ax+w.
These files are not really necessary for any of the deep-learning scripts, which generate the problem on demand.
They are merely provided for better understanding the specific realizations used in the experiments.

## [ista_fista_amp.m](ista_fista_amp.m)

Using the .mat files created by save_problem.py, this octave/matlab script tests the performance of non-learned algorithms ISTA, FISTA, and AMP.

e.g.
```
octave:1> ista_fista_amp
loaded Gaussian A problem
AMP reached NMSE=-35dB at iteration 25
AMP terminal NMSE=-36.7304 dB
FISTA reached NMSE=-35dB at iteration 202
FISTA terminal NMSE=-36.7415 dB
ISTA reached NMSE=-35dB at iteration 3420
ISTA terminal NMSE=-36.7419 dB
```

## [LISTA.py](LISTA.py)

This is an example implementation of LISTA _Learned Iterative Soft Thresholding Algorithm_ by (Gregor&LeCun, 2010 ICML).

## [LAMP.py](LAMP.py)

Example of Learned AMP (LAMP) with a variety of shrinkage functions.

## [LVAMP.py](LVAMP.py)
