# What is this?
This project contains scripts to reproduce experiments from the paper 
[Algorithmic Analysis and Statistical Estimation of SLOPE via Approximate Message Passing](https://arxiv.org/abs/1907.07502)
by Zhiqi Bu, Jason Klusowski, Cynthia Rush, Weijie Su.

# The Problem of Interest
SLOPE is a relatively new convex optimization procedure for high-dimensional linear regression, which includes LASSO as a special case. SLOPE penalizes the solution via the sorted L1 penalty: the larger the rank of the fitted coefficient, the larger the penalty. In this paper, we develop an iterative algorithm, known as approximate message passing (AMP), for SLOPE problem which provably converges to the true minimizer and numerical simulations show that the convergence is surprisingly fast. In addition, AMP allows us to conduct inference on SLOPE minimizer in the asymptotic manner.

# Description of Files

You need to install R-package ['SLOPE'](https://cran.r-project.org/web/packages/SLOPE/index.html) to run the following codes.

## [lambda_to_alpha.py](lambda_to_alpha.R) 

Compute state evolution and calibration between **$$\alpha$$** and **$$\lambda$$** of SLOPE-AMP.

## [AMPfaster](AMPfaster.R)

This is an example implementation of SLOPE-AMP converging much faster than other commonly known iterative algorithms including ISTA and FISTA.
