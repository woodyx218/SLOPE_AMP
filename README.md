# What is this?
This project contains scripts to reproduce experiments from the paper 
[Algorithmic Analysis and Statistical Estimation of SLOPE via Approximate Message Passing](https://arxiv.org/abs/1907.07502)
by Zhiqi Bu, Jason Klusowski, Cynthia Rush, Weijie Su.
and the paper
[Characterizing the SLOPE Trade-off: A Variational Perspective and the Donoho-Tanner Limit](https://arxiv.org/abs/2110.05679)
by Zhiqi Bu, Jason Klusowski, Cynthia Rush, Weijie Su.

# The Problem of Interest
SLOPE is a relatively new convex optimization procedure for high-dimensional linear regression, which includes LASSO as a special case. SLOPE penalizes the solution via the sorted L1 penalty: the larger the rank of the fitted coefficient, the larger the penalty. In this paper, we develop an iterative algorithm, known as approximate message passing (AMP), for SLOPE problem which provably converges to the true minimizer and numerical simulations show that the convergence is surprisingly fast. In addition, AMP allows us to conduct inference on SLOPE minimizer in the asymptotic manner.

# Description of Files

You need to install R-package ['SLOPE'](https://cran.r-project.org/web/packages/SLOPE/index.html) to run the following codes.

## [lambda_to_alpha.R](lambda_to_alpha.R) 

Compute state evolution and calibration between **\alpha** and **\lambda** of SLOPE-AMP. Also include the limiting scalar function in https://arxiv.org/abs/1903.11582

## [AMPfaster.R](AMPfaster.R)

This is an example implementation of SLOPE-AMP converging much faster than other commonly known iterative algorithms including ISTA and FISTA.

## [tradeoff_functions.R](tradeoff_functions.R) 

Compute all quantities used in SLOPE TPP-FDP trade-off. E.g. **q^\star, q_\star, zero-threshold, \epsilon^\star, u_{DT}^\star, t^\star, t_\star**......

## [upper_lower_bounds_plots.R](upper_lower_bounds_plots.R) 

Plots of **q^\star** and **q_\star**, the upper and lower bounds of the true SLOPE TPP-FDP trade-off **q**.


## Citation
```
@article{bu2020algorithmic,
  title={Algorithmic analysis and statistical estimation of SLOPE via approximate message passing},
  author={Bu, Zhiqi and Klusowski, Jason M and Rush, Cynthia and Su, Weijie J},
  journal={IEEE Transactions on Information Theory},
  volume={67},
  number={1},
  pages={506--537},
  year={2020},
  publisher={IEEE}
}

@article{bu2021characterizing,
  title={Characterizing the SLOPE Trade-off: A Variational Perspective and the Donoho-Tanner Limit},
  author={Bu, Zhiqi and Klusowski, Jason and Rush, Cynthia and Su, Weijie J},
  journal={arXiv preprint arXiv:2105.13302},
  year={2021}
}
```
