# CoupledCPF
... which stands for Coupled Conditional Particle Filters

This package accompanies the arXiv report https://arxiv.org/abs/1701.02002
"Smoothing with Couplings of Conditional Particle Filters"
by Pierre E. Jacob, Fredrik Lindsten, Thomas B. Schön

Functions are provided to construct Rhee--Glynn estimators of smoothing functionals.
For comparison, particle filters with fixed-lag and Kalman smoothers are also implemented.
Examples include toy auto-regressive models, the classic nonlinear state space model from Gordon, Salmond & Smith 1993,  and a prey-predator model with an intractable transition density.

Abstract of the arXiv report:

In state space models, smoothing refers to the task of estimating a latent
stochastic process given noisy measurements related to the process. We propose an 
unbiased estimator of smoothing expectations. The lack-of-bias property has
methodological benefits: independent estimators can be generated in parallel,
and confidence intervals can be constructed from the central limit theorem to quantify the approximation error.
To design unbiased estimators, we combine a generic debiasing technique for Markov
chains with a Markov chain Monte Carlo algorithm for smoothing.  The resulting procedure is
widely applicable and we show in numerical experiments that the removal of the
bias comes at a manageable increase in variance.  We establish the validity of
the proposed estimators under mild assumptions. Numerical experiments are
provided on toy models, including a setting of highly-informative observations,
and a realistic Lotka-Volterra model with an intractable transition
density.
