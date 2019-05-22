[![Travis build status](https://travis-ci.org/mrc-ide/drjacoby.svg?branch=master)](https://travis-ci.org/mrc-ide/drjacoby)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/mrc-ide/drjacoby?branch=master&svg=true)](https://ci.appveyor.com/project/mrc-ide/drjacoby)
[![Coverage status](https://codecov.io/gh/mrc-ide/drjacoby/branch/master/graph/badge.svg)](https://codecov.io/github/mrc-ide/drjacoby?branch=master)

# drjacoby

*drjacoby* is designed to run a very simple, but very general and flexible form of Markov chain Monte Carlo (MCMC). The likelihood and prior that go into this MCMC are defined by the user, and can be written in R or C++, with the latter giving significant performance advantages. Currently the main drawback of *drjacoby* is that all parameters must be continuous (i.e. no integer or boolean values), although this is likely to be relaxed in later versions.

The three ingredients that make up this flexible MCMC are:

1. **Reparameterisation**. User-defined parameters are free to take any range (for example [0,infinity], [0,1]). These parameters are then transformed to the [-infinity,infinity] range internally where standard Metropolis-Hastings methods take over. Adjustments for the transformation are made via the Jacobean, which ensures the true posterior distribution is the target.
2. **Metropolis Coupling**. The option exists to run multiple chains at different temperatures, with the "cold chain"" representing our ordinary MCMC chain and the "hot chains" exploring increasingly flat versions of the posterior. The different temperature rungs are linked through Metropolis coupling, which involves swapping values between rungs according to a given acceptance probability. This method is extremely powerful at improving mixing in situations where the posterior is highly correlated and/or multimodal. Note that this is separate from the concept of running multiple chains to diagnose performance, for example through the Gelman-Rubin diagnostic. This can also be done, and can be implemented in parallel over multiple cores.
3. **Tuning*. Some internal parameters of the MCMC, for example the proposal bandwidth, are tuned automatically using the Robbins-Monro method. This should lead to good acceptance rates with minimal user fiddling required.

Together, these methods lead to a flexible MCMC that should provide good results in most simple problems. In some cases this approach can even be preferable to sophisticated methods like Hamiltonian Monte Carlo (HMC), which is excellent when the posterior is highly correlated or oddly shaped, but struggles with highly multimodal distributions.

