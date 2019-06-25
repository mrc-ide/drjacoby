[![Travis build status](https://travis-ci.org/mrc-ide/drjacoby.svg?branch=master)](https://travis-ci.org/mrc-ide/drjacoby)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/mrc-ide/drjacoby?branch=master&svg=true)](https://ci.appveyor.com/project/mrc-ide/drjacoby)
[![Coverage status](https://codecov.io/gh/mrc-ide/drjacoby/branch/master/graph/badge.svg)](https://codecov.io/github/mrc-ide/drjacoby?branch=master)

<br/>
<br/>
<img src="https://raw.githubusercontent.com/mrc-ide/drjacoby/master/R_ignore/images/logo2.png" height="93px" width="300px" />
<br/>

The *drjacoby* package is designed to run a very simple, but very general form of Markov chain Monte Carlo (MCMC). The three ingredients that make up this flexible MCMC are:

1. **Reparameterisation**. User-defined parameters are free to take any range (for example [0,infinity], [0,1]). These parameters are then transformed to the [-infinity,infinity] range internally where standard Metropolis-Hastings methods take over. Adjustments for the transformation are made via the Jacobean, which ensures that the true posterior distribution is explored.
2. **Metropolis Coupling**. The option exists to run multiple chains at different temperatures, with the "cold chain"" representing our ordinary MCMC chain and the "hot chains" exploring increasingly flat versions of the posterior. The different temperature rungs are linked through Metropolis coupling, which involves swapping values between rungs according to a given acceptance probability. This method is extremely powerful at improving mixing in situations where the posterior is highly correlated and/or multimodal.
3. **Dynamic Tuning**. Some internal parameters of the MCMC are tuned automatically, for example the proposal bandwidths are updated using the Robbins-Monro method. This should lead to good acceptance rates with minimal user fiddling required.

Together, these methods lead to a flexible MCMC that should provide good results in most simple problems. In some cases this approach can even be preferable to sophisticated methods like Hamiltonian Monte Carlo (HMC), which is excellent when the posterior is highly correlated or oddly shaped, but struggles with highly multimodal distributions.

The likelihood and prior that go into *drjacoby* are defined by the user, and can be written in R or C++, with the latter leading to significant performance advantages. Currently the main drawback of *drjacoby* is that all parameters must be continuous (i.e. no integer or boolean values), although this is likely to be relaxed in later versions.

After [installing](https://mrc-ide.github.io/drjacoby/articles/installation.html) *drjacoby*, take at look at the first [example application](https://mrc-ide.github.io/drjacoby/articles/example.html).
