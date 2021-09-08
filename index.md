
DrJacoby
Foo
<img src="https://raw.githubusercontent.com/mrc-ide/drjacoby/master/R_ignore/images/logo3.png" height="120px" width="100px"/>
===================================================

[![master build](https://github.com/mrc-ide/drjacoby/workflows/master_build/badge.svg)](https://github.com/mrc-ide/drjacoby/actions)
[![develop build](https://github.com/mrc-ide/drjacoby/workflows/develop_build/badge.svg)](https://github.com/mrc-ide/drjacoby/actions)
[![Coverage status](https://codecov.io/gh/mrc-ide/drjacoby/branch/master/graph/badge.svg)](https://codecov.io/github/mrc-ide/drjacoby?branch=master)


*drjacoby* is a package for running flexible Markov chain Monte Carlo (MCMC) with minimal fiddling required by the user. The likelihood and the priors that go into the model can be written as either R or C++ functions, with the latter typically being much faster to run. Outputs are produced in a standardised format, and can be explored using a range of built in diagnostic plots and statistics.

There are many MCMC programs out there, including more far-reaching programs such as <a href="https://www.mrc-bsu.cam.ac.uk/software/bugs/the-bugs-project-winbugs/">WinBUGS</a>, <a href="http://mcmc-jags.sourceforge.net/">JAGS</a>, <a href="https://cran.r-project.org/web/packages/greta/vignettes/get_started.html">greta</a>, <a href="https://cran.r-project.org/web/packages/rstan/vignettes/rstan.html">STAN</a>, and many more. These programs contain a wide variety of options for specifying complex models, and often run different flavours of MCMC such as Hamiltonian Monte Carlo (HMC). In contrast, *drjacoby* is tailored to a specific type of MCMC problem: those that **mix poorly due to highly correlated and/or multi-modal posteriors**, which can occur in both simple and complex models. If the posterior is particularly peaky then even methods like HMC may fail, meaning (unfortunately) no amount of STAN iterations will get us to the right answer. There are techniques available for dealing with this sort of problem, but they tend to be fairly advanced and tricky to implement. The aim of *drjacoby* is to bring these methods within reach of a "casual" MCMC user, so that reliable results can be obtained without the need to code a custom MCMC program from scratch.

After [installing](https://mrc-ide.github.io/drjacoby/articles/installation.html) *drjacoby*, please take at look at the first [example application](https://mrc-ide.github.io/drjacoby/articles/example.html).
