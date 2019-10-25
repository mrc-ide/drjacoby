[![Travis build status](https://travis-ci.org/mrc-ide/drjacoby.svg?branch=master)](https://travis-ci.org/mrc-ide/drjacoby)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/mrc-ide/drjacoby?branch=master&svg=true)](https://ci.appveyor.com/project/mrc-ide/drjacoby)
[![Coverage status](https://codecov.io/gh/mrc-ide/drjacoby/branch/master/graph/badge.svg)](https://codecov.io/github/mrc-ide/drjacoby?branch=master)

<br/>
<br/>
<img src="https://raw.githubusercontent.com/mrc-ide/drjacoby/master/R_ignore/images/logo2.png" height="93px" width="300px" />
<br/>

*drjacoby* is a package for running flexible Markov chain Monte Carlo (MCMC) with minimal fiddling required by the user. The likelihood and the priors that go into the model can be user-defined as either R functions or C++ functions, with the latter often leading to significant gains in speed. Outputs are produced in a standardised format, and can be explored using a range of diagnostic plots and statistics.

There are many MCMC programs out there, including more far-reaching programs such as <a href="https://www.mrc-bsu.cam.ac.uk/software/bugs/the-bugs-project-winbugs/">WinBUGS</a>, <a href="http://mcmc-jags.sourceforge.net/">JAGS</a>, <a href="https://cran.r-project.org/web/packages/greta/vignettes/get_started.html">greta</a>, <a href="https://cran.r-project.org/web/packages/rstan/vignettes/rstan.html">STAN</a>, and many more. These programs contain a wide range of options for specifying complex models, and often run different flavours of MCMC such as Hamiltonian Monte Carlo (HMC). In contrast, *drjacoby* is tailored to a specific type of MCMC problem: those that **mix poorly due to highly correlated and/or multi-modal posteriors**, which can occur in both simple and complex models. Even modern methods like HMC may fail if the posterior is particularly peaky, so we cannot simply rely on forcing more iterations through existing programs. There are techniques available for dealing with this problem, but they are fairly advanced and tricky to implement. The aim of *drjacoby* is to bring these methods within reach of a casual MCMC user, so that reliable results can be obtained without the need to code an MCMC from scratch.

Currently the main drawback of *drjacoby* is that all parameters must be continuous (i.e. no integer or boolean values), although this is likely to be relaxed in later versions.

After [installing](https://mrc-ide.github.io/drjacoby/articles/installation.html) *drjacoby*, take at look at the first [example application](https://mrc-ide.github.io/drjacoby/articles/example.html).
