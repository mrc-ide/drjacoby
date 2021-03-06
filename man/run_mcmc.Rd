% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{run_mcmc}
\alias{run_mcmc}
\title{Run drjacoby MCMC}
\usage{
run_mcmc(
  data,
  df_params,
  misc = list(),
  loglike,
  logprior,
  burnin = 1000,
  samples = 10000,
  rungs = 1,
  chains = 5,
  coupling_on = TRUE,
  GTI_pow = 1,
  beta_manual = NULL,
  cluster = NULL,
  pb_markdown = FALSE,
  save_data = TRUE,
  silent = FALSE
)
}
\arguments{
\item{data}{a named list of numeric data values. When using C++ versions of
the likelihood and/or prior these values are treated internally as doubles,
so while integer and boolean values can be used, keep in mind that these
will be recast as doubles in the likelihood (i.e. \code{TRUE = 1.0}).}

\item{df_params}{a dataframe of parameters (see \code{?define_params}).}

\item{misc}{optional list object passed to likelihood and prior.}

\item{loglike, logprior}{the log-likelihood and log-prior functions used in
the MCMC. Can either be passed in as R functions, or as character strings
which are compiled in C++ functions.}

\item{burnin}{the number of burn-in iterations.}

\item{samples}{the number of sampling iterations.}

\item{rungs}{the number of temperature rungs used in Metropolis coupling (see
\code{coupling_on}).}

\item{chains}{the number of independent replicates of the MCMC chain to run.
If a \code{cluster} object is defined then these chains are run in
parallel, otherwise they are run in serial.}

\item{coupling_on}{whether to implement Metropolis-coupling over temperature 
rungs.}

\item{GTI_pow}{values in the temperature ladder are raised to this power.
Provides a convenient way of concentrating rungs towards one end of the
temperature scale.}

\item{beta_manual}{option to manually define temperature ladder. These values
are raised to the power \code{GTI_pow}, hence you should use \code{GTI_pow
= 1} if you want to fix powers exactly. If \code{NULL} then an equal
spacing of length \code{rungs} is used between 0 and 1.}

\item{cluster}{option to pass in a cluster environment, allowing chains to be
run in parallel (see package "parallel").}

\item{pb_markdown}{whether to run progress bars in markdown mode, meaning
they are only updated when they reach 100% to avoid large amounts of output
being printed to markdown files.}

\item{save_data}{if \code{TRUE} (the default) the raw input data is stored
for reference in the project output. This allows complete reproducibility
from a project, but may be undesirable when datasets are very large.}

\item{silent}{whether to suppress all console output.}
}
\description{
Run MCMC using defined data object, likelihood function, prior
  function and parameters.
}
