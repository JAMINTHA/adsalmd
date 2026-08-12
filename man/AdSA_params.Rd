\name{AdSA_params}
\alias{AdSA_params}
\title{Default AdSA-ALMD configuration}
\format{A named list of tunable parameters (N, SC_min, SC_tar, Pop_min, Pop_tar, Pop_max, penalty_exp, P0, T0_samples, L, l, eps, alpha, beta, n_runs, n_cores). See R/00_config.R for the current defaults and what each entry controls.}
\usage{AdSA_params}
\description{
Default parameter list used throughout the package and by main.R. Copy
and override individual entries (e.g. \code{params <- AdSA_params; params$SC_min <- 0.6})
rather than editing in place, since this object is shared.
}
