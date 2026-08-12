\name{summarise_runs}
\alias{summarise_runs}
\title{Aggregate summary statistics across multiple runs}
\usage{
summarise_runs(runs_list)
}
\arguments{
\item{runs_list}{List of result objects from adsa_almd().}
}
\value{
data.frame with best, mean and SD of base_fitness.
}
\description{
See R/08_diagnostics.R for details.
}
