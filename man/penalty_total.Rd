\name{penalty_total}
\alias{penalty_total}
\title{Total penalty P(x) = PSC(x) + PPop(x)   [Equation 9]}
\usage{
penalty_total(lma_df, SC_min, Pop_min, Pop_max = Inf, Pop_tar = NULL, r)
}
\arguments{
\item{lma_df}{data.frame from build_lma_df().}
\item{SC_min}{Numeric.}
\item{Pop_min}{Numeric.}
\item{Pop_max}{See Description.}
\item{Pop_tar}{See Description.}
\item{r}{Integer.}
}
\value{
Numeric scalar.
}
\description{
See R/03_objective.R for details.
}
