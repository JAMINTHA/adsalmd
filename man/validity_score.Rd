\name{validity_score}
\alias{validity_score}
\title{Validity score for a single market   [Equation 1]}
\usage{
validity_score(members, W, row_W, col_W, adj, SC_min, Pop_min, Pop_max = Inf)
}
\arguments{
\item{members}{Integer vector.}
\item{W}{Numeric matrix.}
\item{row_W}{Numeric vector.}
\item{col_W}{Numeric vector.}
\item{adj}{Integer/logical matrix.}
\item{SC_min}{Numeric.}
\item{Pop_min}{Numeric.}
\item{Pop_max}{See Description.}
}
\value{
Numeric in [0, 1].
}
\description{
VS(Mi) = min(SC(Mi)/SC_min, 1) * min(Pop(Mi)/Pop_min, 1) * C(Mi)
}
