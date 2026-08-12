\name{compute_scss}
\alias{compute_scss}
\title{Supply-side self-containment for a market   [Equation 10]}
\usage{
compute_scss(members, W, row_W)
}
\arguments{
\item{members}{Integer vector, BGU indices in market.}
\item{W}{Numeric matrix.}
\item{row_W}{Numeric vector.}
}
\value{
Numeric in [0, 1].
}
\description{
SCSS(Ms) = W(Ms, Ms) / W(Ms, G)
}
