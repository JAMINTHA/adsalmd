\name{compute_gci_fast}
\alias{compute_gci_fast}
\title{Fast Global Cohesion Index for large N   [Equation 7]}
\usage{
compute_gci_fast(sol, W, row_W, col_W)
}
\arguments{
\item{sol}{Integer vector of assignments (length N).}
\item{W}{Numeric matrix (N x N).}
\item{row_W}{Numeric vector (length N).}
\item{col_W}{Numeric vector (length N).}
}
\value{
Numeric scalar GCI.
}
\description{
Replaces compute_gci() for N >= 100.  Instead of indexing W[g, M_g] once per BGU (520 separate row-extractions for N=520), this extracts W[Ms, Ms] ONCE per cluster and derives all per-BGU quantities from the sub-matrix row/column sums. Speedup on N=520:  ~15-20x over the naive version.
}
