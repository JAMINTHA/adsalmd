\name{build_results_table}
\alias{build_results_table}
\title{Comparison table matching Table 3 in the paper}
\usage{
build_results_table(results_list, W, adj, row_W, col_W, sc_tol = 0.05)
}
\arguments{
\item{results_list}{Named list: each element is a result from adsa_almd()                     (or an equivalent structure from TTWA/GEA/MSA).}
\item{W}{Numeric matrix.}
\item{adj}{Integer/logical matrix.}
\item{row_W}{Numeric vector.}
\item{col_W}{Numeric vector.}
\item{sc_tol}{Numeric, SC tolerance for misallocation.}
}
\value{
data.frame.
}
\description{
See R/08_diagnostics.R for details.
}
