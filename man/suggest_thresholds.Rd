\name{suggest_thresholds}
\alias{suggest_thresholds}
\title{Suggest SC_min / Pop_min / Pop_max / Pop_tar from a reference partition}
\usage{
suggest_thresholds(sol, W, row_W, col_W, margin = 0.1)
}
\arguments{
\item{sol}{Integer vector, a reference partition.}
\item{W}{Numeric OD matrix (N x N).}
\item{row_W}{Numeric vector, rowSums(W).}
\item{col_W}{Numeric vector, colSums(W).}
\item{margin}{Numeric in [0, 1), how far below the reference   partition's minimum SC/population (and above its maximum   population) to set the suggested thresholds (default 0.10, i.e.   10\%).}
}
\value{
Named list: \\code\{SC_min\}, \\code\{Pop_min\}, \\code\{Pop_max\},   \\code\{Pop_tar\} (the median market population), plus \\code\{lma_df\}   (the per-market table the suggestion was derived from, so you can   see exactly which market drove each number).
}
\description{
If you already know a "good" partition for some data -- a synthetic test case with a known true answer, a previous method's output, or a business rule -- this reads off that partition's own per-market self-containment and population (via \\code\{\\link\{build_lma_df\}\}) and proposes thresholds a small margin below/above the extremes, so real markets in that partition aren't excluded by the very thresholds meant to describe them. This is a starting point, not a substitute for the kind of grid search described in the paper (Table 2) -- treat the output as something to sanity-check and refine, not use blindly.
}
