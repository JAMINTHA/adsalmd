\name{adsa_almd}
\alias{adsa_almd}
\title{Run AdSA-ALMD for a single starting solution   [Algorithm 1]}
\usage{
adsa_almd(init_sol, W, adj, row_W, col_W, SC_min, Pop_min, Pop_max = Inf, Pop_tar = NULL, r, L = 1000, l = 20, alpha = 0.2, beta = 0.2, eps = 0, P0 = 0.8, T0_samples = 100, checkpoint_file = NULL, verbose = TRUE)
}
\arguments{
\item{init_sol}{Integer vector, initial solution (from SARA).}
\item{W}{Numeric matrix, OD matrix (N x N).}
\item{adj}{Integer/logical matrix, adjacency (N x N).}
\item{row_W}{Numeric vector, rowSums(W)  [precomputed].}
\item{col_W}{Numeric vector, colSums(W)  [precomputed].}
\item{SC_min}{Numeric, minimum SC threshold.}
\item{Pop_min}{Numeric, minimum population threshold.}
\item{Pop_max}{See Description.}
\item{Pop_tar}{See Description.}
\item{r}{Integer, penalty exponent (3 or 4).}
\item{L}{Integer, outer trials (default 1000).}
\item{l}{Integer, inner iterations per trial (default 20).}
\item{alpha}{Numeric, ACS increase coefficient.}
\item{beta}{Numeric, ACS decrease coefficient.}
\item{eps}{Numeric, convergence threshold (default 0).}
\item{P0}{Numeric, target acceptance prob for T0 (default 0.8).}
\item{T0_samples}{Integer, samples to estimate T0 (default 100).}
\item{checkpoint_file}{Character path or NULL (default). If supplied,   the current best_sol is written to this CSV after every outer trial.   Leave NULL (the default) for parallel runs, since concurrent workers   writing to the same path will corrupt each other's output.}
\item{verbose}{Logical, print progress every 100 outer trials.}
}
\value{
Named list:   best_sol     — integer vector, best solution found;   best_fitness — numeric, F(best_sol);   gci          — numeric, GCI of best_sol;   base_fitness — numeric, m * GCI (for cross-method comparison);   n_clusters   — integer, number of LMAs;   lma_metrics  — data.frame from build_lma_df();   history      — data.frame, trial-level diagnostics.
}
\description{
See R/07_adsa_almd.R for details.
}
