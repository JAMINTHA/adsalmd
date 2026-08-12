\name{run_adsa_pipeline_parallel}
\alias{run_adsa_pipeline_parallel}
\title{Run run_adsa_pipeline() N times in parallel and keep the best}
\usage{
run_adsa_pipeline_parallel(n_runs, W, adj, row_W, col_W, SC_min, Pop_min, Pop_max = Inf, Pop_tar = NULL, r, L = 1000, l = 20, alpha = 0.2, beta = 0.2, eps = 0, P0 = 0.8, T0_samples = 100, sara_max_it = 5, refine_max_iter = 100, seeds = NULL, n_cores = NULL, pkg = adsalmd)
}
\arguments{
\item{n_runs}{Integer, number of independent pipeline runs.}
\item{W}{See Description.}
\item{adj}{See Description.}
\item{row_W}{See Description.}
\item{col_W}{See Description.}
\item{SC_min}{See Description.}
\item{Pop_min}{See Description.}
\item{Pop_max}{See Description.}
\item{Pop_tar}{See Description.}
\item{r}{See Description.}
\item{L}{See Description.}
\item{l}{See Description.}
\item{alpha}{See Description.}
\item{beta}{See Description.}
\item{eps}{See Description.}
\item{P0}{See Description.}
\item{T0_samples}{See Description.}
\item{sara_max_it}{See Description.}
\item{refine_max_iter}{See Description.}
\item{seeds}{Integer vector of length n_runs, or NULL (default) to   use \\code\{100 + seq_len(n_runs)\}.}
\item{n_cores}{Integer or NULL (default) to use   \\code\{max(1, parallel::detectCores() - 1)\}, capped at n_runs.}
\item{pkg}{Character, package name to load on each worker   (default "adsalmd").}
}
\value{
Named list: \\code\{all_runs\} (list of n_runs results from   run_adsa_pipeline()), \\code\{summary\} (data.frame from   summarise_runs()), \\code\{best_run\} (the single best result).
}
\description{
Starts a PSOCK cluster, runs \\code\{\\link\{run_adsa_pipeline\}\} once per seed (each a full SARA -> AdSA-ALMD -> refinement pass), and returns every run plus the best one by \\code\{base_fitness\} (via \\code\{\\link\{summarise_runs\}\}). Each worker calls \\code\{library(pkg, character.only = TRUE)\} to get access to the package's functions, so the package must already be installed (see the README / main.R for \\code\{install.packages(".", repos = NULL, type = "source")\}).
}
