\name{run_adsa_pipeline}
\alias{run_adsa_pipeline}
\title{Run the full AdSA-ALMD pipeline in one call}
\usage{
run_adsa_pipeline(W, adj, row_W, col_W, SC_min, Pop_min, Pop_max = Inf, Pop_tar = NULL, r, L = 1000, l = 20, alpha = 0.2, beta = 0.2, eps = 0, P0 = 0.8, T0_samples = 100, sara_max_it = 5, refine_max_iter = 100, checkpoint_file = NULL, seed = NULL, verbose = TRUE)
}
\arguments{
\item{W}{Numeric OD matrix (N x N).}
\item{adj}{Integer/logical adjacency matrix (N x N).}
\item{row_W}{Numeric vector, rowSums(W).}
\item{col_W}{Numeric vector, colSums(W).}
\item{SC_min}{Numeric, minimum self-containment threshold.}
\item{Pop_min}{Numeric, minimum population threshold.}
\item{Pop_max}{Numeric, maximum population threshold (default Inf).}
\item{Pop_tar}{Numeric or NULL, target population (default NULL).}
\item{r}{Integer, penalty exponent (3 or 4).}
\item{L}{Integer, AdSA-ALMD outer trials (default 1000).}
\item{l}{Integer, AdSA-ALMD inner iterations per trial (default 20).}
\item{alpha}{Numeric, ACS increase coefficient (default 0.2).}
\item{beta}{Numeric, ACS decrease coefficient (default 0.2).}
\item{eps}{Numeric, convergence threshold (default 0).}
\item{P0}{Numeric, target acceptance probability for T0 (default 0.8).}
\item{T0_samples}{Integer, perturbations used to estimate T0 (default 100).}
\item{sara_max_it}{Integer, passed through to sara_iterated()'s max_it   (default 5).}
\item{refine_max_iter}{Integer, max passes for refine_misallocated_bgus()   (default 100).}
\item{checkpoint_file}{Character path or NULL (default). Forwarded to   adsa_almd(); leave NULL when calling this from   run_adsa_pipeline_parallel().}
\item{seed}{Integer or NULL. If supplied, set.seed(seed) is called   first so the run is reproducible.}
\item{verbose}{Logical, print step-by-step progress (default TRUE).}
}
\value{
Named list; see Details.
}
\description{
Chains the three steps you'd otherwise run by hand: build a SARA initial solution (\\code\{\\link\{sara_iterated\}\}), improve it with simulated annealing (\\code\{\\link\{adsa_almd\}\}), then patch any remaining cohesion-misallocated BGUs (\\code\{\\link\{refine_misallocated_bgus\}\}). Returns the same list structure as \\code\{adsa_almd()\} (elite, best_sol, best_fitness, gci, base_fitness, n_clusters, lma_metrics, history, op_tracker), with \\code\{best_sol\} and its derived metrics (\\code\{gci\}, \\code\{base_fitness\}, \\code\{n_clusters\}, \\code\{lma_metrics\}) updated to reflect the partition \\emph\{after\} refinement, plus the SARA \\code\{init_sol\} that seeded the run.
}
