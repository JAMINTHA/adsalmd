# ============================================================
# FILE: R/11_pipeline.R
# End-to-end AdSA-ALMD pipeline: SARA -> AdSA-ALMD -> refinement.
#
# This wraps the three calls previously made by hand in main.R:
#   sara_iterated()  -> adsa_almd()  -> refine_misallocated_bgus()
# ============================================================

#' Run an expression with `adj` and `params` set as globals
#'
#' operator_1(), operator_6(), operator_8(), operator_10(), operator_11()
#' and the apply_random_operator() dispatcher read a global \code{adj}
#' matrix and a global \code{params} list ($SC_min, $Pop_min, $Pop_max)
#' instead of taking them as arguments. This helper sets both in the
#' calling process's global environment for the duration of a run
#' (restoring whatever was there before, including removing them again
#' if they didn't previously exist), so callers of run_adsa_pipeline()
#' don't have to set up globals themselves. On a parallel worker this
#' runs in that worker's own global environment, so it works the same
#' way under run_adsa_pipeline_parallel().
#'
#' @param adj     Integer/logical matrix (N x N).
#' @param SC_min  Numeric.
#' @param Pop_min Numeric.
#' @param Pop_max Numeric.
#' @param fn      Zero-argument function to run with the globals in place.
#' @return Whatever \code{fn()} returns.
.with_adsa_globals <- function(adj, SC_min, Pop_min, Pop_max, fn) {
  had_params <- exists("params", envir = .GlobalEnv, inherits = FALSE)
  old_params <- if (had_params) get("params", envir = .GlobalEnv) else NULL
  had_adj    <- exists("adj", envir = .GlobalEnv, inherits = FALSE)
  old_adj    <- if (had_adj) get("adj", envir = .GlobalEnv) else NULL

  assign("params", list(SC_min = SC_min, Pop_min = Pop_min, Pop_max = Pop_max),
         envir = .GlobalEnv)
  assign("adj", adj, envir = .GlobalEnv)

  on.exit({
    if (had_params) assign("params", old_params, envir = .GlobalEnv)
    else if (exists("params", envir = .GlobalEnv, inherits = FALSE))
      rm(list = "params", envir = .GlobalEnv)

    if (had_adj) assign("adj", old_adj, envir = .GlobalEnv)
    else if (exists("adj", envir = .GlobalEnv, inherits = FALSE))
      rm(list = "adj", envir = .GlobalEnv)
  }, add = TRUE)

  fn()
}

#' Run the full AdSA-ALMD pipeline in one call
#'
#' Chains the three steps you'd otherwise run by hand: build a SARA
#' initial solution (\code{\link{sara_iterated}}), improve it with
#' simulated annealing (\code{\link{adsa_almd}}), then patch any
#' remaining cohesion-misallocated BGUs
#' (\code{\link{refine_misallocated_bgus}}).
#'
#' Returns the same list structure as \code{adsa_almd()} (elite,
#' best_sol, best_fitness, gci, base_fitness, n_clusters, lma_metrics,
#' history, op_tracker), with \code{best_sol} and its derived metrics
#' (\code{gci}, \code{base_fitness}, \code{n_clusters}, \code{lma_metrics})
#' updated to reflect the partition \emph{after} refinement, plus the
#' SARA \code{init_sol} that seeded the run.
#'
#' @param W       Numeric OD matrix (N x N).
#' @param adj     Integer/logical adjacency matrix (N x N).
#' @param row_W   Numeric vector, rowSums(W).
#' @param col_W   Numeric vector, colSums(W).
#' @param SC_min  Numeric, minimum self-containment threshold.
#' @param Pop_min Numeric, minimum population threshold.
#' @param Pop_max Numeric, maximum population threshold (default Inf).
#' @param Pop_tar Numeric or NULL, target population (default NULL).
#' @param r       Integer, penalty exponent (3 or 4).
#' @param L       Integer, AdSA-ALMD outer trials (default 1000).
#' @param l       Integer, AdSA-ALMD inner iterations per trial (default 20).
#' @param alpha   Numeric, ACS increase coefficient (default 0.2).
#' @param beta    Numeric, ACS decrease coefficient (default 0.2).
#' @param eps     Numeric, convergence threshold (default 0).
#' @param P0      Numeric, target acceptance probability for T0 (default 0.8).
#' @param T0_samples Integer, perturbations used to estimate T0 (default 100).
#' @param sara_max_it Integer, passed through to sara_iterated()'s max_it
#'   (default 5).
#' @param refine_max_iter Integer, max passes for refine_misallocated_bgus()
#'   (default 100).
#' @param checkpoint_file Character path or NULL (default). Forwarded to
#'   adsa_almd(); leave NULL when calling this from
#'   run_adsa_pipeline_parallel().
#' @param seed    Integer or NULL. If supplied, set.seed(seed) is called
#'   first so the run is reproducible.
#' @param verbose Logical, print step-by-step progress (default TRUE).
#' @return Named list; see Details.
#' @export
run_adsa_pipeline <- function(W, adj, row_W, col_W,
                               SC_min, Pop_min, Pop_max = Inf, Pop_tar = NULL,
                               r,
                               L = 1000L, l = 20L,
                               alpha = 0.2, beta = 0.2, eps = 0,
                               P0 = 0.8, T0_samples = 100L,
                               sara_max_it = 5L,
                               refine_max_iter = 100L,
                               checkpoint_file = NULL,
                               seed = NULL,
                               verbose = TRUE) {

  if (!is.null(seed)) set.seed(seed)

  # Force every argument now. If a caller passes e.g. L = params$L, that
  # argument is an unevaluated promise pointing at `params` in the
  # caller's environment. .with_adsa_globals() below overwrites the
  # caller's global `params` (to satisfy operator_1/6/8/10/11's global
  # dependency) with a stripped-down list containing only $SC_min/$Pop_min/
  # $Pop_max -- so if L's promise is still unforced by the time it's
  # finally used, `params$L` resolves against that stripped-down list
  # and silently becomes NULL. Forcing here, before that overwrite,
  # avoids it.
  force(W); force(adj); force(row_W); force(col_W)
  force(SC_min); force(Pop_min); force(Pop_max); force(Pop_tar)
  force(r); force(L); force(l); force(alpha); force(beta); force(eps)
  force(P0); force(T0_samples); force(sara_max_it); force(refine_max_iter)
  force(checkpoint_file); force(verbose)

  .with_adsa_globals(adj, SC_min, Pop_min, Pop_max, function() {
    if (verbose) message("[1/3] SARA initial solution...")
    init_sol <- sara_iterated(
      W = W, adj = adj, row_W = row_W, col_W = col_W,
      SC_min = SC_min, Pop_min = Pop_min, Pop_max = Pop_max,
      max_it = sara_max_it, verbose = verbose
    )
    if (verbose)
      message(sprintf("      done: %d clusters", n_clusters(init_sol)))

    if (verbose) message("[2/3] AdSA-ALMD...")
    run <- adsa_almd(
      init_sol = init_sol, W = W, adj = adj, row_W = row_W, col_W = col_W,
      SC_min = SC_min, Pop_min = Pop_min, Pop_max = Pop_max, Pop_tar = Pop_tar,
      r = r, L = L, l = l, alpha = alpha, beta = beta, eps = eps,
      P0 = P0, T0_samples = T0_samples,
      checkpoint_file = checkpoint_file, verbose = verbose
    )

    if (verbose) message("[3/3] Refining misallocated BGUs...")
    refined_sol <- refine_misallocated_bgus(
      run$best_sol, W, adj, row_W, col_W, max_iter = refine_max_iter
    )

    run$init_sol     <- init_sol
    run$best_sol     <- refined_sol
    run$gci          <- compute_gci(refined_sol, W, row_W, col_W)
    run$base_fitness <- compute_base_fitness(refined_sol, W, row_W, col_W)
    run$n_clusters   <- n_clusters(refined_sol)
    run$lma_metrics  <- build_lma_df(refined_sol, W, row_W, col_W)

    if (verbose) {
      message(sprintf(
        "Pipeline complete | Base fitness = %.3f | GCI = %.3f | m = %d",
        run$base_fitness, run$gci, run$n_clusters
      ))
    }

    run
  })
}

#' Run run_adsa_pipeline() N times in parallel and keep the best
#'
#' Starts a PSOCK cluster, runs \code{\link{run_adsa_pipeline}} once per
#' seed (each a full SARA -> AdSA-ALMD -> refinement pass), and returns
#' every run plus the best one by \code{base_fitness} (via
#' \code{\link{summarise_runs}}).
#'
#' Each worker calls \code{library(pkg, character.only = TRUE)} to get
#' access to the package's functions, so the package must already be
#' installed (see the README / main.R for
#' \code{install.packages(".", repos = NULL, type = "source")}).
#'
#' @param n_runs  Integer, number of independent pipeline runs.
#' @param W,adj,row_W,col_W,SC_min,Pop_min,Pop_max,Pop_tar,r,L,l,alpha,beta,eps,P0,T0_samples,sara_max_it,refine_max_iter
#'   Passed straight through to \code{\link{run_adsa_pipeline}} on every run.
#' @param seeds   Integer vector of length n_runs, or NULL (default) to
#'   use \code{100 + seq_len(n_runs)}.
#' @param n_cores Integer or NULL (default) to use
#'   \code{max(1, parallel::detectCores() - 1)}, capped at n_runs.
#' @param pkg     Character, package name to load on each worker
#'   (default "adsalmd").
#' @return Named list: \code{all_runs} (list of n_runs results from
#'   run_adsa_pipeline()), \code{summary} (data.frame from
#'   summarise_runs()), \code{best_run} (the single best result).
#' @export
run_adsa_pipeline_parallel <- function(n_runs,
                                        W, adj, row_W, col_W,
                                        SC_min, Pop_min, Pop_max = Inf, Pop_tar = NULL,
                                        r,
                                        L = 1000L, l = 20L,
                                        alpha = 0.2, beta = 0.2, eps = 0,
                                        P0 = 0.8, T0_samples = 100L,
                                        sara_max_it = 5L,
                                        refine_max_iter = 100L,
                                        seeds = NULL,
                                        n_cores = NULL,
                                        pkg = "adsalmd") {

  if (is.null(seeds)) seeds <- 100L + seq_len(n_runs)
  stopifnot(length(seeds) == n_runs)

  if (is.null(n_cores)) n_cores <- max(1L, parallel::detectCores() - 1L)
  n_cores <- max(1L, min(n_cores, n_runs))

  # Force all arguments now, in this frame. Left as lazy promises, they'd
  # point back to expressions in the caller's environment (e.g. `adj_test`)
  # which don't exist once the closure below is serialized to a worker.
  force(W); force(adj); force(row_W); force(col_W)
  force(SC_min); force(Pop_min); force(Pop_max); force(Pop_tar)
  force(r); force(L); force(l); force(alpha); force(beta); force(eps)
  force(P0); force(T0_samples); force(sara_max_it); force(refine_max_iter)

  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  parallel::clusterExport(cl, "pkg", envir = environment())
  parallel::clusterEvalQ(cl, library(pkg, character.only = TRUE))

  all_runs <- parallel::parLapply(cl, seeds, function(seed) {
    run_adsa_pipeline(
      W = W, adj = adj, row_W = row_W, col_W = col_W,
      SC_min = SC_min, Pop_min = Pop_min, Pop_max = Pop_max, Pop_tar = Pop_tar,
      r = r, L = L, l = l, alpha = alpha, beta = beta, eps = eps,
      P0 = P0, T0_samples = T0_samples,
      sara_max_it = sara_max_it, refine_max_iter = refine_max_iter,
      checkpoint_file = NULL, seed = seed, verbose = FALSE
    )
  })

  summary_df <- summarise_runs(all_runs)
  best_run   <- all_runs[[summary_df$best_run]]

  list(
    all_runs = all_runs,
    summary  = summary_df,
    best_run = best_run
  )
}
