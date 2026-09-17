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
#' @noRd
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
#' The main entry point. Chains the three stages you would otherwise run by
#' hand — SARA construction, simulated annealing, misallocation refinement —
#' and handles the global-variable setup the operators require. Start here.
#'
#' @details
#' \enumerate{
#'   \item \code{\link{sara_iterated}()} builds a feasible initial partition,
#'     aggregating BGUs until every market is valid and breaking up any market
#'     over \code{Pop_max}.
#'   \item \code{\link{adsa_almd}()} anneals it, perturbing with the eleven
#'     operators and maximising \code{\link{compute_objective}()}.
#'   \item \code{\link{refine_misallocated_bgus}()} moves any remaining
#'     cohesion-misallocated BGU to its best-fitting neighbour.
#' }
#' The stage-2 result is then updated in place: \code{best_sol} becomes the
#' refined partition and \code{gci}, \code{base_fitness}, \code{n_clusters} and
#' \code{lma_metrics} are recomputed from it. \code{best_fitness},
#' \code{history}, \code{elite} and \code{op_tracker} still describe the
#' \emph{pre-refinement} search, so \code{best_fitness} will not in general
#' equal \code{compute_objective(best_sol, ...)}.
#'
#' Stage 3 does not check population bounds (see
#' \code{\link{refine_misallocated_bgus}()}), so re-check
#' \code{lma_metrics} against \code{Pop_min} and \code{Pop_max} before treating
#' the returned partition as feasible.
#'
#' @section Globals and lazy evaluation:
#' Several operators read a global \code{params} list rather than taking
#' thresholds as arguments. This function installs one for the duration of the
#' run and restores whatever was there before — including removing it again if
#' it did not previously exist — so callers need no setup. On a parallel worker
#' this happens in that worker's own global environment, so
#' \code{\link{run_adsa_pipeline_parallel}()} behaves identically.
#'
#' Because that installed \code{params} shadows the caller's, every argument is
#' forced on entry. Without this, calling
#' \code{run_adsa_pipeline(..., L = params$L)} would leave \code{L} as an
#' unevaluated promise that later resolves against the \emph{stripped-down}
#' \code{params} — silently becoming \code{NULL}. Forcing first avoids it, but
#' it is still clearer to pass literal values or a differently-named list.
#'
#' @section Reproducibility:
#' Pass \code{seed} to make a run repeatable; it calls \code{set.seed()} before
#' anything else. Every stage draws from the same stream, so the same seed with
#' the same arguments reproduces the run exactly. Different seeds give
#' materially different partitions — that is the method working as designed,
#' which is why \code{\link{run_adsa_pipeline_parallel}()} exists.
#'
#' @section Runtime:
#' Dominated by stage 2: roughly \code{L * l} objective evaluations plus
#' \code{L} temperature estimates of \code{T0_samples} perturbations each, plus
#' one \code{\link{fix_misallocations}()} sweep per trial. Start with a small
#' \code{L} (say 20) to size a run on your data before committing to the
#' default of 1000.
#'
#' @param W Numeric \code{N x N} origin-destination matrix.
#' @param adj Integer/logical \code{N x N} adjacency matrix.
#' @param row_W Numeric vector, \code{rowSums(W)}.
#' @param col_W Numeric vector, \code{colSums(W)}.
#' @param SC_min Numeric in (0, 1], the minimum self-containment threshold.
#' @param Pop_min Numeric, the minimum market population.
#' @param Pop_max Numeric, the maximum market population (default \code{Inf}).
#' @param Pop_tar Numeric or \code{NULL}, target population (default
#'   \code{NULL}). Accepted and forwarded; not used.
#' @param r Integer, penalty steepness exponent (3 or 4).
#' @param L Integer, annealing outer trials (default 1000).
#' @param l Integer, inner iterations per trial (default 20).
#' @param alpha Numeric, cooling-schedule increase coefficient (default 0.2).
#' @param beta Numeric, cooling-schedule decrease coefficient (default 0.2).
#' @param eps Numeric, convergence threshold (default 0).
#' @param P0 Numeric, target acceptance probability for \code{T0} (default
#'   0.8).
#' @param T0_samples Integer, perturbations per temperature estimate (default
#'   100).
#' @param sara_max_it Integer, forwarded to \code{\link{sara_iterated}()}'s
#'   \code{max_it} (default 5).
#' @param refine_max_iter Integer, forwarded to
#'   \code{\link{refine_misallocated_bgus}()}'s \code{max_iter} (default 100).
#' @param checkpoint_file Character path or \code{NULL} (default). Forwarded to
#'   \code{\link{adsa_almd}()}; leave \code{NULL} under
#'   \code{\link{run_adsa_pipeline_parallel}()}, since concurrent workers would
#'   overwrite each other's file.
#' @param seed Integer or \code{NULL} (default). If supplied,
#'   \code{set.seed(seed)} is called first.
#' @param verbose Logical (default \code{TRUE}); print stage-by-stage progress.
#'   Note \code{\link{adsa_almd}()} emits per-trial messages regardless.
#' @return Named list, the \code{\link{adsa_almd}()} result with two changes:
#'   \describe{
#'     \item{\code{init_sol}}{\emph{added} — the SARA partition that seeded the
#'       run, so you can see what annealing contributed.}
#'     \item{\code{best_sol}}{the partition \emph{after} refinement, with
#'       \code{gci}, \code{base_fitness}, \code{n_clusters} and
#'       \code{lma_metrics} recomputed to match.}
#'     \item{\code{best_fitness}, \code{history}, \code{elite},
#'       \code{op_tracker}}{unchanged from the search; they describe the
#'       pre-refinement partition.}
#'   }
#'
#' @family pipeline
#' @seealso \code{\link{run_adsa_pipeline_parallel}()} for multiple seeds;
#'   \code{\link{sara_iterated}()}, \code{\link{adsa_almd}()} and
#'   \code{\link{refine_misallocated_bgus}()} for the individual stages;
#'   \code{\link{suggest_thresholds}()} for choosing the thresholds;
#'   \code{\link{partitions_match}()} to check a result against a known answer.
#'
#' @examples
#' ## Toy system: six BGUs in a line, forming two 3-BGU markets
#' adj <- matrix(0L, 6, 6); adj[cbind(1:5, 2:6)] <- 1L; adj <- adj + t(adj)
#' W <- matrix(1, 6, 6); W[1:3, 1:3] <- 30; W[4:6, 4:6] <- 30; diag(W) <- 60
#' row_W <- rowSums(W); col_W <- colSums(W)
#'
#' \donttest{
#' res <- suppressMessages(run_adsa_pipeline(
#'   W = W, adj = adj, row_W = row_W, col_W = col_W,
#'   SC_min = 0.55, Pop_min = 200, Pop_max = 500,
#'   r = 3, L = 5L, l = 10L, T0_samples = 20L,
#'   seed = 2026, verbose = FALSE
#' ))
#'
#' res$init_sol      # what SARA produced
#' res$best_sol      # after annealing and refinement
#' res$lma_metrics
#'
#' ## Did it recover the true two-market answer?
#' partitions_match(c(1, 1, 1, 2, 2, 2), res$best_sol)$match
#'
#' ## Re-check feasibility: stage 3 does not
#' with(res$lma_metrics, all(pop >= 200 & pop <= 500))
#'
#' ## The same seed reproduces the run exactly
#' res2 <- suppressMessages(run_adsa_pipeline(
#'   W = W, adj = adj, row_W = row_W, col_W = col_W,
#'   SC_min = 0.55, Pop_min = 200, Pop_max = 500,
#'   r = 3, L = 5L, l = 10L, T0_samples = 20L,
#'   seed = 2026, verbose = FALSE
#' ))
#' identical(res$best_sol, res2$best_sol)
#' }
#'
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

#' Run the pipeline from several seeds in parallel and keep the best
#'
#' Runs \code{\link{run_adsa_pipeline}()} once per seed across a PSOCK cluster
#' and returns every result together with the best one. Because the method is
#' stochastic, this — not a single long run — is how you should produce a final
#' delineation.
#'
#' @details
#' Each worker performs a complete, independent SARA → annealing → refinement
#' pass from its own seed. Runs are compared on \code{base_fitness} (the
#' unpenalised cross-method measure \code{m * GCI}) via
#' \code{\link{summarise_runs}()}, and the winner is returned alongside the
#' full set.
#'
#' Keep the individual runs. The spread across seeds is itself a result: a
#' small \code{sd_fitness} says the search is finding the same structure
#' repeatedly, a large one says the landscape has several competing optima and
#' the best single run may be luck. Compare partitions with
#' \code{\link{partitions_match}()} rather than by fitness alone.
#'
#' More runs than cores is fine — \code{parLapply()} queues them.
#'
#' @section The package must be installed:
#' Each worker calls \code{library(pkg, character.only = TRUE)}, so this will
#' not work against a package merely loaded with \code{devtools::load_all()} or
#' sourced from \code{R/}. Install first:
#' \code{install.packages(".", repos = NULL, type = "source")}.
#'
#' @section Why every argument is forced:
#' Arguments left as lazy promises would point at expressions in the caller's
#' environment (\code{adj_test}, say) that do not exist once the worker closure
#' is serialised, producing an obscure "object not found" from inside a worker.
#' They are therefore all forced in this frame before the cluster starts.
#'
#' Note that \code{W} and \code{adj} are copied to every worker, so peak memory
#' is roughly \code{n_cores} times the size of those matrices.
#'
#' @section Checkpointing is disabled:
#' \code{checkpoint_file} is hard-coded to \code{NULL} for every run, since
#' concurrent workers writing to one path would corrupt each other's output.
#' Use \code{\link{run_adsa_pipeline}()} directly if you need checkpoints.
#'
#' @param n_runs Integer, number of independent pipeline runs.
#' @param W,adj,row_W,col_W,SC_min,Pop_min,Pop_max,Pop_tar,r,L,l,alpha,beta,eps,P0,T0_samples,sara_max_it,refine_max_iter
#'   Passed straight through to \code{\link{run_adsa_pipeline}()} on every run.
#'   See that function for their meanings and defaults.
#' @param seeds Integer vector of length \code{n_runs}, or \code{NULL}
#'   (default) to use \code{100 + seq_len(n_runs)}. Supply your own to
#'   reproduce a previous set of runs.
#' @param n_cores Integer, or \code{NULL} (default) to use
#'   \code{max(1, parallel::detectCores() - 1)}. Capped at \code{n_runs}.
#' @param pkg Character, the package name each worker loads (default
#'   \code{"adsalmd"}).
#' @return Named list:
#'   \describe{
#'     \item{\code{all_runs}}{list of \code{n_runs} results from
#'       \code{\link{run_adsa_pipeline}()}, in seed order.}
#'     \item{\code{summary}}{one-row \code{data.frame} from
#'       \code{\link{summarise_runs}()}: best, mean and SD of base fitness, and
#'       the winning index.}
#'     \item{\code{best_run}}{the single best result, i.e.
#'       \code{all_runs[[summary$best_run]]}.}
#'   }
#'
#' @family pipeline
#' @seealso \code{\link{run_adsa_pipeline}()} for a single run and the full
#'   argument documentation; \code{\link{summarise_runs}()} for the comparison;
#'   \code{\link{build_results_table}()} to tabulate several methods.
#'
#' @examples
#' \dontrun{
#' ## Requires the package to be installed, not just loaded
#' ## install.packages(".", repos = NULL, type = "source")
#'
#' adj <- matrix(0L, 6, 6); adj[cbind(1:5, 2:6)] <- 1L; adj <- adj + t(adj)
#' W <- matrix(1, 6, 6); W[1:3, 1:3] <- 30; W[4:6, 4:6] <- 30; diag(W) <- 60
#' row_W <- rowSums(W); col_W <- colSums(W)
#'
#' out <- run_adsa_pipeline_parallel(
#'   n_runs = 4,
#'   W = W, adj = adj, row_W = row_W, col_W = col_W,
#'   SC_min = 0.55, Pop_min = 200, Pop_max = 500,
#'   r = 3, L = 10L, l = 10L, T0_samples = 20L,
#'   n_cores = 2
#' )
#'
#' out$summary
#' out$best_run$best_sol
#'
#' ## How much do the runs actually disagree?
#' sols <- lapply(out$all_runs, `[[`, "best_sol")
#' vapply(sols, function(s) partitions_match(sols[[1]], s)$match, logical(1))
#' }
#'
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
