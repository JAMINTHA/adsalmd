# ============================================================
# FILE: R/08_diagnostics.R
# Misallocation analysis and results summary (Table 3 in paper).
# ============================================================

#' Is a single BGU cohesion-misallocated?
#'
#' Tests whether a BGU would be better placed in one of its neighbouring
#' markets. This is the quality criterion the paper reports in Table 3: a good
#' delineation should leave few or no BGUs that obviously belong somewhere
#' else, and the count is comparable across methods with different objectives.
#'
#' @details
#' BGU \code{g} is misallocated if there is an adjacent market \code{cl} such
#' that \strong{both} hold:
#' \enumerate{
#'   \item \code{\link{compute_ci}(g, members(cl), ...)} exceeds \code{g}'s
#'     cohesion in its own market (computed against the market \emph{excluding}
#'     \code{g}, so the comparison is like for like);
#'   \item moving \code{g} there would reduce that destination market's
#'     self-containment by no more than \code{sc_tol}.
#' }
#' The function short-circuits on the first market satisfying both, so it tells
#' you \emph{that} a better home exists, not which one — use
#' \code{\link{refine_misallocated_bgus}()} or
#' \code{\link{fix_misallocations}()} to actually move it.
#'
#' The second condition is what stops the test from flagging every BGU on a
#' border. Cohesion alone would happily recommend moves that wreck the
#' receiving market; requiring the destination to stay nearly as self-contained
#' keeps the criterion honest. Note that only the \emph{destination}'s
#' self-containment is checked — the origin market's loss is not.
#'
#' A BGU with no adjacent market other than its own returns \code{FALSE}
#' immediately.
#'
#' @section Choosing sc_tol:
#' \code{sc_tol = 0} means no self-containment loss is acceptable and flags
#' almost nothing; large values approach a pure cohesion test and flag a great
#' deal. The default of 0.05 implements the paper's "without substantial loss
#' in self-containment". Report the tolerance alongside any misallocation count,
#' since the number is meaningless without it.
#'
#' @param g Integer, the BGU index to test.
#' @param sol Integer vector of market assignments (length N).
#' @param W Numeric \code{N x N} origin-destination matrix.
#' @param adj Integer/logical \code{N x N} adjacency matrix.
#' @param row_W Numeric vector, \code{rowSums(W)}.
#' @param col_W Numeric vector, \code{colSums(W)}.
#' @param sc_tol Numeric \eqn{\ge 0}, the maximum tolerated drop in the
#'   destination market's self-containment (default 0.05).
#' @return \code{TRUE} if a better-fitting adjacent market exists,
#'   \code{FALSE} otherwise.
#'
#' @family diagnostics
#' @seealso \code{\link{count_misallocated}()} for the whole-partition count;
#'   \code{\link{refine_misallocated_bgus}()} and
#'   \code{\link{fix_misallocations}()} for the repairs.
#'
#' @examples
#' adj <- matrix(0L, 6, 6); adj[cbind(1:5, 2:6)] <- 1L; adj <- adj + t(adj)
#' W <- matrix(1, 6, 6); W[1:3, 1:3] <- 30; W[4:6, 4:6] <- 30; diag(W) <- 60
#' row_W <- rowSums(W); col_W <- colSums(W)
#'
#' ## BGU 4 is in market 1 but commutes with 5 and 6
#' sol <- c(1, 1, 1, 1, 2, 2)
#' is_bgu_misallocated(4, sol, W, adj, row_W, col_W)
#' is_bgu_misallocated(1, sol, W, adj, row_W, col_W)
#'
#' ## In the correct partition, nothing is misallocated
#' vapply(1:6, is_bgu_misallocated, logical(1),
#'        sol = c(1, 1, 1, 2, 2, 2), W = W, adj = adj,
#'        row_W = row_W, col_W = col_W)
#'
#' ## A zero tolerance flags far less
#' is_bgu_misallocated(4, sol, W, adj, row_W, col_W, sc_tol = 0)
#'
#' @export
is_bgu_misallocated <- function(g, sol, W, adj, row_W, col_W, sc_tol = 0.05) {
  current_cl  <- sol[g]
  current_mbr <- setdiff(cluster_members(current_cl, sol), g)
  current_ci  <- compute_ci(g, current_mbr, W, row_W, col_W)
  current_sc  <- compute_sc(cluster_members(current_cl, sol),
                             W, row_W, col_W)$sc

  adj_cls <- setdiff(unique(sol[which(adj[g, ] > 0L)]), current_cl)
  if (length(adj_cls) == 0L) return(FALSE)

  for (cl in adj_cls) {
    nbr_mbrs <- cluster_members(cl, sol)
    new_ci   <- compute_ci(g, nbr_mbrs, W, row_W, col_W)

    if (new_ci > current_ci) {
      # Would SC of the destination market drop substantially?
      candidate     <- sol
      candidate[g]  <- cl
      new_sc_dest   <- compute_sc(cluster_members(cl, candidate),
                                   W, row_W, col_W)$sc
      orig_sc_dest  <- compute_sc(nbr_mbrs, W, row_W, col_W)$sc

      if ((orig_sc_dest - new_sc_dest) <= sc_tol) return(TRUE)
    }
  }
  FALSE
}

#' Count cohesion-misallocated BGUs in a partition
#'
#' Applies \code{\link{is_bgu_misallocated}()} to every BGU and returns the
#' total. This is the \code{Misallocated_BGU} column of the paper's Table 3 and
#' the headline quality measure that does not depend on any method's own
#' objective — which is what makes it usable for comparing AdSA-ALMD against
#' TTWA, GEA and MSA.
#'
#' @details
#' Read it as a proportion of \code{N}, not as an absolute: 12 misallocated out
#' of 520 BGUs is a good result, out of 20 it is not. Zero is achievable on
#' clean synthetic data but is not the target on real commuting data, where
#' genuinely ambiguous border BGUs exist.
#'
#' Cost is \code{N} calls to \code{\link{is_bgu_misallocated}()}, each of which
#' scores several markets — noticeable on a large map, so compute it once at
#' the end of a run rather than inside a loop.
#'
#' @param sol Integer vector of market assignments (length N).
#' @param W Numeric \code{N x N} origin-destination matrix.
#' @param adj Integer/logical \code{N x N} adjacency matrix.
#' @param row_W Numeric vector, \code{rowSums(W)}.
#' @param col_W Numeric vector, \code{colSums(W)}.
#' @param sc_tol Numeric, tolerated self-containment reduction (default 0.05).
#'   Always report this alongside the count.
#' @return Integer, the number of misallocated BGUs, between 0 and
#'   \code{length(sol)}.
#'
#' @family diagnostics
#' @seealso \code{\link{is_bgu_misallocated}()} for the per-BGU test,
#'   \code{\link{build_results_table}()} which reports this per method,
#'   \code{\link{refine_misallocated_bgus}()} to reduce it.
#'
#' @examples
#' adj <- matrix(0L, 6, 6); adj[cbind(1:5, 2:6)] <- 1L; adj <- adj + t(adj)
#' W <- matrix(1, 6, 6); W[1:3, 1:3] <- 30; W[4:6, 4:6] <- 30; diag(W) <- 60
#' row_W <- rowSums(W); col_W <- colSums(W)
#'
#' count_misallocated(c(1, 1, 1, 2, 2, 2), W, adj, row_W, col_W)  # correct
#' count_misallocated(c(1, 1, 1, 1, 2, 2), W, adj, row_W, col_W)  # BGU 4 adrift
#'
#' ## Sensitivity to the tolerance
#' sol <- c(1, 1, 1, 1, 2, 2)
#' vapply(c(0, 0.05, 0.2, 1), function(tol)
#'          count_misallocated(sol, W, adj, row_W, col_W, sc_tol = tol),
#'        numeric(1))
#'
#' @export
count_misallocated <- function(sol, W, adj, row_W, col_W, sc_tol = 0.05) {
  sum(vapply(seq_along(sol), is_bgu_misallocated,
             logical(1L),
             sol    = sol,
             W      = W,
             adj    = adj,
             row_W  = row_W,
             col_W  = col_W,
             sc_tol = sc_tol))
}

#' Summarise base fitness across independent runs
#'
#' Reduces a list of run results to one row: the best, mean and standard
#' deviation of their base fitness, plus which run was best. This is how the
#' stochastic method's reliability is reported — a single run's score says
#' little without knowing the spread.
#'
#' @details
#' Reads \code{base_fitness} (not \code{best_fitness}) from each element, so
#' the comparison is on the unpenalised cross-method measure \code{m * GCI};
#' see \code{\link{compute_base_fitness}()}. Any object with that component
#' works — results from \code{\link{adsa_almd}()},
#' \code{\link{run_adsa_pipeline}()}, or a hand-built list representing a
#' competing method.
#'
#' \code{best_run} is an index into \code{runs_list}, which is how
#' \code{\link{run_adsa_pipeline_parallel}()} recovers the winning run. A tie
#' resolves to the lowest index.
#'
#' A single-element list gives \code{sd_fitness = NA}, since
#' \code{\link[stats]{sd}()} of one value is undefined.
#'
#' @param runs_list List of run results, each with a \code{base_fitness}
#'   component.
#' @return A one-row \code{data.frame}:
#'   \describe{
#'     \item{\code{best_fitness}}{numeric, the maximum base fitness. Named for
#'       consistency with the run objects; it is a base fitness.}
#'     \item{\code{mean_fitness}}{numeric, the mean across runs.}
#'     \item{\code{sd_fitness}}{numeric, the standard deviation; \code{NA} for
#'       a single run.}
#'     \item{\code{best_run}}{integer, the index of the best run.}
#'   }
#'
#' @family diagnostics
#' @seealso \code{\link{run_adsa_pipeline_parallel}()}, which calls this;
#'   \code{\link{build_results_table}()} for a per-method breakdown.
#'
#' @examples
#' runs <- list(
#'   list(base_fitness = 12.4),
#'   list(base_fitness = 13.1),
#'   list(base_fitness = 12.9)
#' )
#' summarise_runs(runs)
#'
#' ## Recover the winning run
#' s <- summarise_runs(runs)
#' runs[[s$best_run]]
#'
#' @export
summarise_runs <- function(runs_list) {
  bf <- sapply(runs_list, `[[`, "base_fitness")
  data.frame(
    best_fitness = max(bf),
    mean_fitness = mean(bf),
    sd_fitness   = sd(bf),
    best_run     = which.max(bf)
  )
}

#' Cross-method comparison table (paper Table 3)
#'
#' Builds the one-row-per-method summary used to compare AdSA-ALMD against
#' TTWA, GEA, MSA or any other delineation: base fitness, cohesion, mean
#' self-containment on each side, market count, and misallocated BGUs.
#'
#' @details
#' Every method is scored with \emph{this package's} metrics, not with its own
#' objective, which is the point of the table — the numbers are comparable
#' because they are all computed the same way from the same \code{W} and
#' \code{adj}.
#'
#' The \code{names()} of \code{results_list} become the \code{Method} column, so
#' the list must be named. Each element needs four components:
#' \code{best_sol}, \code{base_fitness}, \code{gci} and \code{n_clusters}, plus
#' an \code{lma_metrics} table with \code{scss} and \code{scds} columns. Results
#' from \code{\link{adsa_almd}()} and \code{\link{run_adsa_pipeline}()} already
#' have all of these; for a competing method, construct the same shape (see the
#' examples).
#'
#' \code{SCSS} and \code{SCDS} are \emph{unweighted} means over markets, so a
#' tiny market counts as much as a large one. That is the convention the paper
#' uses; weight by population yourself if you want the alternative reading.
#' Numbers are rounded to three decimals for display, so do not use the
#' returned table for further arithmetic — go back to the run objects.
#'
#' Only \code{Misallocated_BGU} is recomputed here (via
#' \code{\link{count_misallocated}()}, the expensive part); the rest is read
#' straight off each result.
#'
#' @param results_list Named list of run results; see Details for the required
#'   components.
#' @param W Numeric \code{N x N} origin-destination matrix.
#' @param adj Integer/logical \code{N x N} adjacency matrix.
#' @param row_W Numeric vector, \code{rowSums(W)}.
#' @param col_W Numeric vector, \code{colSums(W)}.
#' @param sc_tol Numeric, self-containment tolerance for the misallocation
#'   test (default 0.05).
#' @return A \code{data.frame} with one row per method and columns
#'   \code{Method}, \code{Base_Fitness}, \code{GCI}, \code{SCSS}, \code{SCDS},
#'   \code{m} and \code{Misallocated_BGU}. Rows follow the order of
#'   \code{results_list}.
#'
#' @family diagnostics
#' @seealso \code{\link{count_misallocated}()},
#'   \code{\link{compute_base_fitness}()}, \code{\link{summarise_runs}()}
#'
#' @examples
#' adj <- matrix(0L, 6, 6); adj[cbind(1:5, 2:6)] <- 1L; adj <- adj + t(adj)
#' W <- matrix(1, 6, 6); W[1:3, 1:3] <- 30; W[4:6, 4:6] <- 30; diag(W) <- 60
#' row_W <- rowSums(W); col_W <- colSums(W)
#'
#' ## Wrap any partition into the shape this function expects
#' as_result <- function(sol) list(
#'   best_sol     = sol,
#'   base_fitness = compute_base_fitness(sol, W, row_W, col_W),
#'   gci          = compute_gci(sol, W, row_W, col_W),
#'   n_clusters   = n_clusters(sol),
#'   lma_metrics  = build_lma_df(sol, W, row_W, col_W)
#' )
#'
#' build_results_table(
#'   list(
#'     `AdSA-ALMD` = as_result(c(1, 1, 1, 2, 2, 2)),
#'     Naive       = as_result(c(1, 1, 1, 1, 2, 2)),
#'     Singletons  = as_result(1:6)
#'   ),
#'   W, adj, row_W, col_W
#' )
#'
#' @export
build_results_table <- function(results_list, W, adj, row_W, col_W,
                                 sc_tol = 0.05) {
  do.call(rbind, lapply(names(results_list), function(method) {
    res <- results_list[[method]]
    sol <- res$best_sol
    lma <- res$lma_metrics
    data.frame(
      Method           = method,
      Base_Fitness     = round(res$base_fitness, 3),
      GCI              = round(res$gci,          3),
      SCSS             = round(mean(lma$scss),   3),
      SCDS             = round(mean(lma$scds),   3),
      m                = res$n_clusters,
      Misallocated_BGU = count_misallocated(sol, W, adj, row_W, col_W, sc_tol),
      stringsAsFactors = FALSE
    )
  }))
}
