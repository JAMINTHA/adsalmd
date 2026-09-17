# ============================================================
# FILE: R/12_validation_utils.R
# Helpers for validating the method against known-answer test cases,
# and for choosing SC_min/Pop_min/Pop_max/Pop_tar from a reference
# partition rather than guessing.
# ============================================================

#' Check whether two partitions are identical up to relabelling
#'
#' Two cluster-assignment vectors represent the same grouping if every
#' BGU pair that's co-clustered in one is co-clustered in the other, and
#' vice versa -- regardless of which integer label each cluster happens
#' to have. Use this to check whether an algorithm's output
#' (e.g. \code{result$best_sol}) recovered a known true partition, since
#' cluster labels are arbitrary and won't line up by value.
#'
#' @details
#' Builds the contingency table of \code{sol_a} against \code{sol_b} and
#' declares a match when every row and every column has exactly one non-zero
#' entry — that is, when the two labellings are related by a bijection. This is
#' a strict all-or-nothing test, not a similarity measure: a single BGU in the
#' wrong market makes \code{match} \code{FALSE}.
#'
#' When it fails, the returned table is the diagnosis. A row with two non-zero
#' entries means one true market was split across two predicted ones; a column
#' with two means two true markets were merged. The off-diagonal cell counts
#' tell you how many BGUs are involved, so you can tell a near miss from a
#' structurally different answer.
#'
#' The comparison is symmetric, but the table is labelled \code{true} (rows)
#' and \code{predicted} (columns), so pass the known partition first.
#'
#' For a graded measure rather than a yes/no — the adjusted Rand index, say —
#' use a dedicated clustering-comparison package; this function deliberately
#' has no dependencies.
#'
#' @param sol_a Integer vector, e.g. the true/known partition.
#' @param sol_b Integer vector, e.g. the algorithm's output partition.
#'   Must be the same length as sol_a.
#' @return Named list: \code{match} (logical, TRUE if identical up to
#'   relabelling), \code{table} (the cluster-by-cluster contingency
#'   table -- a clean diagonal-only pattern is what a perfect match
#'   looks like; anything off that diagonal shows exactly which BGUs
#'   / clusters disagree).
#'
#' @family validation
#' @seealso \code{\link{get_ordered_vec}()} for canonical relabelling;
#'   \code{\link{run_adsa_pipeline}()}, whose \code{best_sol} this is meant to
#'   check; \code{\link{count_misallocated}()} for a quality measure that needs
#'   no known answer.
#'
#' @examples
#' truth <- c(1, 1, 1, 2, 2, 2)
#'
#' ## Same grouping, different labels
#' partitions_match(truth, c(7, 7, 7, 3, 3, 3))
#'
#' ## One BGU in the wrong market
#' bad <- partitions_match(truth, c(1, 1, 1, 1, 2, 2))
#' bad$match
#' bad$table
#'
#' ## A market split in two: read it off the row with two entries
#' partitions_match(truth, c(1, 1, 2, 3, 3, 3))$table
#'
#' ## Typical use after a run
#' \dontrun{
#' partitions_match(true_sol, res$best_sol)$match
#' }
#'
#' @export
partitions_match <- function(sol_a, sol_b) {
  stopifnot(length(sol_a) == length(sol_b))
  tab <- table(true = sol_a, predicted = sol_b)
  is_match <- all(rowSums(tab > 0) == 1L) && all(colSums(tab > 0) == 1L)
  list(match = is_match, table = tab)
}

#' Suggest SC_min / Pop_min / Pop_max / Pop_tar from a reference partition
#'
#' If you already know a "good" partition for some data -- a synthetic
#' test case with a known true answer, a previous method's output, or a
#' business rule -- this reads off that partition's own per-market
#' self-containment and population (via \code{\link{build_lma_df}}) and
#' proposes thresholds a small margin below/above the extremes, so real
#' markets in that partition aren't excluded by the very thresholds
#' meant to describe them.
#'
#' This is a starting point, not a substitute for the kind of grid
#' search described in the paper (Table 2) -- treat the output as
#' something to sanity-check and refine, not use blindly.
#'
#' @details
#' From the reference partition's own \code{\link{build_lma_df}()} table:
#' \describe{
#'   \item{\code{SC_min}}{\code{min(sc) * (1 - margin)} — just below the least
#'     self-contained market in the reference.}
#'   \item{\code{Pop_min}}{\code{min(pop) * (1 - margin)} — just below its
#'     smallest market.}
#'   \item{\code{Pop_max}}{\code{max(pop) * (1 + margin)} — just above its
#'     largest.}
#'   \item{\code{Pop_tar}}{\code{median(pop)} — no margin applied, since a
#'     target is not a constraint.}
#' }
#' All floors are clamped at 0.
#'
#' The logic is simply that thresholds meant to \emph{describe} a set of
#' markets should not exclude those markets. Every bound is driven by a single
#' extreme market, so \code{lma_df} is returned with the suggestion: look at it
#' before accepting the numbers, because one unusual market sets each one.
#'
#' @section What this is not:
#' Not a substitute for the grid search described in the paper (Table 2). It
#' reads thresholds off one partition you already believe in, which means a
#' poor reference yields poor thresholds with no warning. Use it to get a run
#' started on new data, then tune.
#'
#' Two failure modes to watch for: a reference containing a degenerate market
#' (a singleton with \code{sc} near 0) drags \code{SC_min} to near zero and
#' makes the constraint toothless; and a wide population spread produces bounds
#' so loose that \code{\link{compute_objective}()}'s hard rejections never
#' fire. In both cases set the offending bound by hand.
#'
#' @param sol   Integer vector, a reference partition.
#' @param W     Numeric OD matrix (N x N).
#' @param row_W Numeric vector, rowSums(W).
#' @param col_W Numeric vector, colSums(W).
#' @param margin Numeric in [0, 1), how far below the reference
#'   partition's minimum SC/population (and above its maximum
#'   population) to set the suggested thresholds (default 0.10, i.e.
#'   10 percent). Use 0 for the reference's exact extremes.
#' @return Named list: \code{SC_min}, \code{Pop_min}, \code{Pop_max},
#'   \code{Pop_tar} (the median market population), plus \code{lma_df}
#'   (the per-market table the suggestion was derived from, so you can
#'   see exactly which market drove each number).
#'
#' @family validation
#' @seealso \code{\link{AdSA_params}} for the study's own values;
#'   \code{\link{build_lma_df}()} for the table this reads;
#'   \code{\link{run_adsa_pipeline}()}, which consumes the thresholds.
#'
#' @examples
#' ## Toy system: six BGUs in a line, forming two 3-BGU markets
#' W <- matrix(1, 6, 6); W[1:3, 1:3] <- 30; W[4:6, 4:6] <- 30; diag(W) <- 60
#' row_W <- rowSums(W); col_W <- colSums(W)
#'
#' th <- suggest_thresholds(c(1, 1, 1, 2, 2, 2), W, row_W, col_W)
#' th[c("SC_min", "Pop_min", "Pop_max", "Pop_tar")]
#'
#' ## Always look at what drove each number
#' th$lma_df
#'
#' ## A tighter margin gives the reference's exact extremes
#' suggest_thresholds(c(1, 1, 1, 2, 2, 2), W, row_W, col_W,
#'                    margin = 0)$SC_min
#'
#' ## A poor reference gives poor thresholds, silently
#' suggest_thresholds(c(1, 1, 1, 1, 1, 2), W, row_W, col_W)$SC_min
#'
#' ## Feed them straight into a run
#' \dontrun{
#' res <- run_adsa_pipeline(W, adj, row_W, col_W,
#'                          SC_min = th$SC_min, Pop_min = th$Pop_min,
#'                          Pop_max = th$Pop_max, Pop_tar = th$Pop_tar,
#'                          r = 3, L = 50L)
#' }
#'
#' @export
suggest_thresholds <- function(sol, W, row_W, col_W, margin = 0.10) {
  lma_df <- build_lma_df(sol, W, row_W, col_W)

  list(
    SC_min  = max(0, min(lma_df$sc)  * (1 - margin)),
    Pop_min = max(0, min(lma_df$pop) * (1 - margin)),
    Pop_max = max(lma_df$pop) * (1 + margin),
    Pop_tar = stats::median(lma_df$pop),
    lma_df  = lma_df
  )
}
