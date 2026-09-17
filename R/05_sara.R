# ============================================================
# FILE: R/05_sara.R
# SARA — Spatially Adaptive Random Aggregation (Algorithm 2)
# and the iterated driver that also respects the population ceiling.
# ============================================================

#' Rank the clusters adjacent to a focal cluster by mutual interaction
#'
#' Scores each neighbouring cluster by the same two-sided squared-flow measure
#' as the cohesion index, but aggregated to cluster level: how much commuting
#' flows between the focal cluster and that neighbour, relative to what each
#' side's marginals could support. Used by \code{\link{sara}()} to pick a merge
#' target.
#'
#' Zero-denominator cases contribute 0, as in \code{\link{compute_ci}()}.
#'
#' @param focal_cl Integer, focal cluster ID.
#' @param adj_cls  Integer vector, IDs of adjacent clusters.
#' @param sol      Integer vector of assignments.
#' @param W        Numeric OD matrix.
#' @param row_W    Numeric vector, rowSums(W).
#' @param col_W    Numeric vector, colSums(W).
#' @return Numeric vector of scores, one per entry of \code{adj_cls}, in the
#'   same order.
#' @noRd
.rank_adjacent_clusters <- function(focal_cl, adj_cls, sol, W, row_W, col_W) {
  focal_mbrs <- cluster_members(focal_cl, sol)

  vapply(adj_cls, function(cl) {
    other_mbrs <- cluster_members(cl, sol)
    w_fo <- sum(W[focal_mbrs, other_mbrs])
    w_of <- sum(W[other_mbrs, focal_mbrs])
    r_f <- sum(row_W[focal_mbrs])
    c_f <- sum(col_W[focal_mbrs])
    r_o <- sum(row_W[other_mbrs])
    c_o <- sum(col_W[other_mbrs])
    t1 <- if (r_f > 0 && c_o > 0)
      w_fo ^ 2 / (r_f * c_o)
    else
      0
    t2 <- if (r_o > 0 && c_f > 0)
      w_of ^ 2 / (r_o * c_f)
    else
      0
    t1 + t2
  }, numeric(1L))
}

#' Merge cl_from into cl_into and relabel consecutively
#'
#' @param sol      Integer vector of assignments.
#' @param cl_from  Integer, cluster that disappears.
#' @param cl_into  Integer, cluster that absorbs it.
#' @return Integer vector relabelled to \code{1 … m-1}.
#' @noRd
.merge_clusters <- function(sol, cl_from, cl_into) {
  sol[sol == cl_from] <- cl_into
  getOrderedVec(sol)
}

#' SARA: Spatially Adaptive Random Aggregation (Algorithm 2)
#'
#' Builds a feasible initial partition by repeated agglomeration: find the
#' least valid market, merge it into one of its strongest neighbours, and
#' repeat until every market satisfies the validity criterion. This is stage 1
#' of the pipeline — the solution it returns is what
#' \code{\link{adsa_almd}()} anneals.
#'
#' @details
#' Starting from \code{init_sol} (usually atomic, one BGU per market), each
#' iteration:
#' \enumerate{
#'   \item scores every current market with \code{\link{validity_score}()};
#'   \item stops if all scores are \code{>= 1} — every market is valid;
#'   \item samples the \strong{focal} market uniformly from the \code{n_low}
#'     lowest-scoring ones;
#'   \item finds the focal market's adjacent markets and ranks them by mutual
#'     interaction, halving the score of any whose merged population would
#'     exceed \code{Pop_max};
#'   \item samples the \strong{target} from the \code{n_high} best-ranked
#'     neighbours and merges focal into target;
#'   \item accepts the merge only if the merged market is contiguous;
#'     otherwise the candidate is discarded and the next iteration tries again.
#' }
#'
#' @section Why it is randomised:
#' Sampling at both choice points (steps 3 and 5) rather than always taking the
#' single worst market and single best neighbour is what makes SARA a
#' \emph{restartable} constructor: different seeds give genuinely different
#' feasible partitions, which is what
#' \code{\link{run_adsa_pipeline_parallel}()} exploits by running several
#' independent pipelines. Set \code{n_low = n_high = 1} for the deterministic
#' greedy variant.
#'
#' @section Termination:
#' Three ways out, and it is worth knowing which one you hit:
#' \itemize{
#'   \item \strong{Success}: every market scores \code{>= 1}.
#'   \item \strong{Stuck}: the market count fails to fall for 100 consecutive
#'     iterations — usually every remaining merge would break contiguity.
#'     Emits a \code{"SARA stuck"} message when \code{verbose = TRUE}.
#'   \item \strong{Capped}: \code{max_iter} iterations elapsed.
#' }
#' In all three cases a relabelled partition is returned, so a result is
#' \strong{not} a guarantee of validity. Check with
#' \code{\link{build_lma_df}()} or \code{\link{validity_score}()} before
#' relying on it.
#'
#' @section Restricting the merge region (active_bgus):
#' \code{active_bgus} confines the work to a subset of the map: only markets
#' containing an active BGU are considered for merging, and neighbours inside
#' the active set are preferred as targets (falling back to any adjacent market
#' when the focal market has no active neighbour). \code{\link{sara_iterated}()}
#' uses this to re-partition just the BGUs it has shattered out of an oversized
#' market, leaving the rest of the map alone. Leave \code{NULL} to work on the
#' whole map.
#'
#' @param init_sol Integer vector, starting partition. \code{seq_len(N)} gives
#'   the atomic start used by \code{\link{sara_iterated}()}.
#' @param W Numeric \code{N x N} origin-destination matrix.
#' @param adj Integer/logical \code{N x N} adjacency matrix.
#' @param row_W Numeric vector, \code{rowSums(W)}.
#' @param col_W Numeric vector, \code{colSums(W)}.
#' @param SC_min Numeric in (0, 1], the minimum self-containment threshold.
#' @param Pop_min Numeric > 0, the minimum market population.
#' @param Pop_max Numeric, the maximum market population. Not a hard
#'   constraint here — merges that would breach it have their score halved
#'   rather than being forbidden, so SARA can still exceed it when there is no
#'   alternative.
#' @param n_low Integer, size of the candidate pool the focal market is sampled
#'   from (default 3). Capped at the number of invalid markets.
#' @param n_high Integer, size of the candidate pool the merge target is
#'   sampled from (default 3).
#' @param max_iter Integer, hard cap on iterations (default 10000).
#' @param active_bgus Integer vector of BGU indices to confine merging to, or
#'   \code{NULL} (default) for the whole map. See the section above.
#' @param verbose Logical (default \code{TRUE}); emit a message if the merge
#'   loop stalls.
#' @return Integer vector of length \code{N}, relabelled to \code{1 … m} by
#'   \code{\link{get_ordered_vec}()}.
#'
#' @family SARA
#' @seealso \code{\link{sara_iterated}()}, the driver you normally call;
#'   \code{\link{validity_score}()} for the criterion being satisfied;
#'   \code{\link{sha}()} for the same idea applied to a BGU subset.
#'
#' @examples
#' ## Toy system: six BGUs in a line, forming two 3-BGU markets
#' adj <- matrix(0L, 6, 6); adj[cbind(1:5, 2:6)] <- 1L; adj <- adj + t(adj)
#' W <- matrix(1, 6, 6); W[1:3, 1:3] <- 30; W[4:6, 4:6] <- 30; diag(W) <- 60
#' row_W <- rowSums(W); col_W <- colSums(W)
#'
#' set.seed(42)
#' sol <- sara(seq_len(6), W, adj, row_W, col_W,
#'             SC_min = 0.55, Pop_min = 200, Pop_max = 1000,
#'             verbose = FALSE)
#' sol
#' build_lma_df(sol, W, row_W, col_W)
#'
#' ## Always confirm the result really is valid — see Termination
#' all(vapply(unique_clusters(sol), function(cl)
#'       validity_score(cluster_members(cl, sol), W, row_W, col_W, adj,
#'                      SC_min = 0.55, Pop_min = 200),
#'     numeric(1)) >= 1)
#'
#' ## Different seeds explore different feasible partitions
#' set.seed(1); a <- sara(seq_len(6), W, adj, row_W, col_W, 0.55, 200, 1000,
#'                        verbose = FALSE)
#' set.seed(7); b <- sara(seq_len(6), W, adj, row_W, col_W, 0.55, 200, 1000,
#'                        verbose = FALSE)
#' partitions_match(a, b)$match
#'
#' @export
sara <- function(init_sol,
                 W,
                 adj,
                 row_W,
                 col_W,
                 SC_min,
                 Pop_min,
                 Pop_max,
                 n_low = 3L,
                 n_high = 3L,
                 max_iter = 10000L,
                 active_bgus = NULL,
                 verbose = TRUE) {
  sol <- getOrderedVec(init_sol)
  iter <- 0L

  if (is.null(active_bgus)) {
    active_bgus <- seq_len(length(sol))
  }

  prev_n <- length(unique(sol[active_bgus]))
  stuck <- 0L

  repeat {
    iter <- iter + 1L
    if (iter > max_iter) {
      break
    }

    clusters <- unique(sol[active_bgus])

    curr_n <- length(clusters)

    if (curr_n == prev_n) {
      stuck <- stuck + 1L
    } else {
      stuck <- 0L
    }

    prev_n <- curr_n

    if (stuck > 100L) {
      if (verbose) message("SARA stuck")
      break
    }

    # Validity scores for every current cluster
    vs <- vapply(clusters, function(cl) {
      validity_score(cluster_members(cl, sol),
                     W,
                     row_W,
                     col_W,
                     adj,
                     SC_min,
                     Pop_min, Pop_max)
    }, numeric(1L))

    # All valid - done
    if (all(vs >= 1))
      break

    # Step 1: sample focal cluster from the n_low least-valid
    n_invalid <- max(1L, sum(vs < 1))
    low_idx <- order(vs)[seq_len(min(n_low, n_invalid))]
    focal_cl <- clusters[sample(low_idx, 1L)]
    focal_mbrs <- cluster_members(focal_cl, sol)

    # Step 2: Find adjacent BGUs and their cluster IDs, preferring other
    # active clusters first so the merge stays confined to the region
    # currently being re-partitioned. Falls back to any adjacent cluster
    # if the focal cluster has no active neighbour.
    adj_bgus <- setdiff(intersect(which(colSums(adj[focal_mbrs, , drop = FALSE]) > 0L), active_bgus), focal_mbrs)
    adj_cls <- setdiff(unique(sol[adj_bgus]), focal_cl)

    if (length(adj_cls) == 0L) {
      adj_bgus <- setdiff(which(colSums(adj[focal_mbrs, , drop = FALSE]) > 0L), focal_mbrs)
      adj_cls <- setdiff(unique(sol[adj_bgus]), focal_cl)
    }

    if (length(adj_cls) == 0L) {
      # No neighbours anywhere in the graph -- a genuinely disconnected
      # fragment. Nothing this function can do; stop gracefully.
      break
    }

    focal_pop <- compute_population(focal_mbrs, row_W)

    adj_pops <- vapply(adj_cls, function(cl) compute_population(cluster_members(cl, sol), row_W), numeric(1L))
    scores <- .rank_adjacent_clusters(focal_cl, adj_cls, sol, W, row_W, col_W)

    for (i in seq_along(adj_cls)) {

      adj_pop <- adj_pops[i]
      merged_pop <- focal_pop + adj_pop

      if (merged_pop > Pop_max) {
        scores[i] <- scores[i] * 0.5
        next
      }

    }

    top_idx <- order(scores)[seq_len(min(n_high, length(adj_cls)))]
    target_cl <- adj_cls[sample(top_idx, 1L)]

    # Step 4: merge focal into target, check contiguity
    candidate <- getOrderedVec(.merge_clusters(sol, focal_cl, target_cl))
    merged_id <- candidate[focal_mbrs[1L]]  # new label of merged cluster
    merged_mbrs <- cluster_members(merged_id, candidate)

    if (is_cluster_contiguous(merged_mbrs, adj)) {
      sol <- candidate
    }
    # If non-contiguous, skip this merge and try again next iteration
  }
  getOrderedVec(sol)
}

#' Score how badly a partition violates the population ceiling
#'
#' Summarises the population-ceiling violations of a partition as two numbers:
#' how many markets are over \code{Pop_max}, and by how much the worst one
#' exceeds it. \code{\link{sara_iterated}()} uses this to decide whether one
#' attempt beat another.
#'
#' @details
#' Deliberately separate from \code{\link{penalty_pop_max_total}()}: that
#' function collapses violations into a single unbounded number dominated by
#' the worst offender, which is the right thing inside a smooth objective but
#' the wrong thing when comparing two infeasible partitions. Keeping excess and
#' count apart lets \code{\link{is_better_max_pop}()} apply them in priority
#' order.
#'
#' Note \code{excess} is measured from the largest market only, not summed over
#' all offenders.
#'
#' @param sol Integer vector of market assignments (length N).
#' @param W Numeric \code{N x N} origin-destination matrix.
#' @param row_W Numeric vector, \code{rowSums(W)}.
#' @param col_W Numeric vector, \code{colSums(W)}.
#' @param Pop_max Numeric, the maximum market population.
#' @return Named list:
#'   \describe{
#'     \item{\code{n_violations}}{integer, how many markets exceed
#'       \code{Pop_max}.}
#'     \item{\code{excess}}{numeric \eqn{\ge 0}, how far the \emph{largest}
#'       market exceeds \code{Pop_max}. Exactly 0 when the partition is
#'       feasible.}
#'   }
#'
#' @family SARA
#' @seealso \code{\link{is_better_max_pop}()} for the comparison,
#'   \code{\link{sara_iterated}()} for the loop that uses both.
#'
#' @examples
#' W <- matrix(1, 6, 6); W[1:3, 1:3] <- 30; W[4:6, 4:6] <- 30; diag(W) <- 60
#' row_W <- rowSums(W); col_W <- colSums(W)
#'
#' ## Two markets of 369 each, against a ceiling of 400: feasible
#' sara_max_pop_score(c(1, 1, 1, 2, 2, 2), W, row_W, col_W, Pop_max = 400)
#'
#' ## One market of 738 against the same ceiling
#' sara_max_pop_score(rep(1, 6), W, row_W, col_W, Pop_max = 400)
#'
#' @export
sara_max_pop_score <- function(sol, W, row_W, col_W, Pop_max) {

  lma_df <- build_lma_df(sol, W, row_W, col_W)

  n_violations <- sum(lma_df$pop > Pop_max)

  max_pop <- max(lma_df$pop)
  excess <- max(max_pop - Pop_max, 0)

  list(
    n_violations = n_violations,
    excess = excess
  )
}

#' Compare two population-ceiling scores
#'
#' Decides whether score \code{a} is better than score \code{b}, using excess
#' as the primary criterion and violation count only as a tie-break.
#'
#' @details
#' \enumerate{
#'   \item If the two \code{excess} values differ, the smaller one wins.
#'   \item Otherwise the smaller \code{n_violations} wins.
#' }
#'
#' Excess comes first because shrinking the single worst market is the move
#' that actually unblocks feasibility: a partition with five markets slightly
#' over the ceiling is closer to done than one with a single market at triple
#' it. The tie-break then prefers fewer offenders among partitions whose worst
#' market is equally bad.
#'
#' Ties on both criteria return \code{FALSE}, so \code{\link{sara_iterated}()}
#' keeps its incumbent rather than churning between equivalent solutions.
#'
#' @param score_a,score_b Lists as returned by
#'   \code{\link{sara_max_pop_score}()}, each with \code{$excess} and
#'   \code{$n_violations}.
#' @return \code{TRUE} if \code{score_a} is strictly better than
#'   \code{score_b}, \code{FALSE} otherwise (including exact ties).
#'
#' @family SARA
#' @seealso \code{\link{sara_max_pop_score}()}
#'
#' @examples
#' a <- list(n_violations = 5L, excess = 100)
#' b <- list(n_violations = 1L, excess = 900)
#' is_better_max_pop(a, b)   # TRUE: smaller worst-market excess wins
#'
#' ## Equal excess: fewer offenders wins
#' is_better_max_pop(list(n_violations = 1L, excess = 100),
#'                   list(n_violations = 4L, excess = 100))
#'
#' ## An exact tie is not an improvement
#' is_better_max_pop(a, a)
#'
#' @export
is_better_max_pop <- function(score_a, score_b) {

  if (score_a$excess != score_b$excess) {
    return(score_a$excess < score_b$excess)
  }

  return(score_a$n_violations <
           score_b$n_violations)
}

#' Iterated SARA: aggregate, then break up markets over the population ceiling
#'
#' The stage-1 entry point. \code{\link{sara}()} on its own satisfies the
#' self-containment and population \emph{floor}, but can leave markets over
#' \code{Pop_max} because it only discourages, never forbids, oversized merges.
#' This driver runs SARA, shatters whichever markets came out too big, and
#' re-aggregates just those BGUs — repeating until the ceiling is met or it
#' runs out of attempts.
#'
#' @details
#' Starting from the atomic partition \code{seq_len(length(row_W))}, each pass:
#' \enumerate{
#'   \item runs \code{\link{sara}()} at the current self-containment threshold,
#'     confined to the currently active BGUs;
#'   \item scores the result with \code{\link{sara_max_pop_score}()} and
#'     \strong{returns immediately} if no market exceeds \code{Pop_max};
#'   \item otherwise records it if it beats the incumbent under
#'     \code{\link{is_better_max_pop}()};
#'   \item explodes every oversized market back into singletons and marks those
#'     BGUs as the active set for the next pass, so the rest of the map is left
#'     untouched.
#' }
#'
#' @section Threshold relaxation:
#' If a pass shatters exactly the same set of BGUs as the previous one, the
#' loop is cycling — re-aggregating the same region to the same oversized
#' market. When that happens the self-containment threshold is multiplied by
#' 0.95, giving \code{\link{sara}()} room to accept splits it previously
#' rejected; any pass that touches a different region resets the threshold to
#' its original value. The outer \code{while} exits if the threshold decays
#' below 0.05.
#'
#' The practical consequence: the returned partition may have been built at a
#' \emph{lower} self-containment threshold than you passed in. The population
#' ceiling takes priority over the self-containment floor.
#'
#' @section What you get back:
#' Either the first fully feasible partition found (early return at step 2), or
#' the best-scoring infeasible one after \code{max_it} passes. There is no flag
#' distinguishing the two — check with \code{\link{sara_max_pop_score}()} if it
#' matters.
#'
#' @param W Numeric \code{N x N} origin-destination matrix.
#' @param adj Integer/logical \code{N x N} adjacency matrix.
#' @param row_W Numeric vector, \code{rowSums(W)}. Its length defines \code{N}.
#' @param col_W Numeric vector, \code{colSums(W)}.
#' @param SC_min Numeric in (0, 1], the starting self-containment threshold.
#'   May be relaxed during the run (see the section above).
#' @param Pop_min Numeric > 0, the minimum market population.
#' @param Pop_max Numeric, the maximum market population. The constraint this
#'   function exists to satisfy.
#' @param max_it Integer, maximum aggregate-and-split passes (default 5).
#' @param verbose Logical (default \code{TRUE}); print the current excess, the
#'   per-market table and the relaxed threshold after each pass.
#' @param temp Integer vector or \code{NULL} (default). A ready-made partition
#'   to use for the first pass instead of calling \code{\link{sara}()} — useful
#'   for resuming from a checkpoint or starting from another method's output.
#'   Consumed once, then ignored.
#' @return Integer vector of length \code{N}, relabelled to \code{1 … m}.
#'
#' @family SARA
#' @seealso \code{\link{sara}()} for the single aggregation pass;
#'   \code{\link{run_adsa_pipeline}()}, which calls this as stage 1;
#'   \code{\link{operator_11}()}, which handles oversized markets during
#'   annealing.
#'
#' @examples
#' ## Toy system: six BGUs in a line, forming two 3-BGU markets
#' adj <- matrix(0L, 6, 6); adj[cbind(1:5, 2:6)] <- 1L; adj <- adj + t(adj)
#' W <- matrix(1, 6, 6); W[1:3, 1:3] <- 30; W[4:6, 4:6] <- 30; diag(W) <- 60
#' row_W <- rowSums(W); col_W <- colSums(W)
#'
#' set.seed(42)
#' init_sol <- sara_iterated(W, adj, row_W, col_W,
#'                           SC_min = 0.55, Pop_min = 200, Pop_max = 500,
#'                           max_it = 5, verbose = FALSE)
#' init_sol
#' build_lma_df(init_sol, W, row_W, col_W)
#'
#' ## Did it actually meet the ceiling, or run out of passes?
#' sara_max_pop_score(init_sol, W, row_W, col_W, Pop_max = 500)
#'
#' ## Seed the first pass from a partition you already have
#' set.seed(42)
#' sara_iterated(W, adj, row_W, col_W, 0.55, 200, 500,
#'               verbose = FALSE, temp = c(1, 1, 1, 2, 2, 2))
#'
#' @export
sara_iterated <- function(W,
                          adj,
                          row_W,
                          col_W,
                          SC_min,
                          Pop_min,
                          Pop_max,
                          max_it = 5L,
                          verbose = TRUE,
                          temp = NULL) {

  start_sol <- getOrderedVec(seq_len(length(row_W)))
  iter <- 0L
  best_sol <- start_sol

  best_score <- list(
    n_violations = Inf,
    excess = Inf
  )

  # Tracks the true best solution across passes (best_sol below is
  # re-aliased to the current attempt every pass for the split logic).
  overall_best_sol <- start_sol
  overall_best_score <- best_score

  min_accept_slf <- 0.5
  Min_self <- SC_min
  init_min_self <- SC_min

  active_bgus <- start_sol
  prev_mbrs <- active_bgus

  while (Min_self > 0.05) {

    iter <- iter + 1L
    if (iter > max_it) {
      break
    }

    if (!is.null(temp)) {
      sol <- temp
      temp <- NULL
    }
    else {
      sol <- sara(
        init_sol = start_sol,
        W = W,
        adj = adj,
        row_W = row_W,
        col_W = col_W,
        SC_min = Min_self,
        Pop_min = Pop_min,
        Pop_max = Pop_max,
        n_low = 5L,
        n_high = 5L,
        max_iter = 10000L,
        active_bgus = active_bgus,
        verbose = verbose
      )
    }

    current_score <- sara_max_pop_score(sol, W, row_W, col_W, Pop_max)

    if (current_score$excess == 0) {
      return(sol)
    }

    if (is_better_max_pop(current_score, overall_best_score)) {
      overall_best_score <- current_score
      overall_best_sol <- sol
    }

    best_score <- current_score
    best_sol <- sol
    lma_df <- build_lma_df(sol, W, row_W, col_W)
    if (verbose) {
      print(current_score$excess)
      print(lma_df)
    }
    big_cls <- lma_df$cluster[lma_df$pop > Pop_max]
    split_bgus <- integer(0)
    cl <- big_cls[1]
    all_mbrs <- c()

    for (cl in big_cls) {
      mbrs <- cluster_members(cl, best_sol)
      all_mbrs <- c(all_mbrs, mbrs)
      split_bgus <- c(split_bgus, mbrs)
      new_ids <- seq(
        max(sol) + 1L,
        max(sol) + length(mbrs))
      sol[mbrs] <- new_ids
    }
    active_bgus <- split_bgus
    if (identical(sort(all_mbrs), sort(prev_mbrs))) {
      Min_self <- Min_self * 0.95
      if (verbose) print(Min_self)
    } else {
      Min_self <- init_min_self
    }

    prev_mbrs <- all_mbrs
    start_sol <- getOrderedVec(sol)

  }

  getOrderedVec(overall_best_sol)
}
