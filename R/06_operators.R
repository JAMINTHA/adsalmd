# ============================================================
# FILE: R/06_operators.R
# Eleven group-based perturbation operators.
# Implements the operators from Martínez-Bernabeu et al. (2012)
# as adapted for AdSA-ALMD.
#
# All operators:
#   - Accept the full parameter set as arguments (no globals).
#   - Return the original sol unchanged if the move is infeasible.
#   - Never break contiguity of any existing cluster.
# ============================================================

# ---- Internal helpers (not exported) -------------------------

#' Border BGUs: members of a cluster that touch another cluster
#'
#' @param members Integer vector of BGU indices in the cluster.
#' @param sol     Integer vector of assignments (unused; kept for signature
#'   symmetry with the other helpers).
#' @param adj     Integer/logical adjacency matrix.
#' @return Integer vector of member indices with at least one neighbour outside
#'   the cluster; \code{integer(0)} for an interior-only or whole-map cluster.
#' @noRd
.border_bgus <- function(members, sol, adj) {
  all_adj      <- which(colSums(adj[members, , drop = FALSE]) > 0L)
  external_adj <- setdiff(all_adj, members)
  if (length(external_adj) == 0L)
    return(integer(0))
  members[vapply(members, function(g)
    any(adj[g, external_adj] > 0L), logical(1L))]
}

#' Cluster IDs adjacent to a given set of members
#'
#' @param members Integer vector of BGU indices.
#' @param sol     Integer vector of assignments.
#' @param adj     Integer/logical adjacency matrix.
#' @return Integer vector of distinct neighbouring cluster labels. May include
#'   the members' own cluster if \code{members} is a strict subset of it.
#' @noRd
.adjacent_clusters <- function(members, sol, adj) {
  adj_bgus <- setdiff(which(colSums(adj[members, , drop = FALSE]) > 0L), members)
  unique(sol[adj_bgus])
}

#' Pairwise cluster interaction score (mutual CI)
#'
#' Cluster-level analogue of \code{\link{compute_ci}()}: two-sided squared
#' flow between clusters A and B, normalised by the relevant marginals.
#' Non-finite results (empty clusters, zero marginals) collapse to 0.
#'
#' @param mbrs_a,mbrs_b Integer vectors of BGU indices.
#' @param W     Numeric OD matrix.
#' @param row_W Numeric vector, rowSums(W).
#' @param col_W Numeric vector, colSums(W).
#' @return Numeric scalar >= 0.
#' @noRd
.cluster_interaction <- function(mbrs_a, mbrs_b, W, row_W, col_W) {
  w_ab <- sum(W[mbrs_a, mbrs_b])
  w_ba <- sum(W[mbrs_b, mbrs_a])
  ra   <- sum(row_W[mbrs_a])
  ca   <- sum(col_W[mbrs_a])
  rb   <- sum(row_W[mbrs_b])
  cb   <- sum(col_W[mbrs_b])
  t1   <- if (rb > 0 && cb > 0)
    w_ab^2 / (ra * cb)
  else 0
  t2   <- if (rb > 0 && ca > 0)
    w_ba^2 / (rb * ca)
  else 0
  result <- t1 + t2

  if(!is.finite(result)) return(0)
  result
}



#' Apply a candidate solution only if all clusters remain contiguous
#'
#' The standard tail call of most operators: accept the move, or silently give
#' back the untouched solution.
#'
#' @param candidate Integer vector, the proposed partition.
#' @param sol       Integer vector, the partition to fall back to.
#' @param adj       Integer/logical adjacency matrix.
#' @return \code{candidate} if every cluster is contiguous, else \code{sol}.
#' @noRd
.accept_if_contiguous <- function(candidate, sol, adj) {
  if (all_clusters_contiguous(candidate, adj))
    candidate
  else
    sol
}

#' Tournament selection over a score vector
#'
#' Picks uniformly at random from the \code{n_way} best entries rather than
#' always taking the single best. This is what keeps the operators stochastic:
#' with \code{n_way = 1} every operator becomes a deterministic greedy move and
#' the annealing search stops exploring.
#'
#' @param values   Numeric vector of scores.
#' @param n_way    Integer, tournament size (default 3); capped at
#'   \code{length(values)}.
#' @param minimize Logical (default TRUE); TRUE selects among the lowest
#'   values, FALSE among the highest.
#' @return Integer index into \code{values}, or NULL if \code{values} is empty.
#' @noRd
.tournament_select <- function(values,
                               n_way = 3L,
                               minimize = TRUE){
  n <- length(values)
  if(n == 0L)
    return(NULL)
  n_way <- min(n_way, n)
  if(minimize){
    idx <- order(values)[seq_len(n_way)]
  } else{
    idx <- order(values, decreasing = TRUE) [seq_len(n_way)]
  }

  sample(idx, 1L)
}

# ---- Operator 1: Merge least-valid cluster into best neighbour ------

#' Operator 1: merge a weakly self-contained market into a neighbour, then
#' re-aggregate
#'
#' Picks a market with low self-containment, merges it into an adjacent market
#' chosen with probability proportional to mutual interaction, then shatters the
#' combined market back into singletons and rebuilds it with
#' \code{\link{sha}()}. The net effect is a wholesale re-partition of two
#' neighbouring markets — the most disruptive move in the operator set after
#' \code{\link{operator_10}()}.
#'
#' @details
#' Step by step:
#' \enumerate{
#'   \item \code{\link{build_lma_df}()}; abort (return \code{sol} unchanged) if
#'     fewer than two markets exist.
#'   \item Tournament-select a focal market from the three lowest \code{sc}
#'     values.
#'   \item Sample one adjacent market with probability proportional to
#'     \code{.cluster_interaction()} (floored at \code{1e-10} so a market with
#'     no measurable interaction can still be chosen).
#'   \item Merge the two, then assign every BGU of the merged market its own
#'     fresh label.
#'   \item Re-aggregate exactly those BGUs with \code{\link{sha}()}
#'     (\code{max_iter = 1000}), which may produce a different number of
#'     markets than the two it started from.
#'   \item Return the result only if every market is contiguous; otherwise
#'     return \code{sol} unchanged.
#' }
#'
#' Because \code{\link{sha}()} re-derives the partition of the merged region
#' from scratch, this operator can both split and join — it is not simply a
#' merge.
#'
#' @section Global variables:
#' Reads \code{params$SC_min} and \code{params$Pop_min} from the calling
#' environment's search path (in practice, the global environment) to configure
#' \code{\link{sha}()}. It will fail with "object 'params' not found" if called
#' directly without one in scope. \code{\link{run_adsa_pipeline}()} sets this
#' up for you; when calling by hand, define
#' \code{params <- list(SC_min = ..., Pop_min = ..., Pop_max = ...)} first.
#'
#' @param sol Integer vector of market assignments (length N).
#' @param W Numeric \code{N x N} origin-destination matrix.
#' @param adj Integer/logical \code{N x N} adjacency matrix.
#' @param row_W Numeric vector, \code{rowSums(W)}.
#' @param col_W Numeric vector, \code{colSums(W)}.
#' @param ... Ignored; present so every operator shares one signature and the
#'   dispatcher can pass thresholds through uniformly.
#' @return Integer vector of length N, relabelled to \code{1 … m}. Equal to
#'   \code{sol} when the move was rejected.
#'
#' @family perturbation operators
#' @seealso \code{\link{apply_random_operator}()} for the dispatcher,
#'   \code{\link{sha}()} for the re-aggregation step,
#'   \code{\link{operator_7}()} for a plain merge with no re-aggregation.
#'
#' @references
#' Martínez-Bernabeu, L., Flórez-Revuelta, F. and Casado-Díaz, J. M. (2012)
#' Grouping genetic operators for the delineation of functional areas based on
#' spatial interaction. \emph{Expert Systems with Applications}.
#'
#' @examples
#' ## Toy system: six BGUs in a line, forming two 3-BGU markets
#' adj <- matrix(0L, 6, 6); adj[cbind(1:5, 2:6)] <- 1L; adj <- adj + t(adj)
#' W <- matrix(1, 6, 6); W[1:3, 1:3] <- 30; W[4:6, 4:6] <- 30; diag(W) <- 60
#' row_W <- rowSums(W); col_W <- colSums(W)
#' sol <- c(1, 1, 1, 2, 2, 2)
#'
#' ## This operator reads a global `params` (see the section above)
#' params <- list(SC_min = 0.55, Pop_min = 100, Pop_max = 1000)
#'
#' set.seed(1)
#' operator_1(sol, W, adj, row_W, col_W)
#'
#' ## A single market has nothing to merge into: returned unchanged
#' identical(operator_1(rep(1, 6), W, adj, row_W, col_W), rep(1, 6))
#'
#' @export
operator_1 <- function(sol, W, adj, row_W, col_W, ...) {
  lma_df   <- build_lma_df(sol, W, row_W, col_W)
  if(nrow(lma_df) < 2L)
    return(sol)

  idx <- .tournament_select(lma_df$sc, n_way = 3L, minimize = TRUE)

  worst_cl <- lma_df$cluster[idx]

  w_mbrs   <- cluster_members(worst_cl, sol)

  adj_cls  <- .adjacent_clusters(w_mbrs, sol, adj)

  if (length(adj_cls) == 0L)
    return(sol)

  scores  <- vapply(adj_cls, function(cl) {
    .cluster_interaction(w_mbrs, cluster_members(cl, sol), W, row_W, col_W)
  }, numeric(1L))

  scores <- pmax(scores, 1e-10)

  best_cl <- adj_cls[sample(length(adj_cls),
                            size = 1L,
                            prob = scores / sum(scores))]

  cand <- sol
  cand[sol == best_cl] <- worst_cl
  cand <- getOrderedVec(cand)

  merged_id <- cand[w_mbrs[1L]]

  merged_mbrs <- cluster_members(merged_id, cand)

  for(i in seq_along(merged_mbrs)){
    new_id <- max(cand) + 1L
    cand[merged_mbrs[i]] <- new_id
  }

  cand <- getOrderedVec(cand)

  cand <- sha(
    bgus = merged_mbrs,
    sol = cand,
    W = W,
    adj = adj,
    row_W = row_W,
    col_W = col_W,
    SC_min = params$SC_min,
    Pop_min = params$Pop_min,
    max_iter = 1000
  )

  if(all_clusters_contiguous(cand, adj)){
    return(cand)
  }

sol
}

# ---- Operator 2: Expand smallest cluster by absorbing one border BGU ----

#' Operator 2: grow an undersized market by one BGU
#'
#' Targets a market with low population and annexes the single adjacent BGU
#' that is most strongly tied to it. The gentlest way the search has of fixing
#' a population-floor violation.
#'
#' @details
#' \enumerate{
#'   \item Tournament-select a focal market from the three smallest \code{pop}
#'     values.
#'   \item Collect every BGU adjacent to it but outside it.
#'   \item Score each candidate with \code{\link{compute_ci}()} against the
#'     focal market's members — how much cohesion it would bring in.
#'   \item Tournament-select from the three highest-scoring candidates and move
#'     it in.
#'   \item Accept only if the result is fully contiguous.
#' }
#'
#' Note that step 4 does not check what the move does to the \emph{donor}
#' market: taking a BGU can push the donor below \code{Pop_min}, or break it in
#' two. The contiguity check catches the second case; the first is left to the
#' objective function, which will simply score the candidate 0 and cause the
#' annealing loop to reject it.
#'
#' Returns \code{sol} unchanged if there are fewer than two markets, or if the
#' focal market has no external neighbours.
#'
#' @inheritParams operator_1
#' @return Integer vector of length N, relabelled to \code{1 … m}. Equal to
#'   \code{sol} when the move was rejected.
#'
#' @family perturbation operators
#' @seealso \code{\link{operator_5}()}, the same annexation aimed at low
#'   self-containment rather than low population;
#'   \code{\link{operator_4}()} for the reverse move.
#'
#' @examples
#' adj <- matrix(0L, 6, 6); adj[cbind(1:5, 2:6)] <- 1L; adj <- adj + t(adj)
#' W <- matrix(1, 6, 6); W[1:3, 1:3] <- 30; W[4:6, 4:6] <- 30; diag(W) <- 60
#' row_W <- rowSums(W); col_W <- colSums(W)
#'
#' ## Market 3 is a lone BGU: the operator grows it
#' sol <- c(1, 1, 2, 2, 2, 3)
#' set.seed(1)
#' operator_2(sol, W, adj, row_W, col_W)
#'
#' @export
operator_2 <- function(sol, W, adj, row_W, col_W, ...) {
  lma_df   <- build_lma_df(sol, W, row_W, col_W)
  if (nrow(lma_df) < 2L)
    return(sol)

  idx <- .tournament_select(lma_df$pop, n_way = 3L, minimize = TRUE)

  small_cl <- lma_df$cluster[idx]
  members  <- cluster_members(small_cl, sol)
  adj_bgus <- setdiff(which(colSums(adj[members, , drop = FALSE]) > 0L), members)
  if (length(adj_bgus) == 0L)
    return(sol)

  scores <- vapply(adj_bgus, function(g) {
    compute_ci(g, members, W, row_W, col_W)
  }, numeric(1L))

  b_idx <- .tournament_select(scores, n_way = 3L, minimize = FALSE)
  best_g <- adj_bgus[b_idx]

  candidate <- getOrderedVec(replace(sol, best_g, small_cl))
  .accept_if_contiguous(candidate, sol, adj)
}



# ---- Operator 3: Dismantle largest cluster — reassign BGUs to neighbours ----

#' Operator 3: shed a border BGU from an oversized market
#'
#' Targets a market with many BGUs and tries to give away one of its border
#' BGUs to whichever neighbour offers it the best cohesion. Works through the
#' border in random order and takes the first move that keeps everything
#' contiguous.
#'
#' @details
#' \enumerate{
#'   \item Tournament-select a focal market from the three largest
#'     \code{n_bgus} counts (note: BGU \emph{count}, not population — see
#'     \code{\link{operator_11}()} for the population-driven equivalent).
#'   \item Compute its border BGUs; abort if it has none or only one member.
#'   \item Walk the border in random order. For each BGU, score its adjacent
#'     markets with \code{\link{compute_ci}()} and tournament-select a target
#'     from the three best.
#'   \item Return the first reassignment that leaves every market contiguous.
#' }
#'
#' The retry loop is what distinguishes this from \code{\link{operator_4}()} and
#' \code{\link{operator_8}()}, which test one candidate and give up. On a
#' tightly connected map that makes it noticeably more likely to produce an
#' actual change.
#'
#' @inheritParams operator_1
#' @return Integer vector of length N, relabelled to \code{1 … m}. Equal to
#'   \code{sol} when no border move was feasible.
#'
#' @family perturbation operators
#' @seealso \code{\link{operator_6}()} to split a large market rather than
#'   erode it; \code{\link{operator_11}()} for population-driven splitting.
#'
#' @examples
#' adj <- matrix(0L, 6, 6); adj[cbind(1:5, 2:6)] <- 1L; adj <- adj + t(adj)
#' W <- matrix(1, 6, 6); W[1:3, 1:3] <- 30; W[4:6, 4:6] <- 30; diag(W) <- 60
#' row_W <- rowSums(W); col_W <- colSums(W)
#'
#' ## Market 1 holds four BGUs; the operator sheds one
#' sol <- c(1, 1, 1, 1, 2, 2)
#' set.seed(3)
#' operator_3(sol, W, adj, row_W, col_W)
#'
#' @export
operator_3 <- function(sol, W, adj, row_W, col_W, ...) {
  lma_df  <- build_lma_df(sol, W, row_W, col_W)
  if(nrow(lma_df) < 2L)
    return(sol)

  idx <- .tournament_select(lma_df$n_bgus, n_way = 3L, minimize = FALSE)

  big_cl  <- lma_df$cluster[idx]
  members <- cluster_members(big_cl, sol)
  if (length(members) <= 1L)
    return(sol)

  borders <- .border_bgus(members, sol, adj)
  if(length(borders) == 0L)
    return(sol)

  for(g in sample(borders)){
    adj_cls <- setdiff(unique(sol[which(adj[g, ] > 0L)]), big_cl)
    if(length(adj_cls) == 0L)
      next

    scores <- vapply(adj_cls, function(cl)
      compute_ci(g, cluster_members(cl, sol), W, row_W, col_W), numeric(1L))
    t_idx <- .tournament_select(scores, n_way = 3L, minimize = FALSE)
    best_cl <- adj_cls[t_idx]

    cand <- getOrderedVec(replace(sol, g, best_cl))

    if(all_clusters_contiguous(cand, adj)) {
      return(cand)
    }

  }

  sol
}






# ---- Operator 4: Remove a border BGU from its cluster, assign to best neighbour ----

#' Operator 4: expel a poorly-fitting BGU from a weakly self-contained market
#'
#' Targets a market with low self-containment, finds one of its border BGUs
#' that contributes \emph{least} cohesion to it, and hands that BGU to its best
#' adjacent market. The complement of \code{\link{operator_5}()}.
#'
#' @details
#' \enumerate{
#'   \item Tournament-select a focal market from the three lowest \code{sc}
#'     values.
#'   \item Score each border BGU by \code{\link{compute_ci}()} against the rest
#'     of the market, and tournament-select from the three \emph{lowest} — the
#'     BGUs the market would miss least.
#'   \item Score that BGU against each adjacent market and move it to the
#'     highest-scoring one.
#'   \item Accept only if the result is fully contiguous.
#' }
#'
#' Note the asymmetry between steps 2 and 3: the BGU to expel is chosen
#' stochastically (tournament), but its destination is chosen greedily
#' (\code{which.max}). The exploration in this operator comes entirely from
#' which BGU moves, not where it lands.
#'
#' Returns \code{sol} unchanged if there are fewer than two markets, the focal
#' market has no border BGUs, or the chosen BGU has no adjacent market to move
#' to.
#'
#' @inheritParams operator_1
#' @return Integer vector of length N, relabelled to \code{1 … m}. Equal to
#'   \code{sol} when the move was rejected.
#'
#' @family perturbation operators
#' @seealso \code{\link{operator_5}()} (pull in rather than push out),
#'   \code{\link{operator_8}()} (the same move chosen entirely at random).
#'
#' @examples
#' adj <- matrix(0L, 6, 6); adj[cbind(1:5, 2:6)] <- 1L; adj <- adj + t(adj)
#' W <- matrix(1, 6, 6); W[1:3, 1:3] <- 30; W[4:6, 4:6] <- 30; diag(W) <- 60
#' row_W <- rowSums(W); col_W <- colSums(W)
#'
#' ## Market 1 wrongly holds BGU 4, which belongs with 5 and 6
#' sol <- c(1, 1, 1, 1, 2, 2)
#' set.seed(2)
#' operator_4(sol, W, adj, row_W, col_W)
#'
#' @export
operator_4 <- function(sol, W, adj, row_W, col_W, ...) {
  lma_df   <- build_lma_df(sol, W, row_W, col_W)
  if(nrow(lma_df) < 2L)
    return(sol)

  idx <- .tournament_select(lma_df$sc, n_way = 3L, minimize = TRUE)
  focal_cl <- lma_df$cluster[idx]
  members  <- cluster_members(focal_cl, sol)
  borders  <- .border_bgus(members, sol, adj)
  if (length(borders) == 0L)
    return(sol)

  b_scores <- vapply(borders, function(g)
    compute_ci(g, setdiff(members, g), W, row_W, col_W), numeric(1L))
  b_idx <- .tournament_select(b_scores, n_way = 3L, minimize = TRUE)

  g       <- borders[b_idx]
  adj_cls <- setdiff(unique(sol[which(adj[g, ] > 0L)]), focal_cl)
  if (length(adj_cls) == 0L)
    return(sol)

  scores  <- vapply(adj_cls, function(cl) {
    compute_ci(g, cluster_members(cl, sol), W, row_W, col_W)
  }, numeric(1L))
  t_idx <- .tournament_select(scores, n_way = 3L, minimize = FALSE)
  best_cl <- adj_cls[which.max(scores)]

  candidate <- getOrderedVec(replace(sol, g, best_cl))
  .accept_if_contiguous(candidate, sol, adj)
}

# ---- Operator 5: Include a border BGU from an adjacent cluster ----

#' Operator 5: annex the best-fitting neighbouring BGU
#'
#' Targets a market with low self-containment and pulls in whichever adjacent
#' BGU would contribute the most cohesion to it. The inverse of
#' \code{\link{operator_4}()}: same focal criterion, opposite direction of
#' travel.
#'
#' @details
#' \enumerate{
#'   \item Tournament-select a focal market from the three lowest \code{sc}
#'     values.
#'   \item Collect every BGU adjacent to it but outside it.
#'   \item Score each with \code{\link{compute_ci}()} against the focal
#'     market's members and tournament-select from the three highest.
#'   \item Accept only if the result is fully contiguous.
#' }
#'
#' Structurally identical to \code{\link{operator_2}()}; the only difference is
#' the focal criterion — \code{sc} here, \code{pop} there. Having both means
#' the search can address a self-containment shortfall and a population
#' shortfall with the same kind of move.
#'
#' As with \code{\link{operator_2}()}, the effect on the donor market is not
#' checked beyond contiguity.
#'
#' @inheritParams operator_1
#' @return Integer vector of length N, relabelled to \code{1 … m}. Equal to
#'   \code{sol} when the move was rejected.
#'
#' @family perturbation operators
#' @seealso \code{\link{operator_2}()}, \code{\link{operator_4}()}
#'
#' @examples
#' adj <- matrix(0L, 6, 6); adj[cbind(1:5, 2:6)] <- 1L; adj <- adj + t(adj)
#' W <- matrix(1, 6, 6); W[1:3, 1:3] <- 30; W[4:6, 4:6] <- 30; diag(W) <- 60
#' row_W <- rowSums(W); col_W <- colSums(W)
#'
#' ## Market 2 is split off from BGU 3, which it should absorb
#' sol <- c(1, 1, 1, 2, 2, 2)
#' set.seed(5)
#' operator_5(sol, W, adj, row_W, col_W)
#'
#' @export
operator_5 <- function(sol, W, adj, row_W, col_W, ...) {
  lma_df   <- build_lma_df(sol, W, row_W, col_W)
  if(nrow(lma_df) < 2L)
    return(sol)

  idx <- .tournament_select(lma_df$sc, n_way = 3L, minimize = TRUE)
  focal_cl <- lma_df$cluster[idx]
  members  <- cluster_members(focal_cl, sol)
  adj_bgus <- setdiff(which(colSums(adj[members, , drop = FALSE]) > 0L), members)
  if (length(adj_bgus) == 0L)
    return(sol)

  scores <- vapply(adj_bgus, function(g) {
    compute_ci(g, members, W, row_W, col_W)
  }, numeric(1L))

  b_idx <- .tournament_select(scores, n_way = 3L, minimize = FALSE)
  best_g <- adj_bgus[b_idx]

  candidate <- getOrderedVec(replace(sol, best_g, focal_cl))
  .accept_if_contiguous(candidate, sol, adj)
}

# ---- Operator 6: Segregation — split a contiguous subset into a new cluster ----

#' Operator 6: segregate the least cohesive fringe of a large market
#'
#' Splits one to three of the worst-fitting BGUs out of a large market into a
#' brand new market. This is the operator that \emph{increases} the market
#' count, and it carries double weight in the dispatcher because raising
#' \code{m} is how the search escapes the over-merged solutions that SARA tends
#' to produce.
#'
#' @details
#' \enumerate{
#'   \item Tournament-select a focal market from the three largest
#'     \code{n_bgus} counts; abort if it has fewer than four members.
#'   \item Score every member by \code{\link{compute_ci}()} against the rest of
#'     the market.
#'   \item Try \code{n_split = 1}, then 2, then 3: take that many
#'     lowest-scoring members as the proposed new market.
#'   \item Accept the first \code{n_split} for which \emph{both} the new market
#'     and the remainder are contiguous, giving the new market the fresh label
#'     \code{max(sol) + 1L}.
#' }
#'
#' Selecting by lowest cohesion rather than by geography is what makes the
#' split meaningful: the BGUs peeled off are the ones whose commuting ties to
#' the market are weakest, which are usually exactly the ones that belong
#' somewhere else.
#'
#' The escalating \code{n_split} matters on a sparse map. A single low-cohesion
#' BGU is often interior to its market, so detaching it alone would break the
#' remainder; taking two or three adjacent ones can leave both halves whole.
#' Unlike the other operators this one does not call
#' \code{.accept_if_contiguous()} — it checks the two affected pieces directly,
#' which is cheaper and sufficient, since no other market is touched.
#'
#' @inheritParams operator_1
#' @return Integer vector of length N, relabelled to \code{1 … m}, with one
#'   more market than \code{sol} on success. Equal to \code{sol} when no split
#'   of size 1–3 kept both pieces contiguous.
#'
#' @family perturbation operators
#' @seealso \code{\link{operator_11}()}, which splits by population rather than
#'   cohesion; \code{\link{operator_3}()}, which erodes a market one BGU at a
#'   time instead of splitting it.
#'
#' @examples
#' adj <- matrix(0L, 6, 6); adj[cbind(1:5, 2:6)] <- 1L; adj <- adj + t(adj)
#' W <- matrix(1, 6, 6); W[1:3, 1:3] <- 30; W[4:6, 4:6] <- 30; diag(W) <- 60
#' row_W <- rowSums(W); col_W <- colSums(W)
#'
#' ## One market of five BGUs, plus a singleton
#' sol <- c(1, 1, 1, 1, 1, 2)
#' set.seed(4)
#' out <- operator_6(sol, W, adj, row_W, col_W)
#' out
#' n_clusters(sol); n_clusters(out)
#'
#' ## Markets with fewer than four members are left alone
#' identical(operator_6(c(1, 1, 1, 2, 2, 2), W, adj, row_W, col_W),
#'           c(1, 1, 1, 2, 2, 2))
#'
#' @export
operator_6 <- function(sol, W, adj, row_W, col_W, ...) {
  lma_df  <- build_lma_df(sol, W, row_W, col_W)
  if(nrow(lma_df) < 2L)
    return(sol)

  idx <- .tournament_select(lma_df$n_bgus, n_way = 3L, minimize = FALSE)
  big_cl  <- lma_df$cluster[idx]
  members <- cluster_members(big_cl, sol)
  if (length(members) < 4L)
    return(sol)

  ci_scores <- vapply(members, function(g){
    M_g <- setdiff(members, g)
    compute_ci(g, M_g, W, row_W, col_W)
  }, numeric(1L))

  for(n_split in c(1L, 2L, 3L)) {
    if(n_split >= length(members))
      next
    subset_g <- members[order(ci_scores)[seq_len(n_split)]]
    remainder <- setdiff(members, subset_g)

    if(is_cluster_contiguous(subset_g, adj) &&
       is_cluster_contiguous(remainder, adj)){
      new_cl <- max(sol) + 1L
      cand <- getOrderedVec(replace(sol, subset_g, new_cl))
      return(cand)
    }
  }
  sol

}

# ---- Operator 7: Annexation — absorb a weakly connected adjacent cluster ----

#' Operator 7: absorb the most weakly connected neighbouring market
#'
#' Picks a market at random and swallows whichever of its neighbours it
#' interacts with \emph{least}. A pure merge: no re-aggregation, no splitting,
#' market count drops by exactly one.
#'
#' @details
#' \enumerate{
#'   \item Sample a focal market uniformly at random.
#'   \item Score its adjacent markets with \code{.cluster_interaction()}.
#'   \item Absorb the \emph{lowest}-scoring one entirely.
#'   \item Accept only if the result is fully contiguous.
#' }
#'
#' Taking the weakest neighbour rather than the strongest is deliberate. The
#' strongest-neighbour merge is what SARA already did during construction, so
#' repeating it here would mostly re-propose moves the search has already
#' settled on. Merging a weak neighbour is a genuinely disruptive move that
#' hands the other operators a large, badly-formed market to re-carve — which
#' is the point of a perturbation.
#'
#' @section Assumption about labels:
#' The focal market is taken as the sampled \emph{index} into the sorted label
#' vector rather than the label at that index. These coincide because every
#' partition in the search has been relabelled to \code{1 … m} by
#' \code{\link{get_ordered_vec}()}. If you call this operator on a partition
#' with gapped or non-consecutive labels (say \code{c(1, 1, 7, 7)}), it will
#' address a market that does not exist and return \code{sol} unchanged.
#' Relabel first.
#'
#' @inheritParams operator_1
#' @return Integer vector of length N, relabelled to \code{1 … m}, with one
#'   fewer market on success. Equal to \code{sol} when the move was rejected.
#'
#' @family perturbation operators
#' @seealso \code{\link{operator_1}()} for a merge followed by re-aggregation;
#'   \code{\link{operator_6}()} for the inverse (splitting) move.
#'
#' @examples
#' adj <- matrix(0L, 6, 6); adj[cbind(1:5, 2:6)] <- 1L; adj <- adj + t(adj)
#' W <- matrix(1, 6, 6); W[1:3, 1:3] <- 30; W[4:6, 4:6] <- 30; diag(W) <- 60
#' row_W <- rowSums(W); col_W <- colSums(W)
#'
#' sol <- c(1, 1, 2, 2, 3, 3)
#' set.seed(6)
#' out <- operator_7(sol, W, adj, row_W, col_W)
#' out
#' n_clusters(sol); n_clusters(out)
#'
#' ## Relabel first if your partition has gapped labels
#' operator_7(get_ordered_vec(c(1, 1, 7, 7, 9, 9)), W, adj, row_W, col_W)
#'
#' @export
operator_7 <- function(sol, W, adj, row_W, col_W, ...) {
  clusters <- unique_clusters(sol)
  if (length(clusters) < 2L)
    return(sol)

  fcl_idx <- sample(length(clusters), 1L)
  focal_cl  <- fcl_idx
  focal_mbrs <- cluster_members(focal_cl, sol)
  adj_cls   <- .adjacent_clusters(focal_mbrs, sol, adj)
  if (length(adj_cls) == 0L)
    return(sol)

  scores   <- vapply(adj_cls, function(cl) {
    .cluster_interaction(focal_mbrs, cluster_members(cl, sol), W, row_W, col_W)
  }, numeric(1L))
  weak_cl  <- adj_cls[which.min(scores)]

  candidate <- getOrderedVec(replace(sol, cluster_members(weak_cl, sol), focal_cl))
  .accept_if_contiguous(candidate, sol, adj)
}

# ---- Operator 8: Random border reassignment ----

#' Operator 8: move one border BGU completely at random
#'
#' The unguided move: a random market, a random border BGU of it, a random
#' adjacent market to send it to. No scoring of any kind.
#'
#' @details
#' \enumerate{
#'   \item Sample a market uniformly from \code{\link{unique_clusters}(sol)}.
#'   \item Sample one of its border BGUs uniformly.
#'   \item Sample one of that BGU's adjacent markets uniformly and move it
#'     there.
#'   \item Accept only if the result is fully contiguous.
#' }
#'
#' Every other operator is greedy in at least one dimension, which means they
#' all propose moves correlated with the current cohesion landscape and can
#' collectively stall in the same basin. This one proposes moves no score would
#' ever suggest, and at high temperature the annealing loop accepts a fair
#' number of them. It is the operator set's source of unbiased exploration, and
#' it is also the cheapest — no \code{\link{build_lma_df}()} call at all.
#'
#' @inheritParams operator_1
#' @return Integer vector of length N, relabelled to \code{1 … m}. Equal to
#'   \code{sol} when the sampled market had no border BGU, the sampled BGU had
#'   no other adjacent market, or the move broke contiguity.
#'
#' @family perturbation operators
#' @seealso \code{\link{operator_4}()} for the scored version of the same move;
#'   \code{\link{operator_9}()} for a random two-way swap.
#'
#' @examples
#' adj <- matrix(0L, 6, 6); adj[cbind(1:5, 2:6)] <- 1L; adj <- adj + t(adj)
#' W <- matrix(1, 6, 6); W[1:3, 1:3] <- 30; W[4:6, 4:6] <- 30; diag(W) <- 60
#' row_W <- rowSums(W); col_W <- colSums(W)
#' sol <- c(1, 1, 1, 2, 2, 2)
#'
#' ## Different seeds give different random moves
#' set.seed(1); operator_8(sol, W, adj, row_W, col_W)
#' set.seed(9); operator_8(sol, W, adj, row_W, col_W)
#'
#' @export
operator_8 <- function(sol, W, adj, row_W, col_W, ...) {
  focal_cl <- sample(unique_clusters(sol), 1L)
  members  <- cluster_members(focal_cl, sol)
  borders  <- .border_bgus(members, sol, adj)
  if (length(borders) == 0L)
    return(sol)

  g       <- sample(borders, 1L)
  adj_cls <- setdiff(unique(sol[which(adj[g, ] > 0L)]), focal_cl)
  if (length(adj_cls) == 0L)
    return(sol)

  target_cl <- sample(adj_cls, 1L)
  candidate <- getOrderedVec(replace(sol, g, target_cl))
  .accept_if_contiguous(candidate, sol, adj)
}

# ---- Operator 9: Exchange border BGUs between two adjacent clusters ----

#' Operator 9: swap border BGUs between two adjacent markets
#'
#' Exchanges one BGU from market A that touches market B with one BGU from B
#' that touches A. The only operator that changes neither the number of markets
#' nor any market's BGU count.
#'
#' @details
#' \enumerate{
#'   \item Sample a market \code{cl1}, then sample one of its neighbours
#'     \code{cl2}.
#'   \item Restrict each market's border BGUs to those actually touching the
#'     other market; abort if either side comes up empty.
#'   \item Sample one BGU from each side and exchange their labels.
#'   \item Accept only if the result is fully contiguous.
#' }
#'
#' Size preservation is what makes this operator useful. Every other move
#' changes at least one market's population, so once the search has found a
#' partition that satisfies \code{Pop_min} and \code{Pop_max} everywhere, most
#' proposals risk pushing some market back outside the bounds and being
#' rejected by \code{\link{compute_objective}()}. A swap of two BGUs of similar
#' size leaves the population profile nearly intact while still reshaping the
#' boundary, so it keeps exploring where the others stall — the classic
#' late-stage refinement move.
#'
#' Both endpoints are sampled uniformly; there is no cohesion scoring.
#'
#' @inheritParams operator_1
#' @return Integer vector of length N, relabelled to \code{1 … m}. Equal to
#'   \code{sol} when fewer than two markets exist, no mutual border was found,
#'   or the swap broke contiguity.
#'
#' @family perturbation operators
#' @seealso \code{\link{operator_8}()}, \code{\link{operator_10}()}
#'
#' @examples
#' adj <- matrix(0L, 6, 6); adj[cbind(1:5, 2:6)] <- 1L; adj <- adj + t(adj)
#' W <- matrix(1, 6, 6); W[1:3, 1:3] <- 30; W[4:6, 4:6] <- 30; diag(W) <- 60
#' row_W <- rowSums(W); col_W <- colSums(W)
#'
#' ## On a line graph a swap almost always breaks contiguity, so the
#' ## solution usually comes back unchanged -- that is the operator
#' ## declining the move, not an error.
#' sol <- c(1, 1, 1, 2, 2, 2)
#' set.seed(8)
#' operator_9(sol, W, adj, row_W, col_W)
#'
#' ## Market sizes are preserved whenever a swap is accepted
#' out <- operator_9(sol, W, adj, row_W, col_W)
#' identical(sort(table(sol)), sort(table(out)))
#'
#' @export
operator_9 <- function(sol, W, adj, row_W, col_W, ...) {
  clusters <- unique_clusters(sol)
  if (length(clusters) < 2L)
    return(sol)

  cl1      <- sample(clusters, 1L)
  adj_cls  <- .adjacent_clusters(cluster_members(cl1, sol), sol, adj)
  if (length(adj_cls) == 0L)
    return(sol)
  cl2      <- sample(adj_cls, 1L)

  mbrs1    <- cluster_members(cl1, sol)
  mbrs2    <- cluster_members(cl2, sol)

  b1_to_2 <- .border_bgus(mbrs1, sol, adj)
  b1_to_2 <- b1_to_2[vapply(b1_to_2, function(g)
    any(adj[g, mbrs2] > 0L), logical(1L))]

  b2_to_1 <- .border_bgus(mbrs2, sol, adj)
  b2_to_1 <- b2_to_1[vapply(b2_to_1, function(g)
    any(adj[g, mbrs1] > 0L), logical(1L))]

  if (length(b1_to_2) == 0L || length(b2_to_1) == 0L)
    return(sol)

  g1 <- sample(b1_to_2, 1L)
  g2 <- sample(b2_to_1, 1L)

  candidate        <- sol
  candidate[g1]    <- cl2
  candidate[g2]    <- cl1
  candidate        <- getOrderedVec(candidate)
  .accept_if_contiguous(candidate, sol, adj)
}

# ---- Operator 10: Transfer a small subset between adjacent clusters ----

# Helper function for operator 10

#' Repair a partition by splitting non-contiguous markets
#'
#' Takes a partition in which some markets may have become disconnected and
#' returns one in which every market is a single connected piece, by giving each
#' extra fragment its own new label. Nothing is merged and no BGU changes
#' neighbourhood — only labels change.
#'
#' @details
#' For each market, a BFS (the same one as
#' \code{\link{is_cluster_contiguous}()}) peels off connected components one at
#' a time. The first component keeps the original label; every subsequent one
#' gets \code{max(sol) + 1L}. A market that is already contiguous, or has one
#' member, is skipped without cost.
#'
#' The market count can therefore only rise, and it rises by exactly the number
#' of extra fragments found.
#'
#' @section Why this exists:
#' Every operator except \code{\link{operator_10}()} proposes a move and
#' \emph{rejects} it if contiguity breaks. That is not viable for
#' \code{\link{operator_10}()}'s crossover branch, which transplants whole
#' markets from a different partition and so almost always leaves fragments
#' behind — rejecting on contiguity would reject nearly every crossover. This
#' function is the alternative: accept the transplant, then repair it. It is
#' the only place in the package where a partition is fixed rather than
#' discarded.
#'
#' \code{max_iter} guards the per-market component loop; reaching it leaves that
#' market partially repaired, so a pathological adjacency matrix degrades
#' rather than hangs.
#'
#' @param sol Integer vector of market assignments (length N).
#' @param adj Integer/logical \code{N x N} adjacency matrix.
#' @param max_iter Integer, cap on component-peeling iterations per market
#'   (default 5000).
#' @return Integer vector of length N, relabelled to \code{1 … m}, in which
#'   \code{\link{all_clusters_contiguous}()} is \code{TRUE}.
#'
#' @family perturbation operators
#' @seealso \code{\link{all_clusters_contiguous}()} to verify the result,
#'   \code{\link{operator_10}()} for the caller.
#'
#' @examples
#' adj <- matrix(0L, 6, 6); adj[cbind(1:5, 2:6)] <- 1L; adj <- adj + t(adj)
#'
#' ## Market 1 is in two pieces: {1} and {3}
#' broken <- c(1, 2, 1, 2, 2, 2)
#' all_clusters_contiguous(broken, adj)
#'
#' fixed <- fix_contiguity(broken, adj)
#' fixed
#' all_clusters_contiguous(fixed, adj)
#'
#' ## An already-contiguous partition is returned unchanged (up to relabelling)
#' fix_contiguity(c(1, 1, 1, 2, 2, 2), adj)
#'
#' @export
fix_contiguity <- function(sol, adj, max_iter = 5000L){
  clusters <- unique_clusters(sol)

  for(cl in clusters){
    members <- cluster_members(cl, sol)
    if(length(members) <= 1L) next
    if(is_cluster_contiguous(members, adj)) next

    remaining <- members
    first <- TRUE
    iter <- 0L

    while (length(remaining) > 0L) {
      iter <- iter + 1L
      if(iter > max_iter){
        break
      }
      seed <- remaining[1L]
      visited <- seed
      frontier <- seed

      while(length(frontier) > 0L){
        candidates <- which(colSums(adj[frontier, , drop = FALSE]) > 0L)
        new_nodes <- setdiff(intersect(candidates, remaining), visited)
        visited <- c(visited, new_nodes)
        frontier <- new_nodes
      }

      if(!first){
        new_cl <- max(sol) + 1L
        sol[visited] <- new_cl
      }

      remaining <- setdiff(remaining, visited)
      first <- FALSE
    }
  }
  getOrderedVec(sol)
}

#' Operator 10: transfer a group of BGUs, or cross over with the elite solution
#'
#' Two operators in one. Without a reference solution it moves a small group of
#' border BGUs between adjacent markets. With one — the search's best solution
#' so far — it performs a \emph{crossover}: it transplants several of that
#' solution's markets wholesale into the current one and repairs the damage.
#' This is the only operator that uses information from outside the current
#' partition.
#'
#' @details
#' \strong{Transfer branch} (\code{best_sol} is \code{NULL} or identical to
#' \code{sol}): sample a market, sample an adjacent market, take one to three of
#' the border BGUs facing it, and move them across as a group. Accept only if
#' the result is fully contiguous. This is essentially
#' \code{\link{operator_8}()} operating on a small block rather than a single
#' BGU.
#'
#' \strong{Crossover branch} (otherwise):
#' \enumerate{
#'   \item Choose \code{r} markets at random from \code{best_sol}, where
#'     \code{r} is between 1 and 6 percent of that solution's market count.
#'   \item Overwrite the current labels of every BGU in those markets with
#'     fresh labels, so that \code{best_sol}'s grouping of those BGUs is
#'     imposed on \code{sol} intact.
#'   \item Repair the fragments this leaves behind in the surrounding markets
#'     with \code{\link{fix_contiguity}()}.
#'   \item Re-aggregate with \code{\link{sha}()} to fold the resulting slivers
#'     back into valid markets.
#'   \item Accept only if everything is contiguous.
#' }
#'
#' @section Cost control:
#' Step 4 is by far the most expensive thing any operator does, so it is
#' throttled. About 70 percent of the time the re-aggregation is restricted to
#' the transplanted BGUs plus their immediate neighbours, with
#' \code{max_iter = 1000}; the rest of the time it runs over the whole map with
#' \code{max_iter = 3000}. This is also why the dispatcher gives operator 10 the
#' lowest weight in the set (0.5 against a baseline of 1).
#'
#' @section Global variables:
#' The crossover branch reads \code{params$SC_min} and \code{params$Pop_min}
#' from the global environment to configure \code{\link{sha}()}; see
#' \code{\link{operator_1}()} for the same caveat. The transfer branch does
#' not, so calling with \code{best_sol = NULL} needs no globals.
#'
#' @param sol Integer vector of market assignments (length N).
#' @param W Numeric \code{N x N} origin-destination matrix.
#' @param adj Integer/logical \code{N x N} adjacency matrix.
#' @param row_W Numeric vector, \code{rowSums(W)}.
#' @param col_W Numeric vector, \code{colSums(W)}.
#' @param best_sol Integer vector or \code{NULL} (default). The partition to
#'   cross over with — in the search, the incumbent best. \code{NULL} or a
#'   partition identical to \code{sol} selects the transfer branch.
#' @param ... Ignored; present for signature compatibility.
#' @return Integer vector of length N, relabelled to \code{1 … m}. Equal to
#'   \code{sol} when the move was rejected.
#'
#' @family perturbation operators
#' @seealso \code{\link{fix_contiguity}()} and \code{\link{sha}()} for the
#'   repair steps; \code{\link{apply_random_operator}()}, which supplies
#'   \code{best_sol}.
#'
#' @references
#' Martínez-Bernabeu, L., Flórez-Revuelta, F. and Casado-Díaz, J. M. (2012)
#' Grouping genetic operators for the delineation of functional areas based on
#' spatial interaction. \emph{Expert Systems with Applications}.
#'
#' @examples
#' adj <- matrix(0L, 6, 6); adj[cbind(1:5, 2:6)] <- 1L; adj <- adj + t(adj)
#' W <- matrix(1, 6, 6); W[1:3, 1:3] <- 30; W[4:6, 4:6] <- 30; diag(W) <- 60
#' row_W <- rowSums(W); col_W <- colSums(W)
#' sol <- c(1, 1, 1, 2, 2, 2)
#'
#' ## Transfer branch: no globals needed
#' set.seed(11)
#' operator_10(sol, W, adj, row_W, col_W, best_sol = NULL)
#'
#' ## Crossover branch against a different incumbent
#' params <- list(SC_min = 0.55, Pop_min = 100, Pop_max = 1000)
#' set.seed(12)
#' operator_10(sol, W, adj, row_W, col_W, best_sol = c(1, 1, 2, 2, 3, 3))
#'
#' @export
operator_10 <- function(sol, W, adj, row_W, col_W, best_sol = NULL, ...){

  low_time_con = TRUE

  if(is.null(best_sol) || identical(sol, best_sol)){
    # simple transfer
    clusters <- unique_clusters(sol)
    if(length(clusters) < 2L) return(sol)
    cl1 <- sample(clusters, 1L)
    adj_cls <- .adjacent_clusters(cluster_members(cl1, sol), sol, adj)
    if(length(adj_cls) == 0L) return(sol)
    cl2 <- sample(adj_cls, 1L)
    borders <- .border_bgus(cluster_members(cl1, sol), sol, adj)
    adj_b <- borders[vapply(borders, function(g) any(adj[g, cluster_members(cl2, sol)] > 0L), logical(1L))]
    if(length(adj_b) == 0L) return(sol)
    move_g <- sample(adj_b, sample(1L:min(3L, length(adj_b)), 1L))
    cand <- getOrderedVec(replace(sol, move_g, cl2))

    return(.accept_if_contiguous(cand, sol, adj))
  }

  # Crossover with sols

  N <- length(sol)
  n_cls_2 <- n_clusters(best_sol)

  r_min <- max(1L, round(n_cls_2*0.01))
  r_max <- max(1L, round(n_cls_2 * 0.06))

  r <- sample(r_min:r_max, 1L)

  cls_2 <- unique_clusters(best_sol)

  selected <- sample(cls_2, min(r, length(cls_2)))

  child <- sol

  offset <- max(child) + 1L

  for(cl in selected){
    mbrs <- cluster_members(cl, best_sol)
    child[mbrs] <- cl + offset
  }

  child <- getOrderedVec(child)

  child <- fix_contiguity(child, adj)

  all_bgus <- seq_len(N)

  affected_bgus <- all_bgus

  max_iter <- 3000L

  if(low_time_con && runif(1L) > 0.3){
    max_iter <- 1000L
    affected_bgus <- unique(c(unlist(lapply(selected, function(cl) cluster_members(cl, best_sol))),
                              unlist(lapply(selected, function(cl){
                                mbrs <- cluster_members(cl, best_sol)
                                setdiff(which(colSums(adj[mbrs,,drop = FALSE]) > 0L), mbrs)
                              }))))
  }

  child <- sha(
    bgus = affected_bgus,
    sol = child,
    W = W,
    adj = adj,
    row_W = row_W,
    col_W = col_W,
    SC_min = params$SC_min,
    Pop_min = params$Pop_min,
    max_iter = max_iter
  )

  if(all_clusters_contiguous(child, adj)){
    return(child)
  }
  sol
}


#' Operator 11: break up a market that exceeds the population ceiling
#'
#' The dedicated repair move for \code{Pop_max} violations. Finds a market over
#' the ceiling, shatters it into singletons, and rebuilds that region from
#' scratch with \code{\link{sha}()} — which will generally produce several
#' smaller markets in its place.
#'
#' @details
#' \enumerate{
#'   \item \code{\link{build_lma_df}()} and select the markets with
#'     \code{pop > Pop_max}. If there are none, return \code{sol} unchanged
#'     immediately — this operator is a no-op on a feasible partition.
#'   \item Tournament-select among the largest offenders (implemented as a
#'     minimising tournament on negated populations).
#'   \item Abort if that market has fewer than three BGUs, since there is no
#'     useful way to split it.
#'   \item Give every one of its BGUs a fresh label, then re-aggregate exactly
#'     those BGUs with \code{\link{sha}()} (\code{max_iter = 500}).
#'   \item Accept only if everything is contiguous.
#' }
#'
#' This is the annealing-time counterpart of what
#' \code{\link{sara_iterated}()} does between passes, and it exists because
#' \code{\link{compute_objective}()} scores an over-ceiling partition 0: without
#' a move that specifically targets the violation, the search would have no
#' gradient to climb back out. It carries a weight of 1.5 in the dispatcher,
#' the second-highest in the set.
#'
#' @section Global variables:
#' Reads \code{params$SC_min} and \code{params$Pop_min} from the global
#' environment; see \code{\link{operator_1}()}. Note that \code{Pop_max} is a
#' proper argument here, but the dispatcher supplies it as
#' \code{params$Pop_max}, so in practice all three come from the same place.
#'
#' @param sol Integer vector of market assignments (length N).
#' @param W Numeric \code{N x N} origin-destination matrix.
#' @param adj Integer/logical \code{N x N} adjacency matrix.
#' @param row_W Numeric vector, \code{rowSums(W)}.
#' @param col_W Numeric vector, \code{colSums(W)}.
#' @param Pop_max Numeric, the population ceiling (default \code{Inf}, which
#'   makes the operator a no-op).
#' @param ... Ignored; present for signature compatibility.
#' @return Integer vector of length N, relabelled to \code{1 … m}, usually with
#'   more markets than \code{sol}. Equal to \code{sol} when no market exceeded
#'   the ceiling or the split was rejected.
#'
#' @family perturbation operators
#' @seealso \code{\link{sara_iterated}()} for the same repair at construction
#'   time; \code{\link{penalty_pop_max_single}()} and
#'   \code{\link{sara_max_pop_score}()} for measuring the violation;
#'   \code{\link{operator_6}()} for cohesion-driven splitting.
#'
#' @examples
#' adj <- matrix(0L, 6, 6); adj[cbind(1:5, 2:6)] <- 1L; adj <- adj + t(adj)
#' W <- matrix(1, 6, 6); W[1:3, 1:3] <- 30; W[4:6, 4:6] <- 30; diag(W) <- 60
#' row_W <- rowSums(W); col_W <- colSums(W)
#' params <- list(SC_min = 0.55, Pop_min = 100, Pop_max = 400)
#'
#' ## One market of all six BGUs (population 738) against a ceiling of 400
#' set.seed(13)
#' out <- operator_11(rep(1, 6), W, adj, row_W, col_W, Pop_max = 400)
#' build_lma_df(out, W, row_W, col_W)
#'
#' ## No violation, no change
#' identical(operator_11(c(1, 1, 1, 2, 2, 2), W, adj, row_W, col_W,
#'                       Pop_max = 400),
#'           c(1, 1, 1, 2, 2, 2))
#'
#' @export
operator_11 <- function(sol, W, adj, row_W, col_W, Pop_max = Inf, ...){

  lma_df <- build_lma_df(sol, W, row_W, col_W)

  big_cls <- lma_df$cluster[lma_df$pop > Pop_max]

  if(length(big_cls) == 0L) return(sol)

  big_pops <- lma_df$pop[match(big_cls, lma_df$cluster)]
  idx <- .tournament_select(-big_pops, n_way = 3L, minimize  = TRUE)

  big_cl <- big_cls[idx]
  big_mbrs <- cluster_members(big_cl, sol)

  if(length(big_mbrs) < 3L) return(sol)

  cand <- sol
  for(i in seq_along(big_mbrs)){
    new_id <- max(cand) + 1L
    cand[big_mbrs[i]] <- new_id
  }

  cand <- getOrderedVec(cand)

  cand <- sha(
    bgus = big_mbrs, sol = cand, W = W, adj = adj, row_W = row_W, col_W = col_W, SC_min = params$SC_min,
    Pop_min = params$Pop_min, max_iter = 500L)

  if(all_clusters_contiguous(cand, adj)){
    return(cand)
  }

  sol
}

#' Create an empty operator-performance tracker
#'
#' Allocates the counter structure that \code{\link{adsa_almd}()} threads
#' through its inner loop to record how each of the eleven operators performed.
#'
#' @details
#' Four integer vectors of length 11, one slot per operator:
#' \itemize{
#'   \item \code{called} — times the dispatcher selected it;
#'   \item \code{unchanged} — times it returned the solution untouched (the
#'     move was infeasible);
#'   \item \code{accepted} — times a changed solution was accepted by the
#'     Metropolis criterion;
#'   \item \code{improved} — times it produced a new global best.
#' }
#' \code{called - unchanged} is the number of real proposals, which is the
#' denominator you usually want when judging an operator.
#'
#' The structure is a plain list, so it survives \code{saveRDS()} and can be
#' pooled across runs by adding the vectors elementwise.
#'
#' @return Named list of four integer vectors of length 11, all zero.
#'
#' @family operator diagnostics
#' @seealso \code{\link{update_op_tracker}()} to record an attempt,
#'   \code{\link{print_op_tracker}()} to display the result.
#'
#' @examples
#' tracker <- init_op_tracker()
#' str(tracker)
#'
#' @export
init_op_tracker <- function(){
  list(
    called = rep(0L, 11L),
    accepted = rep(0L, 11L),
    improved = rep(0L, 11L),
    unchanged = rep(0L, 11L)
  )
}

#' Record one operator attempt in a tracker
#'
#' Updates the counters for a single proposal. Called once per inner-loop
#' iteration by \code{\link{adsa_almd}()}.
#'
#' @details
#' \code{called} always increments. If \code{sol_before} and \code{sol_after}
#' are \code{identical()}, \code{unchanged} increments and the function returns
#' \strong{immediately} — an operator that declined to move is never counted as
#' accepted, even though the annealing loop technically "accepts" the identical
#' solution. Only for a real change are \code{accepted} and \code{improved}
#' considered.
#'
#' \code{F_before} and \code{F_after} are accepted but not currently recorded;
#' they are in the signature so fitness-delta statistics can be added without
#' changing any call site.
#'
#' Trackers are immutable in the R sense: the updated list is returned and must
#' be assigned back.
#'
#' @param tracker List from \code{\link{init_op_tracker}()}.
#' @param op_idx Integer in \code{1:11}, which operator was applied.
#' @param sol_before,sol_after Integer vectors, the partitions before and after
#'   the operator ran. Compared with \code{identical()}.
#' @param F_before,F_after Numeric objective values. Accepted but not used.
#' @param accepted Logical, whether the annealing loop accepted the candidate.
#' @param improved_best Logical, whether the candidate became the new global
#'   best.
#' @return The updated tracker list.
#'
#' @family operator diagnostics
#' @seealso \code{\link{init_op_tracker}()}, \code{\link{print_op_tracker}()}
#'
#' @examples
#' tracker <- init_op_tracker()
#'
#' ## A real, accepted, improving move by operator 6
#' tracker <- update_op_tracker(tracker, op_idx = 6,
#'                              sol_before = c(1, 1, 2, 2),
#'                              sol_after  = c(1, 2, 3, 3),
#'                              F_before = 1.0, F_after = 1.4,
#'                              accepted = TRUE, improved_best = TRUE)
#'
#' ## Operator 9 declined to move: counted as unchanged only
#' tracker <- update_op_tracker(tracker, op_idx = 9,
#'                              sol_before = c(1, 1, 2, 2),
#'                              sol_after  = c(1, 1, 2, 2),
#'                              F_before = 1.4, F_after = 1.4,
#'                              accepted = TRUE, improved_best = FALSE)
#'
#' data.frame(op = 1:11, called = tracker$called,
#'            unchanged = tracker$unchanged, accepted = tracker$accepted)
#'
#' @export
update_op_tracker <- function(tracker,
                              op_idx,
                              sol_before,
                              sol_after,
                              F_before,
                              F_after,
                              accepted,
                              improved_best){
  tracker$called[op_idx] <- tracker$called[op_idx] + 1L

  if(identical(sol_before, sol_after)){
    tracker$unchanged[op_idx] <- tracker$unchanged[op_idx] + 1L

    return(tracker)
  }

  if(accepted){
    tracker$accepted[op_idx] <- tracker$accepted[op_idx] + 1L
  }

  if(improved_best){
    tracker$improved[op_idx] <- tracker$improved[op_idx] + 1L
  }

  tracker
}


#' Print an operator-performance report
#'
#' Writes a formatted table of the tracker counters to the console, one row per
#' operator plus a TOTAL row.
#'
#' @details
#' Columns are Called, Unchanged, Changed (\code{called - unchanged}),
#' Accepted, Improved and Acc\% (\code{100 * accepted / called}).
#'
#' What to look for when tuning:
#' \itemize{
#'   \item A high \strong{Unchanged} share means the operator keeps proposing
#'     infeasible moves on your map — common for \code{\link{operator_9}()} on
#'     sparse adjacency, and for \code{\link{operator_11}()} once the population
#'     ceiling is satisfied.
#'   \item \strong{Acc\%} is computed over \code{Called}, not \code{Changed}, so
#'     an operator that mostly declines to move shows a low acceptance rate for
#'     that reason alone. Compare \code{Accepted / Changed} by hand before
#'     concluding the operator proposes bad moves.
#'   \item \strong{Improved} counts new global bests, which is the column that
#'     actually justifies an operator's weight in
#'     \code{\link{apply_random_operator}()}.
#' }
#'
#' @param tracker List from \code{\link{init_op_tracker}()}, typically
#'   \code{result$op_tracker} after a run.
#' @return Invisibly \code{NULL}; called for the printed output.
#'
#' @family operator diagnostics
#' @seealso \code{\link{init_op_tracker}()}, \code{\link{update_op_tracker}()},
#'   \code{\link{adsa_almd}()}, which returns a tracker as
#'   \code{$op_tracker}.
#'
#' @examples
#' tracker <- init_op_tracker()
#' set.seed(1)
#' tracker$called    <- sample(20:60, 11, replace = TRUE)
#' tracker$unchanged <- pmin(tracker$called, sample(0:20, 11, replace = TRUE))
#' tracker$accepted  <- pmin(tracker$called - tracker$unchanged,
#'                           sample(0:25, 11, replace = TRUE))
#' tracker$improved  <- pmin(tracker$accepted, sample(0:5, 11, replace = TRUE))
#'
#' print_op_tracker(tracker)
#'
#' ## Acceptance rate over real proposals rather than over all calls
#' changed <- tracker$called - tracker$unchanged
#' round(100 * tracker$accepted / pmax(changed, 1), 1)
#'
#' @export
print_op_tracker <- function(tracker){
  cat("\n=== Operator Performance Report ===\n")
  cat(
    sprintf(
      "%-12s %8s %8s %8s %8s %8s %8s \n",
      "Operator",
      "Called",
      "Unchanged",
      "Changed",
      "Accepted",
      "Improved",
      "Acc%"
    )
  )
  cat(paste(rep("-", 68), collapse = ""), "\n")

  for(i in 1:11){
    called <- tracker$called[i]
    unchanged <- tracker$unchanged[i]
    changed <- called - unchanged
    accepted <- tracker$accepted[i]
    improved <- tracker$improved[i]
    acc_pct <- if (called > 0)
      round(100 * accepted / called, 1)
    else
      0

    cat(
      sprintf(
        "%-12s %8s %8s %8s %8s %8s %8s \n",
        paste0("operator_", i),
        called,
        unchanged,
        changed,
        accepted,
        improved,
        acc_pct
      )
    )
  }
  cat(paste(rep("-", 68),collapse = ""), "\n")
  cat(
    sprintf(
      "%-12s %8s %8s %8s %8s %8s %8s \n",
      "TOTAL",
      sum(tracker$called),
      sum(tracker$unchanged),
      sum(tracker$called) - sum(tracker$unchanged),
      sum(tracker$accepted),
      sum(tracker$improved),
      round(100*sum(tracker$accepted) / max(sum(tracker$called), 1), 1)
    )
  )
  cat("======================================\n\n")
}
# ---- Dispatcher ---------------------------------------------------

#' Apply one randomly chosen perturbation operator
#'
#' The single entry point the annealing loop uses to generate a candidate
#' solution. Samples one of the eleven operators according to fixed weights,
#' applies it, and reports back both the new partition and which operator
#' produced it.
#'
#' @details
#' Operators are sampled with these weights:
#' \tabular{lll}{
#'   \strong{Operator} \tab \strong{Weight} \tab \strong{Why} \cr
#'   1–5, 7–9 \tab 1.0 \tab baseline \cr
#'   6 (segregate) \tab 2.0 \tab raising the market count is how the search
#'     escapes SARA's over-merged starting point \cr
#'   10 (crossover) \tab 0.5 \tab by far the most expensive move \cr
#'   11 (split oversized) \tab 1.5 \tab the only way out of a \code{Pop_max}
#'     violation, which scores 0
#' }
#' Weights are normalised internally, so their absolute values do not matter —
#' only their ratios. They are hard-coded in the function body; edit there to
#' retune.
#'
#' Three operators need arguments the others do not, so dispatch is not uniform:
#' \code{\link{operator_10}()} receives \code{best_sol}, and
#' \code{\link{operator_11}()} receives \code{Pop_max = params$Pop_max}. The
#' remaining nine are called through a list with the common signature and
#' \code{...}.
#'
#' A returned partition identical to \code{sol} is normal and expected — it
#' means the selected operator found no feasible move. \code{\link{update_op_tracker}()}
#' distinguishes that case from a genuine proposal.
#'
#' @section Global variables:
#' Reads \code{params$Pop_max} directly, and the operators it dispatches to may
#' read \code{params$SC_min} and \code{params$Pop_min}. A \code{params} list
#' must exist in the global environment; see \code{\link{operator_1}()}.
#'
#' @param sol Integer vector of market assignments (length N).
#' @param W Numeric \code{N x N} origin-destination matrix.
#' @param adj Integer/logical \code{N x N} adjacency matrix.
#' @param row_W Numeric vector, \code{rowSums(W)}.
#' @param col_W Numeric vector, \code{colSums(W)}.
#' @param best_sol Integer vector or \code{NULL} (default), the incumbent best
#'   partition, forwarded to \code{\link{operator_10}()} for crossover.
#' @param ... Passed through to the selected operator.
#' @return Named list:
#'   \describe{
#'     \item{\code{sol}}{integer vector, the perturbed partition. May equal the
#'       input.}
#'     \item{\code{op_idx}}{integer in \code{1:11}, which operator ran.}
#'   }
#'
#' @family perturbation operators
#' @seealso \code{\link{adsa_almd}()} for the search loop,
#'   \code{\link{init_op_tracker}()} for recording what the operators achieve,
#'   \code{\link{estimate_t0}()}, which uses this to sample fitness deltas.
#'
#' @examples
#' adj <- matrix(0L, 6, 6); adj[cbind(1:5, 2:6)] <- 1L; adj <- adj + t(adj)
#' W <- matrix(1, 6, 6); W[1:3, 1:3] <- 30; W[4:6, 4:6] <- 30; diag(W) <- 60
#' row_W <- rowSums(W); col_W <- colSums(W)
#' sol <- c(1, 1, 1, 2, 2, 2)
#'
#' ## The dispatcher and its operators read a global `params`
#' params <- list(SC_min = 0.55, Pop_min = 100, Pop_max = 1000)
#'
#' set.seed(21)
#' res <- apply_random_operator(sol, W, adj, row_W, col_W)
#' res$op_idx
#' res$sol
#'
#' ## How often each operator gets picked, and how often it moves anything
#' set.seed(22)
#' trials <- replicate(30, {
#'   r <- apply_random_operator(sol, W, adj, row_W, col_W)
#'   c(op = r$op_idx, moved = !identical(r$sol, sol))
#' })
#' table(trials["op", ])
#' mean(trials["moved", ])
#'
#' @export
apply_random_operator <- function(sol, W, adj, row_W, col_W, best_sol = NULL, ...) {
  op_fns <- list(
    operator_1,
    operator_2,
    operator_3,
    operator_4,
    operator_5,
    operator_6,
    operator_7,
    operator_8,
    operator_9,
    operator_10,
    operator_11
  )

  weights <- c(1,1,1,1,1,2,1,1,1,0.5,1.5)

  op_idx <- sample(11L, 1L, prob = weights / sum(weights))

  if(op_idx == 10L){
    new_sol <- operator_10(sol, W, adj, row_W, col_W, best_sol = best_sol, ...)
  }else if(op_idx == 11L){
    new_sol <- operator_11(sol, W, adj, row_W, col_W, Pop_max = params$Pop_max, ...)
  }else {
    new_sol <- op_fns[[op_idx]](sol, W, adj, row_W, col_W, ...)
  }
  list(sol = new_sol, op_idx = op_idx)
}
