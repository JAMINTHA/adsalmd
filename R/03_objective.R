# ============================================================
# FILE: R/03_objective.R
# Penalty functions and the AdSA-ALMD objective function.
# All formulas match the paper exactly.
# ============================================================

# ---- Individual penalty terms -------------------------------

#' Self-containment penalty for one market (Equation 13)
#'
#' Charges a market for falling short of the self-containment floor. Returns 0
#' for any market that already clears \code{SC_min}, so the penalty is
#' one-sided: there is no reward for exceeding the threshold.
#'
#' @details
#' \deqn{P_{SC}(M_s) = \begin{cases}
#'   1 - (SC(M_s) / SC_{min})^r & SC(M_s) < SC_{min} \\
#'   0 & \text{otherwise}
#' \end{cases}}
#'
#' The penalty runs from 0 (at \code{sc == SC_min}) to 1 (at \code{sc == 0}).
#' The exponent \code{r} sets its shape, not its range: with \code{r = 3} or
#' \code{4} the curve is flat near the threshold and steep near zero, so a
#' market that just misses the floor is barely charged while one that misses it
#' badly is charged almost the full unit. Raising \code{r} makes the search
#' more tolerant of near misses.
#'
#' Scale matters when reading the result: the penalty is on a 0–1 scale per
#' market and is subtracted from a global cohesion index that grows with
#' \code{N}, so one failing market out of hundreds moves the objective very
#' little. Feasibility at the hard limits is enforced separately by
#' \code{\link{compute_objective}()}, which returns 0 outright.
#'
#' @param sc Numeric, the market's self-containment
#'   (\code{\link{compute_sc}()}\code{$sc}).
#' @param SC_min Numeric in (0, 1], the minimum self-containment threshold.
#' @param r Integer, steepness exponent (the paper uses 3 or 4).
#' @return Numeric scalar in \code{[0, 1]}.
#'
#' @family penalty functions
#' @seealso \code{\link{penalty_sc_total}()} to sum over all markets,
#'   \code{\link{penalty_total}()} for the combined penalty.
#'
#' @examples
#' penalty_sc_single(0.60, SC_min = 0.55, r = 3)   # above the floor: 0
#' penalty_sc_single(0.50, SC_min = 0.55, r = 3)   # just below: small
#' penalty_sc_single(0.10, SC_min = 0.55, r = 3)   # far below: large
#'
#' ## The exponent controls tolerance of near misses
#' sc <- seq(0, 0.55, by = 0.11)
#' data.frame(
#'   sc,
#'   r3 = vapply(sc, penalty_sc_single, numeric(1), SC_min = 0.55, r = 3),
#'   r4 = vapply(sc, penalty_sc_single, numeric(1), SC_min = 0.55, r = 4)
#' )
#'
#' @export
penalty_sc_single <- function(sc, SC_min, r) {
  if (sc >= SC_min) return(0)
  1 - (sc / SC_min)^r
}

#' Population penalty for one undersized market (Equation 14)
#'
#' The population analogue of \code{\link{penalty_sc_single}()}: charges a
#' market for having fewer residents than \code{Pop_min}, and nothing at all
#' once it clears the floor.
#'
#' @details
#' \deqn{P_{Pop}(M_s) = \begin{cases}
#'   1 - (Pop(M_s) / Pop_{min})^r & Pop(M_s) < Pop_{min} \\
#'   0 & \text{otherwise}
#' \end{cases}}
#'
#' Same shape and same 0–1 range as the self-containment penalty, so the two
#' are directly comparable when summed.
#'
#' Note that \code{\link{compute_objective}()} additionally rejects any
#' partition containing an undersized market outright (fitness 0), so in the
#' main search this penalty shapes the landscape only among solutions that are
#' already feasible. It remains useful on its own for diagnosing how close a
#' rejected partition came.
#'
#' @param pop Numeric, the market's population
#'   (\code{\link{compute_population}()}).
#' @param Pop_min Numeric > 0, the minimum population threshold.
#' @param r Integer, steepness exponent.
#' @return Numeric scalar in \code{[0, 1]}.
#'
#' @family penalty functions
#' @seealso \code{\link{penalty_pop_min_total}()},
#'   \code{\link{penalty_pop_max_single}()} for the opposite constraint.
#'
#' @examples
#' penalty_pop_min_single(600, Pop_min = 500, r = 3)   # above the floor: 0
#' penalty_pop_min_single(400, Pop_min = 500, r = 3)
#' penalty_pop_min_single(50,  Pop_min = 500, r = 3)
#'
#' @export
penalty_pop_min_single <- function(pop, Pop_min, r){
  if(pop >= Pop_min) return(0)
  1-(pop/Pop_min) ^ r
}

#' Population penalty for one oversized market
#'
#' Charges a market for exceeding the population ceiling. Unlike the
#' floor penalties this one is unbounded above, and it is deliberately
#' slack: nothing is charged until the market is 50 percent past
#' \code{Pop_max}.
#'
#' @details
#' \deqn{P^{max}_{Pop}(M_s) = \begin{cases}
#'   (Pop(M_s) / Pop_{max} - 1)^r & Pop(M_s) > 1.5\, Pop_{max} \\
#'   0 & \text{otherwise}
#' \end{cases}}
#'
#' Two asymmetries are worth knowing about:
#' \itemize{
#'   \item \strong{The dead zone.} The trigger is \code{1.5 * Pop_max} but the
#'     formula measures excess against \code{Pop_max}, so the penalty jumps
#'     discontinuously from 0 to \code{0.5^r} the moment a market crosses the
#'     trigger. It is a barrier against runaway markets, not a smooth gradient.
#'   \item \strong{No upper bound.} A market at three times \code{Pop_max}
#'     scores \code{2^r}, which with \code{r = 3} already outweighs eight
#'     fully-failing self-containment penalties.
#' }
#'
#' In the main search this term rarely fires, because
#' \code{\link{compute_objective}()} rejects any partition with a market over
#' \code{Pop_max} outright. It is most useful when ranking infeasible
#' intermediates — for instance inside \code{\link{sara_iterated}()}, which has
#' to choose between solutions that all still violate the ceiling.
#'
#' @param pop Numeric, the market's population.
#' @param Pop_max Numeric > 0, the maximum population threshold.
#' @param r Integer, steepness exponent.
#' @return Numeric scalar \eqn{\ge 0}, unbounded above.
#'
#' @family penalty functions
#' @seealso \code{\link{penalty_pop_max_total}()},
#'   \code{\link{operator_11}()}, which exists specifically to break up markets
#'   over \code{Pop_max}.
#'
#' @examples
#' penalty_pop_max_single(3000, Pop_max = 4000, r = 3)   # under: 0
#' penalty_pop_max_single(5000, Pop_max = 4000, r = 3)   # over, inside the
#'                                                       # 1.5x dead zone: 0
#' penalty_pop_max_single(6001, Pop_max = 4000, r = 3)   # past the trigger
#' penalty_pop_max_single(12000, Pop_max = 4000, r = 3)  # runaway market
#'
#' @export
penalty_pop_max_single <- function(pop, Pop_max, r){
  if(pop <= Pop_max*1.5) return(0)
  base_penalty <- (pop/Pop_max - 1)^r

  return(base_penalty)
}

#' Total oversized-population penalty across all markets
#'
#' Sums \code{\link{penalty_pop_max_single}()} over a vector of market
#' populations.
#'
#' @details
#' Because the per-market term is unbounded, this total is dominated by the
#' single largest offender rather than by the count of offenders: one market at
#' four times the ceiling outweighs ten markets just past the trigger. When you
#' care about how \emph{many} markets are over the limit, use
#' \code{\link{sara_max_pop_score}()}, which reports count and excess
#' separately.
#'
#' @param pop_vec Numeric vector of market populations, typically
#'   \code{\link{build_lma_df}(sol, W, row_W, col_W)$pop}.
#' @param Pop_max Numeric > 0, the maximum population threshold.
#' @param r Integer, steepness exponent.
#' @return Numeric scalar \eqn{\ge 0}.
#'
#' @family penalty functions
#' @seealso \code{\link{penalty_pop_max_single}()},
#'   \code{\link{sara_max_pop_score}()}
#'
#' @examples
#' pops <- c(1200, 3800, 6500, 12000)
#' penalty_pop_max_total(pops, Pop_max = 4000, r = 3)
#'
#' ## Dominated by the worst market, not the number of offenders
#' penalty_pop_max_total(rep(6100, 10), Pop_max = 4000, r = 3)
#' penalty_pop_max_total(12000,         Pop_max = 4000, r = 3)
#'
#' @export
penalty_pop_max_total <- function(pop_vec, Pop_max, r){
  sum(vapply(pop_vec, penalty_pop_max_single, numeric(1L), Pop_max = Pop_max, r=r))
}

#' Concentration penalty: discourage one market swallowing the region
#'
#' Charges a partition for unevenness in market sizes, using the
#' Herfindahl-Hirschman index of the population shares. This is what stops the
#' search from drifting towards a single dominant market surrounded by
#' fragments — a failure mode that neither the self-containment nor the
#' per-market population penalties can see, because each individual market can
#' look perfectly healthy.
#'
#' @details
#' With \eqn{p_s = Pop(M_s) / \sum_t Pop(M_t)} the population shares and
#' \eqn{m} the number of markets:
#' \deqn{HHI = \sum_s p_s^2, \qquad
#'       P_{conc} = w \cdot \max\!\left(HHI - \tfrac{1}{m},\ 0\right)}
#'
#' \eqn{1/m} is the HHI of a perfectly even split, so the penalty is the
#' \emph{excess} concentration over what is achievable at that number of
#' markets. It is therefore 0 for an even partition at any \code{m}, and rises
#' towards \eqn{w(1 - 1/m)} as one market takes everything. The
#' \code{max(..., 0)} guard only matters against floating-point noise, since
#' HHI \eqn{\ge 1/m} always.
#'
#' An empty or all-zero population vector returns 0 rather than \code{NaN}.
#'
#' The weight is the term's whole calibration. \code{\link{penalty_total}()}
#' applies it with \code{weight = 2}, which on the Queensland data puts it on
#' roughly the same scale as the summed self-containment penalty. Raising it
#' pushes the search towards evenly sized markets at the cost of cohesion.
#'
#' @param pop_vec Numeric vector of market populations.
#' @param weight Numeric multiplier on the excess concentration (default 1).
#' @return Numeric scalar \eqn{\ge 0}.
#'
#' @family penalty functions
#' @seealso \code{\link{penalty_total}()}, which applies this with weight 2.
#'
#' @examples
#' ## An even split is not penalised at all
#' penalty_concentration(rep(1000, 5))
#'
#' ## One dominant market plus fragments is
#' penalty_concentration(c(9000, 250, 250, 250, 250))
#'
#' ## Scale-free: only shares matter
#' penalty_concentration(c(900, 100)) == penalty_concentration(c(9000, 1000))
#'
#' ## Degenerate inputs are safe
#' penalty_concentration(numeric(0))
#' penalty_concentration(c(0, 0, 0))
#'
#' @export
penalty_concentration <- function(pop_vec, weight = 1.0){
  total <- sum(pop_vec)
  if(total == 0) return(0)

  shares <- pop_vec / total

  hhi <- sum(shares^2)

  m <- length(pop_vec)

  perfect <- 1/m

  excess <- max(hhi - perfect, 0)

  weight * excess
}

# NOTE: penalty_pop_single() was the original Equation 14 implementation.
# It is superseded by penalty_pop_min_single() above and kept here only as a
# record of the original formulation:
#
#   penalty_pop_single <- function(pop, Pop_min, r) {
#     if (pop >= Pop_min) return(0)
#     1 - (pop / Pop_min)^r
#   }

# ---- Aggregate penalties ------------------------------------

#' Total self-containment penalty across all markets (Equation 13, summed)
#'
#' Sums \code{\link{penalty_sc_single}()} over every market in a partition.
#'
#' @details
#' \deqn{P_{SC}(x) = \sum_{M_s \in M} P_{SC}(M_s)}
#'
#' Since each term lies in \code{[0, 1]}, the total is bounded by the number of
#' markets and is readable as "how many markets' worth of self-containment the
#' partition is missing". A total of 2.5 over 40 markets means the shortfall is
#' equivalent to two and a half markets having no self-containment at all.
#'
#' @param sc_vec Numeric vector of market self-containment values, typically
#'   \code{\link{build_lma_df}(sol, W, row_W, col_W)$sc}.
#' @param SC_min Numeric in (0, 1], the minimum self-containment threshold.
#' @param r Integer, steepness exponent.
#' @return Numeric scalar in \code{[0, length(sc_vec)]}.
#'
#' @family penalty functions
#' @seealso \code{\link{penalty_sc_single}()}, \code{\link{penalty_total}()}
#'
#' @examples
#' W <- matrix(1, 6, 6); W[1:3, 1:3] <- 30; W[4:6, 4:6] <- 30; diag(W) <- 60
#' row_W <- rowSums(W); col_W <- colSums(W)
#'
#' lma <- build_lma_df(c(1, 1, 1, 2, 2, 2), W, row_W, col_W)
#' penalty_sc_total(lma$sc, SC_min = 0.55, r = 3)      # a good partition: 0
#'
#' lma_bad <- build_lma_df(c(1, 1, 2, 2, 3, 3), W, row_W, col_W)
#' penalty_sc_total(lma_bad$sc, SC_min = 0.95, r = 3)
#'
#' @export
penalty_sc_total <- function(sc_vec, SC_min, r) {
  sum(vapply(sc_vec, penalty_sc_single, numeric(1L),
             SC_min = SC_min, r = r))
}

#' Total undersized-population penalty across all markets (Equation 14, summed)
#'
#' Sums \code{\link{penalty_pop_min_single}()} over every market in a
#' partition.
#'
#' @details
#' \deqn{P_{Pop}(x) = \sum_{M_s \in M} P_{Pop}(M_s)}
#'
#' Bounded by the number of markets, like \code{\link{penalty_sc_total}()}.
#'
#' Note that \code{\link{penalty_total}()} does \emph{not} include this term —
#' the current objective enforces the population floor as a hard rejection in
#' \code{\link{compute_objective}()} instead. This function remains exported
#' for diagnostics and for building alternative objectives.
#'
#' @param pop_vec Numeric vector of market populations, typically
#'   \code{\link{build_lma_df}(sol, W, row_W, col_W)$pop}.
#' @param Pop_min Numeric > 0, the minimum population threshold.
#' @param r Integer, steepness exponent (default 3).
#' @return Numeric scalar in \code{[0, length(pop_vec)]}.
#'
#' @family penalty functions
#' @seealso \code{\link{penalty_pop_min_single}()},
#'   \code{\link{compute_objective}()} for where the floor is actually enforced.
#'
#' @examples
#' penalty_pop_min_total(c(1200, 800, 300, 90), Pop_min = 500, r = 3)
#' penalty_pop_min_total(c(1200, 800), Pop_min = 500)   # all above the floor
#'
#' @export
penalty_pop_min_total <- function(pop_vec, Pop_min, r = 3L) {
  sum(vapply(pop_vec, penalty_pop_min_single, numeric(1L),
             Pop_min = Pop_min, r = r))
}

#' Combined penalty term of the objective function (Equation 9)
#'
#' The penalty \eqn{P(x)} that \code{\link{compute_objective}()} subtracts from
#' the global cohesion index. Takes the per-market summary table and returns
#' one number.
#'
#' @details
#' As currently implemented the total is
#' \deqn{P(x) = P_{SC}(x) + 2 \cdot P_{conc}(x)}
#' that is, \code{\link{penalty_sc_total}()} on the \code{sc} column plus
#' \code{\link{penalty_concentration}()} on the \code{pop} column with a fixed
#' weight of 2.
#'
#' @section Which arguments are used:
#' \code{lma_df}, \code{SC_min} and \code{r} affect the result.
#' \code{Pop_min}, \code{Pop_max} and \code{Pop_tar} are accepted but
#' \strong{not used}: the population floor and ceiling are enforced as hard
#' rejections inside \code{\link{compute_objective}()} (a violating partition
#' scores 0 outright), so charging a smooth penalty for them as well would be
#' redundant. They stay in the signature because
#' \code{\link{compute_objective}()} forwards its full threshold set here, and
#' because reinstating \code{\link{penalty_pop_min_total}()} or
#' \code{\link{penalty_pop_max_total}()} is then a one-line change. Do not read
#' a non-zero result as evidence that a population constraint fired.
#'
#' The concentration weight of 2 is likewise fixed in the function body rather
#' than exposed; change it there if you are recalibrating the objective.
#'
#' @param lma_df A per-market \code{data.frame} from
#'   \code{\link{build_lma_df}()}; only the \code{sc} and \code{pop} columns are
#'   read.
#' @param SC_min Numeric in (0, 1], the minimum self-containment threshold.
#' @param Pop_min Numeric, the minimum population threshold. Accepted for
#'   signature compatibility; not used (see Details).
#' @param Pop_max Numeric, the maximum population threshold (default
#'   \code{Inf}). Accepted for signature compatibility; not used.
#' @param Pop_tar Numeric or \code{NULL}, target population (default
#'   \code{NULL}). Accepted for signature compatibility; not used.
#' @param r Integer, steepness exponent passed to
#'   \code{\link{penalty_sc_total}()}.
#' @return Numeric scalar \eqn{\ge 0}.
#'
#' @family penalty functions
#' @seealso \code{\link{compute_objective}()}, the only internal caller;
#'   \code{\link{penalty_sc_total}()} and
#'   \code{\link{penalty_concentration}()} for the two live components.
#'
#' @examples
#' W <- matrix(1, 6, 6); W[1:3, 1:3] <- 30; W[4:6, 4:6] <- 30; diag(W) <- 60
#' row_W <- rowSums(W); col_W <- colSums(W)
#' lma <- build_lma_df(c(1, 1, 1, 2, 2, 2), W, row_W, col_W)
#'
#' penalty_total(lma, SC_min = 0.55, Pop_min = 100, r = 3)
#'
#' ## Decompose it: self-containment plus weighted concentration
#' penalty_sc_total(lma$sc, 0.55, 3) + penalty_concentration(lma$pop, 2)
#'
#' ## Pop_min does not change the answer (see the Details section)
#' penalty_total(lma, SC_min = 0.55, Pop_min = 1e6, r = 3)
#'
#' @export
penalty_total <- function(lma_df, SC_min, Pop_min, Pop_max = Inf, Pop_tar = NULL, r) {
  conc_weight <- 2.0
  p_sc <- penalty_sc_total(lma_df$sc, SC_min, r)
  p_conc <- penalty_concentration(lma_df$pop, weight = conc_weight)

  p_sc + p_conc
}

# ---- Objective function -------------------------------------

#' AdSA-ALMD objective function (Equation 5)
#'
#' The single number the whole search maximises: cohesion minus penalties, with
#' every hard constraint applied as an outright rejection to 0. This is what
#' \code{\link{adsa_almd}()} evaluates for every candidate solution, so it is
#' the hottest function in the package.
#'
#' @details
#' \deqn{F(x) = \bigl(GCI(x) - P(x)\bigr) \prod_{i=1}^{m} C(M_i)}
#'
#' where \eqn{C(M_i)} is 1 if market \eqn{i} is spatially contiguous and 0
#' otherwise. The product form makes contiguity a hard constraint: a single
#' broken market zeroes the entire objective.
#'
#' @section Rejection rules:
#' The function returns exactly \code{0} — without computing cohesion at all —
#' when any of these holds. They are checked in this order, cheapest first:
#' \enumerate{
#'   \item any market is non-contiguous
#'     (\code{\link{all_clusters_contiguous}()}), unless
#'     \code{check_contiguity = FALSE};
#'   \item \code{Pop_max} is finite and any market exceeds it;
#'   \item any market falls below \code{Pop_min}.
#' }
#' Otherwise the value is \code{\link{compute_gci}()} minus
#' \code{\link{penalty_total}()}, which can legitimately be negative when the
#' penalties outweigh a weak partition's cohesion.
#'
#' Because rejection and "a genuinely worthless partition" both produce 0, a
#' returned 0 is not diagnostic on its own. To find out \emph{why} a solution
#' scored 0, inspect \code{\link{build_lma_df}()} and
#' \code{\link{all_clusters_contiguous}()} directly.
#'
#' @section The check_contiguity escape hatch:
#' Contiguity checking is a BFS per market and dominates the cost of this
#' function. \code{\link{estimate_t0}()} sets \code{check_contiguity = FALSE}
#' because it only needs the \emph{spread} of fitness differences across random
#' perturbations, not their feasibility — and zeroing out most samples would
#' make that spread meaningless. Leave it \code{TRUE} everywhere else:
#' with it off, non-contiguous partitions receive real, comparable scores and
#' the search will happily accept them.
#'
#' \code{Pop_tar} is accepted and forwarded to \code{\link{penalty_total}()} but
#' does not currently affect the result.
#'
#' @param sol Integer vector of market assignments (length N).
#' @param W Numeric \code{N x N} origin-destination matrix.
#' @param row_W Numeric vector, \code{rowSums(W)}.
#' @param col_W Numeric vector, \code{colSums(W)}.
#' @param SC_min Numeric in (0, 1], the minimum self-containment threshold.
#' @param Pop_min Numeric, the minimum market population. Enforced as a hard
#'   rejection.
#' @param Pop_max Numeric, the maximum market population (default \code{Inf},
#'   i.e. no ceiling). Enforced as a hard rejection when finite.
#' @param Pop_tar Numeric or \code{NULL}, target population (default
#'   \code{NULL}). Forwarded to \code{\link{penalty_total}()}; not used.
#' @param r Integer, penalty steepness exponent (the paper uses 3 or 4).
#' @param adj Integer/logical \code{N x N} adjacency matrix.
#' @param check_contiguity Logical (default \code{TRUE}). Set \code{FALSE} only
#'   for initial-temperature estimation; see the section above.
#' @return Numeric scalar. \code{0} for any rejected partition, otherwise
#'   \code{GCI - P}, which may be negative.
#'
#' @family objective function
#' @seealso \code{\link{compute_gci}()} and \code{\link{penalty_total}()} for
#'   the two components; \code{\link{compute_base_fitness}()} for the
#'   unpenalised cross-method measure; \code{\link{adsa_almd}()} for the search
#'   that maximises this.
#'
#' @examples
#' ## Toy system: six BGUs in a line, forming two 3-BGU markets
#' adj <- matrix(0L, 6, 6); adj[cbind(1:5, 2:6)] <- 1L; adj <- adj + t(adj)
#' W <- matrix(1, 6, 6); W[1:3, 1:3] <- 30; W[4:6, 4:6] <- 30; diag(W) <- 60
#' row_W <- rowSums(W); col_W <- colSums(W)
#'
#' ## The natural partition scores well
#' compute_objective(c(1, 1, 1, 2, 2, 2), W, row_W, col_W,
#'                   SC_min = 0.55, Pop_min = 100, r = 3, adj = adj)
#'
#' ## A non-contiguous market is rejected outright
#' compute_objective(c(1, 2, 1, 2, 2, 2), W, row_W, col_W,
#'                   SC_min = 0.55, Pop_min = 100, r = 3, adj = adj)
#'
#' ## So is a market under the population floor
#' compute_objective(c(1, 1, 1, 2, 2, 2), W, row_W, col_W,
#'                   SC_min = 0.55, Pop_min = 1e6, r = 3, adj = adj)
#'
#' ## Turning off the contiguity check scores the broken partition anyway
#' compute_objective(c(1, 2, 1, 2, 2, 2), W, row_W, col_W,
#'                   SC_min = 0.55, Pop_min = 100, r = 3, adj = adj,
#'                   check_contiguity = FALSE)
#'
#' @export
compute_objective <- function(sol, W, row_W, col_W,
                               SC_min, Pop_min, Pop_max = Inf, Pop_tar = NULL, r, adj,
                               check_contiguity = TRUE) {
  if (check_contiguity && !all_clusters_contiguous(sol, adj)) {
    return(0)
  }
  lma_df <- build_lma_df(sol, W, row_W, col_W)
  if(is.finite(Pop_max) && any(lma_df$pop > Pop_max)){
    return(0)
  }
  if(any(lma_df$pop < Pop_min)){
    return(0)
  }

  gci    <- compute_gci(sol, W, row_W, col_W)

  p      <- penalty_total(lma_df, SC_min, Pop_min, Pop_max, Pop_tar, r)

  gci - p
}
