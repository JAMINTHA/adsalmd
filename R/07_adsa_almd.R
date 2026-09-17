# ============================================================
# FILE: R/07_adsa_almd.R
# AdSA-ALMD main algorithm — Algorithm 1 from the paper.
#
# Key fixes vs original code:
#   - Cooling schedule is ADAPTIVE (Eq. 15), not geometric (T * 0.99)
#   - Temperature is updated INSIDE the inner loop (Algorithm 1 line 22)
#   - MSDA / MSDR track consecutive-streak fitness deviations (Eq. 16–17)
#   - rK and mK reset at the START of every outer trial K
#   - All parameters passed explicitly — no globals
# ============================================================

# ---- T0 estimation ------------------------------------------

#' Estimate the initial annealing temperature (Equation 3)
#'
#' Calibrates the starting temperature to the actual scale of the objective on
#' your data, so that a typical uphill move is accepted with probability
#' \code{P0}. Hard-coding a temperature instead would make the annealing
#' schedule meaningless, since the objective's magnitude depends on \code{N} and
#' on the commuting volumes in \code{W}.
#'
#' @details
#' \deqn{T_0 = \frac{-E[|\Delta F|]}{\ln P_0}}
#'
#' \eqn{E[|\Delta F|]} is estimated by perturbing \code{init_sol}
#' \code{n_samples} times with \code{\link{apply_random_operator}()} and
#' averaging the absolute fitness change. Because \eqn{\ln P_0 < 0} for
#' \eqn{P_0 < 1}, the result is positive. Substituting it into the Metropolis
#' rule \eqn{\exp(\Delta F / T_0)} gives an acceptance probability of about
#' \code{P0} for a move of average badness — a higher \code{P0} therefore means
#' a hotter start and more early exploration.
#'
#' Every perturbation starts from \code{init_sol}, not from a random walk, so
#' the estimate describes the neighbourhood of the starting solution rather
#' than of the whole search space.
#'
#' @section Robustness:
#' Two guards keep the result usable:
#' \itemize{
#'   \item only finite, strictly positive deltas enter the mean, so samples
#'     where the operator declined to move (delta exactly 0) do not drag the
#'     estimate towards zero;
#'   \item if that leaves nothing — every sampled move was a no-op — the mean
#'     falls back to \code{1}. On a small or highly constrained map this
#'     fallback can fire, giving a temperature unrelated to your objective's
#'     scale. If a run behaves as though it is accepting everything or nothing,
#'     check this function's output first.
#' }
#'
#' @section Contiguity is deliberately not checked:
#' Both the baseline and the sampled fitnesses are computed with
#' \code{check_contiguity = FALSE}. Most random perturbations break contiguity
#' somewhere, and each of those would otherwise score exactly 0, so the sampled
#' deltas would degenerate into "distance from \code{F0} to zero" repeated —
#' measuring the objective's level rather than its local variability. Scoring
#' infeasible candidates gives a meaningful spread. The search itself always
#' checks contiguity.
#'
#' @section Global variables:
#' Calls \code{\link{apply_random_operator}()}, which reads a \code{params} list
#' from the global environment. See \code{\link{operator_1}()}.
#'
#' @param init_sol Integer vector, the starting partition (normally SARA's
#'   output).
#' @param W Numeric \code{N x N} origin-destination matrix.
#' @param adj Integer/logical \code{N x N} adjacency matrix.
#' @param row_W Numeric vector, \code{rowSums(W)}.
#' @param col_W Numeric vector, \code{colSums(W)}.
#' @param SC_min Numeric in (0, 1], the minimum self-containment threshold.
#' @param Pop_min Numeric, the minimum market population.
#' @param Pop_max Numeric, the maximum market population.
#' @param Pop_tar Numeric or \code{NULL}, target population. Forwarded to
#'   \code{\link{compute_objective}()}; not used.
#' @param r Integer, penalty steepness exponent.
#' @param P0 Numeric in (0, 1), the target acceptance probability at \code{T0}
#'   (default 0.8). Higher means a hotter start.
#' @param n_samples Integer, number of random perturbations used for the
#'   estimate (default 100). More samples give a steadier temperature at linear
#'   cost.
#' @return Numeric scalar \code{T0 > 0}.
#'
#' @family AdSA-ALMD core
#' @seealso \code{\link{adsa_almd}()}, which calls this once per outer trial;
#'   \code{\link{apply_random_operator}()} for the sampling;
#'   \code{\link{compute_objective}()} for the fitness being differenced.
#'
#' @examples
#' adj <- matrix(0L, 6, 6); adj[cbind(1:5, 2:6)] <- 1L; adj <- adj + t(adj)
#' W <- matrix(1, 6, 6); W[1:3, 1:3] <- 30; W[4:6, 4:6] <- 30; diag(W) <- 60
#' row_W <- rowSums(W); col_W <- colSums(W)
#' params <- list(SC_min = 0.55, Pop_min = 100, Pop_max = 1000)
#'
#' set.seed(31)
#' T0 <- estimate_t0(c(1, 1, 1, 2, 2, 2), W, adj, row_W, col_W,
#'                   SC_min = 0.55, Pop_min = 100, Pop_max = 1000,
#'                   Pop_tar = NULL, r = 3, P0 = 0.8, n_samples = 30)
#' T0
#'
#' ## A hotter start accepts more: raising P0 raises T0
#' set.seed(31)
#' estimate_t0(c(1, 1, 1, 2, 2, 2), W, adj, row_W, col_W,
#'             0.55, 100, 1000, NULL, r = 3, P0 = 0.95, n_samples = 30)
#'
#' @export
estimate_t0 <- function(init_sol,
                        W, 
                        adj, 
                        row_W, 
                        col_W,
                        SC_min, 
                        Pop_min,
                        Pop_max,
                        Pop_tar, 
                        r,
                        P0 = 0.8, 
                        n_samples = 100L) {
# Compute F for the starting point (skip contiguity check for speed)
  F0     <- compute_objective(init_sol, 
                              W, 
                              row_W, 
                              col_W,
                              SC_min, 
                              Pop_min,
                              Pop_max,
                              Pop_tar, 
                              r, 
                              adj,
                              check_contiguity = FALSE)
  
  deltas <- numeric(n_samples)

  for (i in seq_len(n_samples)) {
    cand     <- apply_random_operator(init_sol, W, adj, row_W, col_W,
                                       best_sol = NULL)$sol
    F_cand   <- compute_objective(cand, 
                                  W, 
                                  row_W, 
                                  col_W,
                                  SC_min, 
                                  Pop_min,
                                  Pop_max,
                                  Pop_tar, 
                                  r, 
                                  adj,
                                  check_contiguity = FALSE)
    deltas[i] <- abs(F_cand - F0)
  }

  avg_delta <- mean(deltas[is.finite(deltas) & deltas > 0])
  if (!is.finite(avg_delta) || avg_delta == 0) 
    avg_delta <- 1  # fallback

  -avg_delta / log(P0)
}


# ---- Adaptive Cooling Schedule helpers ----------------------

#' Root mean squared deviation of a streak from the current best (Eq. 16-17)
#'
#' MSD = sqrt( mean((F_j - mu)^2) ) over the fitnesses in one consecutive
#' acceptance or rejection streak.
#'
#' Deviation is measured against mu, the incumbent best, not against the
#' streak's own mean: what the cooling schedule needs to know is how far the
#' search has wandered from the best solution found, not how internally
#' consistent the streak was. An empty streak scores 0, which is what makes
#' the very first iteration of a trial well defined.
#'
#' @param fitness_streak Numeric vector, the consecutive accepted (MSDA) or
#'   rejected (MSDR) fitness values.
#' @param mu Numeric, the current best fitness.
#' @return Numeric scalar >= 0.
#' @noRd
.compute_msd <- function(fitness_streak, mu) {
  if (length(fitness_streak) == 0L) 
    return(0)
  sqrt(mean((fitness_streak - mu)^2))
}

#' Adaptive Cooling Schedule temperature update (Equation 15)
#'
#' Called once per inner iteration (Algorithm 1 line 22) to produce the next
#' temperature. Three effects compose, in order:
#'
#'   1. Baseline geometric cooling: T <- 0.99 * T_k, every step.
#'   2. Streak response. On a long rejection streak (r_k > 100) the temperature
#'      is raised by a factor 1 + min(alpha * MSDR_norm, 0.1) to escape a
#'      basin; on a long acceptance streak (m_k > 100) it is lowered by
#'      max(1 - beta * MSDA_norm, 0.85) to consolidate. MSDA and MSDR are
#'      normalised by max(|mu|, 1) first, so the response is scale-free with
#'      respect to the objective's magnitude.
#'   3. Acceptance-rate correction. Once n_total >= 350, the running acceptance
#'      rate is compared with target_rate: cool 2% extra if it is more than
#'      twice the target, warm 2% if it is less than half.
#'
#' Finally the result is floored at .Machine$double.eps so it stays strictly
#' positive and exp(delta / T) never divides by zero.
#'
#' DEVIATION FROM THE PAPER. Equation 15 as written is additive,
#' T +/- exp(coef * MSD). That form overflows once MSD is large, producing
#' temperature spikes and the high run-to-run variance seen in early testing.
#' The implementation above is multiplicative and clamped: the per-step rise is
#' capped at 10% and the per-step fall at 15%, so no single step can move the
#' temperature by more than a bounded factor. The qualitative behaviour
#' (warm on rejection streaks, cool on acceptance streaks) is preserved; the
#' magnitudes are not those of the paper's formula.
#'
#' Note the rate-correction guard reads `n_total %% 350L`, which is truthy for
#' every n_total that is NOT a multiple of 350 -- so the correction applies on
#' essentially every step past 350 rather than once per 350 steps.
#'
#' @param T_k    Numeric, current temperature.
#' @param r_k    Integer, consecutive rejection count.
#' @param m_k    Integer, consecutive acceptance count.
#' @param MSDA   Numeric, MSD of the accepted-fitness streak.
#' @param MSDR   Numeric, MSD of the rejected-fitness streak.
#' @param alpha  Numeric in (0, 1), increase coefficient.
#' @param beta   Numeric in (0, 1), decrease coefficient.
#' @param mu     Numeric, current best fitness; used only to normalise MSDA/MSDR.
#' @param n_accepted,n_total Integers, running acceptance counts across the
#'   whole run.
#' @param target_rate Numeric, desired long-run acceptance rate (default 0.3).
#' @return Numeric, the updated temperature, strictly positive.
#' @noRd
.update_temperature <- function(T_k,
                                r_k,
                                m_k,
                                MSDA,
                                MSDR,
                                alpha,
                                beta,
                                mu,
                                n_accepted,
                                n_total,
                                target_rate = 0.3) {
  
  T_new <- T_k * 0.99
  
  scale <- max(abs(mu), 1)
  MSDA_norm <- MSDA / scale
  MSDR_norm <- MSDR / scale
  if (r_k > 100L) {
    T_new <- T_new * (1 + min(alpha * MSDR_norm, 0.1))
  } else if (m_k > 100L) {
    T_new <- T_new * max(1 - beta * MSDA_norm, 0.85)
  }
  if (n_total >= 350L && n_total %% 350L) {
    actual_rate <- n_accepted / n_total
    if (actual_rate > target_rate * 2) {
      T_new <- T_new * 0.98
    } else if (actual_rate < target_rate * 0.5) {
      T_new <- T_new * 1.02
    }
  }
  
  T_new <- max(T_new, .Machine$double.eps)
  
  T_new
}

#' Create an empty elite pool
#'
#' Allocates the fixed-size pool of best-scoring solutions that
#' \code{\link{adsa_almd}()} maintains alongside its single incumbent best.
#'
#' @details
#' A plain list of three components: \code{sols} (a list of partitions),
#' \code{fits} (their objective values) and \code{size}. Every slot starts
#' empty, with fitness \code{-Inf} so that the first \code{size} genuine
#' solutions offered to \code{\link{update_elite}()} all get admitted.
#'
#' The pool exists because simulated annealing throws away almost everything it
#' visits: the incumbent best is a single point, and a late uphill excursion can
#' leave the search far from it with no record of the good solutions passed on
#' the way. Keeping the top few lets a run report near-optimal alternatives —
#' often several partitions within a fraction of a percent of each other but
#' with materially different market boundaries, which is exactly what you want
#' to look at when the objective alone cannot separate them.
#'
#' A pool of 5 is what \code{\link{adsa_almd}()} uses. There is no diversity
#' criterion, so a pool can legitimately fill with near-duplicates of the same
#' partition.
#'
#' @param size Integer, number of solutions to retain (default 5).
#' @return Named list:
#'   \describe{
#'     \item{\code{sols}}{list of length \code{size}, all \code{NULL}.}
#'     \item{\code{fits}}{numeric vector of length \code{size}, all
#'       \code{-Inf}.}
#'     \item{\code{size}}{integer, the pool size.}
#'   }
#'
#' @family elite pool
#' @seealso \code{\link{update_elite}()}, \code{\link{get_best_elite}()},
#'   \code{\link{adsa_almd}()}, which returns a pool as \code{$elite}.
#'
#' @examples
#' elite <- init_elite(3)
#' str(elite)
#' elite$fits
#'
#' @export
init_elite <- function(size = 5L) {
  list(
    sols = vector("list",size),
    fits = rep(-Inf,size),
    size = size
  )
}

#' Offer a solution to the elite pool
#'
#' Admits \code{sol} if its fitness beats the pool's current worst entry,
#' displacing that entry. Otherwise the pool is returned untouched.
#'
#' @details
#' The pool is unsorted; the worst slot is located with \code{which.min()} on
#' each call. At \code{size = 5} that is cheaper than maintaining order, and it
#' means the pool's slot positions carry no meaning — use
#' \code{\link{get_best_elite}()} rather than indexing.
#'
#' The comparison is strict (\code{>}), so a solution that exactly ties the
#' worst entry is not admitted. Freshly initialised slots hold \code{-Inf} and
#' are therefore always displaced first.
#'
#' Like the operator tracker, the pool is immutable in the R sense: the updated
#' list is returned and must be assigned back.
#'
#' @param elite List from \code{\link{init_elite}()}.
#' @param sol Integer vector, the candidate partition.
#' @param fitness Numeric, its objective value from
#'   \code{\link{compute_objective}()}.
#' @return The updated pool.
#'
#' @family elite pool
#' @seealso \code{\link{init_elite}()}, \code{\link{get_best_elite}()}
#'
#' @examples
#' elite <- init_elite(3)
#' elite <- update_elite(elite, c(1, 1, 2, 2), 0.8)
#' elite <- update_elite(elite, c(1, 2, 2, 2), 1.2)
#' elite <- update_elite(elite, c(1, 1, 1, 2), 0.5)
#' elite$fits
#'
#' ## A fourth solution displaces the worst entry only if it beats it
#' elite <- update_elite(elite, c(2, 2, 1, 1), 0.6)
#' elite$fits
#' elite <- update_elite(elite, c(2, 1, 1, 1), 0.1)
#' elite$fits
#'
#' @export
update_elite <- function(elite,sol,fitness) {
  worst_idx <- which.min(elite$fits)
  if(fitness > elite$fits[worst_idx]){
    elite$sols [[worst_idx]] <- sol
    elite$fits [worst_idx] <- fitness
  }
  elite
}

#' Best solution in an elite pool
#'
#' Returns the partition with the highest recorded fitness.
#'
#' @details
#' A pool that has never been offered a solution holds \code{NULL} in every
#' slot and \code{-Inf} in every fitness, so this returns \code{NULL} rather
#' than erroring — check for that before using the result. Ties are broken by
#' \code{which.max()}, i.e. the lowest slot index wins.
#'
#' In a completed run this normally equals \code{result$best_sol}, but the two
#' can differ: \code{\link{adsa_almd}()} applies
#' \code{\link{fix_misallocations}()} to its incumbent best after every outer
#' trial, and those repaired partitions are not fed back into the pool.
#'
#' @param elite List from \code{\link{init_elite}()}.
#' @return Integer vector, the best partition in the pool, or \code{NULL} for an
#'   untouched pool.
#'
#' @family elite pool
#' @seealso \code{\link{init_elite}()}, \code{\link{update_elite}()}
#'
#' @examples
#' elite <- init_elite(3)
#' elite <- update_elite(elite, c(1, 1, 2, 2), 0.8)
#' elite <- update_elite(elite, c(1, 2, 2, 2), 1.2)
#' get_best_elite(elite)
#'
#' ## An untouched pool yields NULL
#' is.null(get_best_elite(init_elite(3)))
#'
#' @export
get_best_elite <- function(elite){
  elite$sols[[which.max(elite$fits)]]
}

#' Repair cohesion-misallocated BGUs, respecting the population bounds
#'
#' Sweeps the partition repeatedly, moving any BGU that would be more cohesive
#' in an adjacent market — but only when the move leaves both the donor and the
#' recipient inside \code{[Pop_min, Pop_max]} and keeps every market
#' contiguous. Run by \code{\link{adsa_almd}()} on its incumbent best after
#' every outer trial.
#'
#' @details
#' Each pass walks BGU 1 to N. For a BGU flagged by
#' \code{\link{is_bgu_misallocated}()} it:
#' \enumerate{
#'   \item computes its cohesion in its current market, excluding itself;
#'   \item scores it against every adjacent market with
#'     \code{\link{compute_ci}()} and keeps those that beat the current score;
#'   \item picks the best of those, then refuses the move if the donor would
#'     drop below \code{Pop_min} or the recipient would rise above
#'     \code{Pop_max};
#'   \item applies it only if every market stays contiguous.
#' }
#' Passes repeat until one completes without a single move, or 100 passes have
#' run. Moves take effect immediately, so a later BGU in the same pass is
#' evaluated against the already-updated partition.
#'
#' @section Relationship to refine_misallocated_bgus():
#' \code{\link{refine_misallocated_bgus}()} performs the same sweep \emph{but
#' does not check populations at all}. The difference is deliberate and the two
#' are used at different points:
#' \itemize{
#'   \item this function runs \emph{inside} the annealing loop, where a move
#'     that violates a population bound would make the partition score 0 and be
#'     thrown away — so it declines such moves up front;
#'   \item \code{\link{refine_misallocated_bgus}()} runs once \emph{after} the
#'     search as a final cosmetic pass, prioritising cohesion.
#' }
#' Using the unbounded version inside the loop would repeatedly hand the
#' annealer infeasible partitions; using this one as the final pass would leave
#' misallocations in place.
#'
#' Cost is the reason to care: \code{\link{is_bgu_misallocated}()} is called for
#' every BGU on every pass, and each call scores several candidate markets. On a
#' large map this dominates the per-trial cost of the search.
#'
#' @param best_sol Integer vector, the partition to repair (length N).
#' @param W Numeric \code{N x N} origin-destination matrix.
#' @param adj Integer/logical \code{N x N} adjacency matrix.
#' @param row_W Numeric vector, \code{rowSums(W)}.
#' @param col_W Numeric vector, \code{colSums(W)}.
#' @param Pop_min Numeric, population floor the donor must still satisfy
#'   (default 500).
#' @param Pop_max Numeric, population ceiling the recipient must still satisfy
#'   (default 3500).
#' @return Integer vector of length N. Relabelled by
#'   \code{\link{get_ordered_vec}()} if any move was made; otherwise
#'   \code{best_sol} exactly as supplied.
#'
#' @family AdSA-ALMD core
#' @seealso \code{\link{refine_misallocated_bgus}()} for the unbounded
#'   post-search variant; \code{\link{is_bgu_misallocated}()} for the test;
#'   \code{\link{count_misallocated}()} to measure progress.
#'
#' @examples
#' adj <- matrix(0L, 6, 6); adj[cbind(1:5, 2:6)] <- 1L; adj <- adj + t(adj)
#' W <- matrix(1, 6, 6); W[1:3, 1:3] <- 30; W[4:6, 4:6] <- 30; diag(W) <- 60
#' row_W <- rowSums(W); col_W <- colSums(W)
#'
#' ## BGU 4 sits in the wrong market
#' sol <- c(1, 1, 1, 1, 2, 2)
#' count_misallocated(sol, W, adj, row_W, col_W)
#'
#' ## Loose bounds: the repair goes ahead
#' fixed <- fix_misallocations(sol, W, adj, row_W, col_W,
#'                             Pop_min = 50, Pop_max = 1000)
#' fixed
#' count_misallocated(fixed, W, adj, row_W, col_W)
#'
#' ## Tight bounds block the same move, leaving the partition alone
#' fix_misallocations(sol, W, adj, row_W, col_W,
#'                    Pop_min = 450, Pop_max = 1000)
#'
#' @export
fix_misallocations <- function(best_sol, W, adj, row_W, col_W, Pop_min = 500, Pop_max = 3500) {
  improved <- TRUE
  iter <- 0L
  while (improved && iter < 100L) {
    iter <- iter + 1L
    improved <- FALSE
    for (g in seq_along(best_sol)) {
      if (is_bgu_misallocated(g, best_sol, W, adj, row_W, col_W)) {
        cur_cl <- best_sol[g]
        cur_mbrs <- cluster_members(cur_cl, best_sol)
        adj_cls <- setdiff(
          unique(best_sol[which(adj[g, ] > 0L)]), cur_cl
        )
        if (length(adj_cls) == 0L) next
        ci_cur <- compute_ci(g, setdiff(cur_mbrs, g), W, row_W, col_W)
        scores <- vapply(adj_cls, function(cl)
          compute_ci(g, cluster_members(cl, best_sol), W, row_W, col_W), numeric(1L)
        )
        better <- adj_cls[scores > ci_cur]
        if (length(better) == 0L) next
        best_cl <- adj_cls[which.max(scores[scores > ci_cur])]
        pop_cur_after <- compute_population(setdiff(cur_mbrs, g), row_W)
        pop_tgt_after <- compute_population(cluster_members(best_cl, best_sol), row_W) + compute_population(c(g), row_W)
        if (pop_cur_after < Pop_min) next
        if (pop_tgt_after > Pop_max) next
        cand <- best_sol
        cand[g] <- best_cl
        cand <- getOrderedVec(cand)
        if (all_clusters_contiguous(cand, adj)) {
          best_sol <- cand
          improved <- TRUE
        }
      }
    }
  }
  best_sol
}

# ---- Main algorithm ---------------------------------------------------

#' AdSA-ALMD: adaptive simulated annealing for a single starting solution
#' (Algorithm 1)
#'
#' Stage 2 of the method. Takes a feasible partition — normally
#' \code{\link{sara_iterated}()}'s output — and improves it by simulated
#' annealing under an adaptive cooling schedule, perturbing with the eleven
#' group-based operators and repairing misallocations between trials.
#'
#' @details
#' The search is a nested loop. The \strong{outer loop} runs \code{L} trials;
#' each one resets to the incumbent best, re-estimates the temperature, and
#' ends with a repair pass. The \strong{inner loop} runs \code{l} Metropolis
#' iterations at the current temperature:
#' \enumerate{
#'   \item propose a candidate with \code{\link{apply_random_operator}()},
#'     which receives the incumbent best so \code{\link{operator_10}()} can
#'     cross over with it;
#'   \item score it with \code{\link{compute_objective}()};
#'   \item accept with probability \code{min(1, exp(delta_F / T))} — so any
#'     improvement is taken, and a worsening move is taken with a probability
#'     that falls as the temperature does;
#'   \item extend the acceptance or rejection streak, recompute MSDA or MSDR,
#'     and update the temperature (Algorithm 1 line 22 — inside the inner loop,
#'     not between trials);
#'   \item record the attempt in the operator tracker, and offer any accepted
#'     solution to the elite pool.
#' }
#' After each trial, \code{\link{fix_misallocations}()} is applied to the
#' incumbent best, its fitness is recomputed, one row of diagnostics is
#' appended to the history, and — if \code{checkpoint_file} is set — the current
#' best is written to disk.
#'
#' @section Temperature schedule:
#' \code{\link{estimate_t0}()} is re-run at the start of \emph{every} outer
#' trial, then scaled by \code{(T0 / 100) * 0.5^(K - 1)}. Each trial therefore
#' starts half as hot as the last: trial 1 explores widely, and by trial 20 the
#' schedule is effectively a local search. Combined with the reset to the
#' incumbent best at the top of each trial, this makes the algorithm a sequence
#' of progressively cooler restarts rather than one long anneal. Within a trial
#' the temperature is then adapted by \code{.update_temperature()}.
#'
#' A practical consequence: with a large \code{L} the later trials contribute
#' almost nothing, because \code{0.5^(K-1)} underflows. Most of the work happens
#' in the first few dozen trials, and the convergence check below is what
#' should end the run.
#'
#' @section Convergence and early stopping:
#' From trial 101 onwards (\code{patience = 100}, hard-coded), a trial whose
#' best fitness moved by no more than \code{eps} increments a counter; any real
#' improvement resets it. The run stops when the counter reaches 100. With the
#' default \code{eps = 0} this requires exact equality, so it fires only on a
#' genuinely frozen search.
#'
#' @section Output and progress:
#' Per-trial progress is written with \code{message()} on \emph{every} trial,
#' and the temperature with \code{print()}, regardless of \code{verbose} —
#' which currently gates only the "Converged at trial K" notice. Wrap the call
#' in \code{\link[base]{suppressMessages}()} and
#' \code{\link[utils]{capture.output}()} if you need it quiet;
#' \code{\link{run_adsa_pipeline_parallel}()} runs on workers where this output
#' is discarded anyway.
#'
#' @section Global variables:
#' The operator dispatcher reads a \code{params} list from the global
#' environment. Calling this function directly requires one in scope;
#' \code{\link{run_adsa_pipeline}()} sets it up. See \code{\link{operator_1}()}.
#'
#' @param init_sol Integer vector, the starting partition (normally from
#'   \code{\link{sara_iterated}()}). Relabelled on entry.
#' @param W Numeric \code{N x N} origin-destination matrix.
#' @param adj Integer/logical \code{N x N} adjacency matrix.
#' @param row_W Numeric vector, \code{rowSums(W)}, precomputed.
#' @param col_W Numeric vector, \code{colSums(W)}, precomputed.
#' @param SC_min Numeric in (0, 1], the minimum self-containment threshold.
#' @param Pop_min Numeric, the minimum market population. Enforced as a hard
#'   rejection by the objective.
#' @param Pop_max Numeric, the maximum market population (default \code{Inf}).
#' @param Pop_tar Numeric or \code{NULL}, target population (default
#'   \code{NULL}). Forwarded to the objective; not used.
#' @param r Integer, penalty steepness exponent (3 or 4).
#' @param L Integer, number of outer trials (default 1000). See the temperature
#'   section on why very large values add little.
#' @param l Integer, inner Metropolis iterations per trial (default 20).
#' @param alpha Numeric, cooling-schedule increase coefficient (default 0.2).
#' @param beta Numeric, cooling-schedule decrease coefficient (default 0.2).
#' @param eps Numeric, convergence threshold on the trial-to-trial change in
#'   best fitness (default 0).
#' @param P0 Numeric in (0, 1), target acceptance probability used by
#'   \code{\link{estimate_t0}()} (default 0.8).
#' @param T0_samples Integer, perturbations per temperature estimate (default
#'   100). Paid once per outer trial.
#' @param checkpoint_file Character path or \code{NULL} (default). If supplied,
#'   the current best solution is written to this CSV (columns \code{bgu_id},
#'   \code{cluster}) after every outer trial. Leave \code{NULL} for parallel
#'   runs — concurrent workers writing to one path will corrupt each other's
#'   output.
#' @param verbose Logical (default \code{TRUE}). Currently gates only the
#'   convergence notice; see the output section.
#' @return Named list:
#'   \describe{
#'     \item{\code{elite}}{elite pool (\code{\link{init_elite}()} structure) of
#'       the best solutions seen during the search.}
#'     \item{\code{best_sol}}{integer vector, the best partition found.}
#'     \item{\code{best_fitness}}{numeric, \code{\link{compute_objective}()} of
#'       \code{best_sol}.}
#'     \item{\code{gci}}{numeric, its global cohesion index.}
#'     \item{\code{base_fitness}}{numeric, \code{m * GCI}, for cross-method
#'       comparison.}
#'     \item{\code{n_clusters}}{integer, the number of markets.}
#'     \item{\code{lma_metrics}}{\code{data.frame} from
#'       \code{\link{build_lma_df}()}.}
#'     \item{\code{history}}{\code{data.frame} with one row per completed trial:
#'       \code{trial}, \code{best_fitness}, \code{temp}, \code{n_clusters}.}
#'     \item{\code{op_tracker}}{operator counters; print with
#'       \code{\link{print_op_tracker}()}.}
#'   }
#'
#' @family AdSA-ALMD core
#' @seealso \code{\link{run_adsa_pipeline}()}, which calls this as stage 2 and
#'   handles the globals for you; \code{\link{sara_iterated}()} for stage 1;
#'   \code{\link{refine_misallocated_bgus}()} for stage 3;
#'   \code{\link{compute_objective}()} for what is maximised;
#'   \code{\link{print_op_tracker}()} and \code{\link{summarise_runs}()} for
#'   diagnostics.
#'
#' @examples
#' ## Toy system: six BGUs in a line, forming two 3-BGU markets
#' adj <- matrix(0L, 6, 6); adj[cbind(1:5, 2:6)] <- 1L; adj <- adj + t(adj)
#' W <- matrix(1, 6, 6); W[1:3, 1:3] <- 30; W[4:6, 4:6] <- 30; diag(W) <- 60
#' row_W <- rowSums(W); col_W <- colSums(W)
#'
#' ## The operators read a global `params` (see the section above)
#' params <- list(SC_min = 0.55, Pop_min = 100, Pop_max = 1000)
#'
#' \donttest{
#' set.seed(41)
#' res <- suppressMessages(adsa_almd(
#'   init_sol = c(1, 1, 1, 1, 2, 2),
#'   W = W, adj = adj, row_W = row_W, col_W = col_W,
#'   SC_min = 0.55, Pop_min = 100, Pop_max = 1000,
#'   r = 3, L = 5L, l = 10L, T0_samples = 20L, verbose = FALSE
#' ))
#'
#' res$best_sol
#' res$n_clusters
#' res$lma_metrics
#' res$history
#'
#' ## Which operators actually did the work
#' print_op_tracker(res$op_tracker)
#'
#' ## Near-optimal alternatives the search passed through
#' res$elite$fits
#' }
#'
#' @export
adsa_almd <- function(init_sol,
                      W,
                      adj,
                      row_W,
                      col_W,
                      SC_min,
                      Pop_min,
                      Pop_max = Inf,
                      Pop_tar = NULL,
                      r,
                      L             = 1000L,
                      l             = 20L,
                      alpha         = 0.2,
                      beta          = 0.2,
                      eps           = 0,
                      P0            = 0.8,
                      T0_samples    = 100L,
                      checkpoint_file = NULL,
                      verbose       = TRUE) {
  # ---- Initialise (Algorithm 1 lines 2-6) ----------
  sol <- getOrderedVec(init_sol)
  
  F_cur <- compute_objective(sol, W, row_W, col_W, SC_min, Pop_min, Pop_max, Pop_tar, r, adj)
  
  mu       <- F_cur   # best fitness so far (mu in paper)
  best_sol <- sol      # best solution so far
  
  MSDA <- 0   # initialised once, persist across trials
  MSDR <- 0
  
  history      <- vector("list", L)
  prev_best    <- -Inf
  patience     <- 100L
  max_no_improve <- 100L
  no_improve   <- 0L
  
  op_tracker <- init_op_tracker()
  
  n_accepted_total <- 0L
  n_total_total    <- 0L
  elite <- init_elite(size = 5L)
  elite <- update_elite(elite, sol, F_cur)
  mu    <- F_cur
  
  
  
  # ---- outer loop: K = 1 ...L (Algorithm 1 line 7) -----
  for (K in seq_len(L)) {
    # Reset consecutive-streak counters each trial (Algorithm 1 line 8)
    r_k <- 0L
    m_k <- 0L
    accepted_streak <- numeric(0)
    rejected_streak <- numeric(0)
    
    sol   <- best_sol
    F_cur <- mu
    T_k <- estimate_t0(sol, W, adj, row_W, col_W, SC_min, Pop_min, Pop_max, Pop_tar, r, P0, T0_samples)
    T_k <- (T_k / 100) * (0.5^(K - 1L))
    print(T_k)
    
    # ---- Inner loop: k = 1 ...l (Algorithm 1 line 9) ----
    for (k in seq_len(l)) {
      # Generate candidate (Algorithm 1 line 10)
      op_result <- apply_random_operator(sol, W, adj, row_W, col_W, SC_min, Pop_min, best_sol)
      x_c <- op_result$sol
      op_idx <- op_result$op_idx
      
      F_c <- compute_objective(x_c,
                               W,
                               row_W,
                               col_W,
                               SC_min,
                               Pop_min,
                               Pop_max,
                               Pop_tar,
                               r,
                               adj)
      
      # Acceptance probability (Algorithm 1 lines 11-13)
      delta_F <- F_c - F_cur
      prob    <- min(1, exp(delta_F / T_k))
      
      sol_before <- sol
      accepted   <- FALSE
      improved   <- FALSE
      
      if (runif(1L) < prob) {
        # ---- Accept (Algorithm 1 lines 14-18) ----------
        sol             <- x_c
        F_cur           <- F_c
        m_k             <- m_k + 1L
        r_k             <- 0L
        accepted_streak <- c(accepted_streak, F_c)
        rejected_streak <- numeric(0)
        accepted        <- TRUE
        
        if (F_c > mu) {
          mu       <- F_c
          best_sol <- x_c
          improved <- TRUE
        }
        MSDA <- .compute_msd(accepted_streak, mu)
        
        elite <- update_elite(elite, x_c, F_c)
        
      } else {
        # ---- Reject (Algorithm 1 lines 19 - 21) -----
        r_k             <- r_k + 1L
        m_k             <- 0L
        rejected_streak <- c(rejected_streak, F_c)
        accepted_streak <- numeric(0)
        MSDR            <- .compute_msd(rejected_streak, mu)
      }
      
      # Update tracker
      op_tracker <- update_op_tracker(
        op_tracker,
        op_idx = op_idx,
        sol_before = sol_before,
        sol_after = x_c,
        F_before = F_cur,
        F_after = F_c,
        accepted = accepted,
        improved_best = improved
      )
      
      n_total_total <- n_total_total + 1L
      if (accepted) {
        n_accepted_total <- n_accepted_total + 1L
      }
      
      # Update temperature INSIDE inner loop (Algorithm 1 line 22)
      T_k <- .update_temperature(
        T_k,
        r_k,
        m_k,
        MSDA,
        MSDR,
        alpha,
        beta,
        mu = mu,
        n_accepted = n_accepted_total,
        n_total = n_total_total,
        target_rate = 0.3
      )
      
    } # end inner loop
    
    best_sol <- fix_misallocations(best_sol, W, adj, row_W, col_W, Pop_min, Pop_max)

    if (!is.null(checkpoint_file)) {
      write.csv(
        data.frame(
          bgu_id = seq_along(best_sol),
          cluster = best_sol
        ),
        file = checkpoint_file, row.names = FALSE
      )
    }

    mu <- compute_objective(best_sol,
                            W,
                            row_W,
                            col_W,
                            SC_min,
                            Pop_min,
                            Pop_max,
                            Pop_tar,
                            r,
                            adj)
    
    # Store trial-level diagnostics
    history[[K]] <- data.frame(
      trial         = K,
      best_fitness  = mu,
      temp          = T_k,
      n_clusters    = n_clusters(best_sol),
      stringsAsFactors = FALSE
    )
    
    message(sprintf(
      "Trial %4d | Best F = %8.4f | T = %.6f | m = %d",
      K,
      mu,
      T_k,
      n_clusters(best_sol)
    ))
    
    # Convergence check
    if (K > patience) {
      if (abs(mu - prev_best) <= eps) {
        no_improve <- no_improve + 1L
      } else {
        no_improve <- 0L
      }
      
      if (no_improve >= max_no_improve) {
        if (verbose)
          message(sprintf("Converged at trial %d", K))
        break
      }
    }
    
    prev_best <- mu
  } # end outer loop
  
  # ---- Final metrics ------
  
  history_df <- do.call(rbind, history[!vapply(history, is.null, logical(1L))])
  gci <- compute_gci(best_sol, W, row_W, col_W)
  base_fitness <- compute_base_fitness(best_sol, W, row_W, col_W)
  lma_df <- build_lma_df(best_sol, W, row_W, col_W)
  
  list(
    elite         = elite,
    best_sol      = best_sol,
    best_fitness  = mu,
    gci           = gci,
    base_fitness  = base_fitness,
    n_clusters    = n_clusters(best_sol),
    lma_metrics   = lma_df,
    history       = history_df,
    op_tracker    = op_tracker
  )
}
