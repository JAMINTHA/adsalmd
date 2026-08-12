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

#' Estimate initial temperature T0   [Equation 3]
#'
#' T0 = -E[|Delta F|] / ln(P0)
#'
#' E[|Delta F|] is the mean absolute fitness change across n_samples
#' random perturbations of the initial solution.
#' @param init_sol   Integer vector, starting solution.
#' @param W          Numeric matrix.
#' @param adj        Integer/logical matrix.
#' @param row_W      Numeric vector.
#' @param col_W      Numeric vector.
#' @param SC_min     Numeric.
#' @param Pop_min    Numeric.
#' @param r          Integer, penalty exponent.
#' @param P0         Numeric, target acceptance probability (default 0.8).
#' @param n_samples  Integer, perturbations used (default 100).
#' @return Numeric scalar T0 > 0.
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

#' Mean squared deviation from current best   [Equations 16–17]
#'
#' MSD = sqrt( (1/n) * sum (F_j - mu)^2 )
#'
#' @param fitness_streak Numeric vector, consecutive accepted or rejected fitnesses.
#' @param mu             Numeric, current best fitness value.
#' @return Numeric scalar >= 0.
.compute_msd <- function(fitness_streak, mu) {
  if (length(fitness_streak) == 0L) 
    return(0)
  sqrt(mean((fitness_streak - mu)^2))
}

#' Adaptive Cooling Schedule temperature update   [Equation 15]
#'
#' T_{K+1} = T_K + exp(alpha * MSDR)   if r_K > 0  (consecutive rejects)
#' T_{K+1} = T_K - exp(beta  * MSDA)   if m_K > 0  (consecutive accepts)
#'
#' Temperature is floored at .Machine$double.eps to stay positive.
#'
#' Stability clamps (not in paper formula but required for numerical stability):
#'   - Temperature INCREASE is capped at T_k (cannot more than double per step)
#'   - Temperature DECREASE is capped at 50% of T_k per step
#'   These prevent exp() overflow when MSDA/MSDR are large, which would
#'   cause wild temperature spikes and the high SD seen without clamping.
#'
#' @param T_k   Numeric, current temperature.
#' @param r_k   Integer, consecutive rejection count.
#' @param m_k   Integer, consecutive acceptance count.
#' @param MSDA  Numeric, MSD of accepted-solution fitness values.
#' @param MSDR  Numeric, MSD of rejected-solution fitness values.
#' @param alpha Numeric (0, 1), increase coefficient.
#' @param beta  Numeric (0, 1), decrease coefficient.
#' @return Numeric, updated temperature.
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

init_elite <- function(size = 5L) {
  list(
    sols = vector("list",size),
    fits = rep(-Inf,size),
    size = size
  )
}

update_elite <- function(elite,sol,fitness) {
  worst_idx <- which.min(elite$fits)
  if(fitness > elite$fits[worst_idx]){
    elite$sols [[worst_idx]] <- sol
    elite$fits [worst_idx] <- fitness
  }
  elite
}

get_best_elite <- function(elite){
  elite$sols[[which.max(elite$fits)]]
}

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

#' Run AdSA-ALMD for a single starting solution   [Algorithm 1]
#'
#' @param init_sol   Integer vector, initial solution (from SARA).
#' @param W          Numeric matrix, OD matrix (N x N).
#' @param adj        Integer/logical matrix, adjacency (N x N).
#' @param row_W      Numeric vector, rowSums(W)  [precomputed].
#' @param col_W      Numeric vector, colSums(W)  [precomputed].
#' @param SC_min     Numeric, minimum SC threshold.
#' @param Pop_min    Numeric, minimum population threshold.
#' @param r          Integer, penalty exponent (3 or 4).
#' @param L          Integer, outer trials (default 1000).
#' @param l          Integer, inner iterations per trial (default 20).
#' @param alpha      Numeric, ACS increase coefficient.
#' @param beta       Numeric, ACS decrease coefficient.
#' @param eps        Numeric, convergence threshold (default 0).
#' @param P0         Numeric, target acceptance prob for T0 (default 0.8).
#' @param T0_samples Integer, samples to estimate T0 (default 100).
#' @param checkpoint_file Character path or NULL (default). If supplied,
#'   the current best_sol is written to this CSV after every outer trial.
#'   Leave NULL (the default) for parallel runs, since concurrent workers
#'   writing to the same path will corrupt each other's output.
#' @param verbose    Logical, print progress every 100 outer trials.
#' @return Named list:
#'   best_sol     — integer vector, best solution found;
#'   best_fitness — numeric, F(best_sol);
#'   gci          — numeric, GCI of best_sol;
#'   base_fitness — numeric, m * GCI (for cross-method comparison);
#'   n_clusters   — integer, number of LMAs;
#'   lma_metrics  — data.frame from build_lma_df();
#'   history      — data.frame, trial-level diagnostics.
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
