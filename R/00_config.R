# ============================================================
# FILE: R/00_config.R
# AdSA-ALMD Configuration — All tunable parameters in one place.
# Modify ONLY this file to change algorithm behaviour.
# ============================================================

AdSA_params <- list(

  # ---- Data ----
  N            = 520L,    # BGUs after merging uninhabited SA2s

  # ---- Self-containment thresholds (Table 2 grid search optimum) ----
  SC_min       = 0.55,   # Minimum self-containment threshold
  SC_tar       = 0.75,   # Target  self-containment threshold

  # ---- Population thresholds — 4 % sample of QLD (Table 2) ----
  Pop_min      = 500L,   # Minimum population
  Pop_tar      = 2000L,  # Target  population
  Pop_max      = 4000L,
  
  # ---- Penalty function (paper: r = 3 or 4) ----
  penalty_exp  = 3L,     # Exponent r controlling steepness

  # ---- Initial temperature ----
  P0           = 0.9,    # Target acceptance probability at T0
  T0_samples   = 500L,   # Perturbations used to estimate E[|DeltaF|]

  # ---- AdSA-ALMD inner / outer loops (Algorithm 1) ----
  L            = 1000L,  # Number of outer trials
  l            = 25L,    # Inner iterations per trial
  eps          = 0,      # Convergence threshold epsilon

  # ---- Adaptive Cooling Schedule (Equation 15) ----
  alpha        = 0.2,    # Temperature-increase coefficient (consecutive rejects)
  beta         = 0.2,    # Temperature-decrease coefficient (consecutive accepts)

  # ---- Parallel execution ----
  n_runs       = 10L,    # Independent runs (different random seeds)
  n_cores      = 8L      # CPU cores available in ABS DataLab
)
