# =============================================================================
# FROM ABS CENSUS MICRODATA TO OD AND ADJACENCY MATRICES
# Using exact ABS variable names as they appear in PLIDA / DataLab
#
# Jamintha Samarakoon — QUT / ABS DataLab
#
# INPUTS (as they appear in DataLab):
#   ABSPID   — person identifier
#   PURP     — Place of Usual Residence (coded to SA2 directly)
#   POWP     — Place of Work (coded to Destination Zone / DZN)
#   IFPOWP   — Imputation flag for Place of Work (0 = observed, 1 = imputed)
#   LFSP     — Labour Force Status (filter to employed persons only)
#   AGEP     — Age (filter to working age 15-64)
#
# OUTPUTS:
#   W        — Origin-Destination flow matrix [n_SA2 x n_SA2]
#   A        — Geographic adjacency matrix    [n_SA2 x n_SA2]
# =============================================================================

# ── 0. Packages ───────────────────────────────────────────────────────────────
# install.packages(c("tidyverse", "sf", "spdep"))
library(tidyverse)
library(sf)
library(spdep)


# =============================================================================
# STEP 1: LOAD CENSUS MICRODATA
# Inside the DataLab this would be read from the PLIDA SAS/CSV files.
# Column names match the ABS Census Dictionary variable names exactly.
# =============================================================================

# --- 1a. Simulate data with exact ABS column names ---
# (Replace this block with: census_raw <- read_csv("plida_census_2016.csv")
#  inside the DataLab)

set.seed(42)
n <- 1000

# Simulate 10 Queensland SA2 codes (9-digit format)
sa2_codes <- c(
  "301011001", "301011002", "301011003", "301011004", "301011005",
  "301011006", "301011007", "301011008", "301011009", "301011010"
)

sa2_names <- c(
  "Fortitude Valley", "New Farm",       "Spring Hill",   "Paddington",
  "Milton",           "South Brisbane", "West End",      "Woolloongabba",
  "Kangaroo Point",   "East Brisbane"
)

# Simulate DZN codes (DZNs nest within SA2s — typically 3-5 DZNs per SA2)
# Real DZN codes are 9 digits: first 9 match SA2, then 2 more digits
dzn_sa2_lookup <- tibble(
  DZN_CODE_2016    = c(
    "30101100101", "30101100102",          # DZNs for SA2 301011001
    "30101100201", "30101100202",          # DZNs for SA2 301011002
    "30101100301", "30101100302",          # DZNs for SA2 301011003
    "30101100401",                         # DZNs for SA2 301011004
    "30101100501", "30101100502",          # DZNs for SA2 301011005
    "30101100601",                         # DZNs for SA2 301011006
    "30101100701", "30101100702",          # DZNs for SA2 301011007
    "30101100801",                         # DZNs for SA2 301011008
    "30101100901", "30101100902",          # DZNs for SA2 301011009
    "30101101001"                          # DZNs for SA2 301011010
  ),
  SA2_MAINCODE_2016 = c(
    "301011001", "301011001",
    "301011002", "301011002",
    "301011003", "301011003",
    "301011004",
    "301011005", "301011005",
    "301011006",
    "301011007", "301011007",
    "301011008",
    "301011009", "301011009",
    "301011010"
  )
)

# LFSP categories (ABS Census Dictionary):
#   1 = Employed, worked full-time
#   2 = Employed, worked part-time
#   3 = Employed, away from work
#   4 = Unemployed, looking for full-time work
#   5 = Unemployed, looking for part-time work
#   6 = Not in the labour force
#   @ = Not applicable (under 15)

census_raw <- tibble(
  ABSPID  = sprintf("P%07d", 1:n),
  AGEP    = sample(15:85, n, replace = TRUE),
  LFSP    = sample(c("1","2","3","4","5","6","@"), n, replace = TRUE,
                   prob = c(0.35, 0.20, 0.05, 0.05, 0.03, 0.25, 0.07)),
  PURP    = sample(sa2_codes, n, replace = TRUE,
                   prob = c(0.15,0.12,0.10,0.08,0.10,
                            0.12,0.11,0.08,0.07,0.07)),
  POWP    = sample(dzn_sa2_lookup$DZN_CODE_2016, n, replace = TRUE),
  IFPOWP  = sample(c(0L, 1L), n, replace = TRUE, prob = c(0.889, 0.111))
)

cat("=== RAW CENSUS MICRODATA (exact ABS column names) ===\n")
cat("Columns: ABSPID | AGEP | LFSP | PURP | POWP | IFPOWP\n\n")
print(head(census_raw, 12))
cat(sprintf("\nTotal person records: %d\n\n", nrow(census_raw)))


# =============================================================================
# STEP 2: FILTER TO EMPLOYED WORKING-AGE PERSONS
# LFSP 1, 2, 3 = employed (full-time, part-time, away from work)
# AGEP 15-64   = working age
# IFPOWP == 0  = non-imputed place of work only (optional — depends on project)
# =============================================================================

census_employed <- census_raw |>
  filter(
    AGEP   >= 15 & AGEP <= 64,        # working age
    LFSP   %in% c("1", "2", "3"),     # employed persons only
    IFPOWP == 0L                       # non-imputed POWP records
  )

cat("=== AFTER FILTERING TO EMPLOYED WORKING-AGE PERSONS ===\n")
cat(sprintf(
  "Remaining: %d of %d records (%.1f%%)\n\n",
  nrow(census_employed), nrow(census_raw),
  100 * nrow(census_employed) / nrow(census_raw)
))


# =============================================================================
# STEP 3: CONVERT POWP FROM DZN TO SA2
# POWP is coded to Destination Zone (DZN) in the raw Census data.
# DZNs aggregate to SA2s — the ABS provides DZN_SA2_2016_AUST.csv
# for this conversion (downloadable from ABS website).
#
# Inside DataLab:
#   dzn_sa2_lookup <- read_csv("DZN_SA2_2016_AUST.csv")
# =============================================================================

census_sa2 <- census_employed |>
  left_join(
    dzn_sa2_lookup |> select(DZN_CODE_2016, SA2_MAINCODE_2016),
    by = c("POWP" = "DZN_CODE_2016")
  ) |>
  rename(
    origin_SA2      = PURP,              # SA2 of usual residence
    destination_SA2 = SA2_MAINCODE_2016  # SA2 of place of work (converted)
  ) |>
  filter(!is.na(destination_SA2))        # remove unmatched DZN codes

cat("=== AFTER DZN -> SA2 CONVERSION ===\n")
cat("New columns: origin_SA2 (from PURP) | destination_SA2 (from POWP via DZN)\n\n")
print(
  census_sa2 |>
    select(ABSPID, AGEP, LFSP, PURP = origin_SA2,
           POWP_DZN = POWP, POWP_SA2 = destination_SA2) |>
    head(10)
)


# =============================================================================
# STEP 4: BUILD THE OD FLOW MATRIX
# W[i,j] = number of employed persons residing in SA2 i, working in SA2 j
# =============================================================================

# --- 4a. Long-format flow table ---
od_long <- census_sa2 |>
  group_by(origin_SA2, destination_SA2) |>
  summarise(W_ij = n(), .groups = "drop") |>
  arrange(origin_SA2, destination_SA2)

cat("\n=== OD FLOW TABLE — long format (origin_SA2 | destination_SA2 | W_ij) ===\n")
print(od_long, n = 20)
cat(sprintf(
  "\nNon-zero OD pairs: %d of %d possible (%s sparsity)\n\n",
  nrow(od_long),
  length(sa2_codes)^2,
  scales::percent(1 - nrow(od_long) / length(sa2_codes)^2,
                  accuracy = 0.1)
))

# --- 4b. Wide-format n x n OD matrix ---
W_wide <- od_long |>
  pivot_wider(
    names_from  = destination_SA2,
    values_from = W_ij,
    values_fill = 0L
  ) |>
  arrange(origin_SA2)

# Ensure all SA2s appear as both rows and columns
all_sa2 <- sort(sa2_codes)

W_wide <- W_wide |>
  right_join(tibble(origin_SA2 = all_sa2), by = "origin_SA2") |>
  arrange(origin_SA2)

missing_cols <- setdiff(all_sa2, colnames(W_wide))
for (col in missing_cols) W_wide[[col]] <- 0L
W_wide[is.na(W_wide)] <- 0L

W_matrix <- W_wide |>
  column_to_rownames("origin_SA2") |>
  select(all_of(all_sa2)) |>
  as.matrix()

cat("=== OD MATRIX W [n x n] ===\n")
cat("Rows = origin SA2 (PURP)   |   Cols = destination SA2 (POWP converted)\n")
cat("W[i,j] = commuters residing in SA2 i, working in SA2 j\n\n")
print(W_matrix)

# --- 4c. Self-containment (SC) for each SA2 ---
W_ii  <- diag(W_matrix)                        # internal flows
W_iG  <- rowSums(W_matrix)                     # total outflows (supply)
W_Gi  <- colSums(W_matrix)                     # total inflows  (demand)
SC_SS <- ifelse(W_iG > 0, W_ii / W_iG, NA)    # supply-side SC
SC_DS <- ifelse(W_Gi > 0, W_ii / W_Gi, NA)    # demand-side SC
SC    <- pmin(SC_SS, SC_DS, na.rm = TRUE)      # overall SC = min(SS, DS)

sc_table <- tibble(
  SA2_code  = all_sa2,
  SA2_name  = sa2_names[match(all_sa2, sa2_codes)],
  W_ii      = W_ii,
  W_iG      = W_iG,
  SC_SS     = round(SC_SS, 3),
  SC_DS     = round(SC_DS, 3),
  SC        = round(SC, 3)
)

cat("\n=== SELF-CONTAINMENT TABLE ===\n")
cat("SC_SS = W[i,i] / W[i,G]  (supply-side)\n")
cat("SC_DS = W[i,i] / W[G,i]  (demand-side)\n")
cat("SC    = min(SC_SS, SC_DS)\n\n")
print(sc_table)


# =============================================================================
# STEP 5: BUILD THE GEOGRAPHIC ADJACENCY MATRIX
# A[i,j] = 1 if SA2 i and SA2 j share a boundary, 0 otherwise
#
# In real analysis use ABS SA2 shapefiles:
#   SA2_2016_AUST_GDA2020.shp  (ABS digital boundary files)
# =============================================================================

cat("\n=== STEP 5: GEOGRAPHIC ADJACENCY MATRIX ===\n")

# --- Real shapefile code (run outside DataLab with ABS shapefiles) ---
cat("Real code using ABS SA2 shapefiles:\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat("  sa2_shp <- st_read('SA2_2016_AUST_GDA2020.shp') |>\n")
cat("    filter(STE_NAME16 == 'Queensland') |>\n")
cat("    filter(SA2_MAIN16 %in% your_sa2_codes) |>\n")
cat("    st_make_valid()\n\n")
cat("  # Queen contiguity (shared boundary or point)\n")
cat("  nb <- poly2nb(sa2_shp, queen = TRUE)\n\n")
cat("  # Binary adjacency matrix\n")
cat("  A <- nb2mat(nb, style = 'B', zero.policy = TRUE)\n")
cat("  rownames(A) <- sa2_shp$SA2_MAIN16\n")
cat("  colnames(A) <- sa2_shp$SA2_MAIN16\n")
cat("─────────────────────────────────────────────────────────────────\n\n")

# --- Simulated adjacency for demo (inner Brisbane SA2s) ---
adj_pairs <- list(
  "301011001" = c("301011002","301011003","301011009"),
  "301011002" = c("301011001","301011009","301011010"),
  "301011003" = c("301011001","301011004","301011005"),
  "301011004" = c("301011003","301011005","301011006"),
  "301011005" = c("301011003","301011004","301011006"),
  "301011006" = c("301011004","301011005","301011007"),
  "301011007" = c("301011006","301011008","301011010"),
  "301011008" = c("301011007","301011009","301011010"),
  "301011009" = c("301011001","301011002","301011008"),
  "301011010" = c("301011002","301011007","301011008")
)

A <- matrix(0L, nrow = 10, ncol = 10,
            dimnames = list(all_sa2, all_sa2))

for (i in all_sa2) {
  for (j in adj_pairs[[i]]) {
    A[i, j] <- 1L
    A[j, i] <- 1L
  }
}

cat("=== ADJACENCY MATRIX A [n x n] (simulated from SA2 boundaries) ===\n")
cat("A[i,j] = 1 if SA2 i and SA2 j share a boundary\n\n")
print(A)
cat(sprintf("\nAdjacent pairs: %d  |  Mean neighbours per SA2: %.1f\n\n",
            sum(A) / 2L, mean(rowSums(A))))


# =============================================================================
# STEP 6: EXPORT
# =============================================================================

write_csv(census_sa2 |>
            select(ABSPID, AGEP, LFSP,
                   origin_SA2, POWP_DZN = POWP,
                   destination_SA2, IFPOWP),
          "census_filtered_sa2.csv")

write_csv(od_long,   "W_od_long.csv")

write_csv(
  as.data.frame(W_matrix) |> rownames_to_column("origin_SA2"),
  "W_od_matrix.csv"
)

write_csv(
  as.data.frame(A) |> rownames_to_column("SA2_i"),
  "A_adjacency_matrix.csv"
)

write_csv(sc_table,  "self_containment.csv")

cat("=== OUTPUT FILES ===\n")
cat("  census_filtered_sa2.csv  — filtered person records with SA2 POWP\n")
cat("  W_od_long.csv            — OD flows: origin_SA2 | destination_SA2 | W_ij\n")
cat("  W_od_matrix.csv          — OD matrix W [n x n]\n")
cat("  A_adjacency_matrix.csv   — Adjacency matrix A [n x n]\n")
cat("  self_containment.csv     — SC_SS, SC_DS, SC per SA2\n\n")

cat("=== PIPELINE SUMMARY ===\n")
cat("  PURP (SA2)  →  origin_SA2          [direct, no conversion needed]\n")
cat("  POWP (DZN)  →  join DZN_SA2_lookup →  destination_SA2\n")
cat("  group_by(origin_SA2, destination_SA2) |> count()  →  W[i,j]\n")
cat("  poly2nb(sa2_shapefile)              →  A[i,j]\n")
cat("  W and A are the two inputs to AdSA-ALMD and the Bayesian model\n")


# =============================================================================
# STEP 7: QUALITY ASSURANCE CHECKS
# Confirms that W and A are complete and internally consistent
# before passing them to AdSA-ALMD or the Bayesian model
# =============================================================================

cat("\n")
cat("=============================================================\n")
cat("  STEP 7: QUALITY ASSURANCE CHECKS\n")
cat("=============================================================\n\n")

qa_pass  <- TRUE   # track overall pass/fail
qa_log   <- list() # collect all check results

# ── Helper: print a labelled check result ─────────────────────────────────────
qa_check <- function(label, condition, detail = "") {
  status <- if (condition) "  PASS" else "  FAIL"
  msg    <- paste0(status, "  |  ", label)
  if (nchar(detail) > 0) msg <- paste0(msg, "\n         ", detail)
  cat(msg, "\n")
  if (!condition) qa_pass <<- FALSE
  invisible(condition)
}

cat("── W matrix (OD flows) ──────────────────────────────────────\n\n")

# QA-W1: W is square
qa_check(
  "W is square (n_rows == n_cols)",
  nrow(W_matrix) == ncol(W_matrix),
  sprintf("rows = %d, cols = %d", nrow(W_matrix), ncol(W_matrix))
)

# QA-W2: Row names and column names are identical and match the SA2 reference list
qa_check(
  "W row names match reference SA2 list",
  setequal(rownames(W_matrix), all_sa2),
  sprintf("in W but not reference: %s | in reference but not W: %s",
          paste(setdiff(rownames(W_matrix), all_sa2), collapse = ", ") |>
            (\(x) if (nchar(x) == 0) "none" else x)(),
          paste(setdiff(all_sa2, rownames(W_matrix)), collapse = ", ") |>
            (\(x) if (nchar(x) == 0) "none" else x)())
)

qa_check(
  "W col names match reference SA2 list",
  setequal(colnames(W_matrix), all_sa2),
  sprintf("in W but not reference: %s | in reference but not W: %s",
          paste(setdiff(colnames(W_matrix), all_sa2), collapse = ", ") |>
            (\(x) if (nchar(x) == 0) "none" else x)(),
          paste(setdiff(all_sa2, colnames(W_matrix)), collapse = ", ") |>
            (\(x) if (nchar(x) == 0) "none" else x)())
)

# QA-W3: Every SA2 in reference list appears as an origin (PURP) in raw data
origins_in_data <- unique(census_sa2$origin_SA2)
missing_origins <- setdiff(all_sa2, origins_in_data)
qa_check(
  "All SA2s appear as origin (PURP) in filtered data",
  length(missing_origins) == 0,
  if (length(missing_origins) > 0)
    paste("SA2s with zero residents commuting:", paste(missing_origins, collapse = ", "))
  else "All SA2s have at least one resident commuter"
)

# QA-W4: Every SA2 in reference list appears as a destination (POWP) in raw data
destinations_in_data <- unique(census_sa2$destination_SA2)
missing_destinations <- setdiff(all_sa2, destinations_in_data)
qa_check(
  "All SA2s appear as destination (POWP) in filtered data",
  length(missing_destinations) == 0,
  if (length(missing_destinations) > 0)
    paste("SA2s with zero workers commuting in:", paste(missing_destinations, collapse = ", "))
  else "All SA2s have at least one worker commuting in"
)

# QA-W5: No negative values
qa_check(
  "W contains no negative values",
  all(W_matrix >= 0),
  sprintf("min value = %d", min(W_matrix))
)

# QA-W6: No NA values
qa_check(
  "W contains no NA values",
  !any(is.na(W_matrix)),
  sprintf("NA count = %d", sum(is.na(W_matrix)))
)

# QA-W7: Row sums > 0 for all SA2s (every SA2 has at least one outgoing commuter)
zero_row_sa2 <- rownames(W_matrix)[rowSums(W_matrix) == 0]
qa_check(
  "All SA2 row sums > 0 (every SA2 has at least one resident commuter)",
  length(zero_row_sa2) == 0,
  if (length(zero_row_sa2) > 0)
    paste("SA2s with zero row sum:", paste(zero_row_sa2, collapse = ", "))
  else "All row sums > 0"
)

# QA-W8: Column sums > 0 for all SA2s (every SA2 has at least one incoming commuter)
zero_col_sa2 <- colnames(W_matrix)[colSums(W_matrix) == 0]
qa_check(
  "All SA2 col sums > 0 (every SA2 has at least one worker commuting in)",
  length(zero_col_sa2) == 0,
  if (length(zero_col_sa2) > 0)
    paste("SA2s with zero col sum:", paste(zero_col_sa2, collapse = ", "))
  else "All col sums > 0"
)

# QA-W9: Total flow in W matches total persons in filtered dataset
total_W      <- sum(W_matrix)
total_persons <- nrow(census_sa2)
qa_check(
  "Total flow in W equals total filtered person records",
  total_W == total_persons,
  sprintf("sum(W) = %d | nrow(census_sa2) = %d | difference = %d",
          total_W, total_persons, abs(total_W - total_persons))
)

# QA-W10: No duplicate SA2 codes in rows or columns
qa_check(
  "No duplicate SA2 codes in W rows",
  length(rownames(W_matrix)) == length(unique(rownames(W_matrix))),
  sprintf("duplicates: %s",
          paste(rownames(W_matrix)[duplicated(rownames(W_matrix))],
                collapse = ", "))
)

qa_check(
  "No duplicate SA2 codes in W cols",
  length(colnames(W_matrix)) == length(unique(colnames(W_matrix))),
  sprintf("duplicates: %s",
          paste(colnames(W_matrix)[duplicated(colnames(W_matrix))],
                collapse = ", "))
)

cat("\n── A matrix (adjacency) ─────────────────────────────────────\n\n")

# QA-A1: A is square
qa_check(
  "A is square (n_rows == n_cols)",
  nrow(A) == ncol(A),
  sprintf("rows = %d, cols = %d", nrow(A), ncol(A))
)

# QA-A2: A row and col names match W exactly
qa_check(
  "A row names match W row names",
  identical(sort(rownames(A)), sort(rownames(W_matrix))),
  sprintf("in A not W: %s | in W not A: %s",
          paste(setdiff(rownames(A), rownames(W_matrix)), collapse=", ") |>
            (\(x) if (nchar(x)==0) "none" else x)(),
          paste(setdiff(rownames(W_matrix), rownames(A)), collapse=", ") |>
            (\(x) if (nchar(x)==0) "none" else x)())
)

# QA-A3: A is symmetric
qa_check(
  "A is symmetric (A[i,j] == A[j,i] for all i,j)",
  isSymmetric(A),
  sprintf("max asymmetry: %d", max(abs(A - t(A))))
)

# QA-A4: A diagonal is zero (no self-adjacency)
qa_check(
  "A diagonal is all zeros (no self-adjacency)",
  all(diag(A) == 0),
  sprintf("non-zero diagonal entries: %d", sum(diag(A) != 0))
)

# QA-A5: A contains only 0 and 1
qa_check(
  "A contains only binary values (0 and 1)",
  all(A %in% c(0L, 1L)),
  sprintf("unique values in A: %s", paste(unique(as.vector(A)), collapse = ", "))
)

# QA-A6: Every SA2 has at least one neighbour (no isolated islands)
isolated_sa2 <- rownames(A)[rowSums(A) == 0]
qa_check(
  "No isolated SA2s (every SA2 has at least one neighbour in A)",
  length(isolated_sa2) == 0,
  if (length(isolated_sa2) > 0)
    paste("Isolated SA2s (zero neighbours):", paste(isolated_sa2, collapse = ", "))
  else sprintf("min neighbours = %d, max neighbours = %d",
               min(rowSums(A)), max(rowSums(A)))
)

# QA-A7: No NA values
qa_check(
  "A contains no NA values",
  !any(is.na(A)),
  sprintf("NA count = %d", sum(is.na(A)))
)

cat("\n── Cross-checks: W and A consistency ───────────────────────\n\n")

# QA-C1: W and A have the same dimension
qa_check(
  "W and A have the same dimension",
  nrow(W_matrix) == nrow(A) && ncol(W_matrix) == ncol(A),
  sprintf("W: %dx%d | A: %dx%d",
          nrow(W_matrix), ncol(W_matrix), nrow(A), ncol(A))
)

# QA-C2: W and A have identical SA2 sets
qa_check(
  "W and A cover identical SA2 sets",
  setequal(rownames(W_matrix), rownames(A)),
  sprintf("in W not A: %s | in A not W: %s",
          paste(setdiff(rownames(W_matrix), rownames(A)), collapse=", ") |>
            (\(x) if (nchar(x)==0) "none" else x)(),
          paste(setdiff(rownames(A), rownames(W_matrix)), collapse=", ") |>
            (\(x) if (nchar(x)==0) "none" else x)())
)

# QA-C3: Every adjacent pair in A has a non-zero flow in W or W^T
#         (not a hard requirement but useful to flag for review)
adj_edges     <- which(A == 1, arr.ind = TRUE)
zero_flow_adj <- sum(
  W_matrix[adj_edges] == 0 & t(W_matrix)[adj_edges] == 0
)
qa_check(
  "Adjacent SA2 pairs (A[i,j]=1) with zero flow in both W[i,j] and W[j,i]",
  zero_flow_adj == 0,
  if (zero_flow_adj > 0)
    sprintf("%d adjacent pairs have zero flow in both directions (review recommended)",
            zero_flow_adj)
  else "All adjacent pairs have at least one directional flow"
)

# ── Summary ───────────────────────────────────────────────────────────────────
cat("\n=============================================================\n")
if (qa_pass) {
  cat("  QA RESULT:  ALL CHECKS PASSED\n")
  cat("  W and A are ready for input to AdSA-ALMD\n")
} else {
  cat("  QA RESULT:  ONE OR MORE CHECKS FAILED\n")
  cat("  Review FAIL items above before proceeding\n")
}
cat("=============================================================\n\n")

# Append QA summary to export
qa_summary <- tibble(
  check = c(
    "W is square",
    "W row names match SA2 reference",
    "W col names match SA2 reference",
    "All SA2s appear as origin",
    "All SA2s appear as destination",
    "W no negative values",
    "W no NA values",
    "W all row sums > 0",
    "W all col sums > 0",
    "W total flow == person records",
    "W no duplicate row codes",
    "W no duplicate col codes",
    "A is square",
    "A row names match W",
    "A is symmetric",
    "A diagonal is zero",
    "A binary values only",
    "No isolated SA2s in A",
    "A no NA values",
    "W and A same dimension",
    "W and A same SA2 set",
    "Adjacent pairs have non-zero flow"
  ),
  result = c(
    nrow(W_matrix) == ncol(W_matrix),
    setequal(rownames(W_matrix), all_sa2),
    setequal(colnames(W_matrix), all_sa2),
    length(missing_origins) == 0,
    length(missing_destinations) == 0,
    all(W_matrix >= 0),
    !any(is.na(W_matrix)),
    length(zero_row_sa2) == 0,
    length(zero_col_sa2) == 0,
    total_W == total_persons,
    length(rownames(W_matrix)) == length(unique(rownames(W_matrix))),
    length(colnames(W_matrix)) == length(unique(colnames(W_matrix))),
    nrow(A) == ncol(A),
    identical(sort(rownames(A)), sort(rownames(W_matrix))),
    isSymmetric(A),
    all(diag(A) == 0),
    all(A %in% c(0L, 1L)),
    length(isolated_sa2) == 0,
    !any(is.na(A)),
    nrow(W_matrix) == nrow(A) && ncol(W_matrix) == ncol(A),
    setequal(rownames(W_matrix), rownames(A)),
    zero_flow_adj == 0
  ) |> ifelse("PASS", "FAIL")
)

write_csv(qa_summary, "QA_report.csv")
cat("QA report saved to: QA_report.csv\n")