# =============================================================================
# dev/abs_census_to_od.R
#
# Helpers for turning ABS Census person-microdata (PLIDA / DataLab) into the
# OD flow matrix `W` expected by run_adsa_pipeline(). NOT part of the public
# adsalmd package -- source this file manually in a DataLab session:
#
#   source("dev/abs_census_to_od.R")
#
# Base R only, no extra dependencies.
#
# Reference / worked demo with simulated data: test_census_convert.R
# =============================================================================

#' Filter ABS census microdata to employed, working-age, non-imputed records
#'
#' @param census_raw data.frame/tibble of person records with (by default)
#'   the exact ABS Census Dictionary column names AGEP, LFSP, IFPOWP.
#' @param age_col,age_min,age_max Age column and working-age bounds
#'   (default AGEP, 15-64).
#' @param lfsp_col,lfsp_employed Labour force status column and the codes
#'   treated as "employed" (default LFSP, c("1","2","3") = full-time,
#'   part-time, away from work).
#' @param ifpowp_col Imputation flag for place-of-work column (default
#'   IFPOWP). Only used when `require_non_imputed = TRUE`.
#' @param require_non_imputed Logical (default TRUE); if TRUE, keep only
#'   `ifpowp_col == 0` (observed, not imputed, place of work).
#' @return Filtered tibble, same columns as `census_raw`.
filter_census_employed <- function(census_raw,
                                    age_col = "AGEP", age_min = 15, age_max = 64,
                                    lfsp_col = "LFSP", lfsp_employed = c("1", "2", "3"),
                                    ifpowp_col = "IFPOWP", require_non_imputed = TRUE) {
  stopifnot(all(c(age_col, lfsp_col) %in% names(census_raw)))
  if (require_non_imputed) stopifnot(ifpowp_col %in% names(census_raw))

  keep <- census_raw[[age_col]] >= age_min & census_raw[[age_col]] <= age_max &
    census_raw[[lfsp_col]] %in% lfsp_employed

  if (require_non_imputed) {
    keep <- keep & census_raw[[ifpowp_col]] == 0L
  }

  census_raw[keep, , drop = FALSE]
}

#' Convert a filtered census extract's workplace field from DZN to SA2
#'
#' POWP is coded to Destination Zone (DZN) in raw Census data. This joins a
#' DZN-to-SA2 lookup (the ABS-published `DZN_SA2_<year>_AUST.csv`) to recover
#' the SA2 of the workplace, and renames the residence field for clarity.
#'
#' @param census_filtered Output of \code{\link{filter_census_employed}}
#'   (or equivalent), containing `purp_col` and `powp_col`.
#' @param dzn_sa2_lookup data.frame with a DZN column and an SA2 column
#'   (e.g. the ABS DZN_SA2 concordance).
#' @param purp_col,powp_col Residence (already SA2-coded) and workplace
#'   (DZN-coded) columns in `census_filtered` (default PURP, POWP).
#' @param dzn_lookup_col,sa2_lookup_col DZN and SA2 columns in
#'   `dzn_sa2_lookup` (default DZN_CODE_2016, SA2_MAINCODE_2016).
#' @param origin_name,dest_name Names to give the resulting origin/destination
#'   SA2 columns (default origin_SA2, destination_SA2).
#' @param drop_unmatched Logical (default TRUE); drop rows whose POWP has no
#'   match in the lookup (unmatched DZN codes) rather than keeping NA.
#' @return `census_filtered` with `origin_name` and `dest_name` SA2 columns
#'   added (the original `purp_col`/`powp_col` are left in place).
convert_powp_to_sa2 <- function(census_filtered, dzn_sa2_lookup,
                                 purp_col = "PURP", powp_col = "POWP",
                                 dzn_lookup_col = "DZN_CODE_2016",
                                 sa2_lookup_col = "SA2_MAINCODE_2016",
                                 origin_name = "origin_SA2",
                                 dest_name = "destination_SA2",
                                 drop_unmatched = TRUE) {
  stopifnot(all(c(purp_col, powp_col) %in% names(census_filtered)))
  stopifnot(all(c(dzn_lookup_col, sa2_lookup_col) %in% names(dzn_sa2_lookup)))

  sa2 <- dzn_sa2_lookup[[sa2_lookup_col]][
    match(census_filtered[[powp_col]], dzn_sa2_lookup[[dzn_lookup_col]])
  ]

  out <- census_filtered
  out[[origin_name]] <- out[[purp_col]]
  out[[dest_name]]   <- sa2

  if (drop_unmatched) out <- out[!is.na(out[[dest_name]]), , drop = FALSE]
  out
}

#' Build the square OD flow matrix `W` from origin/destination SA2 columns
#'
#' `W[i, j]` is the count of persons resident in SA2 `i` working in SA2 `j`.
#' Every SA2 in `sa2_codes` appears as both a row and a column, in the order
#' given, even if it has zero flow in the data -- required so `W` and the
#' adjacency matrix from \code{\link{build_adjacency_from_shapefile}} line up.
#'
#' @param census_sa2 data.frame with `origin_col` and `dest_col` SA2 columns
#'   (typically the output of \code{\link{convert_powp_to_sa2}}).
#' @param sa2_codes Character vector, the canonical, ordered set of SA2 codes
#'   for this run. Fixes the row/column order of the returned matrix -- pass
#'   the same vector to \code{\link{build_adjacency_from_shapefile}}.
#' @param origin_col,dest_col Column names in `census_sa2` (default
#'   origin_SA2, destination_SA2).
#' @return Integer `length(sa2_codes) x length(sa2_codes)` matrix, dimnames
#'   set to `sa2_codes` on both margins.
build_od_matrix <- function(census_sa2, sa2_codes,
                             origin_col = "origin_SA2", dest_col = "destination_SA2") {
  stopifnot(all(c(origin_col, dest_col) %in% names(census_sa2)))
  stopifnot(!anyDuplicated(sa2_codes))

  n <- length(sa2_codes)
  W <- matrix(0L, n, n, dimnames = list(sa2_codes, sa2_codes))

  o <- match(census_sa2[[origin_col]], sa2_codes)
  d <- match(census_sa2[[dest_col]], sa2_codes)
  valid <- !is.na(o) & !is.na(d)
  if (any(!valid)) {
    warning(sprintf(
      "%d records reference an SA2 outside `sa2_codes` and were dropped",
      sum(!valid)
    ))
  }

  tab <- table(o[valid], d[valid])
  W[as.integer(rownames(tab)), as.integer(colnames(tab))] <- tab
  storage.mode(W) <- "integer"
  W
}

#' Compute self-containment diagnostics from an OD matrix
#'
#' Useful for choosing \code{SC_min} for \code{\link{run_adsa_pipeline}()}.
#'
#' @param W Square numeric OD matrix with matching row/column dimnames
#'   (e.g. from \code{\link{build_od_matrix}}).
#' @return tibble with one row per SA2: `sa2`, `W_ii` (internal flow),
#'   `W_iG` (row sum / outflow), `W_Gi` (col sum / inflow), `SC_SS`
#'   (supply-side = W_ii / W_iG), `SC_DS` (demand-side = W_ii / W_Gi), and
#'   `SC` = pmin(SC_SS, SC_DS).
compute_self_containment <- function(W) {
  stopifnot(is.matrix(W), nrow(W) == ncol(W))
  stopifnot(identical(rownames(W), colnames(W)))

  W_ii <- diag(W)
  W_iG <- rowSums(W)
  W_Gi <- colSums(W)
  SC_SS <- ifelse(W_iG > 0, W_ii / W_iG, NA_real_)
  SC_DS <- ifelse(W_Gi > 0, W_ii / W_Gi, NA_real_)
  SC    <- pmin(SC_SS, SC_DS)

  data.frame(
    sa2   = rownames(W),
    W_ii  = W_ii,
    W_iG  = W_iG,
    W_Gi  = W_Gi,
    SC_SS = SC_SS,
    SC_DS = SC_DS,
    SC    = SC,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}
