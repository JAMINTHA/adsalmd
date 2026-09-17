# =============================================================================
# dev/validate_lma_inputs.R
#
# QA checks for a (W, adj) pair before handing it to run_adsa_pipeline().
# NOT part of the public adsalmd package -- source this file manually:
#
#   source("dev/validate_lma_inputs.R")
#
# Base R only, no extra dependencies.
# =============================================================================

#' Validate an OD matrix / adjacency matrix pair before running the pipeline
#'
#' Checks structural and consistency properties of `W` and `adj` that
#' \code{\link{run_adsa_pipeline}()} assumes but does not itself verify --
#' matching dimensions and SA2 sets/order, no NA or negative values, `adj`
#' symmetric/binary/zero-diagonal, no isolated units, and (if `census_sa2` is
#' supplied) that `W` accounts for every filtered person record and that
#' every adjacent pair has flow in at least one direction.
#'
#' @param W Square numeric OD matrix with row/column dimnames (e.g. from
#'   `build_od_matrix()`).
#' @param adj Square integer/logical adjacency matrix with row/column
#'   dimnames (e.g. from `build_adjacency_from_shapefile()`).
#' @param census_sa2 Optional data.frame with `origin_col`/`dest_col` SA2
#'   columns -- the filtered, SA2-coded person records `W` was built from.
#'   When supplied, adds checks that every SA2 appears as an origin and a
#'   destination, and that `sum(W)` equals `nrow(census_sa2)`.
#' @param origin_col,dest_col Column names in `census_sa2` (default
#'   origin_SA2, destination_SA2). Ignored if `census_sa2` is NULL.
#' @param verbose Logical (default TRUE); print PASS/FAIL lines as checks run.
#' @return Invisibly, a list with `pass` (logical, TRUE iff every check
#'   passed) and `report` (data.frame of `check`, `result`, `detail`).
validate_W_adj <- function(W, adj, census_sa2 = NULL,
                            origin_col = "origin_SA2", dest_col = "destination_SA2",
                            verbose = TRUE) {
  results <- list()

  add_check <- function(label, condition, detail = "") {
    results[[length(results) + 1]] <<- data.frame(
      check = label, result = isTRUE(condition), detail = detail,
      stringsAsFactors = FALSE
    )
    if (verbose) {
      status <- if (isTRUE(condition)) "  PASS" else "  FAIL"
      msg <- paste0(status, "  |  ", label)
      if (nchar(detail) > 0) msg <- paste0(msg, "\n         ", detail)
      message(msg)
    }
    invisible(condition)
  }

  none_if_empty <- function(x) if (length(x) == 0) "none" else paste(x, collapse = ", ")

  if (verbose) message("-- W matrix (OD flows) --")

  add_check("W is square", nrow(W) == ncol(W),
             sprintf("rows = %d, cols = %d", nrow(W), ncol(W)))
  add_check("W has row and column names", !is.null(rownames(W)) && !is.null(colnames(W)))
  add_check("W row names match column names",
            identical(rownames(W), colnames(W)),
            sprintf("in rows not cols: %s | in cols not rows: %s",
                    none_if_empty(setdiff(rownames(W), colnames(W))),
                    none_if_empty(setdiff(colnames(W), rownames(W)))))
  add_check("W contains no negative values", all(W >= 0, na.rm = TRUE),
             sprintf("min value = %s", suppressWarnings(min(W))))
  add_check("W contains no NA values", !anyNA(W),
             sprintf("NA count = %d", sum(is.na(W))))
  add_check("No duplicate SA2 codes in W rows",
            !anyDuplicated(rownames(W)),
            sprintf("duplicates: %s", none_if_empty(rownames(W)[duplicated(rownames(W))])))

  zero_row <- rownames(W)[rowSums(W) == 0]
  add_check("All W row sums > 0 (every SA2 has an outgoing commuter)",
            length(zero_row) == 0,
            sprintf("zero row sum: %s", none_if_empty(zero_row)))

  zero_col <- colnames(W)[colSums(W) == 0]
  add_check("All W col sums > 0 (every SA2 has an incoming commuter)",
            length(zero_col) == 0,
            sprintf("zero col sum: %s", none_if_empty(zero_col)))

  if (verbose) message("-- adj matrix (adjacency) --")

  add_check("adj is square", nrow(adj) == ncol(adj),
             sprintf("rows = %d, cols = %d", nrow(adj), ncol(adj)))
  add_check("adj is symmetric", isSymmetric(unname(adj)),
             sprintf("max asymmetry = %s", suppressWarnings(max(abs(adj - t(adj))))))
  add_check("adj diagonal is all zero", all(diag(adj) == 0),
             sprintf("non-zero diagonal entries = %d", sum(diag(adj) != 0)))
  add_check("adj contains only 0/1", all(adj %in% c(0, 1)),
             sprintf("unique values = %s", paste(unique(as.vector(adj)), collapse = ", ")))
  add_check("adj contains no NA values", !anyNA(adj),
             sprintf("NA count = %d", sum(is.na(adj))))

  isolated <- rownames(adj)[rowSums(adj) == 0]
  add_check("No isolated units in adj (every unit has >=1 neighbour)",
            length(isolated) == 0,
            sprintf("isolated: %s", none_if_empty(isolated)))

  if (verbose) message("-- Cross-checks: W and adj --")

  add_check("W and adj have the same dimension",
            nrow(W) == nrow(adj) && ncol(W) == ncol(adj),
            sprintf("W: %dx%d | adj: %dx%d", nrow(W), ncol(W), nrow(adj), ncol(adj)))
  add_check("W and adj cover identical SA2 sets",
            setequal(rownames(W), rownames(adj)),
            sprintf("in W not adj: %s | in adj not W: %s",
                    none_if_empty(setdiff(rownames(W), rownames(adj))),
                    none_if_empty(setdiff(rownames(adj), rownames(W)))))
  add_check("W and adj have identical SA2 order",
            identical(rownames(W), rownames(adj)),
            "required so row_W/col_W and adj index the same units -- reorder one to match the other if this fails")

  if (nrow(W) == nrow(adj) && identical(rownames(W), rownames(adj))) {
    adj_edges <- which(adj == 1, arr.ind = TRUE)
    zero_flow_adj <- sum(W[adj_edges] == 0 & t(W)[adj_edges] == 0)
    add_check("Every adjacent pair has flow in at least one direction",
              zero_flow_adj == 0,
              sprintf("%d adjacent pairs have zero flow both ways (review, not necessarily an error)",
                      zero_flow_adj))
  }

  if (!is.null(census_sa2)) {
    if (verbose) message("-- Cross-checks: W vs census_sa2 --")
    stopifnot(all(c(origin_col, dest_col) %in% names(census_sa2)))

    missing_origin <- setdiff(rownames(W), unique(census_sa2[[origin_col]]))
    add_check("Every SA2 appears as an origin in census_sa2",
              length(missing_origin) == 0,
              sprintf("no resident commuters: %s", none_if_empty(missing_origin)))

    missing_dest <- setdiff(colnames(W), unique(census_sa2[[dest_col]]))
    add_check("Every SA2 appears as a destination in census_sa2",
              length(missing_dest) == 0,
              sprintf("no incoming commuters: %s", none_if_empty(missing_dest)))

    add_check("sum(W) equals nrow(census_sa2)",
              sum(W) == nrow(census_sa2),
              sprintf("sum(W) = %d | nrow(census_sa2) = %d", sum(W), nrow(census_sa2)))
  }

  report <- do.call(rbind, results)
  pass <- all(report$result)

  if (verbose) {
    message(if (pass) "QA RESULT: ALL CHECKS PASSED" else "QA RESULT: ONE OR MORE CHECKS FAILED")
  }

  invisible(list(pass = pass, report = report))
}
