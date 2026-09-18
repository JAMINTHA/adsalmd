# =============================================================================
# dev/lma_result_to_abs_correspondence.R
#
# Convert a run_adsa_pipeline() result back into an ABS-style correspondence
# table: SA2 <-> Local Labour Market Area (LMA), in the same column layout as
# an official ABS geography correspondence file (e.g.
# "SA2 (2016) to SA4 (2016) Correspondence"): one row per SA2, its assigned
# higher-level region code/name, and RATIO_FROM_TO / RATIO_TO_FROM.
#
# NOT part of the public adsalmd package -- source manually:
#
#   source("dev/lma_result_to_abs_correspondence.R")
#
# Base R only, no extra dependencies.
# =============================================================================

#' Turn a pipeline partition into an ABS-style SA2-to-LMA correspondence table
#'
#' \code{\link{run_adsa_pipeline}()}'s \code{best_sol} (and \code{init_sol})
#' is a plain integer vector: \code{best_sol[i]} is the LMA label of the SA2
#' at position \code{i} in \code{W}/\code{adj} (i.e. \code{sa2_codes[i]}).
#' This maps that back into the row-per-SA2 layout ABS correspondence files
#' use, with re-numbered, zero-padded LMA codes and, since this is a strict
#' (non-overlapping) partition, \code{RATIO_FROM_TO = 1} for every row.
#'
#' @param sol Integer vector of length \code{length(sa2_codes)}, an
#'   assignment vector from the pipeline (\code{res$best_sol},
#'   \code{res$init_sol}, or any partition in \code{sa2_codes} order).
#' @param sa2_codes Character vector, the SA2 codes in the same order used to
#'   build \code{W}/\code{adj} for this run (i.e. \code{rownames(W)}).
#' @param sa2_names Character vector, same length/order as `sa2_codes`, or
#'   `NULL` (default) to omit the name column.
#' @param population Named or unnamed numeric vector of SA2 population
#'   (aligned to `sa2_codes` by position, or by name if named), used only to
#'   compute \code{RATIO_TO_FROM} (each SA2's population share of its LMA's
#'   total). If `NULL` (default), \code{row_W} (commuter counts) is used as
#'   a proxy if supplied via `row_W`, otherwise every SA2 is weighted equally
#'   within its LMA and a warning is issued.
#' @param row_W Numeric vector aligned to `sa2_codes`, typically
#'   \code{rowSums(W)}. Fallback population proxy when `population` is not
#'   supplied; ignored if `population` is given.
#' @param lma_code_prefix Character, prefix for generated LMA codes (default
#'   `"LMA"`; codes are `sprintf("\%s\%03d", lma_code_prefix, k)` for the
#'   `k`-th LMA in ascending order of `sol`'s original labels).
#' @param sa2_year,lma_year Character, the ASGS year to embed in column
#'   names, matching ABS convention `SA2_MAINCODE_<year>` /
#'   `LMA_CODE_<year>` (default `"2016"` for both).
#' @return data.frame, one row per SA2, columns:
#'   \describe{
#'     \item{\code{SA2_MAINCODE_<sa2_year>}}{the SA2 code.}
#'     \item{\code{SA2_NAME_<sa2_year>}}{if `sa2_names` supplied.}
#'     \item{\code{LMA_CODE_<lma_year>}}{zero-padded generated LMA code.}
#'     \item{\code{LMA_NAME_<lma_year>}}{`"Local Labour Market Area <k>"`.}
#'     \item{\code{RATIO_FROM_TO}}{always 1 (each SA2 belongs wholly to one
#'       LMA).}
#'     \item{\code{RATIO_TO_FROM}}{this SA2's share of its LMA's total
#'       population/weight.}
#'   }
#'   Sorted by `SA2_MAINCODE_<sa2_year>`, matching ABS correspondence file
#'   convention.
build_sa2_lma_correspondence <- function(sol, sa2_codes, sa2_names = NULL,
                                          population = NULL, row_W = NULL,
                                          lma_code_prefix = "LMA",
                                          sa2_year = "2016", lma_year = "2016") {
  stopifnot(length(sol) == length(sa2_codes))
  if (!is.null(sa2_names)) stopifnot(length(sa2_names) == length(sa2_codes))

  weight <- NULL
  if (!is.null(population)) {
    weight <- if (!is.null(names(population))) {
      as.numeric(population[sa2_codes])
    } else {
      stopifnot(length(population) == length(sa2_codes))
      as.numeric(population)
    }
  } else if (!is.null(row_W)) {
    stopifnot(length(row_W) == length(sa2_codes))
    weight <- as.numeric(row_W)
  } else {
    warning(
      "No `population` or `row_W` supplied -- RATIO_TO_FROM will weight ",
      "every SA2 equally within its LMA rather than by population."
    )
    weight <- rep(1, length(sa2_codes))
  }

  labels_sorted <- sort(unique(sol))
  lma_index <- match(sol, labels_sorted)
  lma_code  <- sprintf("%s%03d", lma_code_prefix, lma_index)
  lma_name  <- sprintf("Local Labour Market Area %d", lma_index)

  lma_total_weight <- ave(weight, lma_code, FUN = sum)
  ratio_to_from <- ifelse(lma_total_weight > 0, weight / lma_total_weight, NA_real_)

  out <- data.frame(
    SA2_MAINCODE = sa2_codes,
    stringsAsFactors = FALSE
  )
  names(out) <- sprintf("SA2_MAINCODE_%s", sa2_year)

  if (!is.null(sa2_names)) {
    out[[sprintf("SA2_NAME_%s", sa2_year)]] <- sa2_names
  }

  out[[sprintf("LMA_CODE_%s", lma_year)]] <- lma_code
  out[[sprintf("LMA_NAME_%s", lma_year)]] <- lma_name
  out$RATIO_FROM_TO <- 1
  out$RATIO_TO_FROM <- round(ratio_to_from, 6)

  out[order(out[[sprintf("SA2_MAINCODE_%s", sa2_year)]]), , drop = FALSE]
}

#' Attach per-LMA summary metrics (population, self-containment) to a
#' correspondence table
#'
#' Convenience join of \code{\link{build_sa2_lma_correspondence}()}'s output
#' with \code{run_adsa_pipeline()}'s \code{lma_metrics} (or any
#' \code{build_lma_df()} output), matching on cluster label so you get
#' population/self-containment alongside every SA2's LMA assignment.
#'
#' @param correspondence Output of \code{\link{build_sa2_lma_correspondence}()}.
#' @param sol Same `sol` vector passed to
#'   \code{\link{build_sa2_lma_correspondence}()} (needed to recover the
#'   original cluster label -> generated LMA code mapping).
#' @param lma_metrics data.frame from \code{build_lma_df()} (e.g.
#'   \code{res$lma_metrics}), with a `cluster` column matching `sol`'s labels.
#' @param lma_code_col Name of the LMA code column in `correspondence`
#'   (default `"LMA_CODE_2016"` -- match whatever `lma_year` you used).
#' @return `correspondence` with `pop`, `scss`, `scds`, `sc`, `n_bgus` columns
#'   from `lma_metrics` appended (prefixed `LMA_` to avoid confusion with
#'   per-SA2 figures).
attach_lma_metrics <- function(correspondence, sol, lma_metrics,
                                lma_code_col = "LMA_CODE_2016") {
  stopifnot(lma_code_col %in% names(correspondence))
  labels_sorted <- sort(unique(sol))
  # Recover whatever prefix was actually used to build `correspondence`, in
  # case a non-default lma_code_prefix was passed to
  # build_sa2_lma_correspondence().
  prefix <- sub("[0-9]+$", "", correspondence[[lma_code_col]][1])
  code_by_label <- stats::setNames(
    sprintf("%s%03d", prefix, seq_along(labels_sorted)), labels_sorted
  )

  lma_metrics$.lma_code <- code_by_label[as.character(lma_metrics$cluster)]
  metric_cols <- setdiff(names(lma_metrics), c("cluster", ".lma_code"))
  names(lma_metrics)[match(metric_cols, names(lma_metrics))] <- paste0("LMA_", metric_cols)

  merge(
    correspondence, lma_metrics,
    by.x = lma_code_col, by.y = ".lma_code",
    all.x = TRUE, sort = FALSE
  )
}

#' Join the correspondence table onto SA2 geometry and write a map-ready
#' shapefile
#'
#' \code{\link{build_sa2_lma_correspondence}()} produces a plain data.frame --
#' the ABS correspondence file convention, and what \code{write.csv()}
#' exports. That's not map-ready on its own: nothing ties `LMA_CODE_<year>`
#' back to SA2 *geometry*. This left-joins a correspondence table (or
#' \code{\link{attach_lma_metrics}()}'s output) onto an SA2 boundary
#' shapefile by SA2 code and writes the result out as a shapefile, so it
#' opens directly in QGIS/ArcGIS with the LMA assignment (and any attached
#' metrics) as attribute columns you can symbolise or dissolve by.
#'
#' @param correspondence data.frame, one row per SA2 -- typically
#'   \code{corr_full} / \code{\link{attach_lma_metrics}()}'s output, i.e.
#'   already covering every SA2 in \code{sa2_codes_full} (see
#'   \code{USER_GUIDE.md} step 5.1 -- an SA2 missing from `correspondence`
#'   is simply dropped from the output, not an error).
#' @param shapefile_path Path to the SA2 boundary shapefile to join onto --
#'   the same one \code{sa2_codes_full} was built from (the *full* shapefile
#'   from step 5.1, not the trimmed one from step 2, if you did that step).
#' @param out_path Path to write the joined shapefile to (`.shp`; sidecar
#'   `.dbf`/`.shx`/`.prj` files are written alongside it automatically).
#' @param sa2_code_col Column in `correspondence` holding the SA2 code
#'   (default `"SA2_MAINCODE_2016"`, matching
#'   \code{\link{build_sa2_lma_correspondence}()}'s default `sa2_year`).
#' @param sa2_id_col Column in the shapefile holding the SA2 code (default
#'   `"SA2_MAIN16"`, matching \code{\link{build_adjacency_from_shapefile}()}'s
#'   default).
#' @return Invisibly, the joined `sf` object (also written to `out_path`).
#' @details The shapefile format truncates column names to 10 characters (a
#'   DBF limitation) -- if `correspondence` has several long or similarly-
#'   prefixed columns (e.g. from \code{attach_lma_metrics()}'s `LMA_`
#'   columns), check `names(sf::st_read(out_path))` afterwards to see what
#'   they became, and rename before calling this if that's ambiguous.
export_lma_shapefile <- function(correspondence, shapefile_path, out_path,
                                  sa2_code_col = "SA2_MAINCODE_2016",
                                  sa2_id_col = "SA2_MAIN16") {
  stopifnot(requireNamespace("sf", quietly = TRUE))
  stopifnot(sa2_code_col %in% names(correspondence))

  sa2_shp <- sf::st_read(shapefile_path, quiet = TRUE)
  if (!sa2_id_col %in% names(sa2_shp)) {
    stop(sprintf(
      "`%s` not found in shapefile. Available columns: %s",
      sa2_id_col, paste(names(sa2_shp), collapse = ", ")
    ))
  }

  missing <- setdiff(correspondence[[sa2_code_col]], sa2_shp[[sa2_id_col]])
  if (length(missing) > 0) {
    stop(sprintf(
      "%d SA2 code(s) in `correspondence` not found in shapefile: %s",
      length(missing), paste(missing, collapse = ", ")
    ))
  }

  joined <- merge(
    sa2_shp, correspondence,
    by.x = sa2_id_col, by.y = sa2_code_col,
    all.x = FALSE   # keep only the SA2s actually in `correspondence`
  )

  sf::st_write(joined, out_path, delete_dsn = TRUE, quiet = TRUE)
  invisible(joined)
}
