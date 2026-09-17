# ABS Census → `run_adsa_pipeline()` inputs (internal, not part of the package)

For the full walkthrough — including running the pipeline itself and
converting its output back into an ABS-style correspondence table — see
[`USER_GUIDE.md`](USER_GUIDE.md). This file is the shorter reference for
just the census-to-`W`/`adj` conversion step.

Helper functions for converting ABS Census person-microdata (PLIDA / DataLab)
into the `W` (OD flow matrix) and `adj` (adjacency matrix) arguments required
by `run_adsa_pipeline()`.

**This folder is not part of the installable `adsalmd` package.** It is
excluded via `.Rbuildignore`, none of its functions are exported or
`roxygen`-documented into `man/`. `abs_census_to_od.R` and
`validate_lma_inputs.R` are base R only; `abs_adjacency_from_shapefile.R`
needs `sf` and `spdep` (not declared in `DESCRIPTION` — install them
yourself). `source()` these files manually in a DataLab session:

```r
install.packages(c("sf", "spdep"))  # once, only needed for the adjacency builder

source("dev/abs_census_to_od.R")
source("dev/abs_adjacency_from_shapefile.R")
source("dev/validate_lma_inputs.R")
```

See `test_census_convert.R` (repo root) for a fully worked, self-contained
demo against simulated census data — it shows the same steps inline before
they were factored out into these reusable functions.

## Expected raw input

Person-level microdata with the exact ABS Census Dictionary column names:

| Column   | Meaning                                              |
|----------|-------------------------------------------------------|
| `ABSPID` | Person identifier                                     |
| `AGEP`   | Age                                                    |
| `LFSP`   | Labour force status (`1`/`2`/`3` = employed)           |
| `PURP`   | Place of Usual Residence, coded to SA2 directly        |
| `POWP`   | Place of Work, coded to Destination Zone (DZN)         |
| `IFPOWP` | Imputation flag for POWP (`0` = observed, `1` = imputed) |

Plus two reference datasets you provide:

- A **DZN → SA2 concordance** (ABS-published `DZN_SA2_<year>_AUST.csv`), to
  resolve `POWP`'s DZN code to an SA2.
- An **SA2 boundary shapefile** (ABS ASGS digital boundary file), to derive
  adjacency.

## Pipeline

```r
# 1. Filter to employed, working-age, non-imputed records
census_filtered <- filter_census_employed(census_raw)

# 2. Resolve POWP (DZN) -> destination SA2
census_sa2 <- convert_powp_to_sa2(census_filtered, dzn_sa2_lookup)

# 3. Canonical, ordered SA2 code list for this run -- W and adj must agree
#    on this exact order.
sa2_codes <- sort(unique(c(census_sa2$origin_SA2, census_sa2$destination_SA2)))

# 4. Build the OD matrix
W <- build_od_matrix(census_sa2, sa2_codes)

# 5. (Optional) self-containment diagnostics, to help choose SC_min
sc_table <- compute_self_containment(W)

# 6. Build adjacency from an SA2 boundary shapefile
adj <- build_adjacency_from_shapefile(
  "SA2_2016_AUST_GDA2020.shp", sa2_codes, sa2_id_col = "SA2_MAIN16"
)

# 7. QA before handing off to the pipeline
qa <- validate_W_adj(W, adj, census_sa2 = census_sa2)
stopifnot(qa$pass)

# 8. Feed run_adsa_pipeline()
row_W <- rowSums(W)
col_W <- colSums(W)

res <- run_adsa_pipeline(
  W = W, adj = adj, row_W = row_W, col_W = col_W,
  SC_min = 0.55, Pop_min = 200, Pop_max = 500, r = 3,
  seed = 2026
)
```

## Function reference

- `filter_census_employed()` — age/labour-force/imputation filter (`dev/abs_census_to_od.R`).
- `convert_powp_to_sa2()` — DZN → SA2 join for the workplace field (`dev/abs_census_to_od.R`).
- `build_od_matrix()` — aggregates origin/destination counts into the square `W` matrix, fixed SA2 order (`dev/abs_census_to_od.R`).
- `compute_self_containment()` — per-SA2 `SC_SS`/`SC_DS`/`SC` from `W`, to inform `SC_min` (`dev/abs_census_to_od.R`).
- `build_adjacency_from_shapefile()` — queen/rook contiguity `adj` matrix from an SA2 shapefile via `sf`/`spdep`, reordered to match `W` (`dev/abs_adjacency_from_shapefile.R`).
- `validate_W_adj()` — structural and cross-consistency QA checks on `(W, adj)` before running the pipeline (`dev/validate_lma_inputs.R`).

## Notes

- **SA2 order matters.** `build_od_matrix()` and `build_adjacency_from_shapefile()`
  both take the same `sa2_codes` vector and use it verbatim (not re-sorted)
  as the row/column order — pass the identical vector to both so `W`, `adj`,
  and downstream `row_W`/`col_W` all index the same units in the same order.
  `validate_W_adj()` checks this explicitly.
- `run_adsa_pipeline()` has no direct "population" argument — population is
  read from `row_W`/`col_W` (i.e. it must already be embedded in how you
  construct `W`). If your delineation population differs from commuting
  flow (e.g. resident population rather than commuter counts), that's a
  modelling decision outside the scope of these converters.
- These functions assume 2016-vintage column names
  (`DZN_CODE_2016`/`SA2_MAINCODE_2016`, `SA2_MAIN16`) as defaults, matching
  `test_census_convert.R`. Pass the matching column names explicitly for
  other Census years/ASGS editions (e.g. `SA2_CODE21`).
