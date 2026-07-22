# ============================================================================
# 06_SENSITIVITY_START.R
# ----------------------------------------------------------------------------
# Robustness of conflict-attributable YLL to the coded conflict-ONSET year.
#
# The onset year is a researcher degree of freedom: shifting it moves the
# pre/post boundary, which re-anchors the counterfactual for EVERY estimator at
# once -- the ITS trend origin, the SC i_time and optimization window, the CS
# treatment cohort, and the baseline pre-window. This module re-runs the full
# gap -> YLL chain with conflict_start shifted by each value in `shifts`
# (default -1, 0, +1) and reports how the bottom line moves.
#
# WHY THE DIFFERENCES ARE CLEAN (not noise):
#   run_full_yll() builds the structural CFR/VE draws from a fixed shared seed
#   keyed by COUNTRY (not by timing), and we pass the SAME `seed` for every
#   shift. The common random numbers are therefore identical across shifts, so
#   a difference in the YLL total is attributable to the onset move, not to
#   Monte Carlo variation. (n_sim can be smaller here than the headline run:
#   we want a stable point estimate of the SHIFT, not a tight CI.)
#
# COST: this re-runs gap estimation (incl. SC permutation inference and DiD) and
# the YLL Monte Carlo once PER SHIFT. It is opt-in -- source this file (it only
# defines functions) and call run_start_sensitivity() deliberately, or set
# RUN_START_SENSITIVITY <- TRUE in 00_run_all.R.
#
# Source AFTER 05_figures_and_tables.R (uses theme_nm / nm_palette / METHOD_ORDER
# and the gap + YLL + aggregation functions).
# ============================================================================

# ---------------------------------------------------------------------------
# Shift the onset year and re-derive the dependent pre-window columns.
# ---------------------------------------------------------------------------
# By default ONLY the onset moves (the pre/post boundary); conflict_end is held
# fixed, because the onset -- not the end -- is the arbitrary coding choice in
# question. This means start-1 adds one conflict-year at the front and start+1
# removes the first year, which is the real consequence of mis-dating onset.
# Set shift_end = TRUE for the window-PRESERVING variant (start and end move
# together, holding the post-window length constant), which isolates the
# counterfactual-alignment effect from the accumulation-length effect.
shift_conflict_info <- function(conflict_df, delta, shift_end = FALSE,
                                baseline_window = 3L) {
  out <- conflict_df
  out$conflict_start <- out$conflict_start + as.integer(delta)
  if (shift_end) out$conflict_end <- out$conflict_end + as.integer(delta)
  out$preconflict_year   <- out$conflict_start - 1L
  out$preconflict_window <- purrr::map(out$preconflict_year,
                                       ~ (.x - (baseline_window - 1L)):.x)
  # Drop any country the shift makes degenerate (need at least a 2-year window).
  bad <- out$conflict_start >= out$conflict_end
  if (any(bad)) {
    message("  start-shift ", sprintf("%+d", delta), ": dropping ", sum(bad),
            " country/ies with a degenerate window (",
            paste(out$country[bad], collapse = ", "), ").")
    out <- out[!bad, , drop = FALSE]
  }
  out
}

# ---------------------------------------------------------------------------
# Run the gap -> YLL chain across onset shifts and collect comparable summaries.
# ---------------------------------------------------------------------------
# Returns:
#   $by_method  - global all-disease YLL total per (method, start_shift),
#                 complete-case within shift (same logic as Figure 4 panel C),
#                 plus a per-method % change relative to shift 0.
#   $by_country - per-country YLL total (summed over disease) per
#                 (country, method, start_shift) for localising the sensitivity.
#   $meta       - the call settings.
run_start_sensitivity <- function(coverage_df, covariates,
                                  conflict_df   = conflict_info,
                                  shifts        = c(-2L, -1L, 0L, 1L, 2L),
                                  methods       = ALL_METHODS,
                                  shift_end     = FALSE,
                                  n_sim         = 2000L,
                                  seed          = 42,
                                  framing_focus = if (exists("FRAMING_HEADLINE")) FRAMING_HEADLINE else "campaign_topup",
                                  catchup_focus = 0) {
  stopifnot(0L %in% as.integer(shifts))   # need the reference run to compute % change
  by_method <- list(); by_country <- list(); draws_by_shift <- list()
  
  for (d in as.integer(shifts)) {
    message("\n================  START-SHIFT ", sprintf("%+d yr", d),
            "  ================")
    cinfo_d <- shift_conflict_info(conflict_df, d, shift_end = shift_end)
    gaps_d  <- estimate_all_gaps_plus(coverage_df, cinfo_d, methods = methods)
    yll_d   <- run_full_yll(gaps_d, covariates, conflict_df = cinfo_d,
                            n_sim = n_sim, seed = seed, methods = NULL)
    
    # Retain per-shift draws (headline scenario only) for the COMPOSITION-STABLE
    # fixed-cell aggregation after the loop (.fixed_cell_onset). Without this the
    # per-shift complete-case cell set can change across shifts (a later onset can
    # make a donor-poor country's SC estimable), contaminating the %-change with
    # composition rather than the onset effect (the cells_match_ref=FALSE rows).
    dd <- yll_d$draws %>%
      dplyr::filter(framing == framing_focus, catchup == catchup_focus)
    keep_d <- vapply(dd$d_undisc, function(x) !is.null(x) && length(x) > 0, logical(1))
    draws_by_shift[[as.character(d)]] <- dd[keep_d, , drop = FALSE]
    
    gm <- aggregate_global_by_method_draws(
      yll_d$draws, framing_focus = framing_focus, catchup_focus = catchup_focus)
    if (nrow(gm) > 0)
      by_method[[length(by_method) + 1]] <- gm %>% dplyr::mutate(start_shift = d)
    
    # Per-country totals (sum over disease) under the headline scenario.
    cc <- aggregate_country_yll(yll_d$summary)
    if (nrow(cc) > 0) {
      cc <- cc %>%
        dplyr::filter(framing == framing_focus, catchup == catchup_focus) %>%
        dplyr::group_by(country, method) %>%
        dplyr::summarise(
          yll_undisc_total = sum(yll_undisc_total_mean, na.rm = TRUE),
          yll_disc_total   = sum(yll_disc_total_mean,   na.rm = TRUE),
          .groups = "drop") %>%
        dplyr::mutate(start_shift = d)
      by_country[[length(by_country) + 1]] <- cc
    }
  }
  
  by_method  <- dplyr::bind_rows(by_method)
  by_country <- dplyr::bind_rows(by_country)
  
  # Per-method % change vs the reference (shift 0). Flag composition changes:
  # complete-case cell counts can differ across shifts (a shift may cost SC a
  # donor-poor country), so a moving n_cells means the totals are not on an
  # identical cell set -- read the % change with that in mind.
  if (nrow(by_method) > 0) {
    ref <- by_method %>%
      dplyr::filter(start_shift == 0L) %>%
      dplyr::select(method, ref_undisc = yll_undisc_mean,
                    ref_disc = yll_disc_mean, ref_cells = n_cells)
    by_method <- by_method %>%
      dplyr::left_join(ref, by = "method") %>%
      dplyr::mutate(
        pct_change_undisc = 100 * (yll_undisc_mean - ref_undisc) / ref_undisc,
        pct_change_disc   = 100 * (yll_disc_mean   - ref_disc)   / ref_disc,
        cells_match_ref   = n_cells == ref_cells) %>%
      dplyr::select(-ref_undisc, -ref_disc, -ref_cells)
  }
  
  list(by_method = by_method, by_country = by_country,
       by_method_fixed = .fixed_cell_onset(draws_by_shift, as.integer(shifts)),
       meta = list(shifts = as.integer(shifts), shift_end = shift_end,
                   n_sim = n_sim, seed = seed,
                   framing = framing_focus, catchup = catchup_focus))
}

# ---------------------------------------------------------------------------
# Composition-STABLE onset sensitivity (the primary supplementary table).
# ---------------------------------------------------------------------------
# For EACH method, restrict to the (disease, country) cells the method produced
# at EVERY shift it is present in, then sum the per-cell draw vectors within
# (method, shift) and compute the %-change vs shift 0 on that FIXED cell set.
# This isolates the onset effect from changes in the estimable cell set across
# shifts. Methods absent at some shifts (e.g. DiD at -2, where the pooled CS
# design hits the pre-window data boundary) are evaluated over the shifts where
# they ARE present; n_shifts_used records this. n_cells is constant within a
# method by construction, so cells_match_ref is always TRUE here.
.fixed_cell_onset <- function(draws_by_shift, shifts) {
  if (length(draws_by_shift) == 0) return(tibble::tibble())
  all_d <- dplyr::bind_rows(lapply(names(draws_by_shift), function(k) {
    x <- draws_by_shift[[k]]
    if (nrow(x) == 0) return(NULL)
    x %>% dplyr::mutate(start_shift = as.integer(k))
  }))
  if (is.null(all_d) || nrow(all_d) == 0) return(tibble::tibble())
  out <- list()
  for (m in unique(all_d$method)) {
    dm <- all_d %>% dplyr::filter(method == m)
    sh <- sort(unique(dm$start_shift))
    if (!(0L %in% sh)) next
    fixed <- dm %>% dplyr::distinct(disease, country, start_shift) %>%
      dplyr::count(disease, country, name = "n_sh") %>%
      dplyr::filter(n_sh == length(sh)) %>% dplyr::select(disease, country)
    if (nrow(fixed) == 0) next
    dmf <- dm %>% dplyr::semi_join(fixed, by = c("disease", "country"))
    per_shift <- dmf %>%
      dplyr::group_by(start_shift) %>%
      dplyr::group_modify(~ {
        M_un <- do.call(rbind, .x$d_undisc); M_di <- do.call(rbind, .x$d_disc)
        s_un <- colSums(M_un); s_di <- colSums(M_di)
        tibble::tibble(
          n_cells = nrow(M_un), n_shifts_used = length(sh),
          yll_undisc_mean = mean(s_un),
          yll_undisc_lo = stats::quantile(s_un, .025, names = FALSE),
          yll_undisc_hi = stats::quantile(s_un, .975, names = FALSE),
          yll_disc_mean = mean(s_di),
          yll_disc_lo = stats::quantile(s_di, .025, names = FALSE),
          yll_disc_hi = stats::quantile(s_di, .975, names = FALSE))
      }) %>% dplyr::ungroup() %>% dplyr::mutate(method = m)
    ref <- per_shift %>% dplyr::filter(start_shift == 0L)
    per_shift <- per_shift %>%
      dplyr::mutate(
        pct_change_undisc = 100 * (yll_undisc_mean - ref$yll_undisc_mean) / ref$yll_undisc_mean,
        pct_change_disc   = 100 * (yll_disc_mean   - ref$yll_disc_mean)   / ref$yll_disc_mean,
        cells_match_ref   = TRUE)
    out[[length(out) + 1]] <- per_shift
  }
  dplyr::bind_rows(out)
}

# ---------------------------------------------------------------------------
# Model-agnostic event-time coverage-gap profile (supplementary descriptive).
# ---------------------------------------------------------------------------
# Aligns each country's estimated coverage gap (gap_df; PRIMARY manual-onset
# specification) on event time tau = year - conflict_start, then averages the gap
# by (method, event_time) across countries. Shows the DYNAMIC post-onset coverage
# disruption the YLL rests on. gap_df contains only conflict-window years, so this
# is tau >= 0; the pre-onset trajectory (pre-trend test / pre-onset degradation)
# comes from the CS / HonestDiD event study in 07 and from the onset-shift table.
event_time_gap_profile <- function(gap_df, conflict_df = conflict_info,
                                   target_vaccine = NULL) {
  d <- gap_df
  if (!is.null(target_vaccine)) d <- d %>% dplyr::filter(vaccine == target_vaccine)
  d %>%
    dplyr::inner_join(conflict_df %>% dplyr::select(country, conflict_start),
                      by = "country") %>%
    dplyr::mutate(event_time = year - conflict_start) %>%
    dplyr::group_by(method, event_time) %>%
    dplyr::summarise(
      n_country = dplyr::n_distinct(country),
      mean_gap  = mean(gap, na.rm = TRUE),
      se_gap    = stats::sd(gap, na.rm = TRUE) / sqrt(pmax(dplyr::n(), 1L)),
      .groups   = "drop")
}

# ---------------------------------------------------------------------------
# Figure: onset-shift slopegraph (one line per method across -1 / 0 / +1).
# ---------------------------------------------------------------------------
make_fig_start_sensitivity <- function(sens, discounted = FALSE) {
  bm <- sens$by_method
  if (is.null(bm) || nrow(bm) == 0)
    stop("make_fig_start_sensitivity: empty by_method (no methods produced YLL).")
  ycol  <- if (discounted) "yll_disc_mean" else "yll_undisc_mean"
  lo_c  <- if (discounted) "yll_disc_lo"   else "yll_undisc_lo"
  hi_c  <- if (discounted) "yll_disc_hi"   else "yll_undisc_hi"
  meth_levels <- if (exists("METHOD_ORDER")) METHOD_ORDER else unique(bm$method)
  
  df <- bm %>%
    dplyr::mutate(
      method = factor(method, levels = intersect(meth_levels, unique(method))),
      y = .data[[ycol]], lo = .data[[lo_c]], hi = .data[[hi_c]])
  
  ggplot2::ggplot(df, ggplot2::aes(x = start_shift, y = y,
                                   colour = method, group = method)) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = lo, ymax = hi),
                            alpha = 0.35, linewidth = 0.5) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_x_continuous(breaks = sens$meta$shifts,
                                labels = sprintf("%+d", sens$meta$shifts),
                                name = "Conflict-onset shift (years)") +
    ggplot2::scale_y_continuous(labels = scales::comma_format(),
                                name = paste0("Global conflict-attributable YLL",
                                              if (discounted) " (discounted)" else " (undiscounted)")) +
    { if (exists("nm_palette"))
      ggplot2::scale_colour_manual(values = nm_palette, name = "Estimator")
      else ggplot2::scale_colour_discrete(name = "Estimator") } +
    ggplot2::labs(
      title = "Onset-timing sensitivity | YLL under \u00b11-year shifts in coded conflict start",
      subtitle = paste0("Differences are common-random-number stable (same seed/draws across shifts); ",
                        if (sens$meta$shift_end) "window length held fixed." else "conflict end held fixed."),
      caption = "Lines connect the same estimator across shifts; ranges are draw-level 95% CIs. Comparison-group methods should move less than no-control methods.") +
    { if (exists("theme_nm")) theme_nm() else ggplot2::theme_minimal() }
}