# 06_onset.R | in: coverage_long, covariates_panel, conflict_info, gap_df (shift 0) | out: Table_S2/S2b/S2c/S2d, Figure_S_onset_sensitivity | Methods: Sensitivity (i) onset timing
# Re-runs gap -> YLL with conflict_start shifted (end fixed). Same seeds at every shift, so differences
# are the onset effect, not MC noise. Gaps for all shifts are estimated in parallel (forked) first.

shift_conflict_info <- function(conflict_df, delta, baseline_window = 3L) {
  out <- conflict_df
  out$conflict_start <- out$conflict_start + as.integer(delta)
  out$preconflict_year   <- out$conflict_start - 1L
  out$preconflict_window <- purrr::map(out$preconflict_year, ~ (.x - (baseline_window - 1L)):.x)
  bad <- out$conflict_start >= out$conflict_end          # need >= 2-year window
  if (any(bad)) message("  onset ", sprintf("%+d", delta), ": dropping ", paste(out$country[bad], collapse = ", "))
  out[!bad, , drop = FALSE]
}

run_start_sensitivity <- function(coverage_df, covariates, conflict_df = conflict_info, gap_df0 = NULL,
                                  shifts = c(-2L, -1L, 0L, 1L, 2L), methods = ALL_METHODS,
                                  n_sim = 2000L, seed = SEED_MC, framing_focus = FRAMING_HEADLINE,
                                  catchup_focus = 0, n_cores = N_CORES) {
  stopifnot(0L %in% shifts)
  cinfo <- lapply(shifts, function(d) shift_conflict_info(conflict_df, d)); names(cinfo) <- shifts
  # shift 0 = the main gap_df (identical inputs, deterministic estimators); other shifts:
  # every (shift, estimator) job forked in parallel, then stacked per shift in the fixed order
  todo <- setdiff(as.character(shifts), if (!is.null(gap_df0)) "0")
  jobs <- list(); key <- character(0)
  for (k in todo) {
    jk <- .gap_jobs(coverage_df, cinfo[[k]], methods = methods)
    jobs <- c(jobs, unname(jk)); key <- c(key, paste(k, names(jk), sep = "|"))
  }
  res <- .run_jobs(stats::setNames(jobs, key), n_cores)
  gaps <- list()
  for (k in todo) {
    pk <- res[startsWith(names(res), paste0(k, "|"))]
    names(pk) <- sub("^[^|]*\\|", "", names(pk))
    gaps[[k]] <- .combine_gaps(pk)
  }
  if (!is.null(gap_df0)) gaps[["0"]] <- gap_df0

  by_method <- list(); by_country <- list(); draws_by_shift <- list()
  for (d in as.integer(shifts)) {
    message("  onset shift ", sprintf("%+d", d), ": YLL MC")
    yll_d <- run_full_yll(gaps[[as.character(d)]], covariates, conflict_df = cinfo[[as.character(d)]],
                          n_sim = n_sim, seed = seed)
    dd <- yll_d$draws %>% dplyr::filter(framing == framing_focus, catchup == catchup_focus)
    draws_by_shift[[as.character(d)]] <- .nonempty(dd)
    gm <- aggregate_global_by_method_draws(yll_d$draws, framing_focus = framing_focus, catchup_focus = catchup_focus)
    if (nrow(gm) > 0) by_method[[length(by_method) + 1]] <- gm %>% dplyr::mutate(start_shift = d)
    cc <- aggregate_country_yll(yll_d$summary) %>%
      dplyr::filter(framing == framing_focus, catchup == catchup_focus) %>%
      dplyr::group_by(country, method) %>%
      dplyr::summarise(yll_undisc_total = sum(yll_undisc_total_mean, na.rm = TRUE),
                       yll_disc_total = sum(yll_disc_total_mean, na.rm = TRUE), .groups = "drop") %>%
      dplyr::mutate(start_shift = d)
    by_country[[length(by_country) + 1]] <- cc
  }
  by_method <- dplyr::bind_rows(by_method); by_country <- dplyr::bind_rows(by_country)
  ref <- by_method %>% dplyr::filter(start_shift == 0L) %>%
    dplyr::select(method, ref_undisc = yll_undisc_mean, ref_disc = yll_disc_mean, ref_cells = n_cells)
  by_method <- by_method %>% dplyr::left_join(ref, by = "method") %>%
    dplyr::mutate(pct_change_undisc = 100 * (yll_undisc_mean - ref_undisc) / ref_undisc,
                  pct_change_disc   = 100 * (yll_disc_mean - ref_disc) / ref_disc,
                  cells_match_ref   = n_cells == ref_cells) %>%
    dplyr::select(-ref_undisc, -ref_disc, -ref_cells)
  list(by_method = by_method, by_country = by_country,
       by_method_fixed = .fixed_cell_onset(draws_by_shift),
       meta = list(shifts = as.integer(shifts), n_sim = n_sim, seed = seed))
}

# composition-stable version (S2c): per method, cells present at every shift the method appears in
.fixed_cell_onset <- function(draws_by_shift) {
  all_d <- dplyr::bind_rows(lapply(names(draws_by_shift), function(k) {
    x <- draws_by_shift[[k]]; if (nrow(x) == 0) NULL else x %>% dplyr::mutate(start_shift = as.integer(k)) }))
  out <- list()
  for (m in unique(all_d$method)) {
    dm <- all_d %>% dplyr::filter(method == m)
    sh <- sort(unique(dm$start_shift))
    if (!(0L %in% sh)) next
    fixed <- dm %>% dplyr::distinct(disease, country, start_shift) %>%
      dplyr::count(disease, country, name = "n_sh") %>% dplyr::filter(n_sh == length(sh)) %>%
      dplyr::select(disease, country)
    if (nrow(fixed) == 0) next
    per_shift <- dm %>% dplyr::semi_join(fixed, by = c("disease", "country")) %>%
      dplyr::group_by(start_shift) %>%
      dplyr::group_modify(~ {
        s_un <- colSums(do.call(rbind, .x$d_undisc)); s_di <- colSums(do.call(rbind, .x$d_disc))
        tibble::tibble(n_cells = nrow(.x), n_shifts_used = length(sh),
                       yll_undisc_mean = mean(s_un), yll_undisc_lo = stats::quantile(s_un, .025, names = FALSE),
                       yll_undisc_hi = stats::quantile(s_un, .975, names = FALSE),
                       yll_disc_mean = mean(s_di), yll_disc_lo = stats::quantile(s_di, .025, names = FALSE),
                       yll_disc_hi = stats::quantile(s_di, .975, names = FALSE))
      }) %>% dplyr::ungroup() %>% dplyr::mutate(method = m)
    ref <- per_shift %>% dplyr::filter(start_shift == 0L)
    out[[length(out) + 1]] <- per_shift %>%
      dplyr::mutate(pct_change_undisc = 100 * (yll_undisc_mean - ref$yll_undisc_mean) / ref$yll_undisc_mean,
                    pct_change_disc   = 100 * (yll_disc_mean - ref$yll_disc_mean) / ref$yll_disc_mean,
                    cells_match_ref   = TRUE)
  }
  dplyr::bind_rows(out)
}

# mean gap by (method, event time) on the headline gap_df (S2d)
event_time_gap_profile <- function(gap_df, conflict_df = conflict_info) {
  gap_df %>%
    dplyr::inner_join(conflict_df %>% dplyr::select(country, conflict_start), by = "country") %>%
    dplyr::mutate(event_time = year - conflict_start) %>%
    dplyr::group_by(method, event_time) %>%
    dplyr::summarise(n_country = dplyr::n_distinct(country), mean_gap = mean(gap, na.rm = TRUE),
                     se_gap = stats::sd(gap, na.rm = TRUE) / sqrt(pmax(dplyr::n(), 1L)), .groups = "drop")
}

make_fig_start_sensitivity <- function(sens) {
  df <- sens$by_method %>%
    dplyr::mutate(method = factor(method, levels = intersect(METHOD_ORDER, unique(method))),
                  y = yll_undisc_mean, lo = yll_undisc_lo, hi = yll_undisc_hi)
  ggplot2::ggplot(df, ggplot2::aes(x = start_shift, y = y, colour = method, group = method)) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.35, linewidth = 0.5) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_x_continuous(breaks = sens$meta$shifts,
                                labels = ifelse(sens$meta$shifts == 0, "0", sprintf("%+d", sens$meta$shifts)),
                                name = "Conflict-onset shift (years)") +
    ggplot2::scale_y_continuous(labels = scales::comma_format(),
                                name = "Global conflict-attributable YLL (undiscounted)") +
    ggplot2::scale_colour_manual(values = nm_palette, name = "Estimator") +
    ggplot2::labs(title = sprintf("Onset-timing sensitivity | YLL under shifts of up to ±%d years in coded conflict start",
                                  max(abs(sens$meta$shifts))),
                  subtitle = "Common random numbers across shifts; conflict end held fixed.",
                  caption = "Lines connect the same estimator across shifts; ranges are draw-level 95% CIs.") +
    theme_nm()
}
