
.central_cfr <- function() {
  get1 <- function(param, pop) {
    if (!exists("params")) return(NA_real_)
    r <- params[params$parameter == param &
                  (is.na(pop) | params$population == pop), ]
    if (nrow(r) == 0) NA_real_ else r$mean[1]   # params value column is 'mean'
  }
  tb_treated <- get1("CFR_TB", "Ages 0-4, treated")
  tb_untreat <- get1("CFR_TB", "Ages 0-4")
  p_treat    <- get1("TB_treat_coverage", "Frac paediatric TB treated (conflict)")
  tb_eff <- if (all(is.finite(c(tb_treated, tb_untreat, p_treat))))
    p_treat * tb_treated + (1 - p_treat) * tb_untreat else tb_untreat
  c(Tuberculosis = tb_eff,
    Measles      = get1("CFR_Measles", "Ages <1, community"),
    Pertussis    = get1("CFR_Pertussis", "Ages 0-1"),
    Diphtheria   = get1("CFR_Diphtheria", "Ages 0-4"),
    Tetanus      = get1("CFR_Tetanus", "Post-neonatal/childhood (ltd care)"))
}

implied_attributable_deaths <- function(yll_country,
                                        framing_focus = if (exists("FRAMING_HEADLINE")) FRAMING_HEADLINE else "campaign_topup",
                                        catchup_focus = 0) {
  stopifnot("deaths_total_mean" %in% names(yll_country))
  cfr <- .central_cfr()
  yll_country %>%
    dplyr::filter(framing == framing_focus, catchup == catchup_focus) %>%
    dplyr::group_by(disease, method) %>%
    dplyr::summarise(attributable_deaths = sum(deaths_total_mean, na.rm = TRUE),
                     .groups = "drop") %>%
    dplyr::mutate(
      cfr_central = unname(cfr[disease]),
      implied_cases = ifelse(is.finite(cfr_central) & cfr_central > 0,
                             attributable_deaths / cfr_central, NA_real_))
}


validate_against_benchmark <- function(implied_tab, yll_country,
                                       conflict_df = conflict_info,
                                       out_dir = if (exists("OUT_DIR")) OUT_DIR else ".") {
  bpath <- file.path(out_dir, "external_deaths_benchmark.csv")
  if (!file.exists(bpath)) {
    message("  [validation] no benchmark at ", bpath, ".")
    message("  -> Supply GBD/WHO cause-specific deaths (country,disease,year,deaths) ",
            "for the conflict windows to compute attributable-as-%-of-total. ",
            "Measles+pertussis are ~93% of the YLL total and WILL be queried by reviewers.")
    return(implied_tab)
  }
  bench <- utils::read.csv(bpath, stringsAsFactors = FALSE)
  needed <- c("country", "disease", "year", "deaths")
  if (!all(needed %in% names(bench)))
    stop("external_deaths_benchmark.csv must have columns: ",
         paste(needed, collapse = ", "))
  countries <- unique(yll_country$country)
  bench_win <- bench %>%
    dplyr::inner_join(conflict_df %>% dplyr::select(country, conflict_start, conflict_end),
                      by = "country") %>%
    dplyr::filter(country %in% countries,
                  year >= conflict_start, year <= conflict_end) %>%
    dplyr::group_by(disease) %>%
    dplyr::summarise(benchmark_total_deaths = sum(deaths, na.rm = TRUE),
                     .groups = "drop")
  implied_tab %>%
    dplyr::left_join(bench_win, by = "disease") %>%
    dplyr::mutate(
      attributable_pct_of_total = 100 * attributable_deaths / benchmark_total_deaths,
      flag = dplyr::case_when(
        is.na(benchmark_total_deaths)              ~ "no benchmark",
        attributable_pct_of_total > 100            ~ "IMPLAUSIBLE (>100% of total)",
        attributable_pct_of_total < 0.5            ~ "very small share (recheck)",
        TRUE                                       ~ "plausible"))
}

# ============================================================================
# GBD YLL BENCHMARK (preferred over the deaths benchmark)
# ----------------------------------------------------------------------------

.GBD_YLL_FILES <- c(
  Tuberculosis = "IHME_Tuberculosis_GBD.csv", Measles   = "IHME_Measles_GBD.csv",
  Diphtheria   = "IHME_Diphtheria_GBD.csv",   Pertussis = "IHME_Pertussis_GBD.csv",
  Tetanus      = "IHME_Tetanus_GBD.csv")

# Load GBD total YLL (Number, Both sexes, model's age band) for one disease.
# `alt_files` lets you point at separate YLL downloads if your files aren't the
# multi-measure incidence files (e.g. c(Measles = "IHME_Measles_YLL.csv")).
read_gbd_yll <- function(disease,
                         ihme_dir = if (exists("IHME_DIR")) IHME_DIR else ".",
                         age_band = NULL, alt_files = NULL) {
  hc <- if (exists("harmonise_country")) harmonise_country else function(x) x
  fn <- if (!is.null(alt_files) && disease %in% names(alt_files)) alt_files[[disease]]
  else .GBD_YLL_FILES[[disease]]
  if (is.null(fn) || is.na(fn)) { message("  [yll-bench] no filename for ", disease); return(NULL) }
  path <- file.path(ihme_dir, fn)
  if (!file.exists(path)) {
    message("  [yll-bench] file not found for ", disease, ": ", path); return(NULL) }
  ab <- if (!is.null(age_band)) age_band else
    if (exists("disease_age_band")) disease_age_band(disease) else "<1 year"
  df <- utils::read.csv(path, stringsAsFactors = FALSE)
  need <- c("measure_name","metric_name","sex_name","age_name",
            "location_name","year","val","upper","lower")
  if (!all(need %in% names(df))) {
    message("  [yll-bench] ", disease, ": missing columns (",
            paste(setdiff(need, names(df)), collapse = ", "), "); skipping."); return(NULL) }
  out <- df %>%
    dplyr::filter(measure_name == "YLLs (Years of Life Lost)",
                  metric_name  == "Number",
                  sex_name     == "Both",
                  age_name     == ab) %>%
    dplyr::transmute(country = hc(location_name), disease = .env$disease,
                     year = as.integer(year),
                     gbd_yll = val, gbd_yll_lo = lower, gbd_yll_hi = upper) %>%
    dplyr::group_by(country, disease, year) %>%      # guard duplicate rows
    dplyr::summarise(dplyr::across(c(gbd_yll, gbd_yll_lo, gbd_yll_hi),
                                   ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  if (nrow(out) == 0)
    message("  [yll-bench] ", disease, ": 0 YLL rows after filter ",
            "(check measure='YLLs (Years of Life Lost)', metric='Number', age_name='", ab, "').")
  out
}

load_gbd_yll_benchmark <- function(diseases = names(.GBD_YLL_FILES),
                                   ihme_dir = if (exists("IHME_DIR")) IHME_DIR else ".",
                                   alt_files = NULL) {
  dplyr::bind_rows(lapply(diseases, read_gbd_yll, ihme_dir = ihme_dir, alt_files = alt_files))
}


.apply_gbd_override <- function(gbd_win, gbd_override) {
  if (is.null(gbd_override) || nrow(gbd_override) == 0) return(gbd_win)
  need <- c("country", "disease", "gbd_yll_window")
  if (!all(need %in% names(gbd_override)))
    stop("gbd override needs columns: ", paste(need, collapse = ", "))
  ov <- gbd_override %>%
    dplyr::transmute(country, disease, gbd_yll_window_ov = as.numeric(gbd_yll_window),
                     gbd_source = "override")
  out <- gbd_win %>%
    dplyr::left_join(ov, by = c("country", "disease")) %>%
    dplyr::mutate(
      gbd_yll_window = ifelse(is.finite(gbd_yll_window_ov), gbd_yll_window_ov, gbd_yll_window),
      gbd_source     = ifelse(is.na(gbd_source), "GBD", gbd_source)) %>%
    dplyr::select(-gbd_yll_window_ov)
  n_ov <- sum(out$gbd_source == "override", na.rm = TRUE)
  if (n_ov > 0) message("  [yll-bench] applied ", n_ov,
                        " (country,disease) GBD-denominator override(s).")
  out
}


validate_yll_against_gbd <- function(yll_country, gbd_yll, conflict_df = conflict_info,
                                     framing_focus = if (exists("FRAMING_HEADLINE")) FRAMING_HEADLINE else "campaign_topup",
                                     catchup_focus = 0,
                                     gbd_override = NULL,
                                     over_pct = 100, high_pct = 50) {
  if (is.null(gbd_yll) || nrow(gbd_yll) == 0) return(NULL)
  if (!"yll_undisc_total_mean" %in% names(yll_country)) {
    message("  [yll-bench] yll_country lacks yll_undisc_total_mean; skipping."); return(NULL) }
  gbd_win <- gbd_yll %>%
    dplyr::inner_join(conflict_df %>%
                        dplyr::select(country, conflict_start, conflict_end), by = "country") %>%
    dplyr::filter(year >= conflict_start, year <= conflict_end) %>%
    dplyr::group_by(country, disease) %>%
    dplyr::summarise(gbd_yll_window = sum(gbd_yll, na.rm = TRUE), .groups = "drop")
  gbd_win <- .apply_gbd_override(gbd_win, gbd_override)
  if (!"gbd_source" %in% names(gbd_win)) gbd_win$gbd_source <- "GBD"
  model <- yll_country %>%
    dplyr::filter(framing == framing_focus, catchup == catchup_focus) %>%
    dplyr::select(country, disease, method, model_yll = yll_undisc_total_mean)
  joined <- model %>% dplyr::inner_join(gbd_win, by = c("country", "disease"))
  if (nrow(joined) == 0) {
    message("  [yll-bench] no (country,disease) overlap between model and GBD ",
            "(check harmonise_country() spellings and disease labels)."); return(NULL) }
  .cell_flag <- function(pct) dplyr::case_when(
    !is.finite(pct)   ~ "no GBD denominator",
    pct > over_pct    ~ "EXCEEDS GBD (>100%)",
    pct > high_pct    ~ "high share (>50%)",
    pct < 0.5         ~ "very small share (recheck)",
    TRUE              ~ "plausible")
  by_dm <- joined %>%
    dplyr::group_by(disease, method) %>%
    dplyr::summarise(n_countries      = dplyr::n_distinct(country),
                     model_yll_attrib = sum(model_yll, na.rm = TRUE),
                     gbd_yll_total    = sum(gbd_yll_window, na.rm = TRUE),
                     .groups = "drop") %>%
    dplyr::mutate(attributable_pct_of_gbd = 100 * model_yll_attrib / gbd_yll_total,
                  flag = dplyr::case_when(
                    !is.finite(attributable_pct_of_gbd) ~ "no GBD denominator",
                    attributable_pct_of_gbd > over_pct  ~ "IMPLAUSIBLE (>100% of GBD YLL)",
                    attributable_pct_of_gbd < 0.5       ~ "very small share (recheck)",
                    TRUE                                ~ "plausible"))
  by_country <- joined %>%
    dplyr::mutate(attributable_pct_of_gbd = 100 * model_yll / gbd_yll_window,
                  flag = .cell_flag(attributable_pct_of_gbd)) %>%
    dplyr::arrange(disease, method, dplyr::desc(attributable_pct_of_gbd))
  
  over_attr <- joined %>%
    dplyr::mutate(pct = 100 * model_yll / gbd_yll_window) %>%
    dplyr::group_by(country, disease) %>%
    dplyr::summarise(
      n_methods       = dplyr::n_distinct(method),
      gbd_yll_window  = dplyr::first(gbd_yll_window),
      gbd_source      = dplyr::first(gbd_source),
      pct_min         = min(pct, na.rm = TRUE),
      pct_median      = stats::median(pct, na.rm = TRUE),
      pct_max         = max(pct, na.rm = TRUE),
      .groups = "drop") %>%
    dplyr::mutate(flag = .cell_flag(pct_median)) %>%
    dplyr::filter(pct_max > high_pct) %>%
    dplyr::arrange(dplyr::desc(pct_median))
  list(by_disease_method = by_dm, by_country = by_country,
       over_attribution = over_attr)
}

run_external_validation <- function(yll_country, conflict_df = conflict_info,
                                    gbd_alt_files = NULL,
                                    gbd_yll_override_csv =
                                      if (exists("OUT_DIR"))
                                        file.path(OUT_DIR, "gbd_yll_window_override.csv")
                                    else "gbd_yll_window_override.csv") {
  implied <- implied_attributable_deaths(yll_country)
  out <- validate_against_benchmark(implied, yll_country, conflict_df)
  # SELF-TEST: attributable deaths must be non-negative and finite.
  if (any(implied$attributable_deaths < 0, na.rm = TRUE))
    message("  [validation diag] negative attributable deaths present (check gap signs).")
  if ("attributable_pct_of_total" %in% names(out)) {
    bad <- out %>% dplyr::filter(grepl("IMPLAUSIBLE", flag))
    if (nrow(bad) > 0)
      message("  [validation] IMPLAUSIBLE attributable share for: ",
              paste(unique(bad$disease), collapse = ", "),
              " -> attributable deaths exceed external total; revisit CFR/incidence.")
  }
  if (exists("save_tab")) save_tab(out, "Table_S6_external_validation.csv")
  fig <- tryCatch(.fig_validation(out),
                  error = function(e) { message("  validation fig skipped: ",
                                                conditionMessage(e)); NULL })
  if (!is.null(fig) && exists("save_fig"))
    save_fig(fig, "Figure_S_external_validation", width = 9, height = 5.5)
  
  # ---- GBD YLL benchmark (preferred): attributable YLL as % of GBD total YLL.
  gbd <- tryCatch(load_gbd_yll_benchmark(alt_files = gbd_alt_files),
                  error = function(e) { message("  [yll-bench] load failed: ",
                                                conditionMessage(e)); NULL })
  # Optional per-(country,disease) denominator override for collapsed-surveillance
  # settings (Syria/Yemen). Absent the file, behaviour is unchanged.
  gbd_override <- NULL
  if (!is.null(gbd_yll_override_csv) && file.exists(gbd_yll_override_csv)) {
    gbd_override <- tryCatch(utils::read.csv(gbd_yll_override_csv, stringsAsFactors = FALSE),
                             error = function(e) { message("  [yll-bench] override read failed: ",
                                                           conditionMessage(e)); NULL })
    if (!is.null(gbd_override))
      message("  [yll-bench] using GBD-denominator override: ", gbd_yll_override_csv)
  }
  yll_val <- validate_yll_against_gbd(yll_country, gbd, conflict_df,
                                      gbd_override = gbd_override)
  if (!is.null(yll_val)) {
    if (exists("save_tab")) {
      save_tab(yll_val$by_disease_method, "Table_S6b_yll_validation.csv")
      save_tab(yll_val$by_country,        "Table_S6c_yll_validation_by_country.csv")
      if (!is.null(yll_val$over_attribution) && nrow(yll_val$over_attribution) > 0)
        save_tab(yll_val$over_attribution, "Table_S6d_yll_overattribution.csv")
    }
    # COUNTRY-LEVEL over-attribution diagnostic (the pooled S6b dilutes it). This
    # is the Syria/Yemen pertussis finding -- surface it loudly and by name.
    oa <- yll_val$over_attribution
    if (!is.null(oa)) {
      exceed <- oa %>% dplyr::filter(grepl("EXCEEDS", flag))
      if (nrow(exceed) > 0) {
        message("  [yll-bench] COUNTRY-LEVEL over-attribution (median model YLL > 100% of GBD YLL):")
        for (i in seq_len(nrow(exceed)))
          message(sprintf("      %-12s %-12s median %5.0f%% (range %.0f-%.0f%%, denom=%s)",
                          exceed$country[i], exceed$disease[i], exceed$pct_median[i],
                          exceed$pct_min[i], exceed$pct_max[i], exceed$gbd_source[i]))
        message("    -> in collapsed-surveillance settings GBD under-ascertains the ",
                "denominator; supply gbd_yll_window_override.csv (country,disease,",
                "gbd_yll_window) to benchmark these against an independent source.")
      }
    }
    bad <- yll_val$by_disease_method %>% dplyr::filter(grepl("IMPLAUSIBLE", flag))
    if (nrow(bad) > 0)
      message("  [yll-bench] attributable YLL EXCEEDS GBD total YLL (disease-pooled) for: ",
              paste(unique(paste(bad$disease, bad$method, sep = "/")), collapse = "; "),
              " -> over-attribution or model CFR > GBD implied CFR; inspect.")
    else
      message("  [yll-bench] OK (pooled): every disease x method attributable share < 100% of GBD YLL.")
    # Headline-anchor shares for the log (measles & pertussis dominate the total).
    anchor <- yll_val$by_disease_method %>%
      dplyr::filter(method == "Synthetic Control",
                    disease %in% c("Measles", "Pertussis"))
    if (nrow(anchor) > 0)
      message("  [yll-bench] SC attributable share of GBD YLL: ",
              paste(sprintf("%s %.0f%%", anchor$disease, anchor$attributable_pct_of_gbd),
                    collapse = ", "))
    return(invisible(list(deaths = out, yll = yll_val)))
  }
  message("  [yll-bench] no GBD YLL benchmark produced (IHME files missing/empty; ",
          "pass gbd_alt_files=c(Measles='...') if your YLL data are in separate files).")
  invisible(out)
}

.fig_validation <- function(out) {
  meth_levels <- if (exists("METHOD_ORDER")) METHOD_ORDER else unique(out$method)
  d <- out %>% dplyr::mutate(
    method = factor(method, levels = intersect(meth_levels, unique(method))))
  ggplot2::ggplot(d, ggplot2::aes(x = attributable_deaths, y = disease, colour = method)) +
    ggplot2::geom_point(size = 2.4, alpha = 0.85,
                        position = ggplot2::position_dodge(width = 0.5)) +
    { if (exists("nm_palette"))
      ggplot2::scale_colour_manual(values = nm_palette, name = "Estimator")
      else ggplot2::scale_colour_discrete(name = "Estimator") } +
    ggplot2::scale_x_continuous(labels = scales::comma_format()) +
    ggplot2::labs(
      x = "Implied conflict-attributable deaths (headline scenario)",
      y = NULL,
      title = "Implied attributable deaths by disease, for external benchmarking",
      subtitle = "Compare against GBD/WHO cause-specific deaths; measles & pertussis dominate and warrant explicit validation.",
      caption = "Means summed over countries; supply external_deaths_benchmark.csv for attributable-as-%-of-total.") +
    { if (exists("theme_nm")) theme_nm() else ggplot2::theme_minimal() }
}
