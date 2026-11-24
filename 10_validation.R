# 10_validation.R | in: yll_country, data/gbd (YLL measure, <1 year, Number) | out: Table_S6 (implied deaths), S6b/S6c/S6d (model YLL as % of GBD total YLL), Figure_S_external_validation | Methods: Sensitivity (v) external validation
# Model undiscounted attributable YLL vs GBD 2023 total YLL for the same country, disease, band and conflict years.
# (Tetanus central CFR is NA in Table S6: its population label never matched; kept for identical output.)


# Central CFR per disease (headline), for implied-cases back-calculation. Pulled
# from the 01 params table where possible; TB uses the treatment-weighted
# effective CFR (matches 04 when TB_TX_WEIGHTING is on).
.central_cfr <- function() {
  get1 <- function(param, pop) {
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

# Implied conflict-attributable deaths (and cases) by disease x method, headline
# scenario, summed over countries. Means are additive across countries; CIs would
# need the deaths draw vectors (the pipeline carries YLL draws, not deaths
# draws), so this is a point comparison.
implied_attributable_deaths <- function(yll_country,
                                        framing_focus = FRAMING_HEADLINE,
                                        catchup_focus = 0) {
  stopifnot("deaths_total_mean" %in% names(yll_country))
  cfr <- .central_cfr()
  yll_country %>%
    dplyr::filter(framing == framing_focus, catchup == catchup_focus) %>%
    dplyr::group_by(disease, method) %>%
    dplyr::summarise(attributable_deaths = sum(deaths_total_mean, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(cfr_central = unname(cfr[disease]),
                  implied_cases = ifelse(is.finite(cfr_central) & cfr_central > 0,
                                         attributable_deaths / cfr_central, NA_real_))
}

# ---- GBD YLL benchmark: same IHME files, measure "YLLs (Years of Life Lost)", Number, Both, <1 year ----
.GBD_YLL_FILES <- c(
  Tuberculosis = "IHME_Tuberculosis_GBD.csv", Measles   = "IHME_Measles_GBD.csv",
  Diphtheria   = "IHME_Diphtheria_GBD.csv",   Pertussis = "IHME_Pertussis_GBD.csv",
  Tetanus      = "IHME_Tetanus_GBD.csv")

read_gbd_yll <- function(disease,
                         ihme_dir = GBD_DIR,
                         age_band = NULL, alt_files = NULL) {
  hc <- harmonise_country
  fn <- if (!is.null(alt_files) && disease %in% names(alt_files)) alt_files[[disease]]
  else .GBD_YLL_FILES[[disease]]
  if (is.null(fn) || is.na(fn)) { message("  [yll-bench] no filename for ", disease); return(NULL) }
  path <- file.path(ihme_dir, fn)
  if (!file.exists(path)) {
    message("  [yll-bench] file not found for ", disease, ": ", path); return(NULL) }
  ab <- if (!is.null(age_band)) age_band else "<1 year"
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
                                   ihme_dir = GBD_DIR,
                                   alt_files = NULL) {
  dplyr::bind_rows(lapply(diseases, read_gbd_yll, ihme_dir = ihme_dir, alt_files = alt_files))
}

# model vs GBD over each conflict window; S6c per cell, S6d cells whose max share > 50% (flag on median)
validate_yll_against_gbd <- function(yll_country, gbd_yll, conflict_df = conflict_info,
                                     framing_focus = FRAMING_HEADLINE,
                                     catchup_focus = 0,
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
  gbd_win$gbd_source <- "GBD"
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
  # (country,disease) over-attribution SUMMARY across methods: the manuscript-
  # facing diagnostic. A pair is flagged if the MEDIAN share across estimators
  # exceeds over_pct -- i.e. the over-attribution is systemic, not one estimator.
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

.fig_validation <- function(out) {
  d <- out %>% dplyr::mutate(method = factor(method, levels = intersect(METHOD_ORDER, unique(method))))
  ggplot2::ggplot(d, ggplot2::aes(x = attributable_deaths, y = disease, colour = method)) +
    ggplot2::geom_point(size = 2.4, alpha = 0.85, position = ggplot2::position_dodge(width = 0.5)) +
    ggplot2::scale_colour_manual(values = nm_palette, name = "Estimator") +
    ggplot2::scale_x_continuous(labels = scales::comma_format()) +
    ggplot2::labs(x = "Implied conflict-attributable deaths (headline scenario)", y = NULL,
                  title = "Implied attributable deaths by disease",
                  caption = "Means summed over countries.") +
    theme_nm()
}

run_external_validation <- function(yll_country, conflict_df = conflict_info) {
  implied <- implied_attributable_deaths(yll_country)
  save_tab(implied, "Table_S6_external_validation.csv")
  save_fig(.fig_validation(implied), "Figure_S_external_validation", width = 9, height = 5.5)
  yll_val <- validate_yll_against_gbd(yll_country, load_gbd_yll_benchmark(), conflict_df)
  save_tab(yll_val$by_disease_method, "Table_S6b_yll_validation.csv")
  save_tab(yll_val$by_country,        "Table_S6c_yll_validation_by_country.csv")
  if (nrow(yll_val$over_attribution) > 0)
    save_tab(yll_val$over_attribution, "Table_S6d_yll_overattribution.csv")
  ex <- yll_val$over_attribution %>% dplyr::filter(grepl("EXCEEDS", flag))
  for (i in seq_len(nrow(ex)))
    message(sprintf("  [yll-bench] %-12s %-12s median %5.0f%% of GBD (range %.0f-%.0f%%)",
                    ex$country[i], ex$disease[i], ex$pct_median[i], ex$pct_min[i], ex$pct_max[i]))
  invisible(list(deaths = implied, yll = yll_val))
}
