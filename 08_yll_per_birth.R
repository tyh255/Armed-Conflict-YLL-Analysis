# ============================================================================
# 08_YLL_PER_BIRTH.R   (supplementary; PRIMARY spec unchanged)
# ----------------------------------------------------------------------------
# Population-normalised burden. The absolute YLL total is dominated by Nigeria
# (~70-85% of the complete-case total) simply because Nigeria has by far the
# largest birth cohort -- that answers "where are the most life-years lost" but
# NOT "where is conflict most damaging to vaccination". This module reports YLL
# per 1,000 births over each country's conflict window, which is the INTENSITY of
# the conflict-attributable loss and is the figure that surfaces high-severity /
# smaller-population settings (e.g. Somalia, Syria, Yemen) that the absolute
# total hides. Both views belong in the paper (cf. the team's 2025 systematic
# review: effects concentrate in moderate-to-high-intensity civil wars, a per-
# capita phenomenon).
#
# Denominator: total births over conflict_start:conflict_end, summed from
# covariates_panel$Birth_Cohort (= Birth_Rate * Population_per_1000). Numerator:
# country YLL summed over diseases, per method, headline scenario.
#
# Source AFTER 04/05 (uses yll_country, covariates_panel, conflict_info,
# save_fig/save_tab, nm_palette/theme, FRAMING_HEADLINE).
# ============================================================================

# Total births over each country's conflict window (the per-birth denominator).
conflict_window_births <- function(covariates = covariates_panel,
                                   conflict_df = conflict_info) {
  stopifnot("Birth_Cohort" %in% names(covariates))
  conflict_df %>%
    dplyr::select(country, ISO3, conflict_start, conflict_end) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      births_total = sum(covariates$Birth_Cohort[
        covariates$ISO3 == ISO3 &
          covariates$year >= conflict_start &
          covariates$year <= conflict_end], na.rm = TRUE),
      n_years = conflict_end - conflict_start + 1L) %>%
    dplyr::ungroup()
}

# YLL per 1,000 births by (country, method) for the headline scenario.
yll_per_birth <- function(yll_country, covariates = covariates_panel,
                          conflict_df = conflict_info,
                          framing_focus = if (exists("FRAMING_HEADLINE")) FRAMING_HEADLINE else "campaign_topup",
                          catchup_focus = 0, discount = c("undisc", "disc")) {
  discount <- match.arg(discount)
  ycol <- if (discount == "undisc") "yll_undisc_total_mean" else "yll_disc_total_mean"
  stopifnot(ycol %in% names(yll_country))
  births <- conflict_window_births(covariates, conflict_df)
  
  ctry_method <- yll_country %>%
    dplyr::filter(framing == framing_focus, catchup == catchup_focus) %>%
    dplyr::group_by(country, method) %>%
    dplyr::summarise(yll_total = sum(.data[[ycol]], na.rm = TRUE),
                     .groups = "drop")
  
  ctry_method %>%
    dplyr::left_join(births %>% dplyr::select(country, births_total, n_years),
                     by = "country") %>%
    dplyr::mutate(
      yll_per_1000_births = ifelse(births_total > 0,
                                   1000 * yll_total / births_total, NA_real_),
      discount = discount)
}

# Figure: per-birth YLL ranked by the median across methods (Cleveland dot plot,
# one dot per method), so intensity ordering is visible independent of which
# estimator is used.
make_fig_per_birth <- function(pb_df) {
  if (is.null(pb_df) || nrow(pb_df) == 0)
    stop("make_fig_per_birth: empty per-birth table.")
  ord <- pb_df %>% dplyr::group_by(country) %>%
    dplyr::summarise(med = stats::median(yll_per_1000_births, na.rm = TRUE),
                     .groups = "drop") %>%
    dplyr::arrange(med)
  meth_levels <- if (exists("METHOD_ORDER")) METHOD_ORDER else unique(pb_df$method)
  df <- pb_df %>%
    dplyr::mutate(country = factor(country, levels = ord$country),
                  method  = factor(method, levels = intersect(meth_levels, unique(method))))
  ggplot2::ggplot(df, ggplot2::aes(x = yll_per_1000_births, y = country,
                                   colour = method)) +
    ggplot2::geom_point(size = 2.4, alpha = 0.85) +
    { if (exists("nm_palette"))
      ggplot2::scale_colour_manual(values = nm_palette, name = "Estimator")
      else ggplot2::scale_colour_discrete(name = "Estimator") } +
    ggplot2::scale_x_continuous(labels = scales::comma_format()) +
    ggplot2::labs(
      x = "Conflict-attributable YLL per 1,000 births (over conflict window)",
      y = NULL,
      title = "Intensity of conflict-attributable vaccination loss, population-normalised",
      subtitle = "Per-birth rate separates severity of disruption from population size; contrast with the Nigeria-dominated absolute total.",
      caption = "Headline scenario; one point per estimator. Countries ordered by cross-method median.") +
    { if (exists("theme_nm")) theme_nm() else ggplot2::theme_minimal() }
}

# Driver: build table + figure, with input self-tests.
run_yll_per_birth <- function(yll_country, covariates = covariates_panel,
                              conflict_df = conflict_info) {
  pb <- yll_per_birth(yll_country, covariates, conflict_df)
  # SELF-TEST: every conflict country should have a positive birth denominator;
  # a zero/NA denominator means the WB Birth_Cohort series is missing for that
  # ISO3-year and the rate would be undefined/explosive.
  bad <- pb %>% dplyr::filter(is.na(births_total) | births_total <= 0) %>%
    dplyr::distinct(country)
  if (nrow(bad) > 0)
    message("  [per-birth diag] missing/zero birth denominator for: ",
            paste(bad$country, collapse = ", "),
            " -> per-birth rate is NA there (check WB_Birth Rate / Population).")
  if (exists("save_tab")) save_tab(pb, "Table_S4_yll_per_1000_births.csv")
  fig <- tryCatch(make_fig_per_birth(pb),
                  error = function(e) { message("  per-birth fig skipped: ",
                                                conditionMessage(e)); NULL })
  if (!is.null(fig) && exists("save_fig"))
    save_fig(fig, "Figure_S_yll_per_birth", width = 9, height = 6)
  invisible(pb)
}