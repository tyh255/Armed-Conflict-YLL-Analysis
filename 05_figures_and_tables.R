# ============================================================================
# 05_FIGURES_AND_TABLES.R
# ----------------------------------------------------------------------------
# Publication figures and tables in Nature Medicine style:
#   - Clean classic theme, no chartjunk
#   - Distinct, accessible colour palette
#   - Bold panel titles, clear subtitle/caption
#   - Figure 1: coverage trajectories per country with counterfactuals
#   - Figure 2: forest plot of total YLLs per country-disease faceted by method
#   - Figure 3: catch-up sensitivity (framing x method grid)
#   - Table 1:  country x disease YLL totals at lambda=0
#   - Table 2:  method comparison globals
#   - Table S1: parameter source table
# ============================================================================

# ---- Palettes -------------------------------------------------------------
nm_palette <- c(
  "Baseline"          = "#3B5BA5",   # deep blue
  "ITS"               = "#618B4A",   # muted green
  "Synthetic Control" = "#C44536",   # brick red
  "SC + covariates"   = "#E08A2B",   # amber
  "Augmented SC"      = "#7A5295",   # violet
  "DiD (CS)"          = "#2C7A7B",   # teal
  # baseline-window variants (now first-class methods: 1/3/5-yr)
  "Baseline (1yr)"    = "#27407A",
  "Baseline (3yr)"    = "#3B5BA5",
  "Baseline (5yr)"    = "#8AA0C8"
)
# Method display order (no-control -> control-based; flat -> trend -> match ->
# differencing). Used to order Figure-4 rows and legends.
METHOD_ORDER <- c("Baseline (1yr)", "Baseline (3yr)", "Baseline (5yr)",
                  "ITS", "Synthetic Control",
                  "SC + covariates", "Augmented SC", "DiD (CS)")
nm_disease_palette <- c(
  "Tuberculosis" = "#1F4E79",
  "Measles"      = "#B5651D",
  "Diphtheria"   = "#7A5295",
  "Pertussis"    = "#2E8B57",
  "Tetanus"      = "#8B0000"
)

# --- single source of truth for scenario framings ---------------------------
# Must match CATCHUP_FRAMINGS in 04_yll_monte_carlo.R (line 105).
FRAMING_LEVELS   <- c("routine_recovery", "campaign_topup")
FRAMING_HEADLINE <- "campaign_topup"   # <- set to your chosen headline scenario

FRAMING_LABELS <- c(                   # display names for facets/legends
  routine_recovery = "Routine recovery",
  campaign_topup   = "Campaign top-up"
)

# ---- Theme ----------------------------------------------------------------
theme_nm <- function(base_size = 9) {
  ggplot2::theme_classic(base_size = base_size, base_family = "") +
    ggplot2::theme(
      panel.background  = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background   = ggplot2::element_rect(fill = "white", colour = NA),
      panel.grid.major  = ggplot2::element_line(colour = "grey92", linewidth = 0.3),
      panel.grid.minor  = ggplot2::element_blank(),
      axis.line         = ggplot2::element_line(colour = "black", linewidth = 0.4),
      axis.ticks        = ggplot2::element_line(colour = "black", linewidth = 0.3),
      axis.text         = ggplot2::element_text(colour = "black", size = base_size - 1),
      axis.title        = ggplot2::element_text(colour = "black", size = base_size),
      strip.background  = ggplot2::element_blank(),
      strip.text        = ggplot2::element_text(face = "bold", size = base_size, hjust = 0),
      legend.background = ggplot2::element_rect(fill = "white", colour = NA),
      legend.key        = ggplot2::element_rect(fill = "white", colour = NA),
      legend.title      = ggplot2::element_text(size = base_size - 1, face = "bold"),
      legend.text       = ggplot2::element_text(size = base_size - 1),
      legend.position   = "top",
      plot.title        = ggplot2::element_text(size = base_size + 2, face = "bold", hjust = 0,
                                                margin = ggplot2::margin(b = 4)),
      plot.subtitle     = ggplot2::element_text(size = base_size, hjust = 0,
                                                margin = ggplot2::margin(b = 8)),
      plot.caption      = ggplot2::element_text(size = base_size - 2, hjust = 0,
                                                colour = "grey30",
                                                margin = ggplot2::margin(t = 6)),
      plot.margin       = ggplot2::margin(10, 12, 8, 10)
    )
}

# ---- Formatters -----------------------------------------------------------
fmt_ci <- function(mean_val, lo, hi, digits = 0, big = TRUE) {
  fmt <- function(x) formatC(round(x, digits), format = "f", digits = digits,
                             big.mark = ifelse(big, ",", ""))
  result <- paste0(fmt(mean_val), " (", fmt(lo), "-", fmt(hi), ")")
  result[is.na(mean_val) | is.na(lo) | is.na(hi)] <- "--"
  result
}

# ============================================================================
# fig1_coverage_decline.R  -  REFACTORED FIGURE 1 (drop-in)
# ----------------------------------------------------------------------------
# Replaces make_fig1() in 05_figures_and_tables.R. Source this AFTER 05 (it
# reuses theme_nm() and save_fig()), or paste the bodies over the old make_fig1.
#
# Design goals (Nature, visually striking, decline- and contrast-forward):
#   * ONE figure covering all three antigens (BCG, DTP3, MCV1) instead of three
#     near-identical single-antigen panels.
#   * Panel a - coverage trajectories per country, antigens as bold colour lines;
#     the COVERAGE LOSS (area between each antigen's pre-conflict level and its
#     observed curve, from conflict onset onward) is shaded in the antigen colour
#     so the decline reads instantly. Conflict window shaded; onset marked.
#     Countries ordered worst -> least so the layout itself is a severity gradient.
#   * Panel b - country x antigen heatmap of the percentage-point coverage drop,
#     warm sequential ramp, values printed: the cross-country / cross-antigen
#     contrast at a glance, in the SAME country order as panel a.
#
# Descriptive by design: Figure 1 tells the raw coverage-collapse story;
# counterfactual estimators live in Figures 2/4. Decline is defined model-free
# as (pre-conflict reference level) - (within-conflict trough), in points.
#
# Update 00_run_all.R: replace the three make_fig1()/save_fig() calls with one
# (gap_df is accepted but unused, so the old positional order still works):
#   fig1 <- make_fig1(coverage_long, gap_df, conflict_info)
#   save_fig(fig1, "Figure_1_coverage_decline", width = 11, height = 9)
# ============================================================================

stopifnot(requireNamespace("ggplot2", quietly = TRUE))

# Vivid, colour-blind-safe antigen palette (ColorBrewer Dark2 subset).
VACCINE_LEVELS <- c("BCG", "DTP3", "MCV1")
VACCINE_COLS   <- c(BCG = "#1B9E77", DTP3 = "#7570B3", MCV1 = "#D95F02")
VACCINE_LABS   <- c(BCG = "BCG (TB)", DTP3 = "DTP3 (diphtheria/pertussis/tetanus)",
                    MCV1 = "MCV1 (measles)")

# ---- model-free decline metric: pre-conflict level -> within-conflict trough --
.fig1_declines <- function(coverage_long, conflict_df,
                           vaccines = VACCINE_LEVELS, ref_years = 3L) {
  coverage_long %>%
    dplyr::filter(vaccine %in% vaccines, country %in% conflict_df$country,
                  !is.na(coverage)) %>%
    dplyr::inner_join(dplyr::select(conflict_df, country, conflict_start, conflict_end),
                      by = "country") %>%
    dplyr::group_by(country, vaccine) %>%
    dplyr::summarise(
      pre_level = {
        pre <- coverage[year >= conflict_start[1] - ref_years & year < conflict_start[1]]
        if (length(pre) == 0) coverage[which.min(abs(year - conflict_start[1]))]
        else mean(pre, na.rm = TRUE)
      },
      trough = {
        win <- coverage[year >= conflict_start[1] & year <= conflict_end[1]]
        if (length(win) == 0) NA_real_ else min(win, na.rm = TRUE)
      },
      .groups = "drop") %>%
    dplyr::mutate(decline_pp = pmax(pre_level - trough, 0))
}

# ---- panel a: trajectories with shaded coverage loss -----------------------
.fig1_panel_a <- function(coverage_long, conflict_df, decl, country_order,
                          year_min = 1995L, year_max = 2024L) {
  cov <- coverage_long %>%
    dplyr::filter(vaccine %in% VACCINE_LEVELS, country %in% country_order,
                  year >= year_min, year <= year_max) %>%
    dplyr::mutate(country = factor(country, levels = country_order),
                  vaccine = factor(vaccine, levels = VACCINE_LEVELS))
  shade <- conflict_df %>%
    dplyr::filter(country %in% country_order) %>%
    dplyr::transmute(country = factor(country, levels = country_order),
                     xmin = conflict_start, xmax = conflict_end)
  # loss ribbon: from onset onward, between pmin(obs, pre_level) and pre_level
  ribbon <- cov %>%
    dplyr::left_join(dplyr::select(decl, country, vaccine, pre_level),
                     by = c("country", "vaccine")) %>%
    dplyr::left_join(dplyr::select(conflict_df, country, conflict_start),
                     by = "country") %>%
    dplyr::filter(year >= conflict_start) %>%
    dplyr::mutate(ylo = pmin(coverage, pre_level), yhi = pre_level)
  
  ggplot2::ggplot() +
    ggplot2::geom_rect(data = shade,
                       ggplot2::aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 100),
                       fill = "grey85", alpha = 0.45, inherit.aes = FALSE) +
    ggplot2::geom_segment(data = shade,
                          ggplot2::aes(x = xmin, xend = xmin, y = 0, yend = 100),
                          colour = "grey55", linewidth = 0.3, linetype = "22",
                          inherit.aes = FALSE) +
    ggplot2::geom_ribbon(data = ribbon,
                         ggplot2::aes(x = year, ymin = ylo, ymax = yhi,
                                      fill = vaccine, group = vaccine),
                         alpha = 0.20) +
    ggplot2::geom_line(data = cov,
                       ggplot2::aes(x = year, y = coverage, colour = vaccine),
                       linewidth = 0.75) +
    ggplot2::facet_wrap(~ country, ncol = 5) +
    ggplot2::scale_colour_manual(values = VACCINE_COLS, labels = VACCINE_LABS,
                                 name = NULL, aesthetics = c("colour", "fill")) +
    ggplot2::scale_y_continuous(limits = c(0, 100), expand = c(0, 0),
                                breaks = c(0, 50, 100)) +
    ggplot2::scale_x_continuous(breaks = c(2000, 2012, 2024), expand = c(0.01, 0)) +
    ggplot2::labs(x = NULL, y = "Coverage (%)", tag = "a") +
    theme_nm() +
    ggplot2::theme(legend.position = "top",
                   panel.spacing = ggplot2::unit(0.7, "lines"),
                   plot.tag = ggplot2::element_text(face = "bold", size = 13))
}

# ---- panel b: country x antigen decline heatmap ----------------------------
.fig1_panel_b <- function(decl, country_order) {
  d <- decl %>%
    dplyr::mutate(country = factor(country, levels = country_order),
                  vaccine = factor(vaccine, levels = VACCINE_LEVELS))
  ggplot2::ggplot(d, ggplot2::aes(x = country, y = vaccine, fill = decline_pp)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(is.na(decline_pp), "",
                                                   sprintf("%.0f", decline_pp)),
                                    colour = decline_pp > 28),
                       size = 2.7, fontface = "bold", show.legend = FALSE) +
    ggplot2::scale_fill_gradientn(
      colours = c("#FFF7EC", "#FEE8C8", "#FDBB84", "#FC8D59", "#D7301F", "#7F0000"),
      limits = c(0, NA), name = "Coverage drop\n(percentage points)",
      guide = ggplot2::guide_colourbar(barheight = 0.5, barwidth = 9,
                                       title.position = "top")) +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey25")) +
    ggplot2::scale_y_discrete(limits = rev(VACCINE_LEVELS)) +
    ggplot2::labs(x = NULL, y = NULL, tag = "b") +
    theme_nm() +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      plot.tag = ggplot2::element_text(face = "bold", size = 13))
}

# ---- public: refactored Figure 1 ------------------------------------------
# gap_df kept in the signature for backward compatibility (no longer required;
# Figure 1 is now descriptive). Returns a patchwork object; save with save_fig().
make_fig1 <- function(coverage_long, gap_df = NULL,
                      conflict_df = if (exists("conflict_info")) conflict_info else NULL,
                      vaccines = VACCINE_LEVELS, ref_years = 3L) {
  stopifnot(is.data.frame(conflict_df), "country" %in% names(conflict_df))
  decl <- .fig1_declines(coverage_long, conflict_df, vaccines, ref_years)
  country_order <- decl %>%
    dplyr::group_by(country) %>%
    dplyr::summarise(d = mean(decline_pp, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(d)) %>% dplyr::pull(country)
  
  pa <- .fig1_panel_a(coverage_long, conflict_df, decl, country_order)
  pb <- .fig1_panel_b(decl, country_order)
  
  if (!requireNamespace("patchwork", quietly = TRUE))
    stop("make_fig1 needs {patchwork}; install.packages('patchwork').")
  patchwork::wrap_plots(pa, pb, ncol = 1, heights = c(3.4, 1)) +
    patchwork::plot_annotation(
      title = "Vaccination coverage collapsed across conflict-affected countries",
      subtitle = paste0("Observed coverage (lines) for three antigens; shaded area is the coverage lost ",
                        "below the pre-conflict level from onset onward. Conflict period shaded grey. ",
                        "Countries ordered by mean decline."),
      caption = "Source: WHO/UNICEF WUENIC. Decline = pre-conflict 3-yr mean minus within-conflict trough (percentage points).",
      theme = theme_nm() + ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 13),
        plot.subtitle = ggplot2::element_text(size = 9, colour = "grey25")))
}

# ---- Figure 2: Forest plot of YLLs ---------------------------------------
# NOTE (recommendation #2): undiscounted is the headline metric (GBD-standard
# descriptive burden). 3% discounting is retained as a sensitivity (Table 2 and
# the `discount = "disc"` toggle).
make_fig2 <- function(yll_country_df, lambda_focus = 0,
                      framing_focus = FRAMING_HEADLINE, discount = "undisc") {
  use_mean <- paste0("yll_", discount, "_total_mean")
  use_lo   <- paste0("yll_", discount, "_total_lo")
  use_hi   <- paste0("yll_", discount, "_total_hi")
  
  d <- yll_country_df %>%
    dplyr::filter(catchup == lambda_focus, framing == framing_focus) %>%
    dplyr::mutate(mean_val = .data[[use_mean]],
                  lo       = .data[[use_lo]],
                  hi       = .data[[use_hi]]) %>%
    dplyr::arrange(disease, mean_val) %>%
    dplyr::mutate(country = factor(country, levels = unique(country)))
  
  if (nrow(d) == 0)
    stop(sprintf("make_fig2: no rows for catchup=%s, framing='%s'. Available framings: %s",
                 lambda_focus, framing_focus,
                 paste(unique(yll_country_df$framing), collapse = ", ")), call. = FALSE)
  
  ggplot2::ggplot(d, ggplot2::aes(x = mean_val, y = country, colour = method)) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, xmax = hi), height = 0,
                            position = ggplot2::position_dodge(width = 0.7),
                            linewidth = 0.4) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.7), size = 1.6) +
    ggplot2::facet_wrap(~ disease, scales = "free_x", ncol = 5) +
    ggplot2::scale_colour_manual(values = nm_palette, name = "Method") +
    ggplot2::scale_x_continuous(labels = scales::comma_format(),
                                expand = ggplot2::expansion(mult = c(0.05, 0.12))) +
    ggplot2::labs(
      x = paste0("Conflict-attributable YLLs (",
                 if (discount == "disc") "discounted, 3%" else "undiscounted", ")"),
      y = NULL,
      title = "Figure 2 | Conflict-attributable YLLs by country and disease",
      subtitle = paste0("No catch-up vaccination (lambda=0). 95% CIs from Monte Carlo, ",
                        scales::comma(N_SIM), " draws."),
      caption = "CIs propagate uncertainty in coverage gap, incidence, CFR, and vaccine effectiveness."
    ) +
    theme_nm()
}

# ---- Figure 3: Catch-up sensitivity grid ---------------------------------
make_fig3 <- function(yll_disease_df, discount = "undisc") {
  m_mean <- paste0("yll_", discount, "_mean")
  m_lo   <- paste0("yll_", discount, "_lo")
  m_hi   <- paste0("yll_", discount, "_hi")
  d <- yll_disease_df %>%
    dplyr::mutate(
      framing_label = factor(unname(FRAMING_LABELS[as.character(framing)]),
                             levels = unname(FRAMING_LABELS))
    )
  if (anyNA(d$framing_label))
    stop("make_fig3: framing value(s) not in FRAMING_LABELS: ",
         paste(setdiff(unique(d$framing), names(FRAMING_LABELS)), collapse = ", "),
         call. = FALSE)
  ggplot2::ggplot(d, ggplot2::aes(x = catchup, y = .data[[m_mean]],
                                  colour = disease, fill = disease)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data[[m_lo]], ymax = .data[[m_hi]]),
                         alpha = 0.15, colour = NA) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::geom_point(size = 1.4) +
    ggplot2::facet_grid(framing_label ~ method, scales = "free_y") +
    ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                                breaks = CATCHUP_RATES) +
    ggplot2::scale_y_continuous(labels = scales::comma_format()) +
    ggplot2::scale_colour_manual(values = nm_disease_palette, name = "Disease") +
    ggplot2::scale_fill_manual(values = nm_disease_palette, guide = "none") +
    ggplot2::labs(
      x = "Catch-up intensity (lambda)",
      y = paste0("Total YLLs (",
                 if (discount == "disc") "discounted" else "undiscounted", ")"),
      title = "Figure 3 | Catch-up sensitivity by method and framing",
      subtitle = "Routine recovery: geometric convergence of later cohorts to the counterfactual. Campaign top-up: lambda x reach of the zero-dose cohort.",
      caption = "Each panel sums YLLs across all 10 countries for the disease shown."
    ) +
    theme_nm()
}

# ---- Tables ---------------------------------------------------------------

# F6: distinguish a true ITS null from a FLOORED zero. ITS returns exactly 0 for
# all diseases in Ethiopia/Iraq because the rising pre-conflict trend's
# counterfactual overtakes observed coverage (cf >= obs -> gap <= 0, floored to 0
# YLL in 04). That is a structural limitation, not a true null, and not
# comparable to the donor-based positives. Flag (country,disease) cells where ITS
# was estimated, has a rising pre-trend, and a non-positive central gap.
its_floored_flags <- function(gap_df, vaccine_disease_links) {
  if (!"slope_logit_per_yr" %in% names(gap_df)) return(tibble::tibble())
  gap_df %>%
    dplyr::filter(method == "ITS") %>%
    dplyr::left_join(vaccine_disease_links, by = "vaccine") %>%
    dplyr::group_by(country, disease) %>%
    dplyr::summarise(
      its_floored = any(is.finite(slope_logit_per_yr)) &&
        mean(gap, na.rm = TRUE) <= 0 &&
        mean(slope_logit_per_yr, na.rm = TRUE) > 0,
      .groups = "drop")
}

# F5: headline totals. The SC-anchored headline silently drops Pakistan, Somalia
# and Sri Lanka (SC has no cell there) = 22% of burden under DiD, Pakistan ~2.1M
# alone. Report BOTH (i) the SC 7-country complete-case total, LABELLED as such,
# and (ii) a 10-country total using a per-country best-available estimator
# (prefer SC; fall back to DiD, then Baseline 3yr) so no country is silently 0.
# CIs are draw-based. `yll_draws` is the per disease x country x method object.
assemble_headline_totals <- function(yll_draws,
                                     anchor_priority = c("Synthetic Control",
                                                         "DiD (CS)", "Baseline (3yr)"),
                                     catchup_focus = 0,
                                     framing_focus = FRAMING_HEADLINE) {
  d <- yll_draws %>% dplyr::filter(catchup == catchup_focus, framing == framing_focus)
  keep <- vapply(d$d_undisc, function(x) !is.null(x) && length(x) > 0, logical(1))
  d <- d[keep, , drop = FALSE]
  if (nrow(d) == 0) return(NULL)
  sc <- d %>% dplyr::filter(method == anchor_priority[1])
  sc_tot <- if (nrow(sc) > 0) colSums(do.call(rbind, sc$d_undisc)) else NA_real_
  pick <- d %>% dplyr::filter(method %in% anchor_priority) %>%
    dplyr::mutate(prio = match(method, anchor_priority)) %>%
    dplyr::group_by(disease, country) %>%
    dplyr::slice_min(prio, n = 1, with_ties = FALSE) %>% dplyr::ungroup()
  comp_tot <- colSums(do.call(rbind, pick$d_undisc))
  q <- function(v) if (length(v) == 1 && is.na(v)) c(mean = NA, lo = NA, hi = NA) else
    c(mean = mean(v), lo = stats::quantile(v, .025, names = FALSE),
      hi = stats::quantile(v, .975, names = FALSE))
  list(sc_7country = q(sc_tot), composite_10country = q(comp_tot),
       n_countries_sc = dplyr::n_distinct(sc$country),
       n_countries_comp = dplyr::n_distinct(pick$country),
       anchor_mix = pick %>% dplyr::count(method, name = "n_cells"))
}

make_table1 <- function(yll_country_df, framing_focus = FRAMING_HEADLINE,
                        its_floored = NULL) {
  # Headline = undiscounted (recommendation #2); discounted totals are in Table 2.
  long <- yll_country_df %>%
    dplyr::filter(catchup == 0, framing == framing_focus) %>%
    dplyr::transmute(
      Country = country,
      Disease = disease,
      Method  = method,
      YLLs = fmt_ci(yll_undisc_total_mean, yll_undisc_total_lo,
                    yll_undisc_total_hi, digits = 0))
  # F6: append a dagger to floored ITS cells so a reader does not read them as a
  # true zero burden. Footnote text accompanies Table 1 in the manuscript.
  if (!is.null(its_floored) && nrow(its_floored) > 0) {
    long <- long %>%
      dplyr::left_join(
        its_floored %>% dplyr::transmute(Country = country, Disease = disease, its_floored),
        by = c("Country", "Disease")) %>%
      dplyr::mutate(YLLs = ifelse(Method == "ITS" &
                                    !is.na(its_floored) & its_floored & YLLs != "--",
                                  paste0(YLLs, "\u2020"), YLLs)) %>%
      dplyr::select(-its_floored)
  }
  long %>%
    tidyr::pivot_wider(names_from = Method, values_from = YLLs,
                       values_fill = "--") %>%
    dplyr::arrange(Country, Disease)
}

make_table2 <- function(yll_draws, framing_focus = FRAMING_HEADLINE,
                        catchup_focus = 0) {
  # Draw-level complete-case aggregation (recommendation #4). Each method drops a
  # different set of countries (SC fails without a donor pool; ITS fails with too
  # short a pre-period). Summing each method over whatever it happens to cover
  # makes the totals non-comparable, so we restrict each disease to the country
  # set present under ALL methods, then sum the per-country DRAW VECTORS within
  # method (CRN-correct), taking empirical quantiles. This replaces the prior
  # symmetric-normal, cross-country-independent CI combination.
  agg <- aggregate_method_comparison_draws(yll_draws,
                                           catchup_focus = catchup_focus,
                                           framing_focus = framing_focus)
  if (nrow(agg) == 0) {
    warning("make_table2: no draw rows for catchup=", catchup_focus,
            ", framing='", framing_focus, "'.")
    return(tibble::tibble())
  }
  agg %>%
    dplyr::transmute(
      Disease = disease,
      Method  = method,
      `N countries`     = n_countries,
      YLLs_discounted   = fmt_ci(yll_disc_mean, yll_disc_lo, yll_disc_hi, 0),
      YLLs_undiscounted = fmt_ci(yll_undisc_mean, yll_undisc_lo, yll_undisc_hi, 0)
    ) %>%
    dplyr::arrange(Disease, Method)
}

make_table_s1 <- function() {
  params %>%
    # Drop supplementary-only rows (neonatal-tetanus arm: Seale CFR, maternal-TT
    # VE). They remain in `params` for the TETANUS_ARM='neonatal' sensitivity but
    # belong in the Table S1 note, not the primary table. Helper defined in 01;
    # guard so this still runs if 01's helper is somehow unavailable.
    dplyr::filter(
      if (exists("is_supplementary_row"))
        !is_supplementary_row(parameter, population) else TRUE
    ) %>%
    dplyr::transmute(
      Parameter   = parameter,
      Population  = population,
      `Point estimate` = sprintf("%.4f", mean),
      `95% interval`   = sprintf("%.4f - %.4f", lo, hi),
      Distribution = dist,
      Source       = source
    )
}

# ---- Save helpers ---------------------------------------------------------
save_fig <- function(plot, file, width = 11, height = 6, dpi = 600) {
  for (ext in c("pdf", "png")) {
    path <- file.path(FIG_DIR, paste0(tools::file_path_sans_ext(file), ".", ext))
    ggplot2::ggsave(path, plot = plot, width = width, height = height,
                    dpi = dpi, units = "in", bg = "white")
    message("Saved figure: ", path)
  }
}

save_tab <- function(tbl, file) {
  path <- file.path(TAB_DIR, file)
  utils::write.csv(tbl, path, row.names = FALSE)
  message("Saved table: ", path)
}
# ============================================================================
# Figure 4 | Method triangulation (Nature-style multi-panel comparison)
# ----------------------------------------------------------------------------
# Panel A  Event study: Callaway-Sant'Anna dynamic ATT on coverage by years
#          since conflict onset (negative = conflict-attributable shortfall).
#          The pre-onset coefficients are the parallel-trends falsification:
#          they should sit on zero. (Sun-Abraham line overlaid if supplied.)
# Panel B  Mean conflict-attributable coverage gap (pp) by estimator, with the
#          baseline-window sensitivity (3yr/5yr) shown alongside. This is where
#          the YLL differences originate, since YLL scales ~linearly in the gap.
# Panel C  Total conflict-attributable YLL (undiscounted, no catch-up) by
#          estimator: the bottom-line "do the methods agree?" panel.
#
# Inputs:
#   gap_df            estimate_all_gaps_plus() output (carries attr 'event_study')
#   global_by_method  aggregate_global_by_method_draws() output
#   baseline_windows  gap_baseline_windows() output (optional; Panel B band)
#   es_df             event-study tibble override (else attr(gap_df,'event_study'))
#   sa_df             optional Sun-Abraham event study (did_event_study_sa())
make_fig4_method_comparison <- function(gap_df, global_by_method,
                                        baseline_windows = NULL,
                                        es_df = NULL, sa_df = NULL,
                                        vaccine_focus = "MCV1") {
  if (!requireNamespace("patchwork", quietly = TRUE))
    stop("make_fig4 needs the 'patchwork' package.")
  if (is.null(es_df)) es_df <- attr(gap_df, "event_study")
  
  meth_levels <- METHOD_ORDER
  
  # ---- Panel A: event study -------------------------------------------------
  pA <- NULL
  if (!is.null(es_df) && nrow(es_df) > 0) {
    esd <- es_df
    if ("vaccine" %in% names(esd) && vaccine_focus %in% esd$vaccine)
      esd <- esd %>% dplyr::filter(vaccine == vaccine_focus)
    esd <- esd %>% dplyr::mutate(lo = att - 1.96 * se, hi = att + 1.96 * se,
                                 src = "Callaway-Sant'Anna")
    if (!is.null(sa_df) && nrow(sa_df) > 0) {
      sad <- sa_df
      if ("vaccine" %in% names(sad) && vaccine_focus %in% sad$vaccine)
        sad <- sad %>% dplyr::filter(vaccine == vaccine_focus)
      sad <- sad %>% dplyr::mutate(lo = att - 1.96 * se, hi = att + 1.96 * se,
                                   src = "Sun-Abraham")
      esd <- dplyr::bind_rows(esd, sad)
    }
    pA <- ggplot2::ggplot(esd, ggplot2::aes(x = event_time, y = att,
                                            colour = src, fill = src)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.3) +
      ggplot2::geom_vline(xintercept = -0.5, colour = "grey55",
                          linetype = "dashed", linewidth = 0.3) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
      ggplot2::geom_line(linewidth = 0.55) +
      ggplot2::geom_point(size = 1.3) +
      ggplot2::scale_colour_manual(values = c("Callaway-Sant'Anna" = "#2C7A7B",
                                              "Sun-Abraham" = "#C44536"), name = NULL) +
      ggplot2::scale_fill_manual(values = c("Callaway-Sant'Anna" = "#2C7A7B",
                                            "Sun-Abraham" = "#C44536"), guide = "none") +
      ggplot2::labs(x = "Years since conflict onset",
                    y = paste0("ATT on ", vaccine_focus, " coverage (pp)"),
                    subtitle = "a  Event study (staggered DiD): pre-onset coefficients test parallel trends") +
      theme_nm()
  }
  
  # ---- Panel B: mean coverage gap by method ---------------------------------
  gap_pool <- gap_df %>%
    dplyr::filter(is.finite(gap)) %>%
    dplyr::select(method, country, vaccine, year, gap)
  if (!is.null(baseline_windows) && nrow(baseline_windows) > 0)
    gap_pool <- dplyr::bind_rows(
      gap_pool,
      baseline_windows %>% dplyr::filter(is.finite(gap)) %>%
        dplyr::select(method, country, vaccine, year, gap))
  gapB <- gap_pool %>%
    dplyr::group_by(method) %>%
    dplyr::summarise(n = dplyr::n(), m = mean(gap),
                     se = stats::sd(gap) / sqrt(dplyr::n()), .groups = "drop") %>%
    dplyr::mutate(lo = m - 1.96 * se, hi = m + 1.96 * se,
                  method = factor(method,
                                  levels = rev(c(meth_levels,
                                                 setdiff(unique(method), meth_levels)))))
  pB <- ggplot2::ggplot(gapB, ggplot2::aes(x = m, y = method, colour = method)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.5) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_colour_manual(values = nm_palette, guide = "none") +
    ggplot2::labs(x = "Mean conflict-attributable gap (pp)", y = NULL,
                  subtitle = "b  Coverage gap by estimator (pooled over country-vaccine-year)") +
    theme_nm()
  
  # ---- Panel C: total YLL by method -----------------------------------------
  pC <- NULL
  if (!is.null(global_by_method) && nrow(global_by_method) > 0) {
    gC <- global_by_method %>%
      dplyr::mutate(method = factor(method,
                                    levels = rev(intersect(meth_levels, method))))
    pC <- ggplot2::ggplot(gC, ggplot2::aes(x = yll_undisc_mean, y = method, colour = method)) +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin = yll_undisc_lo, xmax = yll_undisc_hi),
                              height = 0, linewidth = 0.5) +
      ggplot2::geom_point(size = 2) +
      ggplot2::scale_colour_manual(values = nm_palette, guide = "none") +
      ggplot2::scale_x_continuous(labels = scales::comma_format()) +
      ggplot2::labs(x = "Total conflict-attributable YLLs (undiscounted)", y = NULL,
                    subtitle = "c  Bottom-line YLL by estimator (complete-case, draw-level 95% CI)") +
      theme_nm()
  }
  
  # ---- Assemble -------------------------------------------------------------
  bottom <- if (!is.null(pC)) (pB | pC) else pB
  combined <- if (!is.null(pA)) (pA / bottom) + patchwork::plot_layout(heights = c(1, 1)) else bottom
  combined +
    patchwork::plot_annotation(
      title = "Figure 4 | Triangulating conflict-attributable coverage loss across six estimators",
      subtitle = "Agreement across no-control (Baseline, ITS) and comparison-group (SC, SC+cov, augmented SC, DiD) methods is the credibility argument.",
      caption = "Comparison-group estimators net out shocks shared by non-conflict donors (incl. COVID-19); no-control methods do not.",
      theme = theme_nm())
}