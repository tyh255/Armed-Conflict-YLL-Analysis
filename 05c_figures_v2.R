# ============================================================================
# 05c_FIGURES_V2.R   -   Redesigned main figures (Nature spec)
# ----------------------------------------------------------------------------
# Four main figures, rebuilt to the brief and to Nature's artwork rules
# (89 mm single / 183 mm double / 247 mm max height; Helvetica/Arial; panel
# letters 8 pt bold lowercase upright; body text 5-7 pt; colour-blind-safe; no
# red-green encoding; ticks + units on axes; thousands with commas).
#
#   Figure 1  Coverage decline, BROKEN OUT BY ANTIGEN (vaccine x country grid),
#             coverage loss shaded RED; companion drop heatmap with an explicit,
#             switchable decline metric (default = pre-conflict mean -> trough).
#   Figure 2  Cause -> consequence in one 2x2: conflict coverage GAPS by country
#             (a) and method (b); attributable YLL by country (c) and method (d).
#   Figure 3  Triangulation: composite tornado (a); per-method draw-distribution
#             ridgeline (b, colourful); 2-D density contour of the two inferential
#             families draw-by-draw (c) - the CRN-linked triangulation.
#   Figure 4  Per-capita intensity (a); onset-timing falsification (b); external
#             validation heatmap vs GBD with over-attribution flagged (c).
#
# Source AFTER 05 / 05b / 06 / 08 / 09 / 10 (reuses theme_nm, nm_palette,
# VACCINE_*, METHOD_ORDER, FIG4_* and the draw aggregators). Each builder takes
# in-memory objects first and falls back to the saved CSVs, mirroring 05b.
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
})
stopifnot(requireNamespace("ggplot2", quietly = TRUE))

# ---- Nature artwork constants ---------------------------------------------
NAT_W_SINGLE <- 89; NAT_W_HALF <- 120; NAT_W_DOUBLE <- 183; NAT_H_MAX <- 247  # mm
NAT_FONT     <- ""   # "" = device sans (Helvetica/Arial-like); set "Helvetica"
# for final artwork on a cairo device.

# Print-tuned theme: inherits the house theme_nm() look but enforces 7-pt body /
# 8-pt bold lowercase panel tags so text lands in spec at final print size.
theme_nat <- function(base = 7) {
  base_th <- if (exists("theme_nm", mode = "function")) theme_nm(base) else
    ggplot2::theme_classic(base_size = base)
  base_th + ggplot2::theme(
    text            = ggplot2::element_text(family = NAT_FONT),
    plot.title      = ggplot2::element_text(face = "bold", size = base + 1, hjust = 0),
    plot.subtitle   = ggplot2::element_text(size = base - 1, colour = "grey25"),
    plot.caption    = ggplot2::element_text(size = base - 2, colour = "grey35", hjust = 0),
    plot.tag        = ggplot2::element_text(face = "bold", size = 8),  # a, b, c ...
    legend.key.size = ggplot2::unit(3.2, "mm"),
    legend.text     = ggplot2::element_text(size = base - 1),
    legend.title    = ggplot2::element_text(size = base - 1, face = "bold"))
}

# Save at physical mm dimensions; vector PDF (editable text for production) +
# RGB PNG preview at 600 dpi (line-art spec).
save_fig_nat <- function(plot, file, width_mm = NAT_W_DOUBLE, height_mm = 120,
                         dpi = 600, dir = if (exists("FIG_DIR")) FIG_DIR else ".") {
  height_mm <- min(height_mm, NAT_H_MAX)
  stem <- file.path(dir, tools::file_path_sans_ext(file))
  dev_pdf <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  ggplot2::ggsave(paste0(stem, ".pdf"), plot, width = width_mm, height = height_mm,
                  units = "mm", device = dev_pdf)
  ggplot2::ggsave(paste0(stem, ".png"), plot, width = width_mm, height = height_mm,
                  units = "mm", dpi = dpi, bg = "white")
  message("Saved figure: ", stem, ".{pdf,png}  (", width_mm, "x", height_mm, " mm)")
  invisible(stem)
}

.comma  <- function(x) format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE)
.mized  <- function(x) x / 1e6
.fam_of <- function(m) if (exists("FIG4_FAMILY")) FIG4_FAMILY(m) else "All methods"
.fam_col <- if (exists("FIG4_COL")) FIG4_COL else
  c("Design / observed-trend" = "#356B87", "Donor-comparison" = "#C16E2F",
    "Baseline (attenuated)" = "#BFB6A8")

# ============================================================================
# FIGURE 1  -  coverage decline, broken out by antigen
# ============================================================================
# Decline metric is now EXPLICIT and switchable (answers "is the heatmap the
# first-year drop?" -> no, by default it is the maximum drop, pre-conflict mean
# to within-conflict trough). Options:
#   "trough"  pre-conflict mean  -> minimum coverage in the conflict window  (max drop)
#   "onset"   pre-conflict mean  -> coverage in the FIRST conflict year      (on-impact)
#   "endwin"  pre-conflict mean  -> coverage in the LAST conflict year       (sustained)
.f1_decline <- function(coverage_long, conflict_df, vaccines = VACCINE_LEVELS,
                        ref_years = 3L, metric = c("trough", "onset", "endwin")) {
  metric <- match.arg(metric)
  coverage_long %>%
    dplyr::filter(vaccine %in% vaccines, country %in% conflict_df$country, !is.na(coverage)) %>%
    dplyr::inner_join(dplyr::select(conflict_df, country, conflict_start, conflict_end),
                      by = "country") %>%
    dplyr::group_by(country, vaccine) %>%
    dplyr::summarise(
      pre_level = {
        pre <- coverage[year >= conflict_start[1] - ref_years & year < conflict_start[1]]
        if (length(pre) == 0) coverage[which.min(abs(year - conflict_start[1]))] else mean(pre, na.rm = TRUE)
      },
      post_val = {
        cs <- conflict_start[1]; ce <- conflict_end[1]
        win <- coverage[year >= cs & year <= ce]; yrs <- year[year >= cs & year <= ce]
        if (length(win) == 0) NA_real_ else switch(metric,
                                                   trough = min(win, na.rm = TRUE),
                                                   onset  = win[which.min(yrs)],
                                                   endwin = win[which.max(yrs)])
      },
      .groups = "drop") %>%
    dplyr::mutate(drop_pp = pmax(pre_level - post_val, 0))
}

.F1_METRIC_LAB <- c(
  trough = "Maximum coverage drop (pre-conflict mean \u2192 conflict-period trough)",
  onset  = "On-impact coverage drop (pre-conflict mean \u2192 first conflict year)",
  endwin = "End-of-window coverage drop (pre-conflict mean \u2192 last conflict year)")

# panel a: vaccine (rows) x country (cols) small multiples; red loss ribbon.
.f1_panel_a <- function(coverage_long, conflict_df, decl, country_order,
                        year_min = 1995L, year_max = 2024L,
                        loss_col = "#D7301F") {
  cov <- coverage_long %>%
    dplyr::filter(vaccine %in% VACCINE_LEVELS, country %in% country_order,
                  year >= year_min, year <= year_max) %>%
    dplyr::mutate(country = factor(country, levels = country_order),
                  vaccine = factor(vaccine, levels = VACCINE_LEVELS))
  shade <- conflict_df %>% dplyr::filter(country %in% country_order) %>%
    dplyr::transmute(country = factor(country, levels = country_order),
                     xmin = conflict_start, xmax = conflict_end, x0 = conflict_start)
  ribbon <- cov %>%
    dplyr::left_join(dplyr::select(decl, country, vaccine, pre_level), by = c("country", "vaccine")) %>%
    dplyr::left_join(dplyr::select(conflict_df, country, conflict_start), by = "country") %>%
    dplyr::filter(year >= conflict_start) %>%
    dplyr::mutate(ylo = pmin(coverage, pre_level), yhi = pre_level)
  
  # Per-antigen y-range driven by the ACTUALLY OBSERVED coverage across all
  # countries for that antigen (brief: if min observed BCG across the 10 countries
  # is 40, the BCG row should start at 40). We free the y-scale per facet row
  # (vaccine) and inject an invisible blank layer at each row's observed min/max
  # (padded slightly) so every column in a row shares that antigen's data range
  # rather than a forced 0-100 box. Bounds also cover the pre-conflict ribbon top.
  yrng <- cov %>%
    dplyr::left_join(dplyr::select(decl, country, vaccine, pre_level),
                     by = c("country", "vaccine")) %>%
    dplyr::group_by(vaccine) %>%
    dplyr::summarise(
      ymin = min(coverage, na.rm = TRUE),
      ymax = max(c(coverage, pre_level), na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(
      pad  = pmax((ymax - ymin) * 0.04, 1),
      ymin = pmax(floor(ymin - pad), 0),
      ymax = pmin(ceiling(ymax + pad), 100))
  # Two-row blank frame (min, max) per antigen to pin each row's free y-range.
  # vaccine must be the SAME factor type/levels as in `cov` so geom_blank lands
  # in the correct facet rows (a character/factor mismatch can desync faceting).
  blank <- yrng %>%
    tidyr::pivot_longer(c(ymin, ymax), names_to = NULL, values_to = "coverage") %>%
    dplyr::mutate(vaccine = factor(as.character(vaccine), levels = VACCINE_LEVELS),
                  country = factor(country_order[1], levels = country_order),
                  year    = year_min)
  
  ggplot2::ggplot() +
    ggplot2::geom_rect(data = shade, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
                       fill = "grey88", alpha = 0.5, inherit.aes = FALSE) +
    ggplot2::geom_vline(data = shade, ggplot2::aes(xintercept = x0),
                        colour = "grey55", linewidth = 0.25, linetype = "22") +
    ggplot2::geom_blank(data = blank, ggplot2::aes(x = year, y = coverage)) +
    ggplot2::geom_ribbon(data = ribbon, ggplot2::aes(x = year, ymin = ylo, ymax = yhi),
                         fill = loss_col, alpha = 0.30) +
    ggplot2::geom_line(data = cov, ggplot2::aes(x = year, y = coverage, colour = vaccine),
                       linewidth = 0.55) +
    ggplot2::facet_grid(vaccine ~ country, scales = "free_y", labeller = ggplot2::labeller(
      vaccine = function(v) VACCINE_LABS[v])) +
    ggplot2::scale_colour_manual(values = VACCINE_COLS, guide = "none") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.02)),
                                breaks = scales::pretty_breaks(3)) +
    ggplot2::scale_x_continuous(breaks = c(2000, 2020), expand = c(0.02, 0)) +
    ggplot2::labs(x = NULL, y = "Coverage (%)") +
    theme_nat() +
    ggplot2::theme(panel.spacing = ggplot2::unit(0.35, "lines"),
                   strip.text.y = ggplot2::element_text(angle = 0, size = 6),
                   strip.text.x = ggplot2::element_text(size = 6),
                   axis.text = ggplot2::element_text(size = 5))
}

# panel c: number of years of conflict per country (conflict_end - conflict_start),
# ordered to match panel a/b. Descriptive companion requested in review.
.f1_panel_c <- function(conflict_df, country_order) {
  d <- conflict_df %>%
    dplyr::filter(country %in% country_order) %>%
    dplyr::transmute(country = factor(country, levels = country_order),
                     years = conflict_end - conflict_start + 1L)
  ggplot2::ggplot(d, ggplot2::aes(x = years, y = country)) +
    ggplot2::geom_col(width = 0.66, fill = "#4F6D7A") +
    ggplot2::geom_text(ggplot2::aes(label = years), hjust = -0.25,
                       size = 2.1, colour = "grey25") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.18)),
                                breaks = scales::pretty_breaks(4)) +
    ggplot2::labs(x = "Years of conflict", y = NULL,
                  subtitle = "Coded conflict duration (onset \u2192 end, inclusive)") +
    theme_nat()
}

# panel b: country x antigen decline heatmap, warm sequential ramp.
.f1_panel_b <- function(decl, country_order, metric = "trough") {
  d <- decl %>% dplyr::mutate(country = factor(country, levels = country_order),
                              vaccine = factor(vaccine, levels = VACCINE_LEVELS))
  hi_cut <- stats::quantile(d$drop_pp, 0.7, na.rm = TRUE)
  ggplot2::ggplot(d, ggplot2::aes(x = country, y = vaccine, fill = drop_pp)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(is.na(drop_pp), "", sprintf("%.0f", drop_pp)),
                                    colour = drop_pp > hi_cut),
                       size = 2.1, fontface = "bold", show.legend = FALSE) +
    ggplot2::scale_fill_gradientn(
      colours = c("#FFF7EC", "#FEE8C8", "#FDBB84", "#FC8D59", "#D7301F", "#7F0000"),
      limits = c(0, NA), name = "Coverage drop (pp)",
      guide = ggplot2::guide_colourbar(barheight = 0.4, barwidth = 7, title.position = "top")) +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey25")) +
    ggplot2::scale_y_discrete(limits = rev(VACCINE_LEVELS)) +
    ggplot2::labs(x = NULL, y = NULL,
                  subtitle = unname(.F1_METRIC_LAB[metric])) +
    theme_nat() +
    ggplot2::theme(legend.position = "bottom", panel.grid = ggplot2::element_blank(),
                   axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 35, hjust = 1, size = 6))
}

fig1_v2 <- function(coverage_long, conflict_df = if (exists("conflict_info")) conflict_info else NULL,
                    drop_metric = "trough", ref_years = 3L) {
  stopifnot(is.data.frame(conflict_df), "country" %in% names(conflict_df))
  decl <- .f1_decline(coverage_long, conflict_df, ref_years = ref_years, metric = drop_metric)
  country_order <- decl %>% dplyr::group_by(country) %>%
    dplyr::summarise(d = mean(drop_pp, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(d)) %>% dplyr::pull(country)
  # Panel c orders countries by conflict duration (longest at top) for its own
  # readability, independent of panel a/b's decline ordering.
  dur_order <- conflict_df %>%
    dplyr::filter(country %in% country_order) %>%
    dplyr::mutate(years = conflict_end - conflict_start + 1L) %>%
    dplyr::arrange(years) %>% dplyr::pull(country)
  pa <- .f1_panel_a(coverage_long, conflict_df, decl, country_order)
  pb <- .f1_panel_b(decl, country_order, drop_metric)
  pc <- .f1_panel_c(conflict_df, dur_order)
  if (!requireNamespace("patchwork", quietly = TRUE)) stop("fig1_v2 needs {patchwork}.")
  patchwork::wrap_plots(pa, (pb | pc) + patchwork::plot_layout(widths = c(1.5, 1)),
                        ncol = 1, heights = c(2.5, 1)) +
    patchwork::plot_annotation(
      tag_levels = "a",
      title = "Vaccination coverage declined across conflict-affected countries",
      subtitle = paste0("a, Observed coverage by antigen (rows, y-axis scaled to the observed range per antigen) ",
                        "and country (columns, ordered by mean decline); red area is coverage lost below the ",
                        "pre-conflict level from onset; conflict period shaded. b, Country \u00d7 antigen coverage drop. ",
                        "c, Years of conflict per country."),
      caption = "Source: WHO/UNICEF WUENIC. Descriptive (model-free); counterfactual estimators appear in Figs 2-3.",
      theme = theme_nat())
}

# ============================================================================
# FIGURE 2  -  coverage gaps (a, b) -> attributable YLL (c, d)
# ============================================================================
# Coverage GAP (percentage points, counterfactual - observed) is keyed by
# ANTIGEN in gap_df. Panels a/b use difference-in-differences as the common
# estimator (available for all ten countries; the primary in the main text);
# panel d shows every estimator's total to make the convergence explicit.
.f2_post_gap <- function(gap_df, conflict_df, method_keep = NULL) {
  g <- gap_df %>%
    dplyr::inner_join(dplyr::select(conflict_df, country, conflict_start, conflict_end), by = "country") %>%
    dplyr::filter(year >= conflict_start, year <= conflict_end)
  if (!is.null(method_keep)) g <- g %>% dplyr::filter(method %in% method_keep)
  g
}

.f2_gap_by_country <- function(gap_df, conflict_df, method = "DiD (CS)") {
  if (!method %in% unique(gap_df$method)) method <- unique(gap_df$method)[1]
  g <- .f2_post_gap(gap_df, conflict_df, method) %>%
    dplyr::group_by(country, vaccine) %>%
    dplyr::summarise(gap = mean(gap, na.rm = TRUE), .groups = "drop")
  g %>% dplyr::group_by(country) %>%
    dplyr::summarise(gap_mean = mean(gap, na.rm = TRUE),
                     gap_lo = min(gap, na.rm = TRUE), gap_hi = max(gap, na.rm = TRUE),
                     .groups = "drop") %>%
    dplyr::arrange(gap_mean) %>%
    dplyr::mutate(country = factor(country, levels = country), .method = method)
}

.f2_gap_by_method <- function(gap_df, conflict_df) {
  cell <- .f2_post_gap(gap_df, conflict_df) %>%
    dplyr::group_by(method, country, vaccine) %>%
    dplyr::summarise(gap = mean(gap, na.rm = TRUE), .groups = "drop")
  ord <- if (exists("METHOD_ORDER")) METHOD_ORDER else sort(unique(cell$method))
  cell %>% dplyr::group_by(method) %>%
    dplyr::summarise(gap_med = stats::median(gap, na.rm = TRUE),
                     gap_lo = stats::quantile(gap, .25, names = FALSE, na.rm = TRUE),
                     gap_hi = stats::quantile(gap, .75, names = FALSE, na.rm = TRUE),
                     .groups = "drop") %>%
    # ITS rising pre-trend counterfactuals can push the raw cell gap below zero
    # (the documented structural-zero behaviour). The paper treats those as
    # floored at zero (no negative attributable shortfall), so the displayed
    # summary is floored to match - a negative ITS bar would contradict the text.
    dplyr::mutate(gap_med = ifelse(method == "ITS", pmax(gap_med, 0), gap_med),
                  gap_lo  = ifelse(method == "ITS", pmax(gap_lo, 0), gap_lo),
                  gap_hi  = ifelse(method == "ITS", pmax(gap_hi, 0), gap_hi)) %>%
    dplyr::mutate(family = .fam_of(method),
                  method = factor(method, levels = rev(ord[ord %in% method])))
}

# Composite (best-available) per-country YLL with draw-based CI, mirroring the
# headline selection (SC > DiD > Baseline 3yr).
.f2_yll_by_country <- function(yll_draws,
                               anchor_priority = c("Synthetic Control", "DiD (CS)", "Baseline (3yr)"),
                               catchup_focus = 0,
                               framing_focus = if (exists("FRAMING_HEADLINE")) FRAMING_HEADLINE else "campaign_topup") {
  d <- yll_draws %>% dplyr::filter(catchup == catchup_focus, framing == framing_focus)
  keep <- vapply(d$d_undisc, function(x) !is.null(x) && length(x) > 0, logical(1))
  d <- d[keep, , drop = FALSE]
  pick <- d %>% dplyr::filter(method %in% anchor_priority) %>%
    dplyr::mutate(prio = match(method, anchor_priority)) %>%
    dplyr::group_by(disease, country) %>% dplyr::slice_min(prio, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  pick %>% dplyr::group_by(country) %>%
    dplyr::group_modify(~ {
      M <- do.call(rbind, .x$d_undisc); s <- colSums(M)
      tibble::tibble(mean = mean(s), lo = stats::quantile(s, .025, names = FALSE),
                     hi = stats::quantile(s, .975, names = FALSE))
    }) %>% dplyr::ungroup() %>%
    dplyr::arrange(mean) %>% dplyr::mutate(country = factor(country, levels = country))
}

.f2_yll_by_method <- function(yll_draws, catchup_focus = 0,
                              framing_focus = if (exists("FRAMING_HEADLINE")) FRAMING_HEADLINE else "campaign_topup") {
  if (exists("aggregate_global_by_method_draws")) {
    mt <- aggregate_global_by_method_draws(yll_draws, catchup_focus = catchup_focus,
                                           framing_focus = framing_focus) %>%
      dplyr::transmute(method, mean = yll_undisc_mean, lo = yll_undisc_lo, hi = yll_undisc_hi)
  } else {
    d <- yll_draws %>% dplyr::filter(catchup == catchup_focus, framing == framing_focus)
    keep <- vapply(d$d_undisc, function(x) !is.null(x) && length(x) > 0, logical(1)); d <- d[keep, ]
    mt <- d %>% dplyr::group_by(method) %>% dplyr::group_modify(~ {
      M <- do.call(rbind, .x$d_undisc); s <- colSums(M)
      tibble::tibble(mean = mean(s), lo = stats::quantile(s, .025, names = FALSE),
                     hi = stats::quantile(s, .975, names = FALSE)) }) %>% dplyr::ungroup()
  }
  ord <- if (exists("METHOD_ORDER")) METHOD_ORDER else sort(unique(mt$method))
  mt %>% dplyr::mutate(family = .fam_of(method),
                       method = factor(method, levels = rev(ord[ord %in% method])))
}

fig2_v2 <- function(gap_df, yll_draws, conflict_df = if (exists("conflict_info")) conflict_info else NULL,
                    headline_method = "Synthetic Control", gap_method = "DiD (CS)") {
  gc  <- .f2_gap_by_country(gap_df, conflict_df, gap_method)
  gm  <- .f2_gap_by_method(gap_df, conflict_df)
  yc  <- .f2_yll_by_country(yll_draws)
  ym  <- .f2_yll_by_method(yll_draws)
  vcol <- if (exists("VACCINE_COLS")) "#356B87" else "#356B87"
  
  pa <- ggplot2::ggplot(gc, ggplot2::aes(gap_mean, country)) +
    ggplot2::geom_segment(ggplot2::aes(x = gap_lo, xend = gap_hi, yend = country),
                          colour = "grey75", linewidth = 1.3) +
    ggplot2::geom_point(colour = vcol, size = 1.8) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.08))) +
    ggplot2::labs(x = "Mean within-conflict coverage gap (pp)", y = NULL,
                  subtitle = sprintf("Difference-in-differences; whisker = antigen range")) +
    theme_nat()
  
  pb <- ggplot2::ggplot(gm, ggplot2::aes(gap_med, method, colour = family)) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = gap_lo, xmax = gap_hi), height = 0, linewidth = 0.7) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_colour_manual(values = .fam_col, name = NULL) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.04, 0.08))) +
    ggplot2::labs(x = "Within-conflict coverage gap (pp)", y = NULL,
                  subtitle = "Pooled across countries/antigens; whisker = IQR") +
    theme_nat() + ggplot2::theme(legend.position = "top")
  
  pc <- ggplot2::ggplot(yc, ggplot2::aes(.mized(mean), country)) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = .mized(lo), xmax = .mized(hi)), height = 0,
                            colour = "grey70", linewidth = 0.7) +
    ggplot2::geom_point(colour = "#7F2704", size = 1.8) +
    ggplot2::scale_x_continuous(labels = scales::comma, expand = ggplot2::expansion(mult = c(0.02, 0.1))) +
    ggplot2::labs(x = "Attributable YLL (millions, undiscounted)", y = NULL,
                  subtitle = "Composite best-available estimator; 95% CI") +
    theme_nat()
  
  pd <- ggplot2::ggplot(ym, ggplot2::aes(.mized(mean), method, colour = family)) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = .mized(lo), xmax = .mized(hi)), height = 0, linewidth = 0.7) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_colour_manual(values = .fam_col, name = NULL) +
    ggplot2::scale_x_continuous(labels = scales::comma, expand = ggplot2::expansion(mult = c(0.04, 0.08))) +
    ggplot2::labs(x = "7-country total YLL (millions, undiscounted)", y = NULL,
                  subtitle = "Each estimator; 95% CI") +
    theme_nat() + ggplot2::theme(legend.position = "top")
  
  if (!requireNamespace("patchwork", quietly = TRUE)) stop("fig2_v2 needs {patchwork}.")
  (pa | pb) / (pc | pd) +
    patchwork::plot_annotation(
      tag_levels = "a",
      title = "From coverage gaps to years of life lost, by country and estimator",
      subtitle = "Top: conflict-attributable coverage gaps. Bottom: attributable YLL. Left: by country; right: by estimator.",
      theme = theme_nat()) &
    ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 8))
}

# ============================================================================
# FIGURE 3  -  triangulation
# ============================================================================
# Per-method draw-total vectors over a common country set (default: the
# complete-case set shared by all methods) so the distributions and the
# family-vs-family contour are CRN-aligned and comparable.
.f3_method_draw_totals <- function(yll_draws, catchup_focus = 0,
                                   framing_focus = if (exists("FRAMING_HEADLINE")) FRAMING_HEADLINE else "campaign_topup",
                                   common_countries = TRUE) {
  d <- yll_draws %>% dplyr::filter(catchup == catchup_focus, framing == framing_focus)
  keep <- vapply(d$d_undisc, function(x) !is.null(x) && length(x) > 0, logical(1)); d <- d[keep, ]
  if (common_countries) {
    nm <- d %>% dplyr::distinct(method, country) %>% dplyr::count(country, name = "k")
    full <- nm$country[nm$k == max(nm$k)]
    d <- d %>% dplyr::filter(country %in% full)
  }
  d %>% dplyr::group_by(method) %>%
    dplyr::group_map(~ {
      M <- do.call(rbind, .x$d_undisc); tibble::tibble(method = .y$method, draw = seq_len(ncol(M)),
                                                       total = colSums(M)) }) %>%
    dplyr::bind_rows()
}

fig3_tornado <- function(struct_sens = NULL, s5_csv = NULL) {
  s <- if (!is.null(struct_sens)) struct_sens else {
    stopifnot(!is.null(s5_csv)); utils::read.csv(s5_csv, stringsAsFactors = FALSE) }
  s <- s %>% dplyr::filter(axis != "central") %>%
    dplyr::mutate(excluded = grepl("multiplier", scenario, ignore.case = TRUE)) %>%
    dplyr::arrange(pct_change) %>%
    dplyr::mutate(scenario = factor(scenario, levels = scenario),
                  fill = dplyr::case_when(excluded ~ "excluded", pct_change > 0 ~ "up", TRUE ~ "down"))
  ggplot2::ggplot(s, ggplot2::aes(pct_change, scenario, fill = fill)) +
    ggplot2::geom_col(width = 0.66, ggplot2::aes(alpha = excluded)) +
    ggplot2::geom_vline(xintercept = 0, colour = "#1C2B36", linewidth = 0.4) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%+.0f%%", pct_change),
                                    hjust = ifelse(pct_change >= 0, -0.15, 1.12)),
                       size = 2.2, colour = "grey30") +
    ggplot2::scale_fill_manual(values = c(up = "#C16E2F", down = "#356B87", excluded = "#9aa0a6"),
                               guide = "none") +
    ggplot2::scale_alpha_manual(values = c(`TRUE` = 0.55, `FALSE` = 0.9), guide = "none") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.2)) +
    ggplot2::labs(x = "Change from headline (%)", y = NULL,
                  subtitle = "Composite headline; \u2020 measles \u00d72.24 excluded") +
    theme_nat()
}

fig3_density <- function(draw_totals) {
  pal <- if (exists("nm_palette")) nm_palette else NULL
  ord <- if (exists("METHOD_ORDER")) METHOD_ORDER else sort(unique(draw_totals$method))
  dt <- draw_totals %>% dplyr::mutate(method = factor(method, levels = rev(ord[ord %in% method])),
                                      tm = total / 1e6)
  p <- ggplot2::ggplot(dt, ggplot2::aes(x = tm, y = method, fill = method))
  if (requireNamespace("ggridges", quietly = TRUE)) {
    p <- p + ggridges::geom_density_ridges(scale = 1.7, alpha = 0.82, colour = "white",
                                           linewidth = 0.2, rel_min_height = 0.01)
  } else {
    p <- ggplot2::ggplot(dt, ggplot2::aes(x = tm, colour = method, fill = method)) +
      ggplot2::geom_density(alpha = 0.25)
  }
  p + { if (!is.null(pal)) ggplot2::scale_fill_manual(values = pal, guide = "none") } +
    { if (!is.null(pal)) ggplot2::scale_colour_manual(values = pal, guide = "none") } +
    ggplot2::scale_x_continuous(labels = scales::comma) +
    # Cap the visible axis at 15M via the coordinate system (not scale limits) so
    # the full ridge densities are still computed and merely clipped at the view
    # edge - avoids dropping draws or any oob edge-cases with geom_density_ridges.
    ggplot2::coord_cartesian(xlim = c(0, 15)) +
    ggplot2::labs(x = "Total YLL per draw (millions)", y = NULL,
                  subtitle = "Monte Carlo draw distribution by estimator (common country set)") +
    theme_nat()
}

# 2-D density contour: the two inferential families draw-by-draw (CRN-linked).
# Returns NULL (rather than erroring) when a method is absent from the common
# draw set or is too degenerate for a 2-D KDE, so the d/e panels degrade
# gracefully when ITS / SC+covariates lack usable draws.
fig3_contour <- function(draw_totals,
                         x_method = "Synthetic Control", y_method = "DiD (CS)") {
  if (!all(c(x_method, y_method) %in% unique(draw_totals$method))) return(NULL)
  w <- draw_totals %>% dplyr::filter(method %in% c(x_method, y_method)) %>%
    tidyr::pivot_wider(names_from = method, values_from = total)
  if (!all(c(x_method, y_method) %in% names(w))) return(NULL)
  w <- w %>%
    dplyr::rename(x = !!x_method, y = !!y_method) %>%
    dplyr::filter(is.finite(x), is.finite(y)) %>%
    dplyr::mutate(x = x / 1e6, y = y / 1e6)
  # A 2-D kernel density needs spread in both axes and enough points; otherwise
  # fall back gracefully so the panel is skipped, not broken.
  if (nrow(w) < 10 ||
      stats::sd(w$x, na.rm = TRUE) < .Machine$double.eps^0.5 ||
      stats::sd(w$y, na.rm = TRUE) < .Machine$double.eps^0.5) return(NULL)
  rng <- range(c(w$x, w$y), na.rm = TRUE)
  ggplot2::ggplot(w, ggplot2::aes(x, y)) +
    ggplot2::geom_density_2d_filled(contour_var = "ndensity", bins = 8, alpha = 0.9) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "white", linewidth = 0.4, linetype = "22") +
    ggplot2::scale_fill_viridis_d(option = "mako", guide = "none") +
    ggplot2::coord_equal(xlim = rng, ylim = rng) +
    ggplot2::labs(x = sprintf("%s total (millions)", x_method),
                  y = sprintf("%s total (millions)", y_method),
                  subtitle = sprintf("Joint draw density; dashed = equality (%s vs %s)",
                                     y_method, x_method)) +
    theme_nat()
}

fig3_v2 <- function(yll_draws, struct_sens = NULL, s5_csv = NULL,
                    x_method = "Synthetic Control", y_method = "DiD (CS)") {
  dt <- .f3_method_draw_totals(yll_draws)
  pa <- fig3_tornado(struct_sens, s5_csv)
  pb <- fig3_density(dt)
  # Three CRN-linked family-vs-family contours, each anchored on DiD (CS):
  #   c  DiD vs Synthetic Control      (donor lags-only)
  #   d  DiD vs SC + covariates        (donor + covariates)
  #   e  DiD vs ITS                    (design / observed-trend)
  pc <- tryCatch(fig3_contour(dt, x_method = "Synthetic Control", y_method = "DiD (CS)"),
                 error = function(e) { message("  fig3c (DiD vs SC) skipped: ",
                                               conditionMessage(e)); NULL })
  pd <- tryCatch(fig3_contour(dt, x_method = "SC + covariates", y_method = "DiD (CS)"),
                 error = function(e) { message("  fig3d (DiD vs SC+cov) skipped: ",
                                               conditionMessage(e)); NULL })
  pe <- tryCatch(fig3_contour(dt, x_method = "ITS", y_method = "DiD (CS)"),
                 error = function(e) { message("  fig3e (DiD vs ITS) skipped: ",
                                               conditionMessage(e)); NULL })
  if (!requireNamespace("patchwork", quietly = TRUE)) stop("fig3_v2 needs {patchwork}.")
  contour_row <- Filter(Negate(is.null), list(pc, pd, pe))
  if (length(contour_row) == 0) {
    # No contour could be built; fall back to the top row alone.
    (pa | pb) +
      patchwork::plot_annotation(
        tag_levels = "a",
        title = "Triangulation of the attributable-burden estimate",
        subtitle = "a, Structural sensitivity. b, Per-estimator draw distributions.",
        theme = theme_nat()) &
      ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 8))
  } else {
    crow <- patchwork::wrap_plots(contour_row, nrow = 1)
    (pa | pb) / crow + patchwork::plot_layout(heights = c(1, 1.1)) +
      patchwork::plot_annotation(
        tag_levels = "a",
        title = "Triangulation of the attributable-burden estimate",
        subtitle = paste0("a, Structural sensitivity. b, Per-estimator draw distributions. ",
                          "c-e, DiD against each comparison estimator, draw-by-draw ",
                          "(c, Synthetic Control; d, SC + covariates; e, ITS)."),
        theme = theme_nat()) &
      ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 8))
  }
}

# ============================================================================
# FIGURE 4  -  per-capita intensity, onset falsification, external validation
# ============================================================================
fig4_percap <- function(s4 = NULL, s4_csv = NULL, discount_focus = "undisc",
                        best_priority = c("Synthetic Control", "DiD (CS)")) {
  d <- if (!is.null(s4)) s4 else { stopifnot(!is.null(s4_csv))
    utils::read.csv(s4_csv, stringsAsFactors = FALSE) }
  if ("discount" %in% names(d)) d <- d %>% dplyr::filter(discount == discount_focus | is.na(discount))
  # Point = the BEST ESTIMATE per country (Synthetic Control where available,
  # else DiD), mirroring the headline estimator selection - not the cross-method
  # median. Whisker still spans the full estimator range for context.
  rng <- d %>% dplyr::group_by(country) %>%
    dplyr::summarise(lo = min(yll_per_1000_births, na.rm = TRUE),
                     hi = max(yll_per_1000_births, na.rm = TRUE), .groups = "drop")
  best <- d %>% dplyr::filter(method %in% best_priority) %>%
    dplyr::mutate(prio = match(method, best_priority)) %>%
    dplyr::group_by(country) %>% dplyr::slice_min(prio, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(country, best = yll_per_1000_births, best_method = as.character(method))
  # Defensive fallback: if a country carries neither SC nor DiD in the per-birth
  # table, anchor the point on its cross-method median so it is never dropped
  # (which would leave a whisker with no marker). DiD covers all ten countries in
  # practice, so this only guards pathological inputs.
  fb <- d %>% dplyr::group_by(country) %>%
    dplyr::summarise(med = stats::median(yll_per_1000_births, na.rm = TRUE), .groups = "drop")
  s <- rng %>% dplyr::left_join(best, by = "country") %>%
    dplyr::left_join(fb, by = "country") %>%
    dplyr::mutate(best_method = ifelse(is.na(best), "median (no SC/DiD)", best_method),
                  best = ifelse(is.na(best), med, best)) %>%
    dplyr::arrange(best) %>% dplyr::mutate(country = factor(country, levels = country))
  ggplot2::ggplot(s, ggplot2::aes(best, country)) +
    ggplot2::geom_segment(ggplot2::aes(x = lo, xend = hi, yend = country),
                          colour = "grey75", linewidth = 1.2) +
    ggplot2::geom_point(colour = "#7A5295", size = 1.9) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.08))) +
    ggplot2::labs(x = "YLL per 1,000 births", y = NULL,
                  subtitle = "Best estimate (point: SC, else DiD) and estimator range (whisker)") +
    theme_nat()
}

fig4_onset <- function(s2c = NULL, s2c_csv = NULL) {
  d <- if (!is.null(s2c)) s2c else { stopifnot(!is.null(s2c_csv))
    utils::read.csv(s2c_csv, stringsAsFactors = FALSE) }
  pal <- if (exists("nm_palette")) nm_palette else NULL
  ord <- if (exists("METHOD_ORDER")) METHOD_ORDER else sort(unique(d$method))
  d <- d %>% dplyr::mutate(m = yll_undisc_mean / 1e6,
                           method = factor(method, levels = ord[ord %in% unique(method)]))
  # Decluttered: one small panel per estimator (Nature-style faceted small
  # multiples) instead of overplotting all methods on a single axis. A shared
  # free y-scale keeps each estimator's onset response legible.
  ggplot2::ggplot(d, ggplot2::aes(start_shift, m, colour = method, group = method)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3, linetype = "22") +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::geom_point(size = 1.3) +
    ggplot2::facet_wrap(~ method, scales = "free_y") +
    { if (!is.null(pal)) ggplot2::scale_colour_manual(values = pal, guide = "none") } +
    ggplot2::scale_x_continuous(breaks = sort(unique(d$start_shift))) +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::labs(x = "Onset shift (years from coded onset)", y = "Total YLL (millions)",
                  subtitle = "Aggregate maximal at the coded onset (0); one panel per estimator") +
    theme_nat() +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 6),
                   panel.spacing = ggplot2::unit(0.4, "lines"))
}

# external validation: country x disease, median model YLL as % of GBD window,
# diverging at 100% (plausible <-> over GBD), values printed, EXCEEDS flagged.
fig4_validation <- function(s6c = NULL, s6c_csv = NULL, cap = 300) {
  d <- if (!is.null(s6c)) s6c else { stopifnot(!is.null(s6c_csv))
    utils::read.csv(s6c_csv, stringsAsFactors = FALSE) }
  agg <- d %>% dplyr::group_by(country, disease) %>%
    dplyr::summarise(pct = stats::median(attributable_pct_of_gbd, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(pct_c = pmin(pct, cap),
                  exceeds = pct > 100,
                  over_cap = pct > cap,
                  disease = factor(disease),
                  country = factor(country),
                  # Print the CAPPED value (so Syria's huge raw figures don't
                  # overflow the tile and collide with the marker), prefixing ">"
                  # where the true value exceeds the cap. Text colour keys off the
                  # capped fill value so it stays legible on the saturated end.
                  lab = ifelse(is.na(pct), "",
                               ifelse(over_cap, sprintf(">%.0f", cap), sprintf("%.0f", pct_c))),
                  lab_dark = pct_c <= 150)
  ord <- agg %>% dplyr::group_by(country) %>%
    dplyr::summarise(m = stats::median(pct_c, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(m) %>% dplyr::pull(country)
  agg$country <- factor(agg$country, levels = ord)
  ggplot2::ggplot(agg, ggplot2::aes(disease, country, fill = pct_c)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.7) +
    # Exceedance flagged by a small corner tick (top-right of the tile) rather
    # than a centred cross, so it no longer overlaps the printed value.
    ggplot2::geom_point(data = subset(agg, exceeds),
                        ggplot2::aes(x = as.numeric(disease) + 0.34,
                                     y = as.numeric(country) + 0.34),
                        shape = 21, size = 0.9, stroke = 0.3,
                        fill = "white", colour = "grey20", inherit.aes = FALSE) +
    ggplot2::geom_text(ggplot2::aes(label = lab, colour = lab_dark),
                       size = 2, fontface = "bold", show.legend = FALSE) +
    ggplot2::scale_fill_gradientn(
      colours = c("#2166AC", "#92C5DE", "#F7F7F7", "#FDDBC7", "#D6604D", "#67001F"),
      values = scales::rescale(c(0, 50, 100, 150, 220, cap)), limits = c(0, cap),
      oob = scales::squish,
      name = "Model YLL as % of GBD",
      guide = ggplot2::guide_colourbar(barheight = 0.4, barwidth = 7, title.position = "top")) +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "grey25", `FALSE` = "white")) +
    ggplot2::labs(x = NULL, y = NULL,
                  subtitle = paste0("Median across estimators; \u2022 = exceeds GBD (>100%); values capped at ", cap, "% (>", cap, " shown)")) +
    theme_nat() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "bottom",
                   axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
}

fig4_v2 <- function(s4 = NULL, s2c = NULL, s6c = NULL,
                    s4_csv = NULL, s2c_csv = NULL, s6c_csv = NULL) {
  # Onset (S2c) is produced only when RUN_START_SENSITIVITY ran; per-capita (S4)
  # and validation (S6c) only when their modules ran. Build each panel in
  # isolation and assemble whatever is available so one missing input does not
  # sink the figure.
  pa <- tryCatch(fig4_percap(s4, s4_csv),
                 error = function(e) { message("  fig4a (per-capita) skipped: ",
                                               conditionMessage(e)); NULL })
  pb <- tryCatch(fig4_onset(s2c, s2c_csv),
                 error = function(e) { message("  fig4b (onset) skipped: ",
                                               conditionMessage(e)); NULL })
  pc <- tryCatch(fig4_validation(s6c, s6c_csv),
                 error = function(e) { message("  fig4c (validation) skipped: ",
                                               conditionMessage(e)); NULL })
  panels <- Filter(Negate(is.null), list(a = pa, b = pb, c = pc))
  if (length(panels) == 0) stop("fig4_v2: no panels could be built (no S4/S2c/S6c inputs).")
  if (!requireNamespace("patchwork", quietly = TRUE)) stop("fig4_v2 needs {patchwork}.")
  # Panel b is now a faceted small-multiple (one panel per estimator); cramming
  # it into a half-width cell re-clutters it, so give a, b, c their own full-width
  # rows. The per-capita Cleveland plot (a) and heatmap (c) read fine full width.
  fig <- if (!is.null(pa) && !is.null(pb) && !is.null(pc)) {
    pa / pb / pc + patchwork::plot_layout(heights = c(1, 1.1, 1.2))
  } else {
    patchwork::wrap_plots(panels, ncol = 1)
  }
  fig +
    patchwork::plot_annotation(
      tag_levels = "a",
      title = "Intensity, onset falsification, and external validation",
      subtitle = "a, Per-capita burden. b, Onset-timing placebo (one panel per estimator). c, Attributable YLL vs GBD by country and disease.",
      theme = theme_nat()) &
    ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 8))
}

# ---- driver ----------------------------------------------------------------
build_figures_v2 <- function(coverage_long, gap_df, yll_draws,
                             struct_sens = NULL,
                             s4 = NULL, s2c = NULL, s6c = NULL,
                             s5_csv = NULL, s4_csv = NULL, s2c_csv = NULL, s6c_csv = NULL,
                             drop_metric = "trough", height_mm = NULL) {
  # Default the supplementary-table CSV inputs to the saved tables when present,
  # so the driver works whether objects are passed or only the CSVs exist.
  .tab <- function(f) if (exists("TAB_DIR")) file.path(TAB_DIR, f) else f
  if (is.null(s5_csv))  s5_csv  <- .tab("Table_S5_structured_sensitivity.csv")
  if (is.null(s4_csv))  s4_csv  <- .tab("Table_S4_yll_per_1000_births.csv")
  if (is.null(s2c_csv)) s2c_csv <- .tab("Table_S2c_onset_sensitivity_fixed_cells.csv")
  if (is.null(s6c_csv)) s6c_csv <- .tab("Table_S6c_yll_validation_by_country.csv")
  
  f1 <- tryCatch(fig1_v2(coverage_long, drop_metric = drop_metric),
                 error = function(e) { message("  Figure 1 v2 skipped: ", conditionMessage(e)); NULL })
  f2 <- tryCatch(fig2_v2(gap_df, yll_draws),
                 error = function(e) { message("  Figure 2 v2 skipped: ", conditionMessage(e)); NULL })
  f3 <- tryCatch(fig3_v2(yll_draws, struct_sens = struct_sens, s5_csv = s5_csv),
                 error = function(e) { message("  Figure 3 v2 skipped: ", conditionMessage(e)); NULL })
  f4 <- tryCatch(fig4_v2(s4 = s4, s2c = s2c, s6c = s6c,
                         s4_csv = s4_csv, s2c_csv = s2c_csv, s6c_csv = s6c_csv),
                 error = function(e) { message("  Figure 4 v2 skipped: ", conditionMessage(e)); NULL })
  if (exists("save_fig_nat")) {
    if (!is.null(f1)) save_fig_nat(f1, "Figure_1_coverage_decline_v2",    width_mm = NAT_W_DOUBLE, height_mm = 165)
    if (!is.null(f2)) save_fig_nat(f2, "Figure_2_gap_to_yll_v2",          width_mm = NAT_W_DOUBLE, height_mm = 150)
    if (!is.null(f3)) save_fig_nat(f3, "Figure_3_triangulation_v2",       width_mm = NAT_W_DOUBLE, height_mm = 185)
    if (!is.null(f4)) save_fig_nat(f4, "Figure_4_intensity_validation_v2", width_mm = NAT_W_DOUBLE, height_mm = 220)
  }
  invisible(list(fig1 = f1, fig2 = f2, fig3 = f3, fig4 = f4))
}