
stopifnot(requireNamespace("ggplot2", quietly = TRUE))
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
})

# ---- Inferential families (key the palette + grouping off these) ----------
FIG4_DESIGN_METHODS <- c("ITS", "DiD (CS)")
FIG4_DONOR_METHODS  <- c("Synthetic Control", "SC + covariates", "Augmented SC")
# Baseline (1yr) is a pre-conflict-anchored baseline, not a design/observed-trend
# estimator: it carries no observed-trend extrapolation. Group it with the other
# baselines (attenuated family) so the palette/grouping reflect its construction.
FIG4_ATTEN_METHODS  <- c("Baseline (1yr)", "Baseline (3yr)", "Baseline (5yr)")

FIG4_FAMILY <- function(m) dplyr::case_when(
  m %in% FIG4_DESIGN_METHODS ~ "Design / observed-trend",
  m %in% FIG4_DONOR_METHODS  ~ "Donor-comparison",
  TRUE                       ~ "Baseline (attenuated)")

FIG4_COL <- c("Design / observed-trend" = "#356B87",
              "Donor-comparison"        = "#C16E2F",
              "Baseline (attenuated)"   = "#BFB6A8")
FIG4_HEAD   <- "#1C2B36"
FIG4_SHORT  <- c("Baseline (1yr)" = "Baseline 1yr", "Baseline (3yr)" = "Baseline 3yr",
                 "Baseline (5yr)" = "Baseline 5yr", "ITS" = "ITS", "DiD (CS)" = "DiD",
                 "Synthetic Control" = "SC", "SC + covariates" = "SC+cov",
                 "Augmented SC" = "ASCM")

.fig4_theme <- function(base = 9) {
  th <- if (exists("theme_nm", mode = "function")) theme_nm(base) else
    ggplot2::theme_classic(base_size = base)
  th + ggplot2::theme(
    panel.grid.major.y = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(face = "bold", size = base + 3, hjust = 0),
    legend.position = "none")
}

# ---------------------------------------------------------------------------
# Data assembly
# ---------------------------------------------------------------------------
# Draw-level method totals (complete-case across methods). Reuses 03b's
# aggregate_global_by_method_draws() when yll_draws is supplied.
fig4_method_totals <- function(yll_draws = NULL, csv = NULL,
                               catchup_focus = 0, framing_focus = "campaign_topup") {
  if (!is.null(yll_draws) && exists("aggregate_global_by_method_draws")) {
    mt <- aggregate_global_by_method_draws(yll_draws, catchup_focus = catchup_focus,
                                           framing_focus = framing_focus) %>%
      dplyr::transmute(method, mean = yll_undisc_mean,
                       lo = yll_undisc_lo, hi = yll_undisc_hi)
    if (nrow(mt) > 0) return(mt)
  }
  # CSV fallback: sum disease-level point estimates from Table 2 (no valid total CI).
  stopifnot(!is.null(csv))
  t2 <- utils::read.csv(csv, check.names = FALSE, stringsAsFactors = FALSE)
  num <- function(s) as.numeric(gsub(",", "", sub(" .*", "", s)))
  t2$u <- num(t2$YLLs_undiscounted)
  t2 %>% dplyr::group_by(method = Method) %>%
    dplyr::summarise(mean = sum(u), lo = NA_real_, hi = NA_real_, .groups = "drop")
}

# Per disease x method totals with CIs. Prefers the SAME complete-case draw-level
# aggregator Table 2 uses (so panel b matches Table 2 exactly); else a raw
# per-method draw sum; else parse Table 2 CSV.
fig4_disease_method <- function(yll_draws = NULL, csv = NULL,
                                catchup_focus = 0, framing_focus = "campaign_topup",
                                credible_only = TRUE) {
  dm <-
    if (!is.null(yll_draws) && exists("aggregate_method_comparison_draws")) {
      aggregate_method_comparison_draws(yll_draws, catchup_focus = catchup_focus,
                                        framing_focus = framing_focus) %>%
        dplyr::transmute(disease, method, mean = yll_undisc_mean,
                         lo = yll_undisc_lo, hi = yll_undisc_hi)
    } else if (!is.null(yll_draws)) {
      d <- yll_draws %>% dplyr::filter(catchup == catchup_focus, framing == framing_focus)
      keep <- vapply(d$d_undisc, function(x) !is.null(x) && length(x) > 0, logical(1))
      d <- d[keep, , drop = FALSE]
      d %>% dplyr::group_by(disease, method) %>%
        dplyr::group_modify(~ {
          M <- do.call(rbind, .x$d_undisc); s <- colSums(M)
          tibble::tibble(mean = mean(s),
                         lo = stats::quantile(s, .025, names = FALSE),
                         hi = stats::quantile(s, .975, names = FALSE))
        }) %>% dplyr::ungroup()
    } else {
      stopifnot(!is.null(csv))
      t2 <- utils::read.csv(csv, check.names = FALSE, stringsAsFactors = FALSE)
      ci <- function(s) {
        m <- regmatches(s, gregexpr("[0-9,]+", s))
        vapply(m, function(v) as.numeric(gsub(",", "", v)), numeric(3))
      }
      mat <- ci(t2$YLLs_undiscounted)
      tibble::tibble(disease = t2$Disease, method = t2$Method,
                     mean = mat[2, ], lo = mat[1, ], hi = mat[3, ])
    }
  if (credible_only) dm <- dm %>% dplyr::filter(!method %in% FIG4_ATTEN_METHODS)
  dm
}

fig4_sensitivity <- function(struct_sens = NULL, csv = NULL) {
  s <- if (!is.null(struct_sens)) struct_sens else {
    stopifnot(!is.null(csv)); utils::read.csv(csv, stringsAsFactors = FALSE)
  }
  s %>% dplyr::filter(axis != "central") %>%
    dplyr::transmute(scenario, pct_change,
                     excluded = grepl("multiplier", scenario, ignore.case = TRUE)) %>%
    dplyr::arrange(pct_change)
}

# ---------------------------------------------------------------------------
# Panels
# ---------------------------------------------------------------------------
fig4_panel_a <- function(mt, headline = "Synthetic Control") {
  mt <- mt %>% dplyr::mutate(family = FIG4_FAMILY(method),
                             lab = FIG4_SHORT[method],
                             is_head = method == headline) %>%
    dplyr::arrange(mean) %>% dplyr::mutate(rank = dplyr::row_number())
  band <- mt %>% dplyr::filter(method %in% FIG4_DESIGN_METHODS)
  has_band <- nrow(band) > 0
  ggplot2::ggplot(mt, ggplot2::aes(mean / 1e6, rank)) +
    { if (has_band)
      ggplot2::annotate("rect", xmin = min(band$mean) / 1e6, xmax = max(band$mean) / 1e6,
                        ymin = -Inf, ymax = Inf, fill = "grey50", alpha = 0.08) } +
    { if (all(is.finite(mt$lo)))
      ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo / 1e6, xmax = hi / 1e6,
                                           colour = family), height = 0, linewidth = 0.5,
                              alpha = 0.6) } +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = mean / 1e6, yend = rank,
                                       colour = family), linewidth = 0.5, alpha = 0.35) +
    ggplot2::geom_point(ggplot2::aes(colour = family, size = is_head)) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", mean / 1e6)),
                       hjust = -0.35, size = 2.6, colour = "grey35") +
    ggplot2::scale_y_continuous(breaks = mt$rank, labels = mt$lab, expand = c(0.05, 0.4)) +
    ggplot2::scale_x_continuous(limits = c(0, max(mt$mean / 1e6) * 1.18)) +
    ggplot2::scale_colour_manual(values = FIG4_COL) +
    ggplot2::scale_size_manual(values = c(`TRUE` = 3.4, `FALSE` = 2.2), guide = "none") +
    ggplot2::labs(title = "a", x = "7-country total YLL (millions, undiscounted)", y = NULL) +
    .fig4_theme()
}

fig4_panel_b <- function(dm, disease_order = NULL) {
  if (is.null(disease_order))
    disease_order <- dm %>% dplyr::group_by(disease) %>%
      dplyr::summarise(m = max(mean)) %>% dplyr::arrange(m) %>% dplyr::pull(disease)
  # Factor levels must span EVERY method present (not just design+donor); after
  # Baseline (1yr) moved to the attenuated family, restricting levels to
  # design+donor would coerce baseline rows to NA and silently drop them. Use the
  # full method universe so ordering is stable and no row is lost.
  .lvls <- if (exists("METHOD_ORDER"))
    METHOD_ORDER else c(FIG4_DESIGN_METHODS, FIG4_DONOR_METHODS, FIG4_ATTEN_METHODS)
  dm <- dm %>%
    dplyr::mutate(disease = factor(disease, levels = disease_order),
                  family = FIG4_FAMILY(method),
                  method = factor(method, levels = .lvls[.lvls %in% unique(method)]))
  # consensus band per disease = span of design-based point estimates. geom_rect
  # respects the log10 x-transform (geom_tile's linear `width` would not); y via
  # the factor's integer position so the band hugs each disease row.
  band <- dm %>% dplyr::filter(method %in% FIG4_DESIGN_METHODS) %>%
    dplyr::group_by(disease) %>%
    dplyr::summarise(xlo = min(mean), xhi = max(mean), .groups = "drop") %>%
    dplyr::mutate(yc = as.numeric(disease))
  ggplot2::ggplot(dm, ggplot2::aes(mean, disease, colour = family, group = method)) +
    ggplot2::geom_rect(data = band, inherit.aes = FALSE,
                       ggplot2::aes(xmin = xlo, xmax = xhi,
                                    ymin = yc - 0.39, ymax = yc + 0.39),
                       fill = FIG4_COL[["Design / observed-trend"]], alpha = 0.08) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = pmax(lo, 1), xmax = pmax(hi, 1)),
                            position = ggplot2::position_dodge(width = 0.72),
                            height = 0, linewidth = 0.7, alpha = 0.85) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.72), size = 1.8) +
    ggplot2::scale_x_log10(labels = function(x) format(x, big.mark = ",", scientific = FALSE),
                           expand = c(0.04, 0)) +
    ggplot2::scale_colour_manual(values = FIG4_COL, name = NULL) +
    ggplot2::labs(title = "b", x = "7-country total YLL (undiscounted, log scale)", y = NULL) +
    .fig4_theme() +
    ggplot2::theme(legend.position = c(0.18, 0.93),
                   legend.text = ggplot2::element_text(size = 7))
}

fig4_panel_c <- function(sens) {
  sens <- sens %>%
    dplyr::mutate(scenario = factor(scenario, levels = scenario),
                  fill = dplyr::case_when(excluded ~ "excluded",
                                          pct_change > 0 ~ "up", TRUE ~ "down"))
  ggplot2::ggplot(sens, ggplot2::aes(pct_change, scenario, fill = fill)) +
    ggplot2::geom_col(width = 0.66, ggplot2::aes(alpha = excluded)) +
    ggplot2::geom_vline(xintercept = 0, colour = FIG4_HEAD, linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%+.0f%%", pct_change),
                                    hjust = ifelse(pct_change >= 0, -0.15, 1.15)),
                       size = 2.5, colour = "grey30") +
    ggplot2::scale_fill_manual(values = c(up = "#C16E2F", down = "#356B87",
                                          excluded = "#9aa0a6"), guide = "none") +
    ggplot2::scale_alpha_manual(values = c(`TRUE` = 0.55, `FALSE` = 0.9), guide = "none") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.18)) +
    ggplot2::labs(title = "c", x = "Change from headline (%)", y = NULL,
                  caption = "\u2020 measles \u00d72.24 excluded from headline") +
    .fig4_theme()
}

# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
build_fig4_triangulation <- function(yll_draws = NULL, struct_sens = NULL,
                                     table2_csv = NULL, s5_csv = NULL,
                                     headline = "Synthetic Control",
                                     out_pdf = "Fig4_triangulation.pdf",
                                     width = 7.2, height = 8.6) {
  mt   <- fig4_method_totals(yll_draws, csv = table2_csv)
  dm   <- fig4_disease_method(yll_draws, csv = table2_csv)
  sens <- fig4_sensitivity(struct_sens, csv = s5_csv)
  pa <- fig4_panel_a(mt, headline); pb <- fig4_panel_b(dm); pc <- fig4_panel_c(sens)
  
  if (requireNamespace("patchwork", quietly = TRUE)) {
    library(patchwork)
    fig <- (pa | pc) / pb + patchwork::plot_layout(heights = c(1, 1.55))
  } else {
    if (requireNamespace("gridExtra", quietly = TRUE)) {
      top <- gridExtra::arrangeGrob(pa, pc, ncol = 2, widths = c(1.18, 1))
      fig <- gridExtra::arrangeGrob(top, pb, ncol = 1, heights = c(1, 1.55))
    } else stop("install 'patchwork' or 'gridExtra' to assemble Fig 4 panels.")
  }
  ggplot2::ggsave(out_pdf, fig, width = width, height = height, device = grDevices::cairo_pdf)
  ggplot2::ggsave(sub("\\.pdf$", ".png", out_pdf), fig, width = width, height = height, dpi = 300)
  message("Fig 4 written to ", out_pdf)
  invisible(list(panel_a = pa, panel_b = pb, panel_c = pc, mt = mt, dm = dm, sens = sens))
}

# ---- Example -------------------------------------------------------------
if (FALSE) {
  # From the pipeline (valid draw-based CIs):
  build_fig4_triangulation(yll_draws = yll_draws, struct_sens = structured_sensitivity)
  # From saved tables (CSV fallback):
  build_fig4_triangulation(table2_csv = "Table_2_method_comparison.csv",
                           s5_csv = "Table_S5_structured_sensitivity.csv")
}
