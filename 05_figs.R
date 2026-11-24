# 05_figs.R | in: coverage_long, gap_df, yll draws, S2c/S4/S5/S6c tables | out: Tables 1, 2, S1, headline; Figs 1-4; catch-up fig | Methods: Headline composition; figure legends

# ---- palettes, method order, framings, theme, formatters ----
nm_palette <- c(
  "Baseline"          = "#8C8C8C",
  "ITS"               = "#E69F00",   # orange
  "Synthetic Control" = "#D55E00",   # vermillion
  "SC + covariates"   = "#332288",   # indigo (Tol); separates from greys under deuteranopia
  "Augmented SC"      = "#56B4E9",   # sky blue
  "DiD (CS)"          = "#0072B2",   # blue
  # baseline-window variants (first-class methods: 1/3/5-yr)
  "Baseline (1yr)"    = "#4D4D4D",
  "Baseline (3yr)"    = "#8C8C8C",
  "Baseline (5yr)"    = "#BDBDBD"
)
# Method display order (no-control -> control-based; flat -> trend -> match ->
# differencing). Used to order Figure-4 rows and legends.
METHOD_ORDER <- c("Baseline (1yr)", "Baseline (3yr)", "Baseline (5yr)",
                  "ITS", "Synthetic Control",
                  "SC + covariates", "Augmented SC", "DiD (CS)")
nm_disease_palette <- c(          # Okabe-Ito
  "Tuberculosis" = "#0072B2",
  "Measles"      = "#D55E00",
  "Diphtheria"   = "#CC79A7",
  "Pertussis"    = "#56B4E9",
  "Tetanus"      = "#E69F00"
)

FRAMING_LEVELS   <- c("routine_recovery", "campaign_topup")
FRAMING_HEADLINE <- "campaign_topup"   # <- set to your chosen headline scenario

FRAMING_LABELS <- c(                   # display names for facets/legends
  routine_recovery = "Routine recovery",
  campaign_topup   = "Campaign top-up"
)

theme_nm <- function(base_size = 9) {
  ggplot2::theme_classic(base_size = base_size, base_family = "") +
    ggplot2::theme(
      panel.background  = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background   = ggplot2::element_rect(fill = "white", colour = NA),
      panel.grid.major  = ggplot2::element_blank(),   # no gridlines (Nature artwork rule)
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

VACCINE_LEVELS <- c("BCG", "DTP3", "MCV1")

# inferential families (Figs 2-3)
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

# ---- supplementary catch-up figure (lambda grid x framing x method) ----
make_fig_catchup <- function(yll_disease_df, discount = "undisc") {
  m_mean <- paste0("yll_", discount, "_mean")
  m_lo   <- paste0("yll_", discount, "_lo")
  m_hi   <- paste0("yll_", discount, "_hi")
  d <- yll_disease_df %>%
    dplyr::mutate(
      framing_label = factor(unname(FRAMING_LABELS[as.character(framing)]),
                             levels = unname(FRAMING_LABELS))
    )
  if (anyNA(d$framing_label))
    stop("make_fig_catchup: framing value(s) not in FRAMING_LABELS: ",
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

# ---- tables ----
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
# Supplementary figures: vector PDF with live Helvetica text (base pdf device;
# does not need cairo) + 600-dpi PNG.
save_fig <- function(plot, file, width = 11, height = 6, dpi = 600) {
  stem <- path.expand(file.path(FIG_DIR, tools::file_path_sans_ext(file)))
  grDevices::pdf(paste0(stem, ".pdf"), width = width, height = height,
                 family = "Helvetica", useDingbats = FALSE, bg = "white")
  tryCatch(print(plot), finally = grDevices::dev.off())
  message("Saved figure: ", stem, ".pdf")
  ggplot2::ggsave(paste0(stem, ".png"), plot = plot, width = width, height = height,
                  dpi = dpi, units = "in", bg = "white")
  message("Saved figure: ", stem, ".png")
}

save_tab <- function(tbl, file) {
  path <- file.path(TAB_DIR, file)
  utils::write.csv(tbl, path, row.names = FALSE)
  message("Saved table: ", path)
}

# ==== main figures (Nature Portfolio spec: 89/183 mm wide, <= 247 mm tall, Helvetica 5-7 pt) ====
NAT_W_SINGLE <- 89; NAT_W_HALF <- 120; NAT_W_DOUBLE <- 183; NAT_H_MAX <- 247  # mm
NAT_FONT     <- "Helvetica"
NAT_LW   <- function(pt) pt / 2.134          # stroke in pt -> ggplot2 linewidth
NAT_TXT  <- function(pt) pt / ggplot2::.pt   # text in pt   -> geom_text size (mm)

# Okabe-Ito colour-blind-safe palette (Wong 2011, Nat Methods), one source for
# every main figure. Baseline windows are greys (attenuated family). Checked with
# colorspace::deutan/protan/tritan: minimum pairwise CIELAB distance 18.4 in all
# four views (the Okabe-Ito purple fell to 8.7 against grey under deuteranopia,
# so SC + covariates uses Tol indigo instead).
NAT_OI <- c(black = "#000000", orange = "#E69F00", sky = "#56B4E9", green = "#009E73",
            yellow = "#F0E442", blue = "#0072B2", vermillion = "#D55E00", purple = "#CC79A7")
NAT_METHOD_COLS <- c(
  "Baseline (1yr)" = "#4D4D4D", "Baseline (3yr)" = "#8C8C8C", "Baseline (5yr)" = "#BDBDBD",
  "ITS" = "#E69F00", "DiD (CS)" = "#0072B2", "Synthetic Control" = "#D55E00",
  "SC + covariates" = "#332288", "Augmented SC" = "#56B4E9")
NAT_FAMILY_COLS <- c("Design / observed-trend" = "#0072B2", "Donor-comparison" = "#D55E00",
                     "Baseline (attenuated)" = "#8C8C8C")
NAT_LOSS_COL <- "#D55E00"    # coverage-loss shading (Fig 1a)
NAT_UP_COL   <- "#D55E00"; NAT_DOWN_COL <- "#0072B2"; NAT_MUTED_COL <- "#9E9E9E"

# Print theme: classic axes, no gridlines, no in-figure title/subtitle/caption.
theme_nat <- function(base = 7) {
  ggplot2::theme_classic(base_size = base, base_family = NAT_FONT) +
    ggplot2::theme(
      line              = ggplot2::element_line(linewidth = NAT_LW(0.5)),
      axis.line         = ggplot2::element_line(colour = "black", linewidth = NAT_LW(0.5)),
      axis.ticks        = ggplot2::element_line(colour = "black", linewidth = NAT_LW(0.5)),
      axis.ticks.length = ggplot2::unit(0.8, "mm"),
      axis.text         = ggplot2::element_text(size = base - 1, colour = "black"),
      axis.title        = ggplot2::element_text(size = base, colour = "black"),
      panel.grid        = ggplot2::element_blank(),
      panel.background  = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background   = ggplot2::element_rect(fill = "white", colour = NA),
      strip.background  = ggplot2::element_blank(),
      strip.text        = ggplot2::element_text(size = base - 1, colour = "black", hjust = 0,
                                                margin = ggplot2::margin(b = 1.5)),
      legend.background = ggplot2::element_blank(),
      legend.key        = ggplot2::element_blank(),
      legend.key.size   = ggplot2::unit(3, "mm"),
      legend.text       = ggplot2::element_text(size = base - 1),
      legend.title      = ggplot2::element_text(size = base - 1),
      legend.position   = "top",
      legend.margin     = ggplot2::margin(0, 0, 0, 0),
      legend.box.spacing = ggplot2::unit(1, "mm"),
      plot.title        = ggplot2::element_blank(),
      plot.subtitle     = ggplot2::element_blank(),
      plot.caption      = ggplot2::element_blank(),
      plot.tag          = ggplot2::element_text(face = "bold", size = 8, colour = "black"),
      plot.margin       = ggplot2::margin(2, 3, 2, 2, unit = "mm"))
}

# Save at physical mm size: vector PDF with live Helvetica text (base pdf device;
# cairo is not required), plus 600-dpi RGB PNG and LZW TIFF via ragg. Stops if
# the PDF is not written, so a missing vector file can never go unnoticed.
save_fig_nat <- function(plot, file, width_mm = NAT_W_DOUBLE, height_mm = 120,
                         dpi = 600, dir = FIG_DIR) {
  stopifnot(height_mm <= NAT_H_MAX, width_mm <= NAT_W_DOUBLE)
  stem <- path.expand(file.path(dir, tools::file_path_sans_ext(file)))
  grDevices::pdf(paste0(stem, ".pdf"), width = width_mm / 25.4, height = height_mm / 25.4,
                 family = NAT_FONT, useDingbats = FALSE, bg = "white")
  ok <- tryCatch({ print(plot); TRUE }, finally = grDevices::dev.off())
  if (!file.exists(paste0(stem, ".pdf")) || file.size(paste0(stem, ".pdf")) < 1000)
    stop("save_fig_nat: vector PDF was not written for ", file)
  ragg::agg_png(paste0(stem, ".png"), width = width_mm, height = height_mm, units = "mm",
                res = dpi, background = "white")
  tryCatch(print(plot), finally = grDevices::dev.off())
  ragg::agg_tiff(paste0(stem, ".tiff"), width = width_mm, height = height_mm, units = "mm",
                 res = dpi, background = "white", compression = "lzw")
  tryCatch(print(plot), finally = grDevices::dev.off())
  message("Saved figure: ", stem, ".{pdf,png,tiff}  (", width_mm, "x", height_mm, " mm)")
  invisible(stem)
}

.comma  <- function(x) format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE)
.mized  <- function(x) x / 1e6
.fam_of <- function(m) FIG4_FAMILY(m)
.fam_col <- NAT_FAMILY_COLS
.method_col <- NAT_METHOD_COLS
# Wrap long country names onto two lines for narrow facet strips.
.wrap_country <- function(x) sub("^(Burkina|Sri) ", "\\1\n", x)

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
  trough = "Maximum coverage drop (percentage points)",
  onset  = "On-impact coverage drop (percentage points)",
  endwin = "End-of-window coverage drop (percentage points)")

# panel a: vaccine (rows) x country (cols) small multiples; red loss ribbon.
.f1_panel_a <- function(coverage_long, conflict_df, decl, country_order,
                        year_min = 1995L, year_max = 2024L,
                        loss_col = NAT_LOSS_COL) {
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
                       fill = "grey90", inherit.aes = FALSE) +
    ggplot2::geom_vline(data = shade, ggplot2::aes(xintercept = x0),
                        colour = "grey40", linewidth = NAT_LW(0.4), linetype = "22") +
    ggplot2::geom_blank(data = blank, ggplot2::aes(x = year, y = coverage)) +
    ggplot2::geom_ribbon(data = ribbon, ggplot2::aes(x = year, ymin = ylo, ymax = yhi),
                         fill = loss_col, alpha = 0.45) +
    ggplot2::geom_line(data = cov, ggplot2::aes(x = year, y = coverage),
                       colour = "black", linewidth = NAT_LW(0.75)) +
    ggplot2::facet_grid(vaccine ~ country, scales = "free_y", labeller = ggplot2::labeller(
      country = .wrap_country)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.02)),
                                breaks = scales::pretty_breaks(3)) +
    ggplot2::scale_x_continuous(breaks = c(2000, 2020), labels = c("2000", "2020"),
                                expand = c(0.03, 0)) +
    ggplot2::labs(x = "Year", y = "Coverage (%)") +
    theme_nat() +
    ggplot2::theme(panel.spacing.x = ggplot2::unit(1.6, "mm"),
                   panel.spacing.y = ggplot2::unit(2.2, "mm"),
                   strip.text.y = ggplot2::element_text(angle = 0, size = 6, hjust = 0),
                   strip.text.x = ggplot2::element_text(size = 5.5, vjust = 0, lineheight = 0.9),
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
    ggplot2::geom_col(width = 0.66, fill = "grey45") +
    ggplot2::geom_text(ggplot2::aes(label = years), hjust = -0.3,
                       size = NAT_TXT(5.5), colour = "black", family = NAT_FONT) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.18)),
                                breaks = scales::pretty_breaks(4)) +
    ggplot2::labs(x = "Coded conflict duration (years)", y = NULL) +
    theme_nat()
}

# panel b: country x antigen decline heatmap, warm sequential ramp.
.f1_panel_b <- function(decl, country_order, metric = "trough") {
  d <- decl %>% dplyr::mutate(country = factor(country, levels = country_order),
                              vaccine = factor(vaccine, levels = VACCINE_LEVELS))
  hi_cut <- stats::quantile(d$drop_pp, 0.7, na.rm = TRUE)
  ggplot2::ggplot(d, ggplot2::aes(x = country, y = vaccine, fill = drop_pp)) +
    ggplot2::geom_tile(colour = "white", linewidth = NAT_LW(1)) +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(is.na(drop_pp), "", sprintf("%.0f", drop_pp)),
                                    colour = drop_pp > hi_cut),
                       size = NAT_TXT(5.5), family = NAT_FONT, show.legend = FALSE) +
    ggplot2::scale_fill_gradientn(
      colours = c("#FFF7EC", "#FEE8C8", "#FDBB84", "#FC8D59", "#D7301F", "#7F0000"),
      limits = c(0, NA), name = unname(.F1_METRIC_LAB[metric]),
      guide = ggplot2::guide_colourbar(barheight = ggplot2::unit(1.5, "mm"),
                                       barwidth = ggplot2::unit(32, "mm"),
                                       title.position = "top")) +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black")) +
    ggplot2::scale_y_discrete(limits = rev(VACCINE_LEVELS)) +
    ggplot2::labs(x = NULL, y = NULL) +
    theme_nat() +
    ggplot2::theme(legend.position = "bottom",
                   axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 40, hjust = 1, size = 6))
}

fig1_v2 <- function(coverage_long, conflict_df = conflict_info,
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
  patchwork::wrap_plots(pa, (pb | pc) + patchwork::plot_layout(widths = c(1.45, 1)),
                        ncol = 1, heights = c(2.2, 1)) +
    patchwork::plot_annotation(tag_levels = "a", theme = theme_nat())
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
  ord <- METHOD_ORDER
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
                               framing_focus = FRAMING_HEADLINE) {
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
                              framing_focus = FRAMING_HEADLINE) {
  mt <- aggregate_global_by_method_draws(yll_draws, catchup_focus = catchup_focus,
                                         framing_focus = framing_focus) %>%
    dplyr::transmute(method, mean = yll_undisc_mean, lo = yll_undisc_lo, hi = yll_undisc_hi)
  ord <- METHOD_ORDER
  mt %>% dplyr::mutate(family = .fam_of(method),
                       method = factor(method, levels = rev(ord[ord %in% method])))
}

fig2_v2 <- function(gap_df, yll_draws, conflict_df = conflict_info,
                    headline_method = "Synthetic Control", gap_method = "DiD (CS)") {
  gc  <- .f2_gap_by_country(gap_df, conflict_df, gap_method)
  gm  <- .f2_gap_by_method(gap_df, conflict_df)
  yc  <- .f2_yll_by_country(yll_draws)
  ym  <- .f2_yll_by_method(yll_draws)
  
  whisk <- NAT_LW(0.9); pt_sz <- 1.3
  pa <- ggplot2::ggplot(gc, ggplot2::aes(gap_mean, country)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey60", linewidth = NAT_LW(0.4), linetype = "22") +
    ggplot2::geom_segment(ggplot2::aes(x = gap_lo, xend = gap_hi, yend = country),
                          colour = "grey60", linewidth = whisk) +
    ggplot2::geom_point(colour = .fam_col[["Design / observed-trend"]], size = pt_sz) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.03, 0.06))) +
    ggplot2::labs(x = "Mean coverage gap, difference-in-differences\n(percentage points)", y = NULL) +
    theme_nat()
  
  pb <- ggplot2::ggplot(gm, ggplot2::aes(gap_med, method, colour = family)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey60", linewidth = NAT_LW(0.4), linetype = "22") +
    ggplot2::geom_errorbar(orientation = "y", ggplot2::aes(xmin = gap_lo, xmax = gap_hi), width = 0, linewidth = whisk) +
    ggplot2::geom_point(size = pt_sz) +
    ggplot2::scale_colour_manual(values = .fam_col, name = NULL) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.04, 0.06))) +
    ggplot2::labs(x = "Coverage gap pooled across countries and antigens\n(percentage points)", y = NULL) +
    theme_nat()
  
  pc <- ggplot2::ggplot(yc, ggplot2::aes(.mized(mean), country)) +
    ggplot2::geom_errorbar(orientation = "y", ggplot2::aes(xmin = .mized(lo), xmax = .mized(hi)), width = 0,
                            colour = "grey60", linewidth = whisk) +
    ggplot2::geom_point(colour = "black", size = pt_sz) +
    ggplot2::scale_x_continuous(labels = scales::comma, expand = ggplot2::expansion(mult = c(0.02, 0.06))) +
    ggplot2::labs(x = "Attributable YLLs, composite estimator (millions)", y = NULL) +
    theme_nat()
  
  pd <- ggplot2::ggplot(ym, ggplot2::aes(.mized(mean), method, colour = family)) +
    ggplot2::geom_errorbar(orientation = "y", ggplot2::aes(xmin = .mized(lo), xmax = .mized(hi)), width = 0, linewidth = whisk) +
    ggplot2::geom_point(size = pt_sz) +
    ggplot2::scale_colour_manual(values = .fam_col, name = NULL) +
    ggplot2::scale_x_continuous(labels = scales::comma, limits = c(0, NA),
                                expand = ggplot2::expansion(mult = c(0, 0.06))) +
    ggplot2::labs(x = "Seven-country total YLLs (millions)", y = NULL) +
    theme_nat()
  
  (pa | pb) / (pc | pd) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(tag_levels = "a", theme = theme_nat()) &
    ggplot2::theme(legend.position = "bottom")
}

# ============================================================================
# FIGURE 3  -  triangulation
# ============================================================================
# Per-method draw-total vectors over a common country set (default: the
# complete-case set shared by all methods) so the distributions and the
# family-vs-family contour are CRN-aligned and comparable.
.f3_method_draw_totals <- function(yll_draws, catchup_focus = 0,
                                   framing_focus = FRAMING_HEADLINE,
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

# Display labels for the tornado rows (sentence case, plain ASCII for the pdf device).
.F3_SCEN_LAB <- c(
  "discounted (r=0.03)"                       = "Discounting at 3%",
  "catch-up 50%"                              = "Catch-up 50%",
  "catch-up 100%"                             = "Catch-up 100%",
  "routine-recovery framing"                  = "Routine-recovery framing",
  "local national life table"                 = "National life tables",
  "TB untreated CFR (no Tx weight)"           = "TB: untreated CFR for all cases",
  "TB treatment coverage 30%"                 = "TB treatment coverage 30%",
  "TB treatment coverage 60%"                 = "TB treatment coverage 60%",
  "full-year incidence correction (phi = 1)"  = "Full-year incidence correction",
  "no incidence correction"                   = "No incidence correction")

fig3_tornado <- function(struct_sens = NULL, s5_csv = NULL) {
  s <- if (!is.null(struct_sens)) struct_sens else {
    stopifnot(!is.null(s5_csv)); utils::read.csv(s5_csv, stringsAsFactors = FALSE) }
  s <- s %>% dplyr::filter(axis != "central") %>%
    dplyr::mutate(excluded = grepl("multiplier", scenario, ignore.case = TRUE),
                  lab = dplyr::case_when(
                    excluded ~ "Measles CFR multiplier in conflict",
                    scenario %in% names(.F3_SCEN_LAB) ~ unname(.F3_SCEN_LAB[scenario]),
                    TRUE ~ scenario)) %>%
    dplyr::arrange(pct_change) %>%
    dplyr::mutate(lab = factor(lab, levels = lab),
                  fill = dplyr::case_when(excluded ~ "excluded", pct_change > 0 ~ "up", TRUE ~ "down"))
  ggplot2::ggplot(s, ggplot2::aes(pct_change, lab, fill = fill)) +
    ggplot2::geom_col(width = 0.66) +
    ggplot2::geom_vline(xintercept = 0, colour = "black", linewidth = NAT_LW(0.5)) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%+.0f%%", pct_change),
                                    hjust = ifelse(pct_change >= 0, -0.15, 1.15)),
                       size = NAT_TXT(5.5), colour = "black", family = NAT_FONT) +
    ggplot2::scale_fill_manual(values = c(up = NAT_UP_COL, down = NAT_DOWN_COL,
                                          excluded = NAT_MUTED_COL), guide = "none") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.22)) +
    ggplot2::labs(x = "Change from headline estimate (%)", y = NULL) +
    theme_nat()
}

fig3_density <- function(draw_totals) {
  pal <- .method_col
  ord <- METHOD_ORDER
  dt <- draw_totals %>% dplyr::mutate(method = factor(method, levels = rev(ord[ord %in% method])),
                                      tm = total / 1e6)
  # Show the full distributions: upper limit at the 99.5th percentile of all draws.
  xmax <- unname(stats::quantile(dt$tm, 0.995, na.rm = TRUE))
  p <- ggplot2::ggplot(dt, ggplot2::aes(x = tm, y = method, fill = method))
  p <- p + ggridges::geom_density_ridges(scale = 1.5, alpha = 0.85, colour = "white",
                                         linewidth = NAT_LW(0.4), rel_min_height = 0.005)
  p + ggplot2::scale_fill_manual(values = pal, guide = "none") +
    ggplot2::scale_colour_manual(values = pal, guide = "none") +
    ggplot2::scale_x_continuous(labels = scales::comma, expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = c(0.3, 1.6))) +
    ggplot2::coord_cartesian(xlim = c(0, xmax)) +
    ggplot2::labs(x = "Seven-country total YLLs per draw (millions)", y = NULL) +
    theme_nat()
}

# 2-D density contour of two estimators draw-by-draw (CRN-linked). White-to-blue
# sequential fill (lowest band white, so there is no dark background block);
# `rng` fixes a common square range across panels. Returns NULL (rather than
# erroring) when a method is absent or too degenerate for a 2-D KDE.
fig3_contour <- function(draw_totals, x_method = "Synthetic Control", y_method = "DiD (CS)",
                         rng = NULL) {
  if (!all(c(x_method, y_method) %in% unique(draw_totals$method))) return(NULL)
  w <- draw_totals %>% dplyr::filter(method %in% c(x_method, y_method)) %>%
    tidyr::pivot_wider(names_from = method, values_from = total)
  if (!all(c(x_method, y_method) %in% names(w))) return(NULL)
  w <- w %>%
    dplyr::rename(x = !!x_method, y = !!y_method) %>%
    dplyr::filter(is.finite(x), is.finite(y)) %>%
    dplyr::mutate(x = x / 1e6, y = y / 1e6)
  if (nrow(w) < 10 ||
      stats::sd(w$x, na.rm = TRUE) < .Machine$double.eps^0.5 ||
      stats::sd(w$y, na.rm = TRUE) < .Machine$double.eps^0.5) return(NULL)
  if (is.null(rng)) rng <- range(c(w$x, w$y), na.rm = TRUE)
  fills <- c("#FFFFFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#084594")
  ggplot2::ggplot(w, ggplot2::aes(x, y)) +
    ggplot2::geom_density_2d_filled(contour_var = "ndensity", bins = 8) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "black",
                         linewidth = NAT_LW(0.5), linetype = "22") +
    ggplot2::scale_fill_manual(values = fills, guide = "none") +
    ggplot2::scale_x_continuous(expand = c(0, 0)) + ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::coord_equal(xlim = rng, ylim = rng) +
    ggplot2::labs(x = sprintf("%s (million YLLs)", x_method),
                  y = sprintf("%s (million YLLs)", y_method)) +
    theme_nat()
}

fig3_v2 <- function(yll_draws, struct_sens = NULL, s5_csv = NULL,
                    x_method = "Synthetic Control", y_method = "DiD (CS)") {
  dt <- .f3_method_draw_totals(yll_draws)
  pa <- fig3_tornado(struct_sens, s5_csv)
  pb <- fig3_density(dt)
  # Three CRN-linked contours, each with DiD (CS) on the y-axis, on one common
  # square range (0.5th-99.5th percentile of the four estimators involved):
  #   c  DiD vs Synthetic Control;  d  DiD vs SC + covariates;  e  DiD vs ITS
  xs  <- c("Synthetic Control", "SC + covariates", "ITS")
  rng <- unname(stats::quantile(dt$total[dt$method %in% c(xs, "DiD (CS)")] / 1e6,
                                c(0.005, 0.995), na.rm = TRUE))
  crow <- lapply(xs, function(xm) fig3_contour(dt, x_method = xm, y_method = "DiD (CS)", rng = rng))
  crow <- Filter(Negate(is.null), crow)
  # Rows are assembled independently (wrap_elements) so the three square contour
  # panels centre across the full width instead of inheriting panel a's wide
  # label margin; panel letters are therefore set explicitly.
  tag <- function(p, l) p + ggplot2::labs(tag = l)
  top <- patchwork::wrap_elements(full = tag(pa, "a") | tag(pb, "b"))
  if (length(crow) == 0) return(top)
  crow <- Map(tag, crow, letters[2 + seq_along(crow)])
  bottom <- patchwork::wrap_elements(full = patchwork::wrap_plots(crow, nrow = 1))
  top / bottom + patchwork::plot_layout(heights = c(1.2, 1))
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
                          colour = "grey60", linewidth = NAT_LW(0.9)) +
    ggplot2::geom_point(colour = "black", size = 1.3) +
    ggplot2::scale_x_continuous(limits = c(0, NA), expand = ggplot2::expansion(mult = c(0, 0.06))) +
    ggplot2::labs(x = "YLLs per 1,000 births", y = NULL) +
    theme_nat()
}

fig4_onset <- function(s2c = NULL, s2c_csv = NULL) {
  d <- if (!is.null(s2c)) s2c else { stopifnot(!is.null(s2c_csv))
    utils::read.csv(s2c_csv, stringsAsFactors = FALSE) }
  pal <- .method_col
  ord <- METHOD_ORDER
  d <- d %>% dplyr::mutate(m = yll_undisc_mean / 1e6,
                           method = factor(method, levels = ord[ord %in% unique(method)]))
  # Decluttered: one small panel per estimator (Nature-style faceted small
  # multiples) instead of overplotting all methods on a single axis. A shared
  # free y-scale keeps each estimator's onset response legible.
  ggplot2::ggplot(d, ggplot2::aes(start_shift, m, colour = method, group = method)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey50", linewidth = NAT_LW(0.4), linetype = "22") +
    ggplot2::geom_line(linewidth = NAT_LW(0.75)) +
    ggplot2::geom_point(size = 0.9) +
    ggplot2::facet_wrap(~ method, scales = "free_y", nrow = 2) +
    ggplot2::scale_colour_manual(values = pal, guide = "none") +
    ggplot2::scale_x_continuous(breaks = sort(unique(d$start_shift)),
                                labels = function(x) ifelse(x > 0, paste0("+", x), as.character(x))) +
    ggplot2::scale_y_continuous(labels = scales::comma, limits = c(0, NA),
                                expand = ggplot2::expansion(mult = c(0, 0.08))) +
    ggplot2::labs(x = "Shift in coded conflict onset (years)", y = "Total YLLs (millions)") +
    theme_nat() +
    ggplot2::theme(panel.spacing = ggplot2::unit(2.5, "mm"))
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
    ggplot2::geom_tile(colour = "white", linewidth = NAT_LW(1)) +
    # Exceedance flagged by a small corner marker (top-right of the tile) so it
    # does not overlap the printed value.
    ggplot2::geom_point(data = subset(agg, exceeds),
                        ggplot2::aes(x = as.numeric(disease) + 0.36,
                                     y = as.numeric(country) + 0.30),
                        shape = 21, size = 0.7, stroke = 0.3,
                        fill = "white", colour = "black", inherit.aes = FALSE) +
    ggplot2::geom_text(ggplot2::aes(label = lab, colour = lab_dark),
                       size = NAT_TXT(5.5), family = NAT_FONT, show.legend = FALSE) +
    ggplot2::scale_fill_gradientn(
      colours = c("#2166AC", "#92C5DE", "#F7F7F7", "#FDDBC7", "#D6604D", "#67001F"),
      values = scales::rescale(c(0, 50, 100, 150, 220, cap)), limits = c(0, cap),
      oob = scales::squish,
      name = "Model-attributed YLLs as % of GBD total YLLs",
      guide = ggplot2::guide_colourbar(barheight = ggplot2::unit(1.5, "mm"),
                                       barwidth = ggplot2::unit(32, "mm"),
                                       title.position = "top")) +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "white")) +
    ggplot2::labs(x = NULL, y = NULL) +
    theme_nat() +
    ggplot2::theme(legend.position = "bottom",
                   axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(size = 5.5))
}

fig4_v2 <- function(s4 = NULL, s2c = NULL, s6c = NULL, s4_csv = NULL, s2c_csv = NULL, s6c_csv = NULL) {
  # a per-birth burden | b validation heatmap (top row); c onset small multiples (full width)
  pa <- fig4_percap(s4, s4_csv) + ggplot2::labs(tag = "a")
  pb <- fig4_validation(s6c, s6c_csv) + ggplot2::labs(tag = "b")
  pc <- fig4_onset(s2c, s2c_csv) + ggplot2::labs(tag = "c")
  (((pa | pb) + patchwork::plot_layout(widths = c(1, 1.15))) / pc +
     patchwork::plot_layout(heights = c(1.2, 1))) +
    patchwork::plot_annotation(theme = theme_nat())
}

# ---- driver ----------------------------------------------------------------
build_figures_v2 <- function(coverage_long, gap_df, yll_draws,
                             struct_sens = NULL,
                             s4 = NULL, s2c = NULL, s6c = NULL,
                             s5_csv = NULL, s4_csv = NULL, s2c_csv = NULL, s6c_csv = NULL,
                             drop_metric = "trough", height_mm = NULL) {
  # Default the supplementary-table CSV inputs to the saved tables when present,
  # so the driver works whether objects are passed or only the CSVs exist.
  .tab <- function(f) file.path(TAB_DIR, f)
  if (is.null(s5_csv))  s5_csv  <- .tab("Table_S5_structured_sensitivity.csv")
  if (is.null(s4_csv))  s4_csv  <- .tab("Table_S4_yll_per_1000_births.csv")
  if (is.null(s2c_csv)) s2c_csv <- .tab("Table_S2c_onset_sensitivity_fixed_cells.csv")
  if (is.null(s6c_csv)) s6c_csv <- .tab("Table_S6c_yll_validation_by_country.csv")
  
  f1 <- fig1_v2(coverage_long, drop_metric = drop_metric)
  f2 <- fig2_v2(gap_df, yll_draws)
  f3 <- fig3_v2(yll_draws, struct_sens = struct_sens, s5_csv = s5_csv)
  f4 <- fig4_v2(s4 = s4, s2c = s2c, s6c = s6c, s4_csv = s4_csv, s2c_csv = s2c_csv, s6c_csv = s6c_csv)
  save_fig_nat(f1, "Figure_1_coverage_decline_v2",     width_mm = NAT_W_DOUBLE, height_mm = 150)
  save_fig_nat(f2, "Figure_2_gap_to_yll_v2",           width_mm = NAT_W_DOUBLE, height_mm = 115)
  save_fig_nat(f3, "Figure_3_triangulation_v2",        width_mm = NAT_W_DOUBLE, height_mm = 140)
  save_fig_nat(f4, "Figure_4_intensity_validation_v2", width_mm = NAT_W_DOUBLE, height_mm = 165)
  invisible(list(fig1 = f1, fig2 = f2, fig3 = f3, fig4 = f4))
}
