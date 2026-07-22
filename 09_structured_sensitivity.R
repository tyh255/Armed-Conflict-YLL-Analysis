# Reference ex0 (headline life table) and central TB effective CFR.
.REF_EX0 <- if (exists("ref_ex_at_age")) ref_ex_at_age(0) else 88.8718951
.CFR_TB_UNTREATED <- 0.436
.CFR_TB_TREATED   <- 0.019
.cfr_eff_tb <- function(p_treat) p_treat * .CFR_TB_TREATED + (1 - p_treat) * .CFR_TB_UNTREATED
.CFR_TB_CENTRAL <- .cfr_eff_tb(0.45)
# Central conflict measles-CFR multiplier, pulled from the params table
# (Measles_conflict_mult; percentile-matched lognormal central). The realised
# draw MEAN is modestly higher (~2.43) than this median, so using the central
# makes the "measles multiplier ON" tornado bar a mild LOWER bound on that swing.
.MEASLES_MULT_MEAN <- {
  v <- if (exists("params"))
    params$mean[params$parameter == "Measles_conflict_mult"][1] else NA_real_
  if (length(v) == 1 && is.finite(v)) v else 2.24
}

# Per-country life-table rescale factor (local national LE over the conflict
# window vs the reference ex0). Approximates switching LIFE_TABLE_MODE to a
# national life table; < 1 for these countries (national LE ~ 60-65 << 88.87).
.country_le_factor <- function(covariates = covariates_panel,
                               conflict_df = conflict_info) {
  conflict_df %>%
    dplyr::select(country, ISO3, conflict_start, conflict_end) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(le_mean = mean(covariates$Life_Expectancy[
      covariates$ISO3 == ISO3 &
        covariates$year >= conflict_start &
        covariates$year <= conflict_end], na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(factor = ifelse(is.finite(le_mean), le_mean / .REF_EX0, NA_real_)) %>%
    dplyr::select(country, le_factor = factor)
}

# Per-(disease,country) 
.composite_cmap <- function(draws, anchor_priority, catchup_focus = 0,
                            framing_focus = if (exists("FRAMING_HEADLINE")) FRAMING_HEADLINE else "campaign_topup") {
  d <- draws %>% dplyr::filter(catchup == catchup_focus, framing == framing_focus)
  keep <- vapply(d$d_undisc, function(x) !is.null(x) && length(x) > 0, logical(1))
  d <- d[keep, , drop = FALSE]
  if (nrow(d) == 0) return(NULL)
  d %>% dplyr::filter(method %in% anchor_priority) %>%
    dplyr::mutate(prio = match(method, anchor_priority)) %>%
    dplyr::group_by(disease, country) %>%
    dplyr::slice_min(prio, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(disease, country, method)
}

# Core: total over a CELL SET
.tornado_total <- function(draws, anchor_method = NULL, cmap = NULL,
                           disease_factor = function(x) 1,
                           country_factor = NULL,
                           discount = "undisc", catchup_focus = 0,
                           framing_focus = if (exists("FRAMING_HEADLINE")) FRAMING_HEADLINE else "campaign_topup") {
  col <- if (discount == "undisc") "d_undisc" else "d_disc"
  d <- draws %>% dplyr::filter(catchup == catchup_focus, framing == framing_focus)
  d <- if (!is.null(cmap))
    dplyr::inner_join(d, cmap, by = c("disease", "country", "method"))
  else
    dplyr::filter(d, method == anchor_method)
  keep <- vapply(d[[col]], function(x) !is.null(x) && length(x) > 0, logical(1))
  d <- d[keep, , drop = FALSE]
  if (nrow(d) == 0) return(NULL)
  S <- NULL
  for (i in seq_len(nrow(d))) {
    f <- disease_factor(d$disease[i])
    if (!is.null(country_factor)) {
      cf <- country_factor[[d$country[i]]]
      if (is.null(cf) || is.na(cf)) cf <- 1
      f <- f * cf
    }
    v <- f * d[[col]][[i]]
    S <- if (is.null(S)) v else S + v
  }
  c(mean = mean(S),
    lo = stats::quantile(S, 0.025, names = FALSE),
    hi = stats::quantile(S, 0.975, names = FALSE),
    n_cells = nrow(d))
}

run_structured_sensitivity <- function(yll_draws,
                                       anchor_mode = c("composite", "single"),
                                       anchor_method = "Synthetic Control",
                                       anchor_priority = c("Synthetic Control",
                                                           "DiD (CS)", "Baseline (3yr)"),
                                       also_single = TRUE,
                                       out_table = "Table_S5_structured_sensitivity.csv",
                                       out_fig = "Figure_S_structured_sensitivity",
                                       covariates = covariates_panel,
                                       conflict_df = conflict_info) {
  anchor_mode <- match.arg(anchor_mode)
  if (is.null(yll_draws) || nrow(yll_draws) == 0) {
    message("  Tornado: empty draws -> skipping."); return(invisible(NULL))
  }
  le_fac <- .country_le_factor(covariates, conflict_df)
  le_map <- stats::setNames(le_fac$le_factor, le_fac$country)
  
  # ---- choose the cell set: composite headline (default) or single method -----
  cmap <- NULL
  if (anchor_mode == "composite") {
    cmap <- .composite_cmap(yll_draws, anchor_priority)
    if (is.null(cmap) || nrow(cmap) == 0) {
      message("  Tornado: composite cell-map empty; falling back to single anchor '",
              anchor_method, "'."); anchor_mode <- "single"; cmap <- NULL
    }
  }
  if (anchor_mode == "single") {
    if (!anchor_method %in% unique(yll_draws$method)) {
      anchor_method <- unique(yll_draws$method)[1]
      message("  Tornado: anchor method not present; using '", anchor_method, "'.")
    }
    anchor_label <- anchor_method
  } else {
    mix <- cmap %>% dplyr::count(method, name = "n_cells") %>%
      dplyr::arrange(dplyr::desc(n_cells))
    anchor_label <- paste0("Composite (", paste(anchor_priority, collapse = ">"), ")")
    message("  Tornado: composite anchor cell-mix -> ",
            paste(sprintf("%s=%d", mix$method, mix$n_cells), collapse = ", "),
            " (", sum(mix$n_cells), " cells, ",
            dplyr::n_distinct(cmap$country), " countries). This anchors the ",
            "tornado on the SAME total the headline reports.")
  }
  
  # one totalizer bound to the chosen cell set (cmap NULL => single-method path)
  tt <- function(...) .tornado_total(yll_draws, anchor_method = anchor_method,
                                     cmap = cmap, ...)
  
  central <- tt()
  if (is.null(central)) { message("  Tornado: no central cells."); return(invisible(NULL)) }
  central_mean <- central["mean"]
  
  mk <- function(label, axis, res, approx = FALSE) {
    if (is.null(res)) return(NULL)
    tibble::tibble(axis = axis, scenario = label,
                   yll = res["mean"], lo = res["lo"], hi = res["hi"],
                   pct_change = 100 * (res["mean"] - central_mean) / central_mean,
                   approx = approx)
  }
  
  rows <- list(
    mk("central (headline)", "central", central, FALSE),
    # exact-from-draws axes
    mk("discounted (r=0.03)", "discounting", tt(discount = "disc"), FALSE),
    mk("catch-up 50%", "catch_up", tt(catchup_focus = 0.5), FALSE),
    mk("catch-up 100%", "catch_up", tt(catchup_focus = 1.0), FALSE),
    mk("routine-recovery framing", "framing", tt(framing_focus = "routine_recovery"), FALSE),
    # rescale-approx axes
    mk("local national life table", "life_table", tt(country_factor = le_map), TRUE),
    mk(sprintf("measles multiplier ON (x%.2f)", .MEASLES_MULT_MEAN), "measles_mult",
       tt(disease_factor = function(x) if (x == "Measles") .MEASLES_MULT_MEAN else 1), TRUE),
    mk("TB untreated CFR (no Tx weight)", "tb_cfr",
       tt(disease_factor = function(x) if (x == "Tuberculosis") .CFR_TB_UNTREATED / .CFR_TB_CENTRAL else 1), TRUE),
    mk("TB treatment coverage 30%", "tb_cfr",
       tt(disease_factor = function(x) if (x == "Tuberculosis") .cfr_eff_tb(0.30) / .CFR_TB_CENTRAL else 1), TRUE),
    mk("TB treatment coverage 60%", "tb_cfr",
       tt(disease_factor = function(x) if (x == "Tuberculosis") .cfr_eff_tb(0.60) / .CFR_TB_CENTRAL else 1), TRUE)
  )
  tab <- dplyr::bind_rows(rows) %>% dplyr::mutate(anchor_method = anchor_label)
  
  # SELF-TEST: catch-up and discounting must REDUCE the total (monotone, signed)
  .chk <- function(scn, want) {
    r <- tab$pct_change[startsWith(tab$scenario, scn)]
    if (length(r) >= 1 && is.finite(r[1]) &&
        ((want == "down" && r[1] > 1) || (want == "up" && r[1] < -1)))
      message("  [tornado diag] unexpected sign for '", scn, "': ",
              sprintf("%+.1f%%", r[1]), " (expected ", want, ").")
  }
  .chk("catch-up 50%", "down"); .chk("catch-up 100%", "down")
  .chk("discounted", "down"); .chk("local national life table", "down")
  .chk("measles multiplier ON", "up"); .chk("TB untreated CFR", "up")
  
  if (exists("save_tab")) save_tab(tab, out_table)
  fig <- tryCatch(.fig_tornado(tab, central_mean),
                  error = function(e) { message("  tornado fig skipped: ",
                                                conditionMessage(e)); NULL })
  if (!is.null(fig) && exists("save_fig"))
    save_fig(fig, out_fig, width = 9, height = 6)
  message("  Tornado anchor = ", anchor_label,
          "; central undiscounted total = ", format(round(central_mean), big.mark = ","),
          ". NOTE: ITS-horizon axis needs a re-run (CVD_ITS_MAX_HORIZON); see merge_rerun().")
  
  # ---- companion: SC single-method tornado (the old 7-country anchor)
  if (anchor_mode == "composite" && isTRUE(also_single) &&
      anchor_method %in% unique(yll_draws$method)) {
    comp_table <- sub("\\.csv$", "_SCanchor.csv", out_table)        # always distinct
    comp_table <- sub("Table_S5_", "Table_S5b_", comp_table, fixed = TRUE)  # S5 -> S5b
    invisible(tryCatch(
      run_structured_sensitivity(
        yll_draws, anchor_mode = "single", anchor_method = anchor_method,
        anchor_priority = anchor_priority, also_single = FALSE,
        out_table = comp_table,
        out_fig = paste0(out_fig, "_SCanchor"),
        covariates = covariates, conflict_df = conflict_df),
      error = function(e) { message("  S5b SC-anchor companion skipped: ",
                                    conditionMessage(e)); NULL }))
  }
  invisible(tab)
}

# Helper to fold a full-re-run total (e.g. CVD_ITS_MAX_HORIZON=Inf) into the
# tornado table after the fact, keeping the rescale-approx axes honest.
merge_rerun <- function(tornado_tab, scenario, axis, yll_mean, lo = NA, hi = NA) {
  central_mean <- tornado_tab$yll[tornado_tab$scenario == "central (headline)"]
  dplyr::bind_rows(tornado_tab, tibble::tibble(
    axis = axis, scenario = scenario, yll = yll_mean, lo = lo, hi = hi,
    pct_change = 100 * (yll_mean - central_mean) / central_mean,
    approx = FALSE, anchor_method = tornado_tab$anchor_method[1]))
}

.fig_tornado <- function(tab, central_mean) {
  d <- tab %>% dplyr::filter(scenario != "central (headline)") %>%
    dplyr::arrange(pct_change) %>%
    dplyr::mutate(scenario = factor(scenario, levels = scenario))
  ggplot2::ggplot(d, ggplot2::aes(x = pct_change, y = scenario, fill = approx)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey40") +
    ggplot2::scale_fill_manual(values = c("FALSE" = "#2C6E8F", "TRUE" = "#B5651D"),
                               labels = c("exact (from draws)", "rescale approx"),
                               name = NULL) +
    ggplot2::scale_x_continuous(labels = function(x) sprintf("%+d%%", as.integer(x))) +
    ggplot2::labs(
      x = "Change in global conflict-attributable YLL vs headline",
      y = NULL,
      title = "Structured (one-at-a-time) sensitivity of the headline YLL total",
      subtitle = sprintf("Central (headline) undiscounted total = %s. Bars show single-axis perturbations.",
                         format(round(central_mean), big.mark = ",")),
      caption = "Approx axes rescale existing draws by the mean CFR/L factor; ITS-horizon requires a re-run.") +
    { if (exists("theme_nm")) theme_nm() else ggplot2::theme_minimal() }
}
