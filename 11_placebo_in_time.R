
.placebo_conflict_info <- function(conflict_df, baseline_window = 3L) {
  dur <- pmax(conflict_df$conflict_end - conflict_df$conflict_start, 1L)
  conflict_df %>%
    dplyr::mutate(
      .dur             = dur,
      conflict_end     = conflict_start - 1L,
      conflict_start   = conflict_start - .dur,
      preconflict_year = conflict_start - 1L,
      preconflict_window = purrr::map(conflict_start - 1L,
                                      ~ (.x - (baseline_window - 1L)):.x)) %>%
    dplyr::select(-.dur)
}

# Summarise gap magnitude by method over a gap table.
.summarise_gaps_by_method <- function(g, tag) {
  if (is.null(g) || nrow(g) == 0) return(tibble::tibble())
  g %>%
    dplyr::group_by(method) %>%
    dplyr::summarise(
      !!paste0("n_", tag)            := dplyr::n(),
      !!paste0("mean_", tag)         := mean(gap, na.rm = TRUE),
      !!paste0("median_", tag)       := stats::median(gap, na.rm = TRUE),
      !!paste0("mean_abs_", tag)     := mean(abs(gap), na.rm = TRUE),
      !!paste0("frac_gt5pp_", tag)   := mean(abs(gap) > 5, na.rm = TRUE),
      .groups = "drop")
}

# Driver: run placebo gaps for all methods, contrast with real gaps, verdict.
run_placebo_in_time <- function(coverage_df, conflict_df = conflict_info,
                                gap_df_real = NULL,
                                methods = if (exists("ALL_METHODS")) ALL_METHODS else NULL,
                                placebo_tol_pp = 3, ratio_tol = 0.33,
                                seed = if (exists("DEFAULT_SEED")) DEFAULT_SEED else 42) {
  fake <- .placebo_conflict_info(conflict_df)
  message("  Placebo-in-time: estimating gaps on pre-onset windows (this re-runs ",
          "SC/DiD; expect skips where the placebo window pre-dates coverage data).")
  gp <- tryCatch(estimate_all_gaps_plus(coverage_df, fake, methods = methods),
                 error = function(e) { message("    placebo gap step failed: ",
                                               conditionMessage(e)); NULL })
  if (is.null(gp) || nrow(gp) == 0) {
    message("  Placebo-in-time: no placebo gaps produced -> skipping."); return(invisible(NULL))
  }
  placebo <- .summarise_gaps_by_method(gp, "placebo")
  
  # Real gaps for contrast (use the passed gap_df if available, else recompute).
  real_src <- if (!is.null(gap_df_real)) gap_df_real else
    tryCatch(estimate_all_gaps_plus(coverage_df, conflict_df, methods = methods),
             error = function(e) NULL)
  real <- .summarise_gaps_by_method(real_src, "real")
  
  out <- placebo %>% dplyr::left_join(real, by = "method") %>%
    dplyr::mutate(
      placebo_to_real_ratio = ifelse(is.finite(mean_abs_real) & mean_abs_real > 0,
                                     mean_abs_placebo / mean_abs_real, NA_real_),
      verdict = dplyr::case_when(
        abs(mean_placebo) <= placebo_tol_pp &
          (is.na(placebo_to_real_ratio) | placebo_to_real_ratio <= ratio_tol) ~ "PASS (placebo ~ 0)",
        abs(mean_placebo) <= placebo_tol_pp                                    ~ "borderline (small abs, high ratio)",
        TRUE                                                                    ~ "FAIL (systematic placebo gap)"))
  
  # SELF-TEST / DIAGNOSTICS: a placebo mean gap far from 0 is the failure mode.
  failed <- out %>% dplyr::filter(grepl("FAIL", verdict))
  if (nrow(failed) > 0)
    message("  [placebo diag] systematic non-zero placebo gap for: ",
            paste(sprintf("%s (mean %.1fpp)", failed$method, failed$mean_placebo),
                  collapse = ", "),
            " -> that counterfactual is partly trend-driven; weight accordingly.")
  else
    message("  [placebo diag] all estimators return near-zero placebo gaps ",
            "(|mean| <= ", placebo_tol_pp, "pp): falsification passed.")
  
  if (exists("save_tab")) save_tab(out, "Table_S7_placebo_in_time.csv")
  fig <- tryCatch(.fig_placebo(out),
                  error = function(e) { message("  placebo fig skipped: ",
                                                conditionMessage(e)); NULL })
  if (!is.null(fig) && exists("save_fig"))
    save_fig(fig, "Figure_S_placebo_in_time", width = 9, height = 5.5)
  invisible(out)
}

# Figure: placebo vs real mean |gap| by method (placebo bars should be near 0).
.fig_placebo <- function(out) {
  meth_levels <- if (exists("METHOD_ORDER")) METHOD_ORDER else unique(out$method)
  long <- out %>%
    dplyr::select(method, mean_abs_placebo, mean_abs_real) %>%
    tidyr::pivot_longer(c(mean_abs_placebo, mean_abs_real),
                        names_to = "kind", values_to = "mean_abs_gap") %>%
    dplyr::mutate(
      kind = ifelse(kind == "mean_abs_placebo", "Placebo (pre-onset)", "Real (conflict)"),
      method = factor(method, levels = intersect(meth_levels, unique(method))))
  ggplot2::ggplot(long, ggplot2::aes(x = mean_abs_gap, y = method, fill = kind)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7), width = 0.65) +
    ggplot2::scale_fill_manual(values = c("Placebo (pre-onset)" = "#B5651D",
                                          "Real (conflict)" = "#2C6E8F"), name = NULL) +
    ggplot2::labs(
      x = "Mean |coverage gap| (percentage points)", y = NULL,
      title = "In-time placebo falsification: pre-onset gaps should be ~ 0",
      subtitle = "Placebo re-indexes each conflict to an equal-length pre-onset window; large placebo bars indicate a trend-driven counterfactual.",
      caption = "Abadie et al. 2015; Chernozhukov, Wuthrich & Zhu 2021. Real bars shown for contrast.") +
    { if (exists("theme_nm")) theme_nm() else ggplot2::theme_minimal() }
}
