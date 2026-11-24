select_countries <- function(path = file.path(RAW_DIR, "ucdp_acd_v26_1.csv"),
                             years = 2000:2024, min_run = 3L) {
  acd <- utils::read.csv(path, stringsAsFactors = FALSE)
  stopifnot(all(acd$version == 26.1))
  war <- acd %>%
    dplyr::filter(type_of_conflict %in% c(3, 4), intensity_level == 2, year %in% years) %>%
    dplyr::arrange(conflict_id, year) %>%
    dplyr::group_by(conflict_id) %>%
    dplyr::mutate(run = cumsum(c(1L, diff(year) != 1L))) %>%
    dplyr::group_by(conflict_id, run) %>%
    dplyr::summarise(location = dplyr::first(location),
                     gwno_loc = dplyr::first(gwno_loc),
                     war_start = min(year), war_end = max(year), n_years = dplyr::n(),
                     episode_start = as.integer(substr(dplyr::first(start_date2), 1, 4)),
                     side_b = dplyr::first(side_b), .groups = "drop") %>%
    dplyr::filter(n_years >= min_run) %>%
    dplyr::mutate(ISO3 = countrycode::countrycode(as.integer(sub(",.*", "", gwno_loc)), "gwn", "iso3c",
                                                  custom_match = c(`678` = "YEM")),
                  episode_in_window = episode_start >= min(years))
  by_country <- war %>%
    dplyr::group_by(ISO3) %>%
    dplyr::summarise(location = dplyr::first(location),
                     war_runs = paste0(conflict_id, ":", war_start, "-", war_end, " (ep ", episode_start, ")",
                                       collapse = "; "),
                     first_war_year = min(war_start),
                     selected = all(episode_in_window), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(selected), ISO3)
  save_tab(by_country, "Table_country_selection.csv")
  sel <- sort(by_country$ISO3[by_country$selected])
  expected <- sort(c("BFA", "ETH", "IRQ", "MMR", "NGA", "PAK", "SOM", "LKA", "SYR", "YEM"))
  if (!identical(sel, expected))
    stop("Country selection from ACD v26.1 gave ", paste(sel, collapse = " "),
         "; expected ", paste(expected, collapse = " "))
  stopifnot(setequal(sel, conflict_info$ISO3))
  message("Selection: ", length(sel), " countries from ACD v26.1 (", nrow(by_country), " candidates): ",
          paste(sel, collapse = " "))
  invisible(by_country)
}
