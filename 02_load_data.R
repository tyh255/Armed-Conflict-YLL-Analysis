# ---- verify inputs against MANIFEST (sha256) ----
MANIFEST <- utils::read.csv(file.path(ROOT, "data", "MANIFEST.csv"), stringsAsFactors = FALSE)
verify_inputs <- function(man = MANIFEST) {
  for (i in seq_len(nrow(man))) {
    f <- file.path(ROOT, "data", man$file[i])
    if (!file.exists(f)) {
      if (startsWith(man$file[i], "gbd/"))
        stop("GBD input missing: data/", man$file[i], "\n  IHME licence forbids redistribution. ",
             "Download it with the query in data/gbd/README.md and save under that exact name.",
             call. = FALSE)
      stop("Input missing: data/", man$file[i], call. = FALSE)
    }
    h <- digest::digest(file = f, algo = "sha256")
    if (identical(h, man$sha256[i])) next
    if (startsWith(man$file[i], "gbd/"))   # a re-download can differ in row order/formatting only
      warning("sha256 differs from the archived GBD extraction for data/", man$file[i],
              "\n  Results may differ slightly if the GBD values differ; see data/gbd/README.md.", call. = FALSE)
    else stop("sha256 mismatch for data/", man$file[i], "\n  expected ", man$sha256[i], "\n  got      ", h, call. = FALSE)
  }
  message("Inputs checked against MANIFEST sha256: ", nrow(man), " files.")
}
verify_inputs()

read_csv_safe <- function(path) utils::read.csv(path, stringsAsFactors = FALSE)

# ---- country names: GBD/WUENIC spellings -> one key ----
harmonise_country <- function(x) {
  dplyr::case_when(
    x == "American Samoa"                                           ~ "Samoa",
    x == "Côte d'Ivoire"                                       ~ "Cote d'Ivoire",
    x == "occupied Palestinian territory, including east Jerusalem" ~ "State of Palestine",
    x == "Türkiye"                                             ~ "Turkiye",
    x == "United Kingdom of Great Britain and Northern Ireland"     ~ "United Kingdom",
    x == "United States of America"                                 ~ "United States",
    x == "Syrian Arab Republic"                                     ~ "Syria",
    x == "Iran (Islamic Republic of)"                               ~ "Iran",
    x == "Lao People's Democratic Republic"                         ~ "Laos",
    x == "Viet Nam"                                                 ~ "Vietnam",
    x == "Republic of Korea"                                        ~ "South Korea",
    x == "Democratic People's Republic of Korea"                    ~ "North Korea",
    x == "Russian Federation"                                       ~ "Russia",
    x == "United Republic of Tanzania"                              ~ "Tanzania",
    x == "Bolivia (Plurinational State of)"                         ~ "Bolivia",
    x == "Venezuela (Bolivarian Republic of)"                       ~ "Venezuela",
    x == "Republic of Moldova"                                      ~ "Moldova",
    x == "The former Yugoslav Republic of Macedonia"                ~ "North Macedonia",
    x == "Czech Republic"                                           ~ "Czechia",
    x == "Eswatini (Swaziland)"                                     ~ "Eswatini",
    x == "Brunei Darussalam"                                        ~ "Brunei",
    x == "Micronesia (Federated States of)"                         ~ "Micronesia",
    TRUE                                                            ~ x)
}

# ---- World Bank wide -> long ----
to_long_wb <- function(df, value_name) {
  df %>%
    tidyr::pivot_longer(cols = starts_with("X"), names_to = "year", values_to = value_name) %>%
    dplyr::mutate(year = as.integer(stringr::str_extract(year, "\\d{4}"))) %>%
    dplyr::rename(ISO3 = Country.Code) %>%
    dplyr::select(ISO3, year, dplyr::all_of(value_name))
}
WB_births_long <- to_long_wb(read_csv_safe(file.path(RAW_DIR, "wb_birth_rate.csv")), "Birth_Rate")
WB_pop_long    <- to_long_wb(read_csv_safe(file.path(RAW_DIR, "wb_population.csv")), "Population") %>%
  dplyr::mutate(Population_per_1000 = Population / 1000)
WB_life_long   <- to_long_wb(read_csv_safe(file.path(RAW_DIR, "wb_life_expectancy.csv")), "Life_Expectancy")

# ---- GBD 2023 incidence: Incidence / Rate / Both / <1 year; per 100k -> per person ----
read_ihme_incidence <- function(filename, out_prefix, age_band = "<1 year") {
  df <- read_csv_safe(file.path(GBD_DIR, filename))
  out <- df %>%
    dplyr::filter(measure_name == "Incidence", metric_name == "Rate",
                  sex_name == "Both", age_name == age_band) %>%
    dplyr::mutate(country = harmonise_country(location_name), year = as.integer(year),
                  inc_rate = val / 1e5, inc_upper = upper / 1e5, inc_lower = lower / 1e5) %>%
    dplyr::mutate(inc_sigma = dplyr::if_else(               # log-sd from 95% UI
      inc_upper > 0 & inc_lower > 0 & inc_upper > inc_lower,
      (log(inc_upper) - log(inc_lower)) / 3.92, NA_real_)) %>%
    dplyr::select(country, year, inc_rate, inc_upper, inc_lower, inc_sigma) %>%
    dplyr::group_by(country, year) %>%
    dplyr::summarise(inc_rate = mean(inc_rate, na.rm = TRUE), inc_upper = mean(inc_upper, na.rm = TRUE),
                     inc_lower = mean(inc_lower, na.rm = TRUE), inc_sigma = mean(inc_sigma, na.rm = TRUE),
                     .groups = "drop")
  names(out)[3:6] <- paste0(out_prefix, c("_Inc_Rate", "_Inc_Upper", "_Inc_Lower", "_Inc_Sigma"))
  out
}
ihme_tb      <- read_ihme_incidence("IHME_Tuberculosis_GBD.csv", "TB")
ihme_measles <- read_ihme_incidence("IHME_Measles_GBD.csv",      "Measles")
ihme_diphth  <- read_ihme_incidence("IHME_Diphtheria_GBD.csv",   "Diphtheria")
ihme_pert    <- read_ihme_incidence("IHME_Pertussis_GBD.csv",    "Pertussis")
ihme_tet     <- read_ihme_incidence("IHME_Tetanus_GBD.csv",      "Tetanus")   # GBD "Tetanus" excludes neonatal

# ---- WUENIC coverage ----
reshape_coverage <- function(df, label) {
  df %>%
    dplyr::mutate(country = harmonise_country(Location), ISO3 = SpatialDimValueCode,
                  year = suppressWarnings(as.integer(Period)),
                  coverage = suppressWarnings(as.numeric(Value)), vaccine = label) %>%
    dplyr::filter(!is.na(year), !is.na(coverage)) %>%
    dplyr::select(country, ISO3, vaccine, year, coverage)
}
bcg_long  <- reshape_coverage(read_csv_safe(file.path(RAW_DIR, "wuenic_bcg.csv")),  "BCG")
mcv1_long <- reshape_coverage(read_csv_safe(file.path(RAW_DIR, "wuenic_mcv1.csv")), "MCV1")
dpt_long  <- reshape_coverage(read_csv_safe(file.path(RAW_DIR, "wuenic_dtp3.csv")), "DTP3")
iso3_lookup <- bcg_long %>% dplyr::filter(!is.na(ISO3), ISO3 != "") %>% dplyr::distinct(country, ISO3)
backfill_iso3 <- function(df) {
  df %>%
    dplyr::left_join(iso3_lookup, by = "country", suffix = c("", "_lookup")) %>%
    dplyr::mutate(ISO3 = ifelse(is.na(ISO3) | ISO3 == "", ISO3_lookup, ISO3)) %>%
    dplyr::select(-ISO3_lookup)
}
coverage_long <- dplyr::bind_rows(bcg_long, backfill_iso3(mcv1_long), backfill_iso3(dpt_long)) %>%
  dplyr::filter(!is.na(ISO3), ISO3 != "") %>%
  dplyr::arrange(country, vaccine, year)

# ---- country-year covariates ----
covariates_panel <- WB_births_long %>%
  dplyr::full_join(WB_pop_long %>% dplyr::select(ISO3, year, Population, Population_per_1000),
                   by = c("ISO3", "year")) %>%
  dplyr::full_join(WB_life_long, by = c("ISO3", "year")) %>%
  dplyr::arrange(ISO3, year) %>%
  dplyr::mutate(Birth_Cohort = Birth_Rate * Population_per_1000)   # births = CBR x pop/1000
country_iso_lookup <- coverage_long %>% dplyr::distinct(country, ISO3)
covariates_panel <- covariates_panel %>%
  dplyr::left_join(country_iso_lookup, by = "ISO3") %>%
  dplyr::left_join(ihme_tb,      by = c("country", "year")) %>%
  dplyr::left_join(ihme_measles, by = c("country", "year")) %>%
  dplyr::left_join(ihme_diphth,  by = c("country", "year")) %>%
  dplyr::left_join(ihme_pert,    by = c("country", "year")) %>%
  dplyr::left_join(ihme_tet,     by = c("country", "year"))

# linear interpolation within ISO3; rule = 2 carries end values
interpolate_panel <- function(df, vars) {
  df %>%
    dplyr::group_by(ISO3) %>% dplyr::arrange(year) %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(vars), ~ {
      x <- .x
      if (sum(!is.na(x)) < 2) return(x)
      zoo::na.approx(x, x = year, na.rm = FALSE, rule = 2)
    })) %>%
    dplyr::ungroup()
}
incidence_vars <- as.vector(outer(c("TB", "Measles", "Diphtheria", "Pertussis", "Tetanus"),
                                  c("_Inc_Rate", "_Inc_Upper", "_Inc_Lower", "_Inc_Sigma"), paste0) |> t())
covariates_panel <- interpolate_panel(covariates_panel,
                                      c("Birth_Rate", "Population_per_1000", "Birth_Cohort",
                                        "Life_Expectancy", incidence_vars))

# ---- economic / human-development covariates (SC + covariates) ----
WB_gdp_long <- to_long_wb(read_csv_safe(file.path(RAW_DIR, "wb_gdp_pc.csv")), "GDP_pc") %>%
  dplyr::mutate(GDP_pc = suppressWarnings(as.numeric(GDP_pc))) %>%
  dplyr::filter(!is.na(year)) %>%
  dplyr::mutate(log_GDP_pc = ifelse(!is.na(GDP_pc) & GDP_pc > 0, log(GDP_pc), NA_real_))

HDI_INDICATORS <- c(eys = "sch_eys", mys = "sch_mys", gnipc = "gni_pc", hdi = "hdi")
hdi_raw <- read_csv_safe(file.path(RAW_DIR, "undp_hdi.csv"))
present_codes <- intersect(names(HDI_INDICATORS), unique(hdi_raw$indicatorCode))
HDI_long <- hdi_raw %>%
  dplyr::filter(indicatorCode %in% present_codes) %>%
  dplyr::transmute(ISO3 = countryIsoCode, country = harmonise_country(country),
                   year = suppressWarnings(as.integer(year)),
                   col = unname(HDI_INDICATORS[indicatorCode]),
                   value = suppressWarnings(as.numeric(value))) %>%
  dplyr::filter(!is.na(year), !is.na(col)) %>%
  dplyr::group_by(ISO3, country, year, col) %>%
  dplyr::summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = col, values_from = value)

covariates_panel <- covariates_panel %>%
  dplyr::left_join(WB_gdp_long, by = c("ISO3", "year")) %>%
  dplyr::left_join(HDI_long %>% dplyr::select(-country), by = c("ISO3", "year"))
.econ_cols <- intersect(c("GDP_pc", "log_GDP_pc", setdiff(names(HDI_long), c("ISO3", "country", "year"))),
                        names(covariates_panel))
covariates_panel <- interpolate_panel(covariates_panel, .econ_cols)

econ_panel <- covariates_panel %>%
  dplyr::select(dplyr::all_of(c("ISO3", "year", .econ_cols))) %>%
  dplyr::left_join(country_iso_lookup, by = "ISO3") %>%
  dplyr::filter(!is.na(country)) %>%
  dplyr::distinct(ISO3, country, year, .keep_all = TRUE)
ECON_PREDICTORS <- intersect(c("log_GDP_pc", "sch_eys", "sch_mys", "hdi"), names(econ_panel))

stopifnot(all(conflict_info$ISO3 %in% covariates_panel$ISO3))
message("Data: coverage rows ", nrow(coverage_long), " | covariate rows ", nrow(covariates_panel),
        " | SC predictors {", paste(ECON_PREDICTORS, collapse = ", "), "}")
