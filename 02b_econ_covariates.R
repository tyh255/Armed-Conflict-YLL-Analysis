stopifnot(exists("covariates_panel"), exists("harmonise_country"),
          exists("to_long_wb"), exists("YLL_DATA_DIR"))

GDP_PC_FILE <- Sys.getenv("CVD_GDP_FILE",
                          unset = file.path(YLL_DATA_DIR, "WB GDPpC.csv"))
HDI_FILE    <- Sys.getenv("CVD_HDI_FILE",
                          unset = file.path(YLL_DATA_DIR, "hdi_data.csv"))

# HDR indicatorCodes to pull, mapped to tidy column names. Extra codes that are
# absent from the file are skipped with a message (HDR releases vary).
HDI_INDICATORS <- c(
  eys   = "sch_eys",     # Expected years of schooling
  mys   = "sch_mys",     # Mean years of schooling
  gnipc = "gni_pc",      # GNI per capita (2017 PPP $)
  hdi   = "hdi"          # Human Development Index (composite)
)

# ---- 1. World Bank GDP per capita -----------------------------------------
WB_gdp_long <- NULL
if (file.exists(GDP_PC_FILE)) {
  gdp_raw <- utils::read.csv(GDP_PC_FILE, stringsAsFactors = FALSE,
                             check.names = TRUE)
  # to_long_wb() keys on Country.Code + columns beginning "X" (read.csv prefixes
  # bare numeric year headers with "X"). Matches the other WB loaders in 02.
  WB_gdp_long <- to_long_wb(gdp_raw, "GDP_pc") %>%
    dplyr::mutate(GDP_pc = suppressWarnings(as.numeric(GDP_pc))) %>%
    dplyr::filter(!is.na(year)) %>%
    dplyr::mutate(log_GDP_pc = ifelse(!is.na(GDP_pc) & GDP_pc > 0,
                                      log(GDP_pc), NA_real_))
  message("Econ covariates: loaded GDP per capita for ",
          dplyr::n_distinct(WB_gdp_long$ISO3), " ISO3 codes from ",
          basename(GDP_PC_FILE), ".")
} else {
  warning("GDP per capita file not found at\n  ", GDP_PC_FILE,
          "\n  -> GDP_pc/log_GDP_pc will be NA; SC(+cov)/augSC/DiD covariate ",
          "matching on GDP is disabled.")
}

# ---- 2. UNDP HDR human-development series ----------------------------------
HDI_long <- NULL
if (file.exists(HDI_FILE)) {
  hdi_raw <- utils::read.csv(HDI_FILE, stringsAsFactors = FALSE)
  need <- c("countryIsoCode", "indicatorCode", "year", "value")
  miss <- setdiff(need, names(hdi_raw))
  if (length(miss) > 0) {
    warning("hdi_data.csv missing columns: ", paste(miss, collapse = ", "),
            " -> HDI covariates disabled.")
  } else {
    present_codes <- intersect(names(HDI_INDICATORS), unique(hdi_raw$indicatorCode))
    absent_codes  <- setdiff(names(HDI_INDICATORS), present_codes)
    if (length(absent_codes) > 0)
      message("Econ covariates: HDR indicatorCodes not in file (skipped): ",
              paste(absent_codes, collapse = ", "),
              " -- adjust HDI_INDICATORS if your release names them differently.")
    HDI_long <- hdi_raw %>%
      dplyr::filter(indicatorCode %in% present_codes) %>%
      dplyr::transmute(
        ISO3    = countryIsoCode,
        country = harmonise_country(if ("country" %in% names(hdi_raw)) country else countryIsoCode),
        year    = suppressWarnings(as.integer(year)),
        col     = unname(HDI_INDICATORS[indicatorCode]),
        value   = suppressWarnings(as.numeric(value))
      ) %>%
      dplyr::filter(!is.na(year), !is.na(col)) %>%
      dplyr::group_by(ISO3, country, year, col) %>%
      dplyr::summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = col, values_from = value)
    message("Econ covariates: loaded HDR indicators {",
            paste(present_codes, collapse = ", "), "} for ",
            dplyr::n_distinct(HDI_long$ISO3), " ISO3 codes from ",
            basename(HDI_FILE), ".")
  }
} else {
  warning("HDI file not found at\n  ", HDI_FILE,
          "\n  -> schooling/HDI covariates will be NA.")
}

# ---- 3. Merge into covariates_panel (by ISO3, year) ------------------------
.econ_cols <- character(0)
if (!is.null(WB_gdp_long)) {
  covariates_panel <- covariates_panel %>%
    dplyr::left_join(WB_gdp_long, by = c("ISO3", "year"))
  .econ_cols <- c(.econ_cols, "GDP_pc", "log_GDP_pc")
}
if (!is.null(HDI_long)) {
  hdi_cols <- setdiff(names(HDI_long), c("ISO3", "country", "year"))
  covariates_panel <- covariates_panel %>%
    dplyr::left_join(HDI_long %>% dplyr::select(-country),
                     by = c("ISO3", "year"))
  .econ_cols <- c(.econ_cols, hdi_cols)
}

# Interpolate within ISO3 (same rule=2 LOCF-at-ends convention as 02). Only
# touch columns that exist.
if (length(.econ_cols) > 0 && exists("interpolate_panel")) {
  .econ_present <- intersect(.econ_cols, names(covariates_panel))
  covariates_panel <- interpolate_panel(covariates_panel, .econ_present)
}

# ---- 4. Tidy econ_panel for SC/DiD predictors ------------------------------
# Keyed by (ISO3, country, year). country uses the coverage-panel spelling so
# joins against coverage_long (the SC unit identifier) succeed.
if (length(.econ_cols) > 0) {
  .keep <- c("ISO3", "year", intersect(.econ_cols, names(covariates_panel)))
  econ_panel <- covariates_panel %>%
    dplyr::select(dplyr::all_of(.keep)) %>%
    dplyr::left_join(country_iso_lookup, by = "ISO3") %>%
    dplyr::filter(!is.na(country)) %>%
    dplyr::distinct(ISO3, country, year, .keep_all = TRUE)
  ECON_PREDICTORS <- intersect(c("log_GDP_pc", "sch_eys", "sch_mys", "hdi"),
                               names(econ_panel))
  message("Econ panel ready: ", nrow(econ_panel), " ISO3-year rows; ",
          "candidate SC/DiD predictors = {",
          paste(ECON_PREDICTORS, collapse = ", "), "}.")
} else {
  econ_panel <- NULL
  ECON_PREDICTORS <- character(0)
  warning("No economic covariates loaded; econ_panel is NULL and the ",
          "covariate-adjusted estimators in 03b will be skipped.")
}
