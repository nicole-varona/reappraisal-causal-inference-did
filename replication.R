# =============================================================================
# replication.R
# =============================================================================
# Heavy pre-computation script for the reappraisal of Schaff (2024).
#
# This script estimates all DCDH specifications and event-study figures
# reported in MY457_reappraisal_62628.qmd. Outputs are saved as .rds files
# in output/. The .qmd document loads these files at render time.
#
# REPRODUCIBILITY NOTE:
# All estimation blocks below are commented out by default. The .rds files
# in output/ are already present in the repository, so the .qmd renders
# without re-running this script. To regenerate any specific .rds, uncomment
# the corresponding block. To regenerate everything, uncomment all blocks
# (total computation ~20 hours on a 16 GB machine).
#
# IMPORTANT: The "+ Controls" panels of Figure 3 (Panels A.II, B.II, C.II)
# and the equivalent panels of Figures 4–6 require more memory than a 16 GB
# machine permits. Their estimation code is preserved (commented out) for
# execution on higher-memory hardware but was NOT executed for the rendered
# HTML. This limitation is discussed in the qmd document.
#
# STRUCTURE:
#   0. Setup (libraries, paths, helper functions, data)
#   1. Table 1 - Computed inside the .qmd (lightweight OLS)
#   2. Table 2 - Computed inside the .qmd (manual within-FE)
#   3. Table 3 - DCDH individual DiD (9 specs)
#   4. Figure 3 - Event-study (3 Baseline + 3 "+ Controls" not executed)
#   5. Table 4 - Mechanism tests (8 specs)
#   6. Figure 4 - Event-study clerks (Baseline)
#   7. Table 5 - Triple-difference (7 specs)
#   8. Figure 5 - Event-study triple-diff (2 Baseline)
#   9. Figure 6 - Event-study clerks during the war (Baseline)
#
# Author: 62628 | MY457 Causal Inference, WT 2026, LSE.
# =============================================================================

# =============================================================================
# 0. SETUP
# =============================================================================

## ---- 0.1 Packages -------------------------------------------------------- ##

# macOS Tahoe (26.x) fix: rgl does not find OpenGL. Set NULL mode before
# loading DIDmultiplegtDYN (which depends transitively on rgl).
options(rgl.useNULL = TRUE)
Sys.setenv(RGL_USE_NULL = TRUE)

# To install (uncomment on first run):
# install.packages(c(
#   "haven", "dplyr", "tidyr", "fixest",
#   "DIDmultiplegtDYN", "polars", "DIDmultiplegt",
#   "ggplot2", "modelsummary", "kableExtra",
#   "patchwork", "broom", "purrr", "here"
# ))

library(haven)
library(dplyr)
library(tidyr)
library(fixest)
library(DIDmultiplegtDYN)
library(polars)              # backend of DIDmultiplegtDYN (required)
library(ggplot2)
library(modelsummary)
library(kableExtra)
library(patchwork)
library(broom)
library(purrr)
library(here)

## ---- 0.2 Paths ----------------------------------------------------------- ##

data_dir <- here::here()
out_dir  <- here::here("output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

## ---- 0.3 Global options and seed ----------------------------------------- ##

set.seed(123456789)            # Same seed as Schaff for bootstrap
options(digits = 4, scipen = 999)

# Unified plot theme for all figures
theme_schaff <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title       = element_text(face = "plain", size = 11),
    axis.title       = element_text(size = 10)
  )

## ---- 0.4 Load data ------------------------------------------------------- ##

df_noe <- read_dta(file.path(data_dir, "data", "Noerdlingen_panel_replication_JEH.dta"))

# Thirty Years' War subset (Table 5, Figure 5, Figure 6)
df_30yw <- df_noe %>%
  filter(units_30yw == 1, year >= 1603, year <= 1646)

## ---- 0.5 Auxiliary variables --------------------------------------------- ##

# The .do file uses occupation_* as a wildcard; build it explicitly in R.
occ_vars <- grep("^occupation_", names(df_noe), value = TRUE)

# Occupational controls for Table 4 and Figure 4 (excluding the clerk cluster
# 670-720 when the treatment is cityadmin_f, to avoid collinearity).
occ_vars_no_clerks <- setdiff(occ_vars,
                              c("occupation_670", "occupation_680",
                                "occupation_690", "occupation_700",
                                "occupation_710", "occupation_720"))

## ---- 0.6 Wrapper for the DCDH estimator ---------------------------------- ##

# Wraps did_multiplegt_dyn with sensible defaults matching the paper's
# Stata specification: bootstrap with 100 replications, clustered SEs.
run_dcdh <- function(data, y, d,
                     controls    = NULL,
                     effects     = 0,
                     placebo     = 0,
                     time_var    = "periodID",
                     id_var      = "ID",
                     nboot       = 100,
                     seed        = 123456789) {
  
  vars_used <- c(y, d, time_var, id_var, controls)
  if (!all(vars_used %in% names(data))) {
    miss <- vars_used[!vars_used %in% names(data)]
    stop("Variables not found: ", paste(miss, collapse = ", "))
  }
  
  d_clean <- data[complete.cases(data[, vars_used]), ]
  
  args <- list(
    df         = d_clean,
    outcome    = y,
    group      = id_var,
    time       = time_var,
    treatment  = d,
    effects    = max(effects, 1),
    cluster    = id_var,
    bootstrap  = list(nboot, seed),
    graph_off  = TRUE
  )
  
  if (placebo > 0) args$placebo <- placebo
  if (!is.null(controls) && length(controls) > 0) args$controls <- controls
  
  set.seed(seed)
  do.call(did_multiplegt_dyn, args)
}

## ---- 0.7 Extract ATT summary from a did_multiplegt_dyn object ------------ ##

extract_dcdh_summary <- function(res, label_y, label_spec) {
  eff <- res$results$Effects
  tibble(
    Outcome       = label_y,
    Specification = label_spec,
    ATT           = round(eff[1, "Estimate"], 3),
    SE            = round(eff[1, "SE"],       3),
    Lo_95_CI      = round(eff[1, "LB CI"],    3),
    Up_95_CI      = round(eff[1, "UB CI"],    3),
    N             = eff[1, "N"]
  )
}

## ---- 0.8 Event-study plot ------------------------------------------------ ##

# Plots placebo and effect coefficients from a did_multiplegt_dyn result.
# Placebos are indexed k = -2, -3, ..., (k = -1 is omitted reference).
# Effects are indexed k = 0, 1, 2, ...
plot_event_dcdh <- function(res, title = "", y_label = "Treatment effect",
                            y_limits = NULL) {
  
  pl <- as.data.frame(res$results$Placebos)
  pl$period <- -(seq_len(nrow(pl)) + 1)
  
  ef <- as.data.frame(res$results$Effects)
  ef$period <- seq_len(nrow(ef)) - 1
  
  names(pl)[names(pl) == "Estimate"] <- "est"
  names(ef)[names(ef) == "Estimate"] <- "est"
  names(pl)[names(pl) == "LB CI"]    <- "lb"
  names(ef)[names(ef) == "LB CI"]    <- "lb"
  names(pl)[names(pl) == "UB CI"]    <- "ub"
  names(ef)[names(ef) == "UB CI"]    <- "ub"
  
  df_plot <- bind_rows(
    pl %>% select(period, est, lb, ub),
    ef %>% select(period, est, lb, ub)
  )
  
  p <- ggplot(df_plot, aes(x = period, y = est)) +
    geom_hline(yintercept = 0, color = "grey50", linewidth = 0.3) +
    geom_vline(xintercept = -0.5, linetype = "dashed", color = "red",
               linewidth = 0.3) +
    geom_errorbar(aes(ymin = lb, ymax = ub), width = 0.1, color = "navy") +
    geom_point(color = "navy", size = 2) +
    labs(title = title, x = "Periods since the event", y = y_label) +
    theme_schaff
  
  if (!is.null(y_limits)) p <- p + ylim(y_limits)
  
  p
}

# =============================================================================
# 1. TABLE 1 - PRE-TREATMENT BALANCE
# =============================================================================
# Note: Table 1 is fully computed inside the .qmd document (it uses only OLS
# univariate regressions and runs in seconds). The block below is preserved
# as a reference implementation; it does not feed into the rendered HTML.
# -----------------------------------------------------------------------------

# df_pre <- df_noe %>%
#   group_by(ID) %>%
#   mutate(councilmember_tinvar = as.integer(any(councilmemberterm == 1, na.rm = TRUE))) %>%
#   ungroup() %>%
#   filter(councilmemberterm == 0)
#
# vars_balance <- c(
#   "logwealth", "wealthpercentile", "top5percent",
#   "ntaxpayments", "womentaxpayer", "cityadmin_f",
#   "occupation_500", "occupation_640",
#   "occupation_900", "occupation_100"
# )
#
# get_balance_row <- function(var) {
#   m <- lm(as.formula(paste(var, "~ councilmember_tinvar")), data = df_pre)
#   tibble(
#     variable = var,
#     diff     = coef(m)["councilmember_tinvar"],
#     se       = sqrt(vcov(m)["councilmember_tinvar", "councilmember_tinvar"]),
#     p_value  = summary(m)$coefficients["councilmember_tinvar", "Pr(>|t|)"],
#     mean_all = mean(df_pre[[var]], na.rm = TRUE)
#   )
# }
#
# tab1 <- purrr::map_dfr(vars_balance, get_balance_row)


# =============================================================================
# 2. TABLE 2 - CITY-LEVEL TWFE
# =============================================================================
# Note: Table 2 is fully computed inside the .qmd document. The within-FE
# estimator is implemented manually to match Stata's xtreg, fe vce(cluster)
# convention exactly (including the within-R^2 statistic, which fixest does
# not report by default). The block below is preserved as reference.
# -----------------------------------------------------------------------------

# df_city <- read_dta(file.path(data_dir, "data", "City_panel_replication_JEH.dta"))
#
# controls_city <- c("conflict25", "logpopulation", "epidemic",
#                    "reformationintro", "loguniversitydist")
#
# run_city_reg <- function(outcome, with_controls) {
#   rhs <- if (with_controls) {
#     paste(c("electionparticipation", controls_city), collapse = " + ")
#   } else {
#     "electionparticipation"
#   }
#   fml <- as.formula(paste0(outcome, " ~ ", rhs, " | IDlocality + yearineq25"))
#   feols(fml, data = df_city, cluster = ~ IDlocality)
# }


# =============================================================================
# 3. TABLE 3 - INDIVIDUAL DiD (DCDH)
# =============================================================================
# 9 specifications of the individual DiD with 100-replication clustered
# bootstrap. Each spec saves its .rds before clearing memory; if R crashes
# on spec K, the process can resume from spec K+1 without re-running the
# previous ones.
#
# All blocks are commented out by default. The .rds files in output/ are
# already present from the original run. To regenerate any block, uncomment
# the relevant lines.
# -----------------------------------------------------------------------------

## ---- 3.1 Estimate each spec ---------------------------------------------- ##

# Spec 1/9 (dCdH_simpleDiD_0.rds): logwealth, Only FE
# res_tmp <- run_dcdh(df_noe, y = "logwealth", d = "councilmemberterm",
#                     controls = NULL, effects = 1, nboot = 100)
# saveRDS(res_tmp, file.path(out_dir, "dCdH_simpleDiD_0.rds"))
# rm(res_tmp); gc()

# Spec 2/9 (dCdH_simpleDiD_1.rds): logwealth, Baseline
# res_tmp <- run_dcdh(df_noe, y = "logwealth", d = "councilmemberterm",
#                     controls = c("ntaxpayments", "ntaxpayments2"),
#                     effects = 1, nboot = 100)
# saveRDS(res_tmp, file.path(out_dir, "dCdH_simpleDiD_1.rds"))
# rm(res_tmp); gc()

# Spec 3/9 (dCdH_simpleDiD_2.rds): logwealth, Baseline + Controls
# res_tmp <- run_dcdh(df_noe, y = "logwealth", d = "councilmemberterm",
#                     controls = c("ntaxpayments", "ntaxpayments2", occ_vars, "womentaxpayer"),
#                     effects = 1, nboot = 100)
# saveRDS(res_tmp, file.path(out_dir, "dCdH_simpleDiD_2.rds"))
# rm(res_tmp); gc()

# Spec 4/9 (dCdH_simpleDiD_3.rds): wealthpercentile, Baseline
# res_tmp <- run_dcdh(df_noe, y = "wealthpercentile", d = "councilmemberterm",
#                     controls = c("ntaxpayments", "ntaxpayments2"),
#                     effects = 1, nboot = 100)
# saveRDS(res_tmp, file.path(out_dir, "dCdH_simpleDiD_3.rds"))
# rm(res_tmp); gc()

# Spec 5/9 (dCdH_simpleDiD_4.rds): wealthpercentile, Baseline + Controls
# res_tmp <- run_dcdh(df_noe, y = "wealthpercentile", d = "councilmemberterm",
#                     controls = c("ntaxpayments", "ntaxpayments2", occ_vars, "womentaxpayer"),
#                     effects = 1, nboot = 100)
# saveRDS(res_tmp, file.path(out_dir, "dCdH_simpleDiD_4.rds"))
# rm(res_tmp); gc()

# Spec 6/9 (dCdH_simpleDiD_5.rds): top5percent, Baseline
# res_tmp <- run_dcdh(df_noe, y = "top5percent", d = "councilmemberterm",
#                     controls = c("ntaxpayments", "ntaxpayments2"),
#                     effects = 1, nboot = 100)
# saveRDS(res_tmp, file.path(out_dir, "dCdH_simpleDiD_5.rds"))
# rm(res_tmp); gc()

# Spec 7/9 (dCdH_simpleDiD_6.rds): top5percent, Baseline + Controls
# res_tmp <- run_dcdh(df_noe, y = "top5percent", d = "councilmemberterm",
#                     controls = c("ntaxpayments", "ntaxpayments2", occ_vars, "womentaxpayer"),
#                     effects = 1, nboot = 100)
# saveRDS(res_tmp, file.path(out_dir, "dCdH_simpleDiD_6.rds"))
# rm(res_tmp); gc()

# Spec 8/9 (dCdH_simpleDiD_7.rds): logwealth, Baseline (mayor sub-sample)
# res_tmp <- run_dcdh(df_noe, y = "logwealth", d = "mayorterm",
#                     controls = c("ntaxpayments", "ntaxpayments2"),
#                     effects = 1, nboot = 100)
# saveRDS(res_tmp, file.path(out_dir, "dCdH_simpleDiD_7.rds"))
# rm(res_tmp); gc()

# Spec 9/9 (dCdH_simpleDiD_8.rds): logwealth, Baseline + Controls (mayor sub-sample)
# res_tmp <- run_dcdh(df_noe, y = "logwealth", d = "mayorterm",
#                     controls = c("ntaxpayments", "ntaxpayments2", occ_vars, "womentaxpayer"),
#                     effects = 1, nboot = 100)
# saveRDS(res_tmp, file.path(out_dir, "dCdH_simpleDiD_8.rds"))
# rm(res_tmp); gc()


# =============================================================================
# 4. FIGURE 3 - EVENT-STUDY (FLEXIBLE DiD)
# =============================================================================
# 3 Baseline panels (logwealth, wealthpercentile, top5percent). The
# "+ Controls" panels (A.II, B.II, C.II) require more memory than a 16 GB
# machine permits, owing to the 61-dimensional occupational control set
# combined with the 100-replication clustered bootstrap. The estimation code
# for those panels is preserved here for execution on higher-memory hardware
# but was not executed for the rendered HTML.
# -----------------------------------------------------------------------------

## ---- 4.1 Panel A.I: logwealth, Baseline ---------------------------------- ##

# res_tmp <- run_dcdh(df_noe, y = "logwealth", d = "councilmemberterm",
#                     controls = c("ntaxpayments", "ntaxpayments2"),
#                     effects = 5, placebo = 3)
# saveRDS(res_tmp, file.path(out_dir, "fig3_logwealth_base.rds"))
# rm(res_tmp); gc()

## ---- 4.2 Panel A.II: logwealth, Baseline + Controls (NOT EXECUTED) ------ ##

# WARNING: Exceeds 16 GB RAM. Preserved for higher-memory hardware.
# res_tmp <- run_dcdh(df_noe, y = "logwealth", d = "councilmemberterm",
#                     controls = c("ntaxpayments", "ntaxpayments2", occ_vars, "womentaxpayer"),
#                     effects = 5, placebo = 3)
# saveRDS(res_tmp, file.path(out_dir, "fig3_logwealth_ctrl.rds"))
# rm(res_tmp); gc()

## ---- 4.3 Panel B.I: wealthpercentile, Baseline -------------------------- ##

# res_tmp <- run_dcdh(df_noe, y = "wealthpercentile", d = "councilmemberterm",
#                     controls = c("ntaxpayments", "ntaxpayments2"),
#                     effects = 5, placebo = 3)
# saveRDS(res_tmp, file.path(out_dir, "fig3_pct_base.rds"))
# rm(res_tmp); gc()

## ---- 4.4 Panel B.II: wealthpercentile, Baseline + Controls (NOT EXECUTED) ##

# WARNING: Exceeds 16 GB RAM. Preserved for higher-memory hardware.
# res_tmp <- run_dcdh(df_noe, y = "wealthpercentile", d = "councilmemberterm",
#                     controls = c("ntaxpayments", "ntaxpayments2", occ_vars, "womentaxpayer"),
#                     effects = 5, placebo = 3)
# saveRDS(res_tmp, file.path(out_dir, "fig3_pct_ctrl.rds"))
# rm(res_tmp); gc()

## ---- 4.5 Panel C.I: top5percent, Baseline ------------------------------- ##

# res_tmp <- run_dcdh(df_noe, y = "top5percent", d = "councilmemberterm",
#                     controls = c("ntaxpayments", "ntaxpayments2"),
#                     effects = 5, placebo = 3)
# saveRDS(res_tmp, file.path(out_dir, "fig3_top5_base.rds"))
# rm(res_tmp); gc()

## ---- 4.6 Panel C.II: top5percent, Baseline + Controls (NOT EXECUTED) ---- ##

# WARNING: Exceeds 16 GB RAM. Preserved for higher-memory hardware.
# res_tmp <- run_dcdh(df_noe, y = "top5percent", d = "councilmemberterm",
#                     controls = c("ntaxpayments", "ntaxpayments2", occ_vars, "womentaxpayer"),
#                     effects = 5, placebo = 3)
# saveRDS(res_tmp, file.path(out_dir, "fig3_top5_ctrl.rds"))
# rm(res_tmp); gc()

# =============================================================================
# 5. TABLE 4 - MECHANISM TESTS
# =============================================================================
# 8 specifications testing the candidate channels of magistrate enrichment:
# dynastic inheritance, business opportunities (merchants), and patronage
# (city clerks). Outcome: ln-wealth.
#
# All blocks commented out by default. The .rds files in output/ are already
# present from the original run.
# -----------------------------------------------------------------------------

## ---- 5.1 Estimate each spec ---------------------------------------------- ##

# Spec 1/8 (dCdH_t4_spec_1.rds): counciltermXcouncilinherit, Baseline
# res_tmp <- run_dcdh(df_noe, y = "logwealth", d = "counciltermXcouncilinherit",
#                     controls = c("councilmemberterm", "ntaxpayments", "ntaxpayments2"),
#                     effects = 1, nboot = 100)
# saveRDS(res_tmp, file.path(out_dir, "dCdH_t4_spec_1.rds"))
# rm(res_tmp); gc()

# Spec 2/8 (dCdH_t4_spec_2.rds): counciltermXcouncilinherit, Baseline + Controls
# res_tmp <- run_dcdh(df_noe, y = "logwealth", d = "counciltermXcouncilinherit",
#                     controls = c("councilmemberterm", occ_vars, "womentaxpayer",
#                                  "ntaxpayments", "ntaxpayments2"),
#                     effects = 1, nboot = 100)
# saveRDS(res_tmp, file.path(out_dir, "dCdH_t4_spec_2.rds"))
# rm(res_tmp); gc()

# Spec 3/8 (dCdH_t4_spec_3.rds): councilmemberterm + inherit, Baseline
# res_tmp <- run_dcdh(df_noe, y = "logwealth", d = "councilmemberterm",
#                     controls = c("counciltermXcouncilinherit", "ntaxpayments", "ntaxpayments2"),
#                     effects = 1, nboot = 100)
# saveRDS(res_tmp, file.path(out_dir, "dCdH_t4_spec_3.rds"))
# rm(res_tmp); gc()

# Spec 4/8 (dCdH_t4_spec_4.rds): councilmemberterm + inherit, Baseline + Controls
# res_tmp <- run_dcdh(df_noe, y = "logwealth", d = "councilmemberterm",
#                     controls = c("counciltermXcouncilinherit", occ_vars, "womentaxpayer",
#                                  "ntaxpayments", "ntaxpayments2"),
#                     effects = 1, nboot = 100)
# saveRDS(res_tmp, file.path(out_dir, "dCdH_t4_spec_4.rds"))
# rm(res_tmp); gc()

# Spec 5/8 (dCdH_t4_spec_5.rds): counciltermXmerchant, Baseline
# res_tmp <- run_dcdh(df_noe, y = "logwealth", d = "counciltermXmerchant",
#                     controls = c("councilmemberterm", "occupation_500",
#                                  "ntaxpayments", "ntaxpayments2"),
#                     effects = 1, nboot = 100)
# saveRDS(res_tmp, file.path(out_dir, "dCdH_t4_spec_5.rds"))
# rm(res_tmp); gc()

# Spec 6/8 (dCdH_t4_spec_6.rds): counciltermXmerchant, Baseline + Controls
# res_tmp <- run_dcdh(df_noe, y = "logwealth", d = "counciltermXmerchant",
#                     controls = c("councilmemberterm", occ_vars, "womentaxpayer",
#                                  "ntaxpayments", "ntaxpayments2"),
#                     effects = 1, nboot = 100)
# saveRDS(res_tmp, file.path(out_dir, "dCdH_t4_spec_6.rds"))
# rm(res_tmp); gc()

# Spec 7/8 (dCdH_t4_spec_7.rds): cityadmin_f, Baseline
# res_tmp <- run_dcdh(df_noe, y = "logwealth", d = "cityadmin_f",
#                     controls = c("ntaxpayments", "ntaxpayments2"),
#                     effects = 1, nboot = 100)
# saveRDS(res_tmp, file.path(out_dir, "dCdH_t4_spec_7.rds"))
# rm(res_tmp); gc()

# Spec 8/8 (dCdH_t4_spec_8.rds): cityadmin_f, Baseline + Controls
# res_tmp <- run_dcdh(df_noe, y = "logwealth", d = "cityadmin_f",
#                     controls = c(occ_vars_no_clerks, "womentaxpayer",
#                                  "ntaxpayments", "ntaxpayments2"),
#                     effects = 1, nboot = 100)
# saveRDS(res_tmp, file.path(out_dir, "dCdH_t4_spec_8.rds"))
# rm(res_tmp); gc()


# =============================================================================
# 6. FIGURE 4 - EVENT-STUDY CITY CLERKS
# =============================================================================
# Baseline panel of the clerk treatment effect dynamics. The "+ Controls"
# panel was not executed (exceeds 16 GB RAM).
# -----------------------------------------------------------------------------

## ---- 6.1 Panel Baseline -------------------------------------------------- ##

# res_tmp <- run_dcdh(df_noe, y = "logwealth", d = "cityadmin_f",
#                     controls = c("ntaxpayments", "ntaxpayments2"),
#                     effects = 5, placebo = 3)
# saveRDS(res_tmp, file.path(out_dir, "fig4_base.rds"))
# rm(res_tmp); gc()

## ---- 6.2 Panel Baseline + Controls (NOT EXECUTED) ----------------------- ##

# WARNING: Exceeds 16 GB RAM. Preserved for higher-memory hardware.
# res_tmp <- run_dcdh(df_noe, y = "logwealth", d = "cityadmin_f",
#                     controls = c(occ_vars_no_clerks, "womentaxpayer",
#                                  "ntaxpayments", "ntaxpayments2"),
#                     effects = 5, placebo = 3)
# saveRDS(res_tmp, file.path(out_dir, "fig4_ctrl.rds"))
# rm(res_tmp); gc()


# =============================================================================
# 7. TABLE 5 - TRIPLE-DIFFERENCE (THIRTY YEARS' WAR)
# =============================================================================
# 7 specifications of the triple-difference Council × Post × War on the
# restricted sample 1603-1646. Tests whether wealth gains from council
# membership were amplified during the war.
# -----------------------------------------------------------------------------

## ---- 7.1 Estimate each spec ---------------------------------------------- ##

# Spec 1/7 (dCdH_t5_spec_1.rds): logwealth, Only FE
# res_tmp <- run_dcdh(df_30yw, y = "logwealth", d = "counciltermXpost1618",
#                     controls = c("councilmemberterm"),
#                     effects = 1, nboot = 100)
# saveRDS(res_tmp, file.path(out_dir, "dCdH_t5_spec_1.rds"))
# rm(res_tmp); gc()

# Spec 2/7 (dCdH_t5_spec_2.rds): logwealth, Baseline
# res_tmp <- run_dcdh(df_30yw, y = "logwealth", d = "counciltermXpost1618",
#                     controls = c("councilmemberterm", "ntaxpayments", "ntaxpayments2"),
#                     effects = 1, nboot = 100)
# saveRDS(res_tmp, file.path(out_dir, "dCdH_t5_spec_2.rds"))
# rm(res_tmp); gc()

# Spec 3/7 (dCdH_t5_spec_3.rds): logwealth, Baseline + Controls
# res_tmp <- run_dcdh(df_30yw, y = "logwealth", d = "counciltermXpost1618",
#                     controls = c(occ_vars, "womentaxpayer", "councilmemberterm",
#                                  "ntaxpayments", "ntaxpayments2"),
#                     effects = 1, nboot = 100)
# saveRDS(res_tmp, file.path(out_dir, "dCdH_t5_spec_3.rds"))
# rm(res_tmp); gc()

# Spec 4/7 (dCdH_t5_spec_4.rds): wealthpercentile, Baseline
# res_tmp <- run_dcdh(df_30yw, y = "wealthpercentile", d = "counciltermXpost1618",
#                     controls = c("councilmemberterm", "ntaxpayments", "ntaxpayments2"),
#                     effects = 1, nboot = 100)
# saveRDS(res_tmp, file.path(out_dir, "dCdH_t5_spec_4.rds"))
# rm(res_tmp); gc()

# Spec 5/7 (dCdH_t5_spec_5.rds): wealthpercentile, Baseline + Controls
# res_tmp <- run_dcdh(df_30yw, y = "wealthpercentile", d = "counciltermXpost1618",
#                     controls = c(occ_vars, "womentaxpayer", "councilmemberterm",
#                                  "ntaxpayments", "ntaxpayments2"),
#                     effects = 1, nboot = 100)
# saveRDS(res_tmp, file.path(out_dir, "dCdH_t5_spec_5.rds"))
# rm(res_tmp); gc()

# Spec 6/7 (dCdH_t5_spec_6.rds): top5percent, Baseline
# res_tmp <- run_dcdh(df_30yw, y = "top5percent", d = "counciltermXpost1618",
#                     controls = c("councilmemberterm", "ntaxpayments", "ntaxpayments2"),
#                     effects = 1, nboot = 100)
# saveRDS(res_tmp, file.path(out_dir, "dCdH_t5_spec_6.rds"))
# rm(res_tmp); gc()

# Spec 7/7 (dCdH_t5_spec_7.rds): top5percent, Baseline + Controls
# res_tmp <- run_dcdh(df_30yw, y = "top5percent", d = "counciltermXpost1618",
#                     controls = c(occ_vars, "womentaxpayer", "councilmemberterm",
#                                  "ntaxpayments", "ntaxpayments2"),
#                     effects = 1, nboot = 100)
# saveRDS(res_tmp, file.path(out_dir, "dCdH_t5_spec_7.rds"))
# rm(res_tmp); gc()


# =============================================================================
# 8. FIGURE 5 - EVENT-STUDY TRIPLE-DIFFERENCE
# =============================================================================
# 2 Baseline panels (logwealth, wealthpercentile) of the Council × Post × War
# interaction. "+ Controls" panels not executed (exceed 16 GB RAM).
# -----------------------------------------------------------------------------

## ---- 8.1 Panel A.I: logwealth, Baseline --------------------------------- ##

# res_tmp <- run_dcdh(df_30yw, y = "logwealth", d = "counciltermXpost1618",
#                     controls = c("councilmemberterm", "ntaxpayments", "ntaxpayments2"),
#                     effects = 5, placebo = 3)
# saveRDS(res_tmp, file.path(out_dir, "fig5_logw_base.rds"))
# rm(res_tmp); gc()

## ---- 8.2 Panel A.II: logwealth, Baseline + Controls (NOT EXECUTED) ------ ##

# WARNING: Exceeds 16 GB RAM. Preserved for higher-memory hardware.
# res_tmp <- run_dcdh(df_30yw, y = "logwealth", d = "counciltermXpost1618",
#                     controls = c(occ_vars, "womentaxpayer", "councilmemberterm",
#                                  "ntaxpayments", "ntaxpayments2"),
#                     effects = 5, placebo = 3)
# saveRDS(res_tmp, file.path(out_dir, "fig5_logw_ctrl.rds"))
# rm(res_tmp); gc()

## ---- 8.3 Panel B.I: wealthpercentile, Baseline -------------------------- ##

# res_tmp <- run_dcdh(df_30yw, y = "wealthpercentile", d = "counciltermXpost1618",
#                     controls = c("councilmemberterm", "ntaxpayments", "ntaxpayments2"),
#                     effects = 5, placebo = 3)
# saveRDS(res_tmp, file.path(out_dir, "fig5_pct_base.rds"))
# rm(res_tmp); gc()

## ---- 8.4 Panel B.II: wealthpercentile, Baseline + Controls (NOT EXECUTED) ##

# WARNING: Exceeds 16 GB RAM. Preserved for higher-memory hardware.
# res_tmp <- run_dcdh(df_30yw, y = "wealthpercentile", d = "counciltermXpost1618",
#                     controls = c(occ_vars, "womentaxpayer", "councilmemberterm",
#                                  "ntaxpayments", "ntaxpayments2"),
#                     effects = 5, placebo = 3)
# saveRDS(res_tmp, file.path(out_dir, "fig5_pct_ctrl.rds"))
# rm(res_tmp); gc()


# =============================================================================
# 9. FIGURE 6 - EVENT-STUDY CITY CLERKS DURING THE WAR
# =============================================================================
# Baseline panel of the clerk × war interaction. "+ Controls" panel not
# executed (exceeds 16 GB RAM).
# -----------------------------------------------------------------------------

## ---- 9.1 Panel Baseline -------------------------------------------------- ##

# res_tmp <- run_dcdh(df_30yw, y = "logwealth", d = "cityadmin_fXpost1618",
#                     controls = c("cityadmin_f", "ntaxpayments", "ntaxpayments2"),
#                     effects = 4, placebo = 2)
# saveRDS(res_tmp, file.path(out_dir, "fig6_base.rds"))
# rm(res_tmp); gc()

## ---- 9.2 Panel Baseline + Controls (NOT EXECUTED) ----------------------- ##

# WARNING: Exceeds 16 GB RAM. Preserved for higher-memory hardware.
# res_tmp <- run_dcdh(df_30yw, y = "logwealth", d = "cityadmin_fXpost1618",
#                     controls = c("cityadmin_f", occ_vars_no_clerks,
#                                  "womentaxpayer", "ntaxpayments", "ntaxpayments2"),
#                     effects = 4, placebo = 2)
# saveRDS(res_tmp, file.path(out_dir, "fig6_ctrl.rds"))
# rm(res_tmp); gc()


# =============================================================================
# END OF SCRIPT
# =============================================================================