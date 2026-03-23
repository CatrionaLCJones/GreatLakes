# ============================================================================
# Statistical Analysis and Plotting
# Lake Erie Central Basin Ecosystem Metabolism Study
#
# Produces all figures and statistical results reported in the manuscript
# and supplementary information.
#
# Figures produced:
#   Fig 3   - Temporal trends in GPP, R, and NEP by position and year
#   Fig 4   - Environmental variables (cDOM, Chla, Kd, SM, Temp, Shortwave)
#             by month and position
#   Fig 5   - GPP vs environmental drivers (scatter + mixed-effects fits)
#   Fig 6   - Cumulative carbon fixation vs release by site and year
#   Fig 7   - Carbon fixation time series with CyAN HABs maps
#   Fig 8   - Carbon release time series with CyAN HABs maps
#   Fig S1  - Full vs Submerged dataset comparison (code not recovered from
#             provided scripts; see Table S2 for model statistics)
#   Fig S2  - Respiration vs environmental drivers (scatter + mixed-effects
#             fits; supplementary companion to Fig 5)
#
# Statistical results in manuscript / supplementary:
#   Results section: position * month LME models for GPP, R, NEP
#   Table S2:        Full vs Submerged comparison LME models
#   Table S3:        Environmental driver LME models (GPP, R, NEP ~ each
#                    driver +/- position interaction)
# ============================================================================

# --- Libraries ---------------------------------------------------------------
library(ggplot2)
library(ggpubr)
library(patchwork)
library(dplyr)
library(lubridate)
library(tidyr)
library(reshape2)
library(viridis)
library(lme4)
library(lmerTest)
library(emmeans)
library(performance)
library(car)
library(broom)
library(gridExtra)
library(ggeffects)
library(terra)
library(fields)
library(png)
library(doBy)

mean_na <- function(x) mean(x, na.rm = TRUE)
se      <- function(x) sqrt(var(x, na.rm = TRUE) / sum(!is.na(x)))

# ============================================================================
# PART 1: Load and prepare data
# ============================================================================

# Primary compiled metabolism dataset (produced by the LakeMetabolizer
# modelling step, not included here; see Winslow et al. 2016).
Metab <- readRDS("lakeerie_dataset_20250916.rds")

# --- Merge daily satellite water-quality data --------------------------------
# Load interpolated satellite products (produced by
# satellite_data_download_assembly.R).
chla_interp <- readRDS("/Users/jdh/lakeerie/data/lake_erie_chla_satellite.rds")
cdom_interp <- readRDS("/Users/jdh/lakeerie/data/lake_erie_cdom_satellite.rds")
kd_interp   <- readRDS("/Users/jdh/lakeerie/data/lake_erie_kd_satellite.rds")
sm_interp   <- readRDS("/Users/jdh/lakeerie/data/lake_erie_sm_satellite.rds")

Metab$date <- as.Date(Metab$date)

# Aggregate satellite data to daily site means
chla_daily <- summaryBy(
  chlor_a_nn + chlor_a ~ site + date,
  data = chla_interp[!is.na(chla_interp$chlor_a_nn), ],
  FUN  = mean_na
)
cdom_daily <- summaryBy(
  cdom + nn_mean ~ site + date,
  data = cdom_interp[!is.na(cdom_interp$nn_mean_interpolated), ],
  FUN  = mean_na
)
kd_daily <- summaryBy(
  kd + kd_nn_mean ~ site + date,
  data = kd_interp[!is.na(kd_interp$kd_nn_mean), ],
  FUN  = mean_na
)
sm_daily <- summaryBy(
  sm + sm_nn_mean ~ site + date,
  data = sm_interp[!is.na(sm_interp$sm_nn_mean), ],
  FUN  = mean_na
)

# Merge satellite products into the metabolism dataset
Metab <- Metab %>%
  left_join(chla_daily[, c("date", "site", "chlor_a_nn.mean_na", "chlor_a.mean_na")],
            by = c("site", "date")) %>%
  left_join(cdom_daily[, c("date", "site", "cdom.mean_na", "nn_mean.mean_na")],
            by = c("site", "date")) %>%
  left_join(kd_daily[, c("date", "site", "kd.mean_na", "kd_nn_mean.mean_na")],
            by = c("site", "date")) %>%
  left_join(sm_daily[, c("date", "site", "sm.mean_na", "sm_nn_mean.mean_na")],
            by = c("site", "date"))

# --- Derived variables -------------------------------------------------------
Metab <- Metab %>%
  mutate(
    month_name = month(date, label = TRUE, abbr = TRUE),
    year       = year(date)
  )

# Site classification
offshore_sites  <- c("CHRP2", "CHRP4", "CHRP9")
nearshore_sites <- c("CHRP1", "CHRP3", "CHRP5", "CHRP7", "CHRP8")

Metab <- Metab %>%
  mutate(
    position = case_when(
      site %in% nearshore_sites ~ "Nearshore",
      site %in% offshore_sites  ~ "Offshore",
      TRUE ~ NA_character_
    )
  )

# Log-modulus transformation: sign(x) * log10(1 + |x|)
# Applied to GPP, R, and NEP for mixed-effects model fitting
# (Whittaker et al. 2005; preserves sign of metabolism values).
log_modulus <- function(x, base = 10) {
  sign(x) * log(1 + abs(x), base)
}

Metab$R_logmod   <- log_modulus(Metab$R)
Metab$GPP_logmod <- log_modulus(Metab$GPP)
Metab$NEP_logmod <- log_modulus(Metab$NEP)

# Spatial (depth-integrated) estimates: volumetric rate * mixing depth (m)
Metab$R_sp   <- Metab$R   * Metab$mean_mix
Metab$GPP_sp <- Metab$GPP * Metab$mean_mix
Metab$NEP_sp <- Metab$NEP * Metab$mean_mix

Metab$R_sp_logmod   <- log_modulus(Metab$R_sp)
Metab$GPP_sp_logmod <- log_modulus(Metab$GPP_sp)
Metab$NEP_sp_logmod <- log_modulus(Metab$NEP_sp)

# Restrict statistical models to May-September due to limited October data
# (nearshore n = 6, offshore n = 15 in October).
Metab_oct <- Metab %>% filter(month_name != "Oct")

# Helper: extract key statistics from an lmer model
extr_lmer <- function(x) {
  lmer_formula <- as.character(x@call)[2]
  pvalue  <- anova(x)[[6]]
  numDF   <- anova(x)[[3]]
  denDF   <- anova(x)[[4]]
  Fvalue  <- anova(x)[[5]]
  cond_r2 <- r2_nakagawa(x)[[1]]
  marg_r2 <- r2_nakagawa(x)[[2]]
  print(summary(x))
  print(anova(x))
  data.frame(lmer_formula, pvalue, numDF, denDF, Fvalue, cond_r2, marg_r2)
}

# ============================================================================
# PART 2: Overall temporal and spatial trends in ecosystem metabolism
# (Results: "Overall temporal and spatial trends"; Fig 3)
# Models: GPP_logmod / R_logmod / NEP_logmod ~ position * month_name + (1|site)
# Restricted to May-September (n = 2,291 observations)
# ============================================================================

Position_month_GPP_model <- lmer(GPP_logmod ~ position * month_name + (1|site),
                                 data = Metab_oct)
anova(Position_month_GPP_model)
GPP_pos_emmeans <- emmeans(Position_month_GPP_model, ~ position * month_name)
contrast(GPP_pos_emmeans, "tukey")
extr_lmer(Position_month_GPP_model)

Position_month_R_model <- lmer(R_logmod ~ position * month_name + (1|site),
                                data = Metab_oct)
anova(Position_month_R_model)
R_pos_emmeans <- emmeans(Position_month_R_model, ~ position * month_name)
contrast(R_pos_emmeans, "tukey")
extr_lmer(Position_month_R_model)

Position_month_NEP_model <- lmer(NEP_logmod ~ position * month_name + (1|site),
                                  data = Metab_oct)
anova(Position_month_NEP_model)
NEP_pos_emmeans <- emmeans(Position_month_NEP_model, ~ position * month_name)
contrast(NEP_pos_emmeans, "tukey")
extr_lmer(Position_month_NEP_model)

# ============================================================================
# PART 3: Figure 3 - Monthly mean GPP, R, and NEP by position and year
# ============================================================================

# --- Prepare plotting data ---------------------------------------------------
top_plot_data <- Metab %>%
  select(year, month_name, position, GPP, R) %>%
  pivot_longer(cols = c(GPP, R), names_to = "variable", values_to = "value") %>%
  mutate(variable = recode(variable, GPP = "GPP_vol", R = "R_vol")) %>%
  group_by(year, month_name, position, variable) %>%
  summarise(
    mean_value = mean(value,  na.rm = TRUE),
    se_value   = sd(value,   na.rm = TRUE) / sqrt(n()),
    .groups    = "drop"
  )

bottom_plot_data <- Metab %>%
  select(year, month_name, position, NEP) %>%
  group_by(year, month_name, position) %>%
  summarise(
    mean_value = mean(NEP, na.rm = TRUE),
    se_value   = sd(NEP,  na.rm = TRUE) / sqrt(n()),
    .groups    = "drop"
  )

# --- Top panel: GPP and R ----------------------------------------------------
top_plot <- top_plot_data %>%
  ggplot(aes(x = month_name, y = mean_value,
             color = position, shape = variable)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_line(aes(group = interaction(position, variable)), alpha = 0.7) +
  geom_errorbar(aes(ymin = mean_value - se_value,
                    ymax = mean_value + se_value),
                width = 0.2, alpha = 0.7) +
  geom_point(size = 2.5) +
  facet_wrap(~ year, ncol = 1) +
  scale_color_manual(values = c(Nearshore = "#7b3294", Offshore = "#fdae61")) +
  scale_shape_manual(
    name   = "Metabolism\ncomponent",
    values = c(GPP_vol = 17, R_vol = 16),
    labels = c(GPP_vol = "GPP", R_vol = "R")
  ) +
  labs(x = "",
       y = expression("Ecosystem Metabolism (mg-O"[2]~"L"^-1~"day"^-1*")")) +
  theme_bw() +
  theme(panel.grid       = element_blank(),
        strip.background = element_rect(fill = "white"),
        legend.position  = "right",
        axis.text.x      = element_blank())

# --- Bottom panel: NEP -------------------------------------------------------
bottom_plot <- bottom_plot_data %>%
  ggplot(aes(x = month_name, y = mean_value, color = position)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_line(aes(group = position), alpha = 0.7) +
  geom_errorbar(aes(ymin = mean_value - se_value,
                    ymax = mean_value + se_value),
                width = 0.2, alpha = 0.7) +
  geom_point(size = 2.5, shape = 15) +
  facet_wrap(~ year, ncol = 1) +
  scale_color_manual(values = c(Nearshore = "#7b3294", Offshore = "#fdae61")) +
  labs(x = "",
       y = expression("NEP (mg-O"[2]~"L"^-1~"day"^-1*")")) +
  theme_bw() +
  theme(panel.grid       = element_blank(),
        strip.background = element_rect(fill = "white"),
        legend.position  = "right")

# Combine panels and export (width = 600, height = 800)
ggpubr::ggarrange(top_plot, bottom_plot, ncol = 1)

# ============================================================================
# PART 4: Figure 4 - Environmental variables by month and position
# ============================================================================

monthly_summary <- Metab %>%
  group_by(month_name, position) %>%
  summarise(
    mean_cdom  = mean_na(nn_mean.mean_na),
    se_cdom    = se(nn_mean.mean_na),
    mean_chla  = mean_na(chlor_a_nn.mean_na),
    se_chla    = se(chlor_a_nn.mean_na),
    mean_kd    = mean_na(kd_nn_mean.mean_na),
    se_kd      = se(kd_nn_mean.mean_na),
    mean_sm    = mean_na(sm_nn_mean.mean_na),
    se_sm      = se(sm_nn_mean.mean_na),
    mean_temp  = mean_na(mean_ts),
    se_temp    = se(mean_ts),
    mean_insol = mean_na(mean_swnet),
    se_insol   = se(mean_swnet),
    .groups    = "drop"
  ) %>%
  filter(month_name %in% c("May", "Jun", "Jul", "Aug", "Sep", "Oct"))

clr <- c(Nearshore = "#7b3294", Offshore = "#fdae61")

p_cdom <- ggplot(monthly_summary,
                  aes(x = month_name, y = mean_cdom,
                      color = position, group = position)) +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_cdom - se_cdom,
                    ymax = mean_cdom + se_cdom), width = 0.2) +
  scale_color_manual(values = clr) +
  labs(x = NULL, y = expression(paste("cDOM absorption (m"^-1*")")),
       color = "Position", tag = "A") +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none",
        axis.text.x = element_blank(),
        plot.tag = element_text(size = 12, face = "bold"))

p_chla <- ggplot(monthly_summary,
                  aes(x = month_name, y = mean_chla,
                      color = position, group = position)) +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_chla - se_chla,
                    ymax = mean_chla + se_chla), width = 0.2) +
  scale_color_manual(values = clr) +
  labs(x = NULL,
       y = expression(paste("Chlorophyll-", italic("a"), " (\u03bcg L"^-1*")")),
       color = "Position", tag = "B") +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none",
        axis.text.x = element_blank(),
        plot.tag = element_text(size = 12, face = "bold"))

p_kd <- ggplot(monthly_summary,
                aes(x = month_name, y = mean_kd,
                    color = position, group = position)) +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_kd - se_kd, ymax = mean_kd + se_kd),
                width = 0.2) +
  scale_color_manual(values = clr) +
  labs(x = NULL, y = expression(paste(italic("K")[d], " (m"^-1*")")),
       color = "Position", tag = "C") +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none",
        axis.text.x = element_blank(),
        plot.tag = element_text(size = 12, face = "bold"))

p_sm <- ggplot(monthly_summary,
                aes(x = month_name, y = mean_sm,
                    color = position, group = position)) +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_sm - se_sm, ymax = mean_sm + se_sm),
                width = 0.2) +
  scale_color_manual(values = clr) +
  labs(x = NULL, y = expression(paste("Suspended Minerals (mg L"^-1*")")),
       color = "Position", tag = "D") +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none",
        axis.text.x = element_blank(),
        plot.tag = element_text(size = 12, face = "bold"))

p_temp <- ggplot(monthly_summary,
                  aes(x = month_name, y = mean_temp,
                      color = position, group = position)) +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_temp - se_temp,
                    ymax = mean_temp + se_temp), width = 0.2) +
  scale_color_manual(values = clr) +
  labs(x = "Month", y = expression("Temperature (\u00b0C)"),
       color = "Position", tag = "E") +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none",
        plot.tag = element_text(size = 12, face = "bold"))

p_insol <- ggplot(monthly_summary,
                   aes(x = month_name, y = mean_insol,
                       color = position, group = position)) +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_insol - se_insol,
                    ymax = mean_insol + se_insol), width = 0.2) +
  scale_color_manual(values = clr) +
  labs(x = "Month", y = expression("Net Shortwave Radiation (W m"^-2*")"),
       color = "Position", tag = "F") +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none",
        plot.tag = element_text(size = 12, face = "bold"))

ggpubr::ggarrange(p_cdom, p_chla, p_kd, p_sm, p_temp, p_insol,
                  ncol = 2, nrow = 3, common.legend = TRUE)

# Simple LM: nearshore vs offshore differences in mean environmental conditions
# (values cited in Results section)
anova(lm(chlor_a_nn.mean_na ~ position, data = Metab))
anova(lm(nn_mean.mean_na    ~ position, data = Metab))
anova(lm(kd_nn_mean.mean_na ~ position, data = Metab))
anova(lm(sm_nn_mean.mean_na ~ position, data = Metab))

# ============================================================================
# PART 5: Environmental driver models (Table S3) and Figure 5 / Figure S2
# Models: GPP_logmod / R_logmod / NEP_logmod ~ log10(driver) [* position]
# Site as random intercept; interaction included initially, dropped if n.s.
# ============================================================================

# --- GPP models (Table S3, Fig 5) -------------------------------------------
# Temperature x Position interaction
gpp_temp_lmer <- lmer(GPP_logmod ~ log10(mean_ts) * position + (1|site),
                       data = Metab_oct)
extr_lmer(gpp_temp_lmer)

# cDOM + Position (no interaction; Table S3)
gpp_cdom_lmer <- lmer(GPP_logmod ~ log10(nn_mean.mean_na) + position + (1|site),
                       data = Metab_oct)
extr_lmer(gpp_cdom_lmer)

# Suspended minerals x Position interaction
gpp_sm_lmer <- lmer(GPP_logmod ~ log10(sm_nn_mean.mean_na) * position + (1|site),
                     data = Metab_oct)
extr_lmer(gpp_sm_lmer)

# Kd x Position interaction
gpp_kd_lmer <- lmer(GPP_logmod ~ log10(kd_nn_mean.mean_na) * position + (1|site),
                     data = Metab_oct)
extr_lmer(gpp_kd_lmer)

# Chlorophyll-a (no position effect; Table S3)
gpp_chla_lmer <- lmer(GPP_logmod ~ log10(chlor_a_nn.mean_na) + (1|site),
                       data = Metab_oct)
extr_lmer(gpp_chla_lmer)

# Net shortwave radiation (no significant effect; Table S3)
gpp_sw_lmer <- lmer(GPP_logmod ~ log10(mean_swnet) + (1|site),
                     data = Metab_oct)
extr_lmer(gpp_sw_lmer)

# --- Respiration models (Table S3, Fig S2) -----------------------------------
r_temp_lmer <- lmer(R_logmod ~ log10(mean_ts) * position + (1|site),
                     data = Metab_oct)
extr_lmer(r_temp_lmer)

r_cdom_lmer <- lmer(R_logmod ~ log10(nn_mean.mean_na) * position + (1|site),
                     data = Metab_oct)
extr_lmer(r_cdom_lmer)

r_sm_lmer <- lmer(R_logmod ~ log10(sm_nn_mean.mean_na) * position + (1|site),
                   data = Metab_oct)
extr_lmer(r_sm_lmer)

r_kd_lmer <- lmer(R_logmod ~ log10(kd_nn_mean.mean_na) * position + (1|site),
                   data = Metab_oct)
extr_lmer(r_kd_lmer)

r_chla_lmer <- lmer(R_logmod ~ log10(chlor_a_nn.mean_na) + (1|site),
                     data = Metab_oct)
extr_lmer(r_chla_lmer)

r_sw_lmer <- lmer(R_logmod ~ log10(mean_swnet) + (1|site),
                   data = Metab_oct)
extr_lmer(r_sw_lmer)

# --- NEP models (Table S3) ---------------------------------------------------
nep_temp_lmer <- lmer(NEP_logmod ~ log10(mean_ts)          + (1|site), data = Metab_oct)
nep_cdom_lmer <- lmer(NEP_logmod ~ log10(nn_mean.mean_na)   + (1|site), data = Metab_oct)
nep_sm_lmer   <- lmer(NEP_logmod ~ log10(sm_nn_mean.mean_na)+ (1|site), data = Metab_oct)
nep_kd_lmer   <- lmer(NEP_logmod ~ log10(kd_nn_mean.mean_na)+ (1|site), data = Metab_oct)
nep_chla_lmer <- lmer(NEP_logmod ~ log10(chlor_a_nn.mean_na)+ (1|site), data = Metab_oct)
nep_sw_lmer   <- lmer(NEP_logmod ~ log10(mean_swnet)         + (1|site), data = Metab_oct)

extr_lmer(nep_temp_lmer); extr_lmer(nep_cdom_lmer)
extr_lmer(nep_sm_lmer);   extr_lmer(nep_kd_lmer)
extr_lmer(nep_chla_lmer); extr_lmer(nep_sw_lmer)

# --- Figure 5: GPP vs environmental drivers ----------------------------------
# Monthly means used as plotting data; mixed-effects model predictions overlaid.
# Prediction bands are marginal (fixed effects only; no random effects).

cmc_month <- Metab %>%
  group_by(site, year, month_name, position) %>%
  summarise(across(c(GPP, R, NEP,
                     chlor_a_nn.mean_na, nn_mean.mean_na,
                     kd_nn_mean.mean_na, sm_nn_mean.mean_na,
                     mean_ts, mean_swnet),
                   mean_na),
            .groups = "drop")

# Helper: build prediction data and SE for a fitted lmer
make_pred <- function(model, x_var, x_range, group_var = NULL,
                      group_levels = NULL) {
  if (!is.null(group_var)) {
    pred <- expand.grid(
      x      = seq(x_range[1], x_range[2], length.out = 100),
      group  = group_levels
    )
    names(pred) <- c(x_var, group_var)
    mm <- model.matrix(reformulate(paste0("log10(", x_var, ") *", group_var)),
                       data = pred)
  } else {
    pred <- data.frame(x = seq(x_range[1], x_range[2], length.out = 100))
    names(pred) <- x_var
    mm <- model.matrix(reformulate(paste0("log10(", x_var, ")")), data = pred)
  }
  pred$fit <- as.vector(mm %*% fixef(model))
  pred$se  <- sqrt(diag(mm %*% vcov(model) %*% t(mm)))
  pred$lwr <- pred$fit - 1.96 * pred$se
  pred$upr <- pred$fit + 1.96 * pred$se
  pred
}

# cDOM (panel a) - no interaction, separate position offsets
pred_cdom <- make_pred(gpp_cdom_lmer, "nn_mean.mean_na",
                       range(cmc_month$nn_mean.mean_na, na.rm = TRUE),
                       group_var = "position",
                       group_levels = c("Nearshore", "Offshore"))

cdom_gpp_plot <- ggplot(cmc_month,
                         aes(x = nn_mean.mean_na, y = GPP,
                             color = position)) +
  theme_classic() + theme(legend.position = "none") +
  scale_x_log10() +
  geom_ribbon(data = pred_cdom,
              aes(x = nn_mean.mean_na, ymin = lwr, ymax = upr, fill = position),
              alpha = 0.2, color = NA, inherit.aes = FALSE) +
  geom_line(data = pred_cdom,
            aes(x = nn_mean.mean_na, y = fit, color = position),
            linewidth = 1, inherit.aes = FALSE) +
  geom_point(size = 3) +
  scale_color_manual(values = clr) +
  scale_fill_manual(values  = clr) +
  xlab(expression("cDOM Absorbance (m"^-1*")")) +
  ylab(expression("GPP (g-m"^-2~d^-1*")"))

# Suspended minerals (panel b) - interaction with position
pred_sm <- make_pred(gpp_sm_lmer, "sm_nn_mean.mean_na",
                     range(cmc_month$sm_nn_mean.mean_na, na.rm = TRUE),
                     group_var = "position",
                     group_levels = c("Nearshore", "Offshore"))

sm_gpp_plot <- ggplot(cmc_month,
                       aes(x = sm_nn_mean.mean_na, y = GPP,
                           color = position)) +
  theme_classic() + theme(legend.position = "none") +
  scale_x_log10() +
  geom_ribbon(data = pred_sm,
              aes(x = sm_nn_mean.mean_na, ymin = lwr, ymax = upr,
                  fill = position),
              alpha = 0.2, color = NA, inherit.aes = FALSE) +
  geom_line(data = pred_sm,
            aes(x = sm_nn_mean.mean_na, y = fit, color = position),
            linewidth = 1, inherit.aes = FALSE) +
  geom_point(size = 3) +
  scale_color_manual(values = clr) +
  scale_fill_manual(values  = clr) +
  xlab(expression("Suspended Matter (mg L"^-1*")")) +
  ylab(expression("GPP (g-m"^-2~d^-1*")"))

# Kd (panel c) - interaction with position
pred_kd <- make_pred(gpp_kd_lmer, "kd_nn_mean.mean_na",
                     range(cmc_month$kd_nn_mean.mean_na, na.rm = TRUE),
                     group_var = "position",
                     group_levels = c("Nearshore", "Offshore"))

kd_gpp_plot <- ggplot(cmc_month,
                       aes(x = kd_nn_mean.mean_na, y = GPP,
                           color = position)) +
  theme_classic() + theme(legend.position = "none") +
  scale_x_log10() +
  geom_ribbon(data = pred_kd,
              aes(x = kd_nn_mean.mean_na, ymin = lwr, ymax = upr,
                  fill = position),
              alpha = 0.2, color = NA, inherit.aes = FALSE) +
  geom_line(data = pred_kd,
            aes(x = kd_nn_mean.mean_na, y = fit, color = position),
            linewidth = 1, inherit.aes = FALSE) +
  geom_point(size = 3) +
  scale_color_manual(values = clr) +
  scale_fill_manual(values  = clr) +
  xlab(expression("K"[d]~"(m"^-1*")")) +
  ylab(expression("GPP (g-m"^-2~d^-1*")"))

# Chlorophyll-a (panel d) - no position effect; single pooled line
pred_chla <- make_pred(gpp_chla_lmer, "chlor_a_nn.mean_na",
                        range(cmc_month$chlor_a_nn.mean_na, na.rm = TRUE))

chla_gpp_plot <- ggplot(cmc_month,
                         aes(x = chlor_a_nn.mean_na, y = GPP,
                             color = position)) +
  theme_classic() + theme(legend.position = "none") +
  scale_x_log10() +
  geom_ribbon(data = pred_chla,
              aes(x = chlor_a_nn.mean_na, ymin = lwr, ymax = upr),
              alpha = 0.2, fill = "gray50", color = NA, inherit.aes = FALSE) +
  geom_line(data = pred_chla,
            aes(x = chlor_a_nn.mean_na, y = fit),
            linewidth = 1, color = "black", inherit.aes = FALSE) +
  geom_point(size = 3) +
  scale_color_manual(values = clr) +
  xlab(expression("Chlorophyll "*italic(a)~"(\u03bcg L"^-1*")")) +
  ylab(expression("GPP (g-m"^-2~d^-1*")"))

# Temperature (panel e) - interaction with position
pred_temp <- make_pred(gpp_temp_lmer, "mean_ts",
                       range(cmc_month$mean_ts, na.rm = TRUE),
                       group_var = "position",
                       group_levels = c("Nearshore", "Offshore"))

temp_gpp_plot <- ggplot(cmc_month,
                         aes(x = mean_ts, y = GPP, color = position)) +
  theme_classic() +
  scale_x_log10() +
  geom_ribbon(data = pred_temp,
              aes(x = mean_ts, ymin = lwr, ymax = upr, fill = position),
              alpha = 0.2, color = NA, inherit.aes = FALSE) +
  geom_line(data = pred_temp,
            aes(x = mean_ts, y = fit, color = position),
            linewidth = 1, inherit.aes = FALSE) +
  geom_point(size = 3) +
  scale_color_manual(values = clr) +
  scale_fill_manual(values  = clr) +
  xlab("Temperature (\u00b0C)") +
  ylab(expression("GPP (g-m"^-2~d^-1*")"))

# Figure 5 combined (save at ~9 x 12 in, 300 dpi)
combined_gpp_plot <- (cdom_gpp_plot + sm_gpp_plot) /
                     (kd_gpp_plot   + chla_gpp_plot) /
                     (temp_gpp_plot + plot_spacer()) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 14))

ggsave("~/gpp_plots_combined.png", combined_gpp_plot,
       width = 9, height = 12, dpi = 300)

# --- Figure S2: Respiration vs environmental drivers -------------------------
# Structured identically to Fig 5; code abbreviated for clarity.

pred_r_cdom <- make_pred(r_cdom_lmer, "nn_mean.mean_na",
                         range(cmc_month$nn_mean.mean_na, na.rm = TRUE),
                         "position", c("Nearshore", "Offshore"))
cdom_r_plot <- ggplot(cmc_month,
                       aes(x = nn_mean.mean_na, y = R, color = position)) +
  theme_classic() + theme(legend.position = "none") +
  scale_x_log10() +
  geom_ribbon(data = pred_r_cdom,
              aes(x = nn_mean.mean_na, ymin = lwr, ymax = upr, fill = position),
              alpha = 0.2, color = NA, inherit.aes = FALSE) +
  geom_line(data = pred_r_cdom,
            aes(x = nn_mean.mean_na, y = fit, color = position),
            linewidth = 1, inherit.aes = FALSE) +
  geom_point(size = 3) +
  scale_color_manual(values = clr) + scale_fill_manual(values = clr) +
  xlab(expression("cDOM Absorbance (m"^-1*")")) +
  ylab(expression("R (g-m"^-2~d^-1*")"))

pred_r_sm <- make_pred(r_sm_lmer, "sm_nn_mean.mean_na",
                       range(cmc_month$sm_nn_mean.mean_na, na.rm = TRUE),
                       "position", c("Nearshore", "Offshore"))
sm_r_plot <- ggplot(cmc_month,
                     aes(x = sm_nn_mean.mean_na, y = R, color = position)) +
  theme_classic() + theme(legend.position = "none") +
  scale_x_log10() +
  geom_ribbon(data = pred_r_sm,
              aes(x = sm_nn_mean.mean_na, ymin = lwr, ymax = upr, fill = position),
              alpha = 0.2, color = NA, inherit.aes = FALSE) +
  geom_line(data = pred_r_sm,
            aes(x = sm_nn_mean.mean_na, y = fit, color = position),
            linewidth = 1, inherit.aes = FALSE) +
  geom_point(size = 3) +
  scale_color_manual(values = clr) + scale_fill_manual(values = clr) +
  xlab(expression("Suspended Matter (mg L"^-1*")")) +
  ylab(expression("R (g-m"^-2~d^-1*")"))

pred_r_kd <- make_pred(r_kd_lmer, "kd_nn_mean.mean_na",
                       range(cmc_month$kd_nn_mean.mean_na, na.rm = TRUE),
                       "position", c("Nearshore", "Offshore"))
kd_r_plot <- ggplot(cmc_month,
                     aes(x = kd_nn_mean.mean_na, y = R, color = position)) +
  theme_classic() + theme(legend.position = "none") +
  scale_x_log10() +
  geom_ribbon(data = pred_r_kd,
              aes(x = kd_nn_mean.mean_na, ymin = lwr, ymax = upr, fill = position),
              alpha = 0.2, color = NA, inherit.aes = FALSE) +
  geom_line(data = pred_r_kd,
            aes(x = kd_nn_mean.mean_na, y = fit, color = position),
            linewidth = 1, inherit.aes = FALSE) +
  geom_point(size = 3) +
  scale_color_manual(values = clr) + scale_fill_manual(values = clr) +
  xlab(expression("K"[d]~"(m"^-1*")")) +
  ylab(expression("R (g-m"^-2~d^-1*")"))

pred_r_chla <- make_pred(r_chla_lmer, "chlor_a_nn.mean_na",
                          range(cmc_month$chlor_a_nn.mean_na, na.rm = TRUE))
chla_r_plot <- ggplot(cmc_month,
                       aes(x = chlor_a_nn.mean_na, y = R, color = position)) +
  theme_classic() + theme(legend.position = "none") +
  scale_x_log10() +
  geom_ribbon(data = pred_r_chla,
              aes(x = chlor_a_nn.mean_na, ymin = lwr, ymax = upr),
              alpha = 0.2, fill = "gray50", color = NA, inherit.aes = FALSE) +
  geom_line(data = pred_r_chla,
            aes(x = chlor_a_nn.mean_na, y = fit),
            linewidth = 1, color = "black", inherit.aes = FALSE) +
  geom_point(size = 3) +
  scale_color_manual(values = clr) +
  xlab(expression("Chlorophyll "*italic(a)~"(\u03bcg L"^-1*")")) +
  ylab(expression("R (g-m"^-2~d^-1*")"))

pred_r_temp <- make_pred(r_temp_lmer, "mean_ts",
                          range(cmc_month$mean_ts, na.rm = TRUE),
                          "position", c("Nearshore", "Offshore"))
temp_r_plot <- ggplot(cmc_month,
                       aes(x = mean_ts, y = R, color = position)) +
  theme_classic() +
  scale_x_log10() +
  geom_ribbon(data = pred_r_temp,
              aes(x = mean_ts, ymin = lwr, ymax = upr, fill = position),
              alpha = 0.2, color = NA, inherit.aes = FALSE) +
  geom_line(data = pred_r_temp,
            aes(x = mean_ts, y = fit, color = position),
            linewidth = 1, inherit.aes = FALSE) +
  geom_point(size = 3) +
  scale_color_manual(values = clr) + scale_fill_manual(values = clr) +
  xlab("Temperature (\u00b0C)") +
  ylab(expression("R (g-m"^-2~d^-1*")"))

# Figure S2 combined
combined_r_plot <- (cdom_r_plot + sm_r_plot) /
                   (kd_r_plot   + chla_r_plot) /
                   (temp_r_plot + plot_spacer()) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 14))

ggsave("~/respiration_plots_combined.png", combined_r_plot,
       width = 9, height = 12, dpi = 300)

# ============================================================================
# PART 6: Carbon flux calculations and bootstrap uncertainty propagation
# (Results: "Net carbon flux"; Fig 6)
# ============================================================================

# Calculate carbon fixation and release across all combinations of
# photosynthetic quotient (PQ: 0.8-1.6) and respiratory quotient
# (RQ: 0.5-2.0) from the literature (Methods; n = 9 PQ x 4 RQ = 36 combos).
pq_values <- seq(0.8, 1.6, by = 0.1)
rq_values <- seq(0.5, 2.0, by = 0.5)

carbon_flux <- Metab %>%
  select(date, GPP_sp, R_sp, site, position, year) %>%
  crossing(pq = pq_values, rq = rq_values) %>%
  mutate(
    c_fixed    = (GPP_sp    * (1 / pq)) * (12 / 32),
    c_released = (abs(R_sp) * rq)       * (12 / 32),
    Julian     = yday(date)
  )

# Daily mean and SD across PQ/RQ combinations (propagated uncertainty)
carbon_sum <- carbon_flux %>%
  group_by(site, year, Julian, position) %>%
  summarise(
    mean_C_fix = mean(c_fixed,    na.rm = TRUE),
    sd_C_fix   = sd(c_fixed,     na.rm = TRUE),
    mean_C_rel = mean(c_released, na.rm = TRUE),
    sd_C_rel   = sd(c_released,  na.rm = TRUE),
    .groups    = "drop"
  )

# Bootstrap: n = 1,000 samples per site-day from normal distributions
# parameterised by the mean and SD from the PQ/RQ combinations above.
generate_bootstrap_samples <- function(data, n_samples = 1000) {
  results_list <- vector("list", nrow(data))
  for (i in seq_len(nrow(data))) {
    results_list[[i]] <- data.frame(
      site     = data$site[i],
      year     = data$year[i],
      Julian   = data$Julian[i],
      sample_id = seq_len(n_samples),
      c_fix    = rnorm(n_samples, data$mean_C_fix[i], data$sd_C_fix[i]),
      c_rel    = rnorm(n_samples, data$mean_C_rel[i], data$sd_C_rel[i])
    )
  }
  bind_rows(results_list)
}

bootstrap_samples <- generate_bootstrap_samples(carbon_sum)

# Bootstrap summary: mean and SD per site-day
bootstrap_summary <- bootstrap_samples %>%
  group_by(site, year, Julian) %>%
  summarise(
    mean_c_fix_bootstrap = mean(c_fix),
    sd_c_fix_bootstrap   = sd(c_fix),
    mean_c_rel_bootstrap = mean(c_rel),
    sd_c_rel_bootstrap   = sd(c_rel),
    .groups = "drop"
  ) %>%
  mutate(date = as.Date(Julian - 1, origin = paste0(year, "-01-01")))

# Add coordinates and position
coords <- Metab %>% group_by(site) %>%
  distinct(site, Latitude, Longitude)

Carbon_flux <- bootstrap_summary %>%
  left_join(coords, by = "site") %>%
  mutate(position = case_when(
    site %in% nearshore_sites ~ "Nearshore",
    site %in% offshore_sites  ~ "Offshore",
    TRUE ~ NA_character_
  )) %>%
  filter(mean_c_fix_bootstrap > 0, mean_c_rel_bootstrap > 0)

write.csv(Carbon_flux, "Carbon_flux.csv", row.names = FALSE)

# --- Correlation between fixation and release --------------------------------
carbon_cor <- cor.test(Carbon_flux$mean_c_fix_bootstrap,
                        Carbon_flux$mean_c_rel_bootstrap,
                        method = "pearson")
carbon_cor

# Stratified Pearson correlations by year and position (cited in Results)
stratified_cors <- Carbon_flux %>%
  group_by(year, position) %>%
  summarise(
    cor = cor(mean_c_fix_bootstrap, mean_c_rel_bootstrap),
    .groups = "drop"
  )

# --- Fixation:release ratios (Table 2) ---------------------------------------
# Geometric means with 95% CI computed on log-transformed ratios
fix_rel_ratios <- Carbon_flux %>%
  group_by(site, position) %>%
  summarise(
    n          = n(),
    mean_log_pr = mean(log(mean_c_fix_bootstrap / mean_c_rel_bootstrap),
                       na.rm = TRUE),
    se_log_pr  = sd(log(mean_c_fix_bootstrap / mean_c_rel_bootstrap),
                    na.rm = TRUE) /
      sqrt(sum(!is.na(mean_c_fix_bootstrap / mean_c_rel_bootstrap))),
    mean_pr    = exp(mean_log_pr),
    CI_pr_lower = exp(mean_log_pr - 1.96 * se_log_pr),
    CI_pr_upper = exp(mean_log_pr + 1.96 * se_log_pr),
    .groups    = "drop"
  )

# Annual cumulative ratios (Table 3)
year_end_values <- Carbon_flux %>%
  group_by(site, year, position) %>%
  summarise(
    final_cum_cfix = sum(mean_c_fix_bootstrap),
    final_cum_crel = sum(mean_c_rel_bootstrap),
    cfix_error     = sqrt(sum(sd_c_fix_bootstrap^2)),
    crel_error     = sqrt(sum(sd_c_rel_bootstrap^2)),
    .groups = "drop"
  )

year_end_annual_ratios <- year_end_values %>%
  mutate(log_ratio = log(final_cum_cfix / final_cum_crel)) %>%
  group_by(position, year) %>%
  summarise(
    n           = n(),
    mean_log_pr  = mean(log_ratio, na.rm = TRUE),
    se_log_pr    = sd(log_ratio, na.rm = TRUE) / sqrt(n()),
    mean_pr      = exp(mean_log_pr),
    CI_pr_lower  = exp(mean_log_pr - 1.96 * se_log_pr),
    CI_pr_upper  = exp(mean_log_pr + 1.96 * se_log_pr),
    .groups      = "drop"
  )

# ============================================================================
# PART 7: Figure 6 - Cumulative C fixation vs C release by site and year
# ============================================================================

ggplot(year_end_values,
       aes(x = final_cum_crel, y = final_cum_cfix,
           color = site, shape = position)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "gray50") +
  geom_errorbar(aes(ymin = final_cum_cfix - cfix_error,
                    ymax = final_cum_cfix + cfix_error),
                width = 1.5, alpha = 0.8, linewidth = 0.8) +
  geom_errorbarh(aes(xmin = final_cum_crel - crel_error,
                     xmax = final_cum_crel + crel_error),
                 height = 1.5, alpha = 0.8, linewidth = 0.8) +
  geom_point(size = 3) +
  facet_wrap(~ year, ncol = 1) +
  scale_color_viridis_d() +
  scale_shape_manual(values = c(Nearshore = 15, Offshore = 8)) +
  theme_bw() +
  theme(panel.grid       = element_blank(),
        aspect.ratio     = 0.7,
        legend.position  = "right") +
  labs(x     = expression("Cumulative C release (g-C m"^-2*")"),
       y     = expression("Cumulative C fixation (g-C m"^-2*")"),
       color = "Site",
       shape = "Position")
# Export: width = 450, height = 1350

# ============================================================================
# PART 8: Figures 7 & 8 - Daily carbon flux time series with CyAN HABs maps
# ============================================================================
# CyAN imagery (NASA Ocean Biology Distributed Active Archive Center)
# is read as GeoTIFF files saved locally (one file per target date).
# The five days centred on each target date are averaged to smooth
# short-term variability.

# Target dates are the peak fixation and release events visible in the
# daily time series for 2017, 2018, and 2019.
original_dates <- c(
  "2017-06-09", "2017-07-19", "2017-08-19", "2017-09-16",
  "2018-08-04", "2018-08-24", "2018-09-19",
  "2019-06-25", "2019-07-10", "2019-09-15", "2019-09-24"
)

expand_dates <- function(date) {
  date <- as.Date(date)
  data.frame(
    expanded_date = seq(date - 2, date + 2, by = "day"),
    central_date  = date,
    year          = year(date),
    Julian        = yday(date)
  )
}

date_mapping <- do.call(rbind, lapply(original_dates, expand_dates))

# Compute mean C flux (fixation - release) and join to date mapping
Carbon_flux <- Carbon_flux %>%
  mutate(mean_c_flux_bootstrap =
           mean_c_fix_bootstrap - mean_c_rel_bootstrap)

NOAA_maps <- Carbon_flux %>%
  inner_join(date_mapping %>% select(-year, -Julian),
             by = c("date" = "expanded_date"))

NOAA_maps_avg <- NOAA_maps %>%
  group_by(central_date, site) %>%
  summarise(
    across(where(is.numeric) & !c(Latitude, Longitude),
           ~ mean(.x, na.rm = TRUE)),
    Latitude  = first(Latitude),
    Longitude = first(Longitude),
    .groups   = "drop"
  )

# Map TIFF filenames to dates
get_tiff_path <- function(date) {
  paste0(format(date, "%Y-%m-%d"), "_NASA_HABs.tif")
}

# Helper for safe PNG creation
safe_create_png <- function(filename, width, height, res, plot_fn) {
  tryCatch({
    png(filename, width = width, height = height, res = res)
    plot_fn()
    dev.off()
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    stop(paste("Error creating", filename, ":", e$message))
  })
}

# --- Create one panel per target date ----------------------------------------
# The value_var argument controls which column is plotted as point colour:
#   "mean_c_fix_bootstrap"  → Figure 7 (carbon fixation)
#   "mean_c_rel_bootstrap"  → Figure 8 (carbon release)
create_cflux_maps <- function(value_var,
                               value_range = c(-4, 1),
                               out_file    = "cflux_main_plots.png") {
  unique_dates <- sort(unique(NOAA_maps_avg$central_date))

  # Compute consistent map extent across all dates
  first_tiff   <- terra::rast(get_tiff_path(unique_dates[1]))
  all_sf       <- terra::vect(NOAA_maps_avg,
                               geom = c("Longitude", "Latitude"),
                               crs  = "EPSG:4326")
  all_trans    <- terra::project(all_sf, terra::crs(first_tiff))
  all_coords   <- terra::crds(all_trans)
  buffer       <- 50000   # 50 km
  map_ext      <- terra::ext(
    range(all_coords[, 1])[1] - buffer,
    range(all_coords[, 1])[2] + buffer,
    range(all_coords[, 2])[1] - buffer,
    range(all_coords[, 2])[2] + buffer
  )

  plot_fn <- function() {
    par(mfrow = c(3, 4), mar = c(2, 2, 3, 2))

    for (current_date in unique_dates) {
      current_data <- NOAA_maps_avg[
        NOAA_maps_avg$central_date == current_date &
          !is.na(NOAA_maps_avg[[value_var]]), ]
      tiff_path <- get_tiff_path(current_date)
      if (!file.exists(tiff_path)) { warning(paste("Not found:", tiff_path)); next }

      tryCatch({
        base_map   <- terra::rast(tiff_path)
        pts_sf     <- terra::vect(current_data, geom = c("Longitude", "Latitude"),
                                   crs = "EPSG:4326")
        pts_trans  <- terra::project(pts_sf, terra::crs(base_map))
        current_data$x <- terra::crds(pts_trans)[, 1]
        current_data$y <- terra::crds(pts_trans)[, 2]
        cropped    <- terra::crop(base_map, map_ext)

        terra::plot(cropped, axes = FALSE, legend = FALSE)

        plot_vals <- pmin(pmax(current_data[[value_var]],
                               value_range[1]), value_range[2])
        col_idx   <- findInterval(plot_vals,
                                   seq(value_range[1], value_range[2],
                                       length.out = 100),
                                   all.inside = TRUE)
        point_cols <- viridis::viridis(100)[col_idx]
        points(current_data$x, current_data$y,
               pch = 21, bg = point_cols, col = "black", cex = 8)
        text(current_data$x, current_data$y,
             labels = current_data$site, col = "white", font = 2, cex = 3.5)
        title(paste0(year(current_date), " JD ", yday(current_date)),
              cex.main = 3.5, line = 0.8)
      }, error = function(e) {
        cat("Skipping", as.character(current_date), ":", e$message, "\n")
      })
    }
    # Fill remaining grid cells
    remaining <- 12 - length(unique_dates)
    if (remaining > 0) for (i in seq_len(remaining)) plot.new()
  }

  safe_create_png(out_file, 4400, 1800, 120, plot_fn)

  # Standalone legend
  legend_fn <- function() {
    par(mar = c(5, 0, 8, 2))
    plot.new()
    fields::image.plot(
      legend.only  = TRUE, zlim = value_range,
      col          = viridis::viridis(99),
      legend.shrink = 0.9, legend.width = 1.5, legend.mar = 6,
      legend.args  = list(text = "Mean C Flux\n(g-O2/m2)",
                          side = 3, line = 1, cex = 1.8),
      axis.args    = list(cex.axis = 1.5, las = 1)
    )
  }
  safe_create_png(gsub("main_plots", "legend", out_file),
                  600, 1800, 120, legend_fn)
}

# Figure 7: C fixation maps
create_cflux_maps("mean_c_fix_bootstrap",
                   value_range = c(-4, 1),
                   out_file    = "cfixation_main_plots.png")

# Figure 8: C release maps
create_cflux_maps("mean_c_rel_bootstrap",
                   value_range = c(-4, 1),
                   out_file    = "crelease_main_plots.png")

# --- Time-series ribbon plots used in Figs 7 & 8 ----------------------------
viridis_site_cols <- c(
  CHRP1 = "#440154ff", CHRP2 = "#453781FF", CHRP3 = "#33638DFF",
  CHRP4 = "#1F968BFF", CHRP5 = "#3CBB75FF", CHRP7 = "#73D055FF",
  CHRP8 = "#95D840FF", CHRP9 = "#FDE725FF"
)

# Carbon fixation time series (Fig 7)
ggplot(Carbon_flux, aes(x = Julian, color = site, fill = site)) +
  geom_ribbon(aes(ymin = pmax(0, mean_c_fix_bootstrap - sd_c_fix_bootstrap),
                  ymax = mean_c_fix_bootstrap + sd_c_fix_bootstrap),
              alpha = 0.2) +
  geom_line(aes(y = mean_c_fix_bootstrap), linewidth = 0.8) +
  facet_wrap(~ year, scales = "free", ncol = 2) +
  scale_color_manual(values = viridis_site_cols) +
  scale_fill_manual(values  = viridis_site_cols) +
  theme_bw(base_size = 14) +
  labs(x = "Julian Day",
       y = expression(paste("Carbon Fixation (g C ", m^-2, " ", d^-1, ")"))) +
  theme(legend.position  = "right",
        panel.grid       = element_blank(),
        strip.background = element_rect(fill = "white"))

# Carbon release time series (Fig 8)
ggplot(Carbon_flux, aes(x = Julian, color = site, fill = site)) +
  geom_ribbon(aes(ymin = pmax(0, mean_c_rel_bootstrap - sd_c_rel_bootstrap),
                  ymax = mean_c_rel_bootstrap + sd_c_rel_bootstrap),
              alpha = 0.2) +
  geom_line(aes(y = mean_c_rel_bootstrap), linewidth = 0.8) +
  facet_wrap(~ year, scales = "free", ncol = 2) +
  scale_color_manual(values = viridis_site_cols) +
  scale_fill_manual(values  = viridis_site_cols) +
  theme_bw(base_size = 14) +
  labs(x = "Julian Day",
       y = expression(paste("Carbon Release (g C ", m^-2, " ", d^-1, ")"))) +
  theme(legend.position  = "right",
        panel.grid       = element_blank(),
        strip.background = element_rect(fill = "white"))
# Export: width = 1600, height = 783

# NOTE ON FIGURE S1 (Full vs Submerged dataset comparison; Table S2):
# The linear mixed-effects models comparing GPP, R, and NEP between the Full
# (shallow/surface) and Submerged (deep) datasets are summarised in Table S2
# of the supplementary information. The code that generated these models and
# Fig S1 was not recovered from the provided scripts and is therefore not
# included here. The models used site and date as random effects, with the
# Full/Submerged grouping as a fixed effect; model statistics are given in
# Table S2.
