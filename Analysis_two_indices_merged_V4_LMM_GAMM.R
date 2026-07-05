# ============================================================
# Analysis_two_indices_merged_V3.R
# Original LOESS visualisation  +
# Section B: Linear Mixed-Effects Models (LMM)
# Section C: Generalized Additive Mixed Models (GAMM)
# Section D: Summary tables + confirmatory charts
# ============================================================

library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)      # plot_grid + get_legend
library(lme4)         # lmer()
library(lmerTest)     # p-values for LMM via Satterthwaite
library(mgcv)         # gam()
library(ggeffects)    # ggpredict() — marginal predictions for LMM
library(knitr)        # kable() for console tables

# setwd("C:/Users/santamaria/Documents/Spatial_analysis/Rstudio")
setwd("E:/SDM_project")
# ── colour palette (shared throughout) ───────────────────────
cat_cols <- c("Reserve"    = "#1B9E77",
              "No Reserve" = "#D95F02")

# ============================================================
# 1.  DATA LOADING & HARMONISATION   (unchanged from V2)
# ============================================================

load_country <- function(ndvi_file, evi_file,
                         out_label = "outReserve",
                         drop_cat  = NULL,
                         country   = "CO") {

  read_index <- function(file, index_name) {
    d  <- read.csv(file)
    d2 <- pivot_longer(d, cols = colnames(d)[4:29], names_to = "Year")
    d2$category <- ifelse(d2$nombre != out_label, "Reserve", "No Reserve")
    d2$index    <- index_name
    d2
  }

  dat <- rbind(read_index(ndvi_file, "NDVI"),
               read_index(evi_file,  "EVI"))

  dat$Year     <- as.numeric(gsub("X", "", dat$Year))
  dat$category <- as.factor(dat$category)
  dat$country  <- country

  if (!is.null(drop_cat) && "categoria" %in% names(dat))
    dat <- dat[dat$categoria != drop_cat, ]

  dat
}

dat_CO <- load_country("data/points_CO_NDVI.csv",    "data/points_CO_EVI.csv",
                       out_label = "outReserve",  country = "Colombia (CO)")
dat_AR <- load_country("data/points_AR_PD_NDVI.csv", "data/points_AR_PD_EVI.csv",
                       out_label = "outReserves",
                       drop_cat  = "areas-importantes-para-la-consevacion-de-las-aves",
                       country   = "Argentina (AR)")
dat_PY <- load_country("data/points_PY_NDVI.csv",    "data/points_PY_EVI.csv",
                       out_label = "outReserve",
                       drop_cat  = "areas-importantes-para-la-consevacion-de-las-aves",
                       country   = "Paraguay (PY)")

common_cols <- Reduce(intersect, list(names(dat_CO), names(dat_AR), names(dat_PY)))
dat_all     <- rbind(dat_CO[, common_cols],
                     dat_AR[, common_cols],
                     dat_PY[, common_cols])

dat_all$country  <- factor(dat_all$country,
                           levels = c("Colombia (CO)", "Argentina (AR)", "Paraguay (PY)"))
dat_all$category <- factor(dat_all$category, levels = c("Reserve", "No Reserve"))

# Centre Year to help model convergence (Year_c = 0 at 2000)
dat_all$Year_c <- dat_all$Year - min(dat_all$Year)

# ============================================================
# 2.  ORIGINAL LOESS VISUALISATION   (unchanged from V2)
# ============================================================

base_theme <- theme_bw(base_size = 10) +
  theme(
    axis.text.x      = element_text(angle = 45, vjust = 0.5),
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold"),
    legend.position  = "none"
  )

make_plot <- function(dat_country, country_label, show_x_axis = FALSE) {

  hline_overall <- dat_country %>%
    group_by(index) %>%
    summarise(mean_val = mean(value, na.rm = TRUE), .groups = "drop")

  hline_cat <- dat_country %>%
    group_by(index, category) %>%
    summarise(mean_val = mean(value, na.rm = TRUE), .groups = "drop")

  p <- ggplot(dat_country, aes(x = Year, y = value, color = category)) +
    geom_hline(data = hline_overall, aes(yintercept = mean_val),
               color = "grey30", linetype = "dashed",
               linewidth = 0.65, alpha = 0.85, inherit.aes = FALSE) +
    geom_hline(data = hline_cat, aes(yintercept = mean_val, color = category),
               linetype = "dashed", linewidth = 0.55, alpha = 0.80,
               inherit.aes = FALSE) +
    stat_summary(fun = mean, geom = "line",
                 aes(group = category), linewidth = 0.4) +
    geom_smooth(method = "loess", level = 0.90, se = TRUE, linewidth = 1.1) +
    facet_wrap(~ index, scales = "free_y", nrow = 1) +
    scale_color_manual(values = cat_cols, name = "Category") +
    labs(y = "Index value", title = country_label,
         x = if (show_x_axis) "Year" else NULL) +
    base_theme +
    theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5))

  if (!show_x_axis)
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  p
}

p_CO <- make_plot(dat_all[dat_all$country == "Colombia (CO)", ],  "Colombia (CO)",  FALSE)
p_AR <- make_plot(dat_all[dat_all$country == "Argentina (AR)", ], "Argentina (AR)", FALSE)
p_PY <- make_plot(dat_all[dat_all$country == "Paraguay (PY)", ],  "Paraguay (PY)",  TRUE)

legend_plot    <- p_CO + theme(legend.position = "right",
                               legend.title = element_text(face = "bold"))
shared_legend  <- get_legend(legend_plot)
panel_col      <- plot_grid(p_CO, p_AR, p_PY, ncol = 1,
                            align = "v", axis = "lr",
                            rel_heights = c(1, 1, 1.15))
combined       <- plot_grid(panel_col, shared_legend, ncol = 2, rel_widths = c(1, 0.15))

png("index_reserves_combined.png", width = 8, height = 10, units = "in", res = 600)
print(combined); dev.off()
message("LOESS combined figure saved.")

# ============================================================
# SECTION B — LINEAR MIXED-EFFECTS MODELS (LMM)
# ============================================================
# Model per country × index:
#   value ~ category * Year_c + (1 | nombre)
#
#   Fixed effects
#     category          : baseline difference Reserve vs No Reserve
#     Year_c            : average trend over time (across both categories)
#     category:Year_c   : KEY TERM — do the two categories diverge/converge?
#
#   Random effect
#     (1 | nombre)      : each sampling site has its own intercept,
#                         accounting for repeated measurements over 26 years
# ============================================================

countries <- levels(dat_all$country)
indices   <- c("NDVI", "EVI")

# ── B1. Fit all 6 LMMs and store them ───────────────────────
lmm_list <- list()

for (ctry in countries) {
  for (idx in indices) {
    key <- paste(ctry, idx, sep = " | ")
    sub <- dat_all %>% filter(country == ctry, index == idx)

    lmm_list[[key]] <- lmer(
      value ~ category * Year_c + (1 | nombre),
      data    = sub,
      REML    = TRUE,
      control = lmerControl(optimizer = "bobyqa")
    )
  }
}

# ── B2. Extract fixed-effect table for all models ───────────
lmm_summary_df <- bind_rows(lapply(names(lmm_list), function(key) {
  m   <- lmm_list[[key]]
  cf  <- as.data.frame(coef(summary(m)))
  cf$term    <- rownames(cf)
  cf$model   <- key
  # Parse country / index back out
  parts <- strsplit(key, " \\| ")[[1]]
  cf$country <- parts[1]
  cf$index   <- parts[2]
  cf
})) %>%
  rename(estimate = Estimate, se = `Std. Error`,
         df = df, t_val = `t value`, p_val = `Pr(>|t|)`) %>%
  mutate(
    sig = case_when(
      p_val < 0.001 ~ "***",
      p_val < 0.01  ~ "**",
      p_val < 0.05  ~ "*",
      p_val < 0.10  ~ ".",
      TRUE          ~ "n.s."
    )
  ) %>%
  select(country, index, term, estimate, se, df, t_val, p_val, sig)

# Print to console
cat("\n========================================================\n")
cat("  LMM FIXED-EFFECTS SUMMARY\n")
cat("========================================================\n")
print(kable(lmm_summary_df, digits = 4, format = "simple"))

# AFTER
gamm_summary_df_display <- gamm_summary_df %>%
  mutate(p_val = ifelse(p_val < 0.0001,
                        formatC(p_val, format = "e", digits = 2),
                        formatC(p_val, format = "f", digits = 4)))
print(kable(gamm_summary_df_display, format = "simple"))

# Save to CSV
write.csv(lmm_summary_df, "LMM_fixed_effects_summary.csv", row.names = FALSE)
message("LMM summary CSV saved.")

# ── B3. Marginal predictions from each LMM (for plotting) ───
lmm_preds <- bind_rows(lapply(names(lmm_list), function(key) {
  m     <- lmm_list[[key]]
  parts <- strsplit(key, " \\| ")[[1]]

  # ggpredict averages over random effects — population-level prediction
  preds <- ggpredict(m, terms = c("Year_c [all]", "category"))
  preds <- as.data.frame(preds)
  names(preds)[names(preds) == "x"]     <- "Year_c"
  names(preds)[names(preds) == "group"] <- "category"

  # Back-transform Year_c to real Year
  yr_min <- min(dat_all$Year)
  preds$Year    <- preds$Year_c + yr_min
  preds$country <- parts[1]
  preds$index   <- parts[2]
  preds
}))

lmm_preds$country  <- factor(lmm_preds$country,  levels = levels(dat_all$country))
lmm_preds$category <- factor(lmm_preds$category, levels = c("Reserve", "No Reserve"))

# ── B4. Raw annual means (for background points) ─────────────
raw_means <- dat_all %>%
  group_by(country, index, category, Year) %>%
  summarise(mean_val = mean(value, na.rm = TRUE), .groups = "drop")

# ── B5. Interaction significance labels for the plot ─────────
interaction_labels <- lmm_summary_df %>%
  filter(grepl("category.*Year|Year.*category", term, ignore.case = TRUE)) %>%
  mutate(label = paste0("Interaction p ", sig,
                        "\n(p = ", formatC(p_val, digits = 3, format = "f"), ")"))

# ── B6. LMM CHART ────────────────────────────────────────────
# One facet per country × index showing:
#   • thin lines = annual raw means per category
#   • thick ribbons = LMM marginal fit ± 95 % CI
#   • annotation = interaction p-value

p_lmm <- ggplot() +

  # Raw annual means as thin background lines
  geom_line(data = raw_means,
            aes(x = Year, y = mean_val, color = category),
            linewidth = 0.35, alpha = 0.45, linetype = "solid") +

  # LMM marginal fit ribbon
  geom_ribbon(data = lmm_preds,
              aes(x = Year, ymin = conf.low, ymax = conf.high, fill = category),
              alpha = 0.18) +

  # LMM marginal fit line
  geom_line(data = lmm_preds,
            aes(x = Year, y = predicted, color = category),
            linewidth = 1.2) +

  # Interaction p-value annotation
  geom_text(data = interaction_labels,
            aes(label = label),
            x = -Inf, y = Inf,
            hjust = -0.08, vjust = 1.4,
            size = 2.8, fontface = "italic", color = "grey20",
            inherit.aes = FALSE) +

  facet_grid(country ~ index, scales = "free_y") +
  scale_color_manual(values = cat_cols, name = "Category") +
  scale_fill_manual(values  = cat_cols, name = "Category") +
  scale_x_continuous(breaks = seq(1998, 2024, by = 4)) +
  labs(
    title    = "LMM marginal fits — Reserve vs No Reserve",
    subtitle = "Thick line = population-level prediction ± 95 % CI  |  Thin line = annual raw means",
    x        = "Year",
    y        = "Index value"
  ) +
  theme_bw(base_size = 10) +
  theme(
    strip.background  = element_rect(fill = "grey90"),
    strip.text        = element_text(face = "bold", size = 9),
    legend.position   = "bottom",
    legend.title      = element_text(face = "bold"),
    axis.text.x       = element_text(angle = 45, vjust = 0.5),
    plot.title        = element_text(face = "bold", hjust = 0.5),
    plot.subtitle     = element_text(size = 8, hjust = 0.5, color = "grey40"),
    panel.spacing     = unit(0.7, "lines")
  )

png("LMM_marginal_fits.png", width = 9, height = 10, units = "in", res = 600)
print(p_lmm); dev.off()
message("LMM chart saved.")

# ── B7. Coefficient plot — interaction term only ─────────────
# Visualises effect size and uncertainty of the key term across all 6 models

interact_coef <- lmm_summary_df %>%
  filter(grepl("category.*Year|Year.*category", term, ignore.case = TRUE)) %>%
  mutate(
    ci_lo = estimate - 1.96 * se,
    ci_hi = estimate + 1.96 * se,
    label = paste(country, index, sep = "\n")
  )

p_coef <- ggplot(interact_coef,
                 aes(x = estimate, y = reorder(label, estimate),
                     color = sig == "n.s.")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.25, linewidth = 0.8) +
  geom_point(size = 3.5) +
  scale_color_manual(values = c("FALSE" = "#1B9E77", "TRUE" = "grey60"),
                     labels = c("FALSE" = "Significant (p < 0.05)",
                                "TRUE"  = "Not significant"),
                     name = NULL) +
  labs(
    title    = "LMM — category × Year interaction coefficients",
    subtitle = "Positive = Reserve trend steeper than No Reserve over time",
    x        = "Interaction coefficient ± 95 % CI",
    y        = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title     = element_text(face = "bold", hjust = 0.5),
    plot.subtitle  = element_text(size = 8, hjust = 0.5, color = "grey40"),
    legend.position = "bottom",
    panel.grid.major.y = element_line(color = "grey92")
  )

png("LMM_interaction_coefplot.png", width = 7, height = 5, units = "in", res = 600)
print(p_coef); dev.off()
message("LMM coefficient plot saved.")

# ── B8. LMM residual diagnostics ─────────────────────────────

png("LMM_residual_diagnostics.png", width = 10, height = 8, units = "in", res = 600)
par(mfrow = c(3, 4), mar = c(4, 4, 3, 1))
for (key in names(lmm_list)) {
  m   <- lmm_list[[key]]
  res <- resid(m)
  fit <- fitted(m)
  plot(fit, res,
       main  = key,
       xlab  = "Fitted", ylab  = "Residuals",
       pch   = 16, cex = 0.5, col = alpha("steelblue", 0.4))
  abline(h = 0, col = "red", lty = 2)
}
dev.off()

png("LMM_qqplots.png", width = 10, height = 8, units = "in", res = 600)
par(mfrow = c(3, 4), mar = c(4, 4, 3, 1))
for (key in names(lmm_list)) {
  m   <- lmm_list[[key]]
  qqnorm(resid(m), main = key, pch = 16, cex = 0.5, col = alpha("steelblue", 0.5))
  qqline(resid(m), col = "red", lty = 2)
}
dev.off()
message("LMM diagnostic plots saved.")

# ============================================================
# SECTION C — GENERALIZED ADDITIVE MIXED MODELS (GAMM)
# ============================================================
# Model per country × index:
#   value ~ category + s(Year_c, by = category, k = 10) + s(nombre, bs = "re")
#
#   s(Year_c, by = category)  : separate smooth per category — this is the
#                               formal confirmatory equivalent of the LOESS
#   s(nombre, bs = "re")      : random intercept per site (re = random effect)
#   k = 10                    : max 10 basis functions, enough for 26 years
# ============================================================

gamm_list <- list()

for (ctry in countries) {
  for (idx in indices) {
    key <- paste(ctry, idx, sep = " | ")
    sub <- dat_all %>% filter(country == ctry, index == idx)

    # nombre must be a factor for the random-effect smooth
    sub$nombre_f <- factor(sub$nombre)

    gamm_list[[key]] <- gam(
      value ~ category +
        s(Year_c, by = category, k = 10) +
        s(nombre_f, bs = "re"),
      data   = sub,
      method = "REML",
      family = gaussian()
    )
  }
}

# ── C1. Smooth significance table ────────────────────────────
gamm_summary_df <- bind_rows(lapply(names(gamm_list), function(key) {
  m     <- gamm_list[[key]]
  parts <- strsplit(key, " \\| ")[[1]]
  sm    <- summary(m)$s.table
  df_sm <- as.data.frame(sm)
  df_sm$smooth  <- rownames(sm)
  df_sm$country <- parts[1]
  df_sm$index   <- parts[2]
  df_sm
})) %>%
  rename(edf = edf, ref_df = Ref.df, F_val = F, p_val = `p-value`) %>%
  filter(!grepl("nombre", smooth)) %>%   # drop site random-effect rows
  mutate(
    sig = case_when(
      p_val < 0.001 ~ "***",
      p_val < 0.01  ~ "**",
      p_val < 0.05  ~ "*",
      p_val < 0.10  ~ ".",
      TRUE          ~ "n.s."
    )
  ) %>%
  select(country, index, smooth, edf, ref_df, F_val, p_val, sig)

cat("\n========================================================\n")
cat("  GAMM SMOOTH-TERM SIGNIFICANCE TABLE\n")
cat("========================================================\n")
print(kable(gamm_summary_df, digits = 4, format = "simple"))

gamm_summary_csv <- gamm_summary_df
write.csv(gamm_summary_csv, "GAMM_smooth_significance.csv", row.names = FALSE)
write.csv(gamm_summary_df, "GAMM_smooth_significance.csv", row.names = FALSE)
message("GAMM summary CSV saved.")

# ── C2. Predicted smooth curves per model ────────────────────
gamm_preds <- bind_rows(lapply(names(gamm_list), function(key) {
  m     <- gamm_list[[key]]
  parts <- strsplit(key, " \\| ")[[1]]
  ctry  <- parts[1]
  idx   <- parts[2]

  sub <- dat_all %>% filter(country == ctry, index == idx)
  sub$nombre_f <- factor(sub$nombre)

  yr_seq <- seq(min(sub$Year_c), max(sub$Year_c), length.out = 100)
  yr_min <- min(dat_all$Year)

  preds_list <- lapply(c("Reserve", "No Reserve"), function(cat_lbl) {
    new_d <- data.frame(
      Year_c   = yr_seq,
      category = factor(cat_lbl, levels = levels(sub$category)),
      nombre_f = factor(levels(sub$nombre_f)[1],   # reference site — excluded by exclude
                        levels = levels(sub$nombre_f))
    )

    # exclude = "s(nombre_f)" drops the random intercept → population-level prediction
    pred <- predict(m, newdata = new_d, se.fit = TRUE,
                    exclude = "s(nombre_f)")
    data.frame(
      Year     = yr_seq + yr_min,
      predicted = pred$fit,
      se        = pred$se.fit,
      conf.low  = pred$fit - 1.96 * pred$se.fit,
      conf.high = pred$fit + 1.96 * pred$se.fit,
      category  = cat_lbl,
      country   = ctry,
      index     = idx
    )
  })
  bind_rows(preds_list)
}))

gamm_preds$country  <- factor(gamm_preds$country,  levels = levels(dat_all$country))
gamm_preds$category <- factor(gamm_preds$category, levels = c("Reserve", "No Reserve"))

# ── C3. Smooth significance labels for GAMM chart ────────────
gamm_sig_labels <- gamm_summary_df %>%
  mutate(
    category = case_when(
      grepl("No Reserve", smooth) ~ "No Reserve",
      grepl("Reserve",    smooth) ~ "Reserve",
      TRUE                        ~ "Reserve"
    ),
    label = paste0(category, ": edf=", round(edf, 1), " ", sig)
  ) %>%
  group_by(country, index) %>%
  summarise(label = paste(label, collapse = "\n"), .groups = "drop")

# ── C4. GAMM CHART ───────────────────────────────────────────
p_gamm <- ggplot() +

  geom_line(data = raw_means,
            aes(x = Year, y = mean_val, color = category),
            linewidth = 0.35, alpha = 0.4, linetype = "solid") +

  geom_ribbon(data = gamm_preds,
              aes(x = Year, ymin = conf.low, ymax = conf.high, fill = category),
              alpha = 0.18) +

  geom_line(data = gamm_preds,
            aes(x = Year, y = predicted, color = category),
            linewidth = 1.2) +

  geom_text(data = gamm_sig_labels,
            aes(label = label),
            x = -Inf, y = Inf,
            hjust = -0.06, vjust = 1.35,
            size = 2.5, fontface = "italic", color = "grey20",
            inherit.aes = FALSE) +

  facet_grid(country ~ index, scales = "free_y") +
  scale_color_manual(values = cat_cols, name = "Category") +
  scale_fill_manual(values  = cat_cols, name = "Category") +
  scale_x_continuous(breaks = seq(1998, 2024, by = 4)) +
  labs(
    title    = "GAMM smooth fits — Reserve vs No Reserve",
    subtitle = "Thick line = population-level GAM smooth ± 95 % CI  |  Thin line = annual raw means",
    x        = "Year",
    y        = "Index value"
  ) +
  theme_bw(base_size = 10) +
  theme(
    strip.background   = element_rect(fill = "grey90"),
    strip.text         = element_text(face = "bold", size = 9),
    legend.position    = "bottom",
    legend.title       = element_text(face = "bold"),
    axis.text.x        = element_text(angle = 45, vjust = 0.5),
    plot.title         = element_text(face = "bold", hjust = 0.5),
    plot.subtitle      = element_text(size = 8, hjust = 0.5, color = "grey40"),
    panel.spacing      = unit(0.7, "lines")
  )

png("GAMM_smooth_fits.png", width = 9, height = 10, units = "in", res = 600)
print(p_gamm); dev.off()
message("GAMM chart saved.")
# ── C5. Difference smooths — Reserve minus No Reserve ────────
# A positive difference means Reserve > No Reserve at that year.
# Where the 90% CI does not cross zero, the difference is significant.
# No gratia dependency — uses predict() directly on the gam object.

diff_preds <- bind_rows(lapply(names(gamm_list), function(key) {
  m     <- gamm_list[[key]]
  parts <- strsplit(key, " \\| ")[[1]]
  ctry  <- parts[1]
  idx   <- parts[2]

  sub      <- dat_all %>% filter(country == ctry, index == idx)
  sub$nombre_f <- factor(sub$nombre)
  yr_min   <- min(dat_all$Year)
  yr_seq   <- seq(min(sub$Year_c), max(sub$Year_c), length.out = 200)

  make_nd <- function(cat_lbl) {
    data.frame(
      Year_c   = yr_seq,
      category = factor(cat_lbl, levels = levels(sub$category)),
      nombre_f = factor(levels(sub$nombre_f)[1], levels = levels(sub$nombre_f))
    )
  }

  p_res  <- predict(m, newdata = make_nd("Reserve"),
                    se.fit = TRUE, exclude = "s(nombre_f)",
                    newdata.guaranteed = TRUE)
  p_nres <- predict(m, newdata = make_nd("No Reserve"),
                    se.fit = TRUE, exclude = "s(nombre_f)",
                    newdata.guaranteed = TRUE)

  diff_fit <- p_res$fit - p_nres$fit
  diff_se  <- sqrt(p_res$se.fit^2 + p_nres$se.fit^2)

  # 90% CI (z = 1.645) — more sensitive than 95% for detecting sig. periods
  z <- 1.645

  data.frame(
    Year     = yr_seq + yr_min,
    diff     = diff_fit,
    diff_lo  = diff_fit - z * diff_se,
    diff_hi  = diff_fit + z * diff_se,
    sig_band = (diff_fit - z * diff_se > 0) |
               (diff_fit + z * diff_se < 0),
    country  = ctry,
    index    = idx
  )
}))

diff_preds$country <- factor(diff_preds$country, levels = levels(dat_all$country))

p_diff <- ggplot(diff_preds, aes(x = Year, y = diff)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_ribbon(aes(ymin = diff_lo, ymax = diff_hi,
                  fill = sig_band), alpha = 0.30) +
  geom_line(aes(color = sig_band), linewidth = 1.0) +
  facet_grid(country ~ index, scales = "free_y") +
  scale_fill_manual(values  = c("TRUE" = "#1B9E77", "FALSE" = "grey70"),
                    labels  = c("TRUE" = "Sig. difference", "FALSE" = "Not significant"),
                    name    = NULL) +
  scale_color_manual(values = c("TRUE" = "#1B9E77", "FALSE" = "grey50"),
                     labels = c("TRUE" = "Sig. difference", "FALSE" = "Not significant"),
                     name   = NULL) +
  scale_x_continuous(breaks = seq(1998, 2024, by = 4)) +
  labs(
    title    = "GAMM difference smooth — Reserve minus No Reserve",
    subtitle = "Green = periods where the 90% CI excludes zero (significant difference)",
    x        = "Year",
    y        = "Difference in index value"
  ) +
  theme_bw(base_size = 10) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold", size = 9),
    legend.position  = "bottom",
    axis.text.x      = element_text(angle = 45, vjust = 0.5),
    plot.title       = element_text(face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(size = 8, hjust = 0.5, color = "grey40"),
    panel.spacing    = unit(0.7, "lines")
  )

png("GAMM_difference_smooths.png", width = 9, height = 10, units = "in", res = 600)
print(p_diff); dev.off()
message("GAMM difference-smooth chart saved.")


# ============================================================
# SECTION D — COMBINED SIGNIFICANCE SUMMARY CHART
# ============================================================
# Heatmap of p-values: rows = country × index,
# columns = LMM interaction + GAMM Reserve smooth + GAMM No Reserve smooth

lmm_heat <- lmm_summary_df %>%
  filter(grepl("category.*Year|Year.*category", term, ignore.case = TRUE)) %>%
  mutate(test = "LMM\ncategory×Year") %>%
  select(country, index, test, p_val, sig)

gamm_heat <- gamm_summary_df %>%
  mutate(
    category = case_when(
      grepl("No Reserve", smooth) ~ "No Reserve",
      TRUE                        ~ "Reserve"
    ),
    test = paste0("GAMM smooth\n", category)
  ) %>%
  select(country, index, test, p_val, sig)

heat_df <- bind_rows(lmm_heat, gamm_heat) %>%
  mutate(
    country = factor(country, levels = levels(dat_all$country)),
    test    = factor(test, levels = c("LMM\ncategory×Year",
                                      "GAMM smooth\nReserve",
                                      "GAMM smooth\nNo Reserve")),
    neg_log_p = pmin(-log10(p_val), 4)   # cap at 4 for colour range
  )

p_heat <- ggplot(heat_df,
                 aes(x = test, y = interaction(index, country, sep = "\n"),
                     fill = neg_log_p)) +
  geom_tile(color = "white", linewidth = 0.6) +
  # AFTER
  geom_text(aes(label = paste0(sig, "\n",
                               ifelse(p_val < 0.0001,
                                      paste0("p<0.0001"),
                                      paste0("p=", formatC(p_val, digits = 3, format = "f"))))),
    size = 2.8, color = "white", fontface = "bold") +
  scale_fill_gradient(low = "grey85", high = "#1B6E5E",
                      name = expression(-log[10](p)),
                      limits = c(0, 4),
                      breaks = c(0, 1, 2, 3, 4),
                      labels = c("1", "0.1", "0.01", "0.001", "≤0.0001")) +
  labs(
    title    = "Significance summary — LMM & GAMM",
    subtitle = "Darker = more significant  |  *** p<0.001  ** p<0.01  * p<0.05",
    x        = NULL,
    y        = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title      = element_text(face = "bold", hjust = 0.5),
    plot.subtitle   = element_text(size = 8,  hjust = 0.5, color = "grey40"),
    axis.text.x     = element_text(size = 9),
    axis.text.y     = element_text(size = 9),
    legend.position = "right",
    panel.grid      = element_blank()
  )

print(p_heat)
png("Significance_heatmap.png", width = 8, height = 6, units = "in", res = 600)
print(p_heat); dev.off()
message("Significance heatmap saved.")

# ============================================================
# SECTION E — CONSOLE SUMMARY
# ============================================================
cat("\n")
cat("========================================================\n")
cat("  FILES SAVED\n")
cat("========================================================\n")
cat("  index_reserves_combined.png   — LOESS visualisation\n")
cat("  LMM_fixed_effects_summary.csv — full LMM coefficients\n")
cat("  LMM_marginal_fits.png         — LMM predicted trends\n")
cat("  LMM_interaction_coefplot.png  — effect sizes ± CI\n")
cat("  LMM_residual_diagnostics.png  — residual vs fitted\n")
cat("  LMM_qqplots.png               — normality of residuals\n")
cat("  GAMM_smooth_significance.csv  — GAMM smooth p-values\n")
cat("  GAMM_smooth_fits.png          — GAMM predicted curves\n")
cat("  GAMM_difference_smooths.png   — Reserve - No Reserve\n")
cat("  Significance_heatmap.png      — combined p-value tile\n")
cat("========================================================\n")
cat("\nINTERPRETATION GUIDE\n")
cat("--------------------------------------------------------\n")
cat("LMM category:Year_c term\n")
cat("  Significant (p<0.05) → the TWO TRENDS differ in SLOPE\n")
cat("  Positive coef        → Reserve is gaining faster over time\n")
cat("  Negative coef        → No Reserve gaining faster\n\n")
cat("GAMM smooth term (per category)\n")
cat("  edf >> 1             → non-linear temporal pattern\n")
cat("  Significant smooth   → the trend is NOT flat for that category\n\n")
cat("GAMM difference smooth chart\n")
cat("  Green segments       → years where Reserve ≠ No Reserve (sig.)\n")
cat("  Positive difference  → Reserve > No Reserve at that year\n")
cat("========================================================\n")
