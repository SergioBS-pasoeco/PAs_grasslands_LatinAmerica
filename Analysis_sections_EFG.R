# ============================================================
# SECTIONS E, F, G — Advanced confirmatory analyses
# Append this block to Analysis_two_indices_merged_V4_LMM_GAMM.R
#
# E: Segmented regression / changepoint analysis  (segmented)
# F: Functional Data Analysis                     (fda + fdapace)
# G: Bayesian GAM                                 (brms + tidybayes)
# ============================================================

# ── New packages needed ──────────────────────────────────────
# Install once if not present:
# install.packages(c("segmented", "fda", "fdapace",
#                    "brms", "tidybayes", "patchwork"))
#
# brms requires a C++ toolchain (Rtools on Windows) and
# Stan backend:  install.packages("cmdstanr", repos = "https://mc-stan.org/r-packages/")
#                cmdstanr::install_cmdstan()

library(segmented)   # breakpoint / changepoint in linear models
library(fda)         # functional data objects & FPCA
library(fdapace)     # functional PCA with FPCA()
library(brms)        # Bayesian regression via Stan
library(tidybayes)   # tidy Stan posterior draws
library(patchwork)   # combine ggplots with + / /

# ── Annual means (needed by all three sections) ──────────────
# One mean per country × index × category × Year
annual_means <- dat_all %>%
  group_by(country, index, category, Year) %>%
  summarise(mean_val = mean(value, na.rm = TRUE),
            .groups  = "drop")

yr_min <- min(dat_all$Year)   # 1998 typically

# ============================================================
# SECTION E — SEGMENTED REGRESSION (CHANGEPOINT ANALYSIS)
# ============================================================
# Strategy:
#   1. Fit a simple linear model:  mean_val ~ Year
#   2. Pass it to segmented() to detect 1 breakpoint
#   3. Do this separately for Reserve and No Reserve,
#      per country × index  (24 models total)
#   4. Collect breakpoint year + CI, slopes before/after,
#      and the Davies test p-value (tests H0: no breakpoint)
#   5. Plot fitted piecewise lines + breakpoint markers
# ============================================================

seg_results <- list()


for (ctry in levels(dat_all$country)) {
  for (idx in c("NDVI", "EVI")) {
    for (cat_lbl in c("Reserve", "No Reserve")) {

      key <- paste(ctry, idx, cat_lbl, sep = " | ")

      sub <- annual_means %>%
        filter(country == ctry, index == idx,
               category == cat_lbl) %>%
        arrange(Year)

      if (nrow(sub) < 6) next   # need enough points

      lm_base <- lm(mean_val ~ Year, data = sub)

      seg_fit <- tryCatch(
        segmented(lm_base,
                  seg.Z  = ~ Year,
                  psi = list(Year = median(sub$Year)), 
                  control = seg.control(it.max = 500,
                                        tol    = 1e-6,
                                        display = FALSE)),
        error = function(e) NULL
      )

      if (is.null(seg_fit)) {
        seg_results[[key]] <- data.frame(
          country  = ctry, index = idx, category = cat_lbl,
          bp_year  = NA, bp_lo = NA, bp_hi = NA,
          slope_b4 = NA, slope_af = NA,
          davies_p = NA, converged = FALSE
        )
        next
      }

      # Davies test: H0 = no breakpoint
      dav_p <- tryCatch(
        davies.test(lm_base, seg.Z = ~ Year)$p.value,
        error = function(e) NA_real_
      )

      # Breakpoint estimate and 95% CI — guarded against outdistanced BPs
      bp <- tryCatch(seg_fit$psi[, "Est."], error = function(e) NA_real_)
      
      bp_ci_bp <- tryCatch({
        ci <- confint(seg_fit)
        ci[grep("psi", rownames(ci)), , drop = FALSE]
      }, error = function(e) matrix(NA_real_, nrow = 1, ncol = 2))
      
      bp_lo <- tryCatch(bp_ci_bp[1, 1], error = function(e) NA_real_)
      bp_hi <- tryCatch(bp_ci_bp[1, 2], error = function(e) NA_real_)
      
      sl <- tryCatch(slope(seg_fit)$Year[, "Est."], error = function(e) c(NA_real_, NA_real_))
      
      seg_results[[key]] <- data.frame(
        country   = ctry,
        index     = idx,
        category  = cat_lbl,
        bp_year   = bp[1],
        bp_lo     = bp_lo,
        bp_hi     = bp_hi,
        slope_b4  = sl[1],
        slope_af  = sl[2],
        davies_p  = dav_p,
        converged = !is.na(bp[1])
      )
    }
  }
}

seg_df <- bind_rows(seg_results)
seg_df$country  <- factor(seg_df$country,  levels = levels(dat_all$country))
seg_df$category <- factor(seg_df$category, levels = c("Reserve", "No Reserve"))
seg_df$sig_bp   <- ifelse(!is.na(seg_df$davies_p) & seg_df$davies_p < 0.05,
                          "Significant breakpoint", "No breakpoint")
seg_df

write.csv(seg_df, "Segmented_changepoint_results.csv", row.names = FALSE)

cat("\n========================================================\n")
cat("  SEGMENTED REGRESSION — BREAKPOINT SUMMARY\n")
cat("========================================================\n")
print(kable(seg_df %>%
              select(country, index, category,
                     bp_year, bp_lo, bp_hi,
                     slope_b4, slope_af, davies_p, sig_bp),
            digits = 3, format = "simple"))

# ── E1. Generate piecewise fitted lines for plotting ─────────
seg_fitted <- bind_rows(lapply(names(seg_results), function(key) {
  res <- seg_results[[key]]
  if (!res$converged) return(NULL)

  parts   <- strsplit(key, " \\| ")[[1]]
  ctry    <- parts[1]
  idx     <- parts[2]
  cat_lbl <- parts[3]

  sub <- annual_means %>%
    filter(country == ctry, index == idx, category == cat_lbl) %>%
    arrange(Year)

  yr_seq <- seq(min(sub$Year), max(sub$Year), length.out = 200)

  # Piecewise fit manually from slopes + breakpoint
  bp   <- res$bp_year
  sl_b <- res$slope_b4
  sl_a <- res$slope_af
  icpt <- sub$mean_val[1] - sl_b * sub$Year[1]   # approx intercept

  fitted_val <- ifelse(
    yr_seq <= bp,
    icpt + sl_b * yr_seq,
    icpt + sl_b * bp + sl_a * (yr_seq - bp)
  )

  data.frame(
    Year     = yr_seq,
    fitted   = fitted_val,
    country  = ctry,
    index    = idx,
    category = cat_lbl
  )
}))

seg_fitted$country  <- factor(seg_fitted$country,  levels = levels(dat_all$country))
seg_fitted$category <- factor(seg_fitted$category, levels = c("Reserve", "No Reserve"))

# ── E2. Changepoint chart ─────────────────────────────────────
p_seg <- ggplot() +

  # Raw annual means
  geom_point(data = annual_means,
             aes(x = Year, y = mean_val, color = category),
             size = 1.2, alpha = 0.55) +

  # Piecewise fitted lines
  geom_line(data = seg_fitted,
            aes(x = Year, y = fitted, color = category),
            linewidth = 1.0) +

  # Breakpoint vertical lines — only significant ones
  geom_vline(data = seg_df %>% filter(sig_bp == "Significant breakpoint"),
             aes(xintercept = bp_year, color = category),
             linetype = "dashed", linewidth = 0.7, alpha = 0.8) +

  # Breakpoint CI band
  geom_rect(data = seg_df %>% filter(sig_bp == "Significant breakpoint"),
            aes(xmin = bp_lo, xmax = bp_hi,
                ymin = -Inf, ymax = Inf,
                fill = category),
            alpha = 0.07, inherit.aes = FALSE) +

  # Davies p annotation
  geom_text(data = seg_df %>%
              filter(converged) %>%
              mutate(label = paste0(category, ": bp=",
                                    round(bp_year, 1),
                                    "\np=", formatC(davies_p,
                                                    digits = 3,
                                                    format = "f"))) %>%
              group_by(country, index) %>%
              summarise(label = paste(label, collapse = "\n"),
                        .groups = "drop"),
            aes(label = label),
            x = -Inf, y = Inf,
            hjust = -0.05, vjust = 1.4,
            size = 2.3, fontface = "italic",
            color = "grey25", inherit.aes = FALSE) +

  facet_grid(country ~ index, scales = "free_y") +
  scale_color_manual(values = cat_cols, name = "Category") +
  scale_fill_manual(values  = cat_cols, name = "Category") +
  scale_x_continuous(breaks = seq(1998, 2024, by = 4)) +
  labs(
    title    = "Segmented regression — changepoint analysis",
    subtitle = "Dashed line + shaded band = breakpoint year ± 95% CI (Davies p < 0.05 only)",
    x        = "Year",
    y        = "Index value (annual mean)"
  ) +
  theme_bw(base_size = 10) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold", size = 9),
    legend.position  = "bottom",
    legend.title     = element_text(face = "bold"),
    axis.text.x      = element_text(angle = 45, vjust = 0.5),
    plot.title       = element_text(face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(size = 8, hjust = 0.5,
                                    color = "grey40"),
    panel.spacing    = unit(0.7, "lines")
  )
p_seg

png("Segmented_changepoint_chart.png",
    width = 9, height = 10, units = "in", res = 600)
print(p_seg); dev.off()
message("Segmented regression chart saved.")

# ── E3. Slope-change summary chart ───────────────────────────
# Shows before/after slope per model as a dot + arrow,
# making it easy to see where the trajectory accelerated/reversed

slope_long <- seg_df %>%
  filter(converged, !is.na(davies_p), davies_p < 0.05) %>%
  select(country, index, category, slope_b4, slope_af, bp_year) %>%
  pivot_longer(cols = c(slope_b4, slope_af),
               names_to  = "period",
               values_to = "slope") %>%
  mutate(
    period  = ifelse(period == "slope_b4", "Before breakpoint",
                     "After breakpoint"),
    period  = factor(period,
                     levels = c("Before breakpoint",
                                "After breakpoint")),
    label   = paste(country, index, sep = "\n")
  )

if (nrow(slope_long) > 0) {

  p_slope <- ggplot(slope_long,
                    aes(x = period, y = slope,
                        color = category,
                        group = interaction(category, label))) +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "grey60") +
    geom_line(linewidth = 0.8, alpha = 0.7) +
    geom_point(size = 3) +
    facet_wrap(~ label, scales = "free_y") +
    scale_color_manual(values = cat_cols, name = "Category") +
    labs(
      title    = "Slope change before vs after breakpoint",
      subtitle = "Only models with significant Davies test (p < 0.05)",
      x        = NULL,
      y        = "Annual slope (index units / year)"
    ) +
    theme_bw(base_size = 10) +
    theme(
      strip.background = element_rect(fill = "grey90"),
      strip.text       = element_text(face = "bold", size = 9),
      legend.position  = "bottom",
      legend.title     = element_text(face = "bold"),
      plot.title       = element_text(face = "bold", hjust = 0.5),
      plot.subtitle    = element_text(size = 8, hjust = 0.5,
                                      color = "grey40")
    )

  png("Segmented_slope_change.png",
      width = 8, height = 6, units = "in", res = 600)
  print(p_slope); dev.off()
  message("Slope change chart saved.")

} else {
  message("No significant breakpoints found — slope chart skipped.")
}

# ============================================================
# SECTION F — FUNCTIONAL DATA ANALYSIS (FDA)
# ============================================================
# Strategy:
#   Each *site* (nombre) contributes a time series of 26 annual
#   mean values → this is treated as a functional observation.
#
#   Steps:
#   1. Pivot to wide: rows = sites, cols = years
#   2. Smooth each site's time series into a functional data
#      object using B-spline basis (fda package)
#   3. Run Functional PCA (FPCA) to find main modes of
#      temporal variation across sites
#   4. Compare PC scores between Reserve and No Reserve
#      with a Wilcoxon test (non-parametric, robust)
#   5. Plot mean functions per category + first 2 FPCs
# ============================================================

# ── F1. Build site-level annual time series ──────────────────
site_annual <- dat_all %>%
  group_by(country, index, category, nombre, Year) %>%
  summarise(val = mean(value, na.rm = TRUE), .groups = "drop")

fda_results <- list()

for (ctry in levels(dat_all$country)) {
  for (idx in c("NDVI", "EVI")) {

    key <- paste(ctry, idx, sep = " | ")

    sub <- site_annual %>%
      filter(country == ctry, index == idx) %>%
      arrange(nombre, Year)

    # Wide matrix: rows = sites, cols = years
    wide <- sub %>%
      pivot_wider(id_cols    = c(nombre, category),
                  names_from = Year,
                  values_from = val) %>%
      arrange(nombre)

    yr_cols <- sort(unique(sub$Year))
    mat     <- as.matrix(wide[, as.character(yr_cols)])

    # Drop sites with any NA (incomplete time series)
    complete_rows <- complete.cases(mat)
    mat      <- mat[complete_rows, ]
    cats     <- wide$category[complete_rows]
    sites    <- wide$nombre[complete_rows]

    if (nrow(mat) < 4) {
      message("Skipping FDA for ", key, " — too few complete sites")
      next
    }

    # ── B-spline basis: 10 basis functions over 26 years ──
    t_pts  <- yr_cols
    n_basis <- min(10, length(t_pts) - 2)
    basis   <- create.bspline.basis(rangeval = range(t_pts),
                                    nbasis   = n_basis)
    fd_obj  <- smooth.basis(t_pts, t(mat), basis)$fd

    # ── Functional PCA ────────────────────────────────────
    fpca_res <- pca.fd(fd_obj, nharm = 3)

    # PC scores per site
    scores_df <- as.data.frame(fpca_res$scores)
    names(scores_df) <- paste0("PC", 1:3)
    scores_df$category <- cats
    scores_df$nombre   <- sites
    scores_df$country  <- ctry
    scores_df$index    <- idx

    # Variance explained
    var_exp <- fpca_res$varprop * 100

    # Wilcoxon test on PC1 scores: Reserve vs No Reserve
    pc1_res  <- scores_df$PC1[scores_df$category == "Reserve"]
    pc1_nres <- scores_df$PC1[scores_df$category == "No Reserve"]

    wx_p <- tryCatch(
      wilcox.test(pc1_res, pc1_nres)$p.value,
      error = function(e) NA_real_
    )

    fda_results[[key]] <- list(
      fd_obj    = fd_obj,
      fpca_res  = fpca_res,
      scores_df = scores_df,
      var_exp   = var_exp,
      wx_p      = wx_p,
      yr_cols   = yr_cols,
      cats      = cats,
      ctry      = ctry,
      idx       = idx
    )
  }
}

# ── F2. Mean function plot per category ───────────────────────
# Evaluate mean functional curve for Reserve and No Reserve

fda_mean_df <- bind_rows(lapply(names(fda_results), function(key) {
  res  <- fda_results[[key]]
  ctry <- res$ctry
  idx  <- res$idx

  yr_fine <- seq(min(res$yr_cols), max(res$yr_cols), length.out = 200)
  cats_u  <- c("Reserve", "No Reserve")

  bind_rows(lapply(cats_u, function(cat_lbl) {
    idx_cat <- which(res$cats == cat_lbl)
    if (length(idx_cat) == 0) return(NULL)

    # Mean of evaluated functional objects for this category
    fd_cat  <- res$fd_obj[idx_cat]
    eval_mat <- eval.fd(yr_fine, fd_cat)   # rows=time, cols=sites
    mean_curve <- rowMeans(eval_mat, na.rm = TRUE)
    se_curve   <- apply(eval_mat, 1, sd, na.rm = TRUE) /
                  sqrt(ncol(eval_mat))

    data.frame(
      Year     = yr_fine,
      mean_val = mean_curve,
      se       = se_curve,
      lo       = mean_curve - 1.96 * se_curve,
      hi       = mean_curve + 1.96 * se_curve,
      category = cat_lbl,
      country  = ctry,
      index    = idx
    )
  }))
}))

fda_mean_df$country  <- factor(fda_mean_df$country,
                                levels = levels(dat_all$country))
fda_mean_df$category <- factor(fda_mean_df$category,
                                levels = c("Reserve", "No Reserve"))

# Wilcoxon summary table
wx_summary <- bind_rows(lapply(names(fda_results), function(key) {
  res <- fda_results[[key]]
  data.frame(
    country  = res$ctry,
    index    = res$idx,
    PC1_var  = round(res$var_exp[1], 1),
    PC2_var  = round(res$var_exp[2], 1),
    wx_PC1_p = res$wx_p,
    sig      = case_when(
      is.na(res$wx_p)  ~ "n.s.",
      res$wx_p < 0.001 ~ "***",
      res$wx_p < 0.01  ~ "**",
      res$wx_p < 0.05  ~ "*",
      TRUE             ~ "n.s."
    )
  )
}))

wx_summary$country <- factor(wx_summary$country,
                              levels = levels(dat_all$country))

cat("\n========================================================\n")
cat("  FDA — FPCA VARIANCE EXPLAINED & WILCOXON PC1 TEST\n")
cat("========================================================\n")
print(kable(wx_summary, digits = 4, format = "simple"))
write.csv(wx_summary, "FDA_FPCA_wilcoxon_summary.csv", row.names = FALSE)

# ── F3. Mean function chart ───────────────────────────────────
p_fda_mean <- ggplot(fda_mean_df,
                     aes(x = Year, y = mean_val,
                         color = category, fill = category)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15,
              color = NA) +
  geom_line(linewidth = 1.1) +

  # Wilcoxon annotation
  geom_text(data = wx_summary %>%
              mutate(label = paste0("Wilcoxon PC1\np=",
                                    formatC(wx_PC1_p,
                                            digits = 3,
                                            format = "f"),
                                    " ", sig)),
            aes(label = label),
            x = -Inf, y = Inf,
            hjust = -0.05, vjust = 1.4,
            size = 2.5, fontface = "italic",
            color = "grey25", inherit.aes = FALSE) +

  facet_grid(country ~ index, scales = "free_y") +
  scale_color_manual(values = cat_cols, name = "Category") +
  scale_fill_manual(values  = cat_cols, name = "Category") +
  scale_x_continuous(breaks = seq(1998, 2024, by = 4)) +
  labs(
    title    = "FDA — mean functional curves per category",
    subtitle = "Smoothed B-spline mean ± 95% CI across sites  |  Wilcoxon test on PC1 scores",
    x        = "Year",
    y        = "Index value"
  ) +
  theme_bw(base_size = 10) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold", size = 9),
    legend.position  = "bottom",
    legend.title     = element_text(face = "bold"),
    axis.text.x      = element_text(angle = 45, vjust = 0.5),
    plot.title       = element_text(face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(size = 8, hjust = 0.5,
                                    color = "grey40"),
    panel.spacing    = unit(0.7, "lines")
  )

png("FDA_mean_functions.png",
    width = 9, height = 10, units = "in", res = 600)
print(p_fda_mean); dev.off()
message("FDA mean function chart saved.")

# ── F4. PC1 score distribution chart ─────────────────────────
scores_all <- bind_rows(lapply(fda_results, `[[`, "scores_df"))
scores_all$country  <- factor(scores_all$country,
                               levels = levels(dat_all$country))
scores_all$category <- factor(scores_all$category,
                               levels = c("Reserve", "No Reserve"))

p_pc1 <- ggplot(scores_all,
                aes(x = category, y = PC1,
                    fill = category, color = category)) +
  geom_violin(alpha = 0.4, trim = TRUE) +
  geom_boxplot(width = 0.15, alpha = 0.6,
               outlier.shape = NA, color = "grey30") +
  geom_jitter(width = 0.07, size = 0.9, alpha = 0.5) +
  facet_grid(country ~ index, scales = "free_y") +
  scale_fill_manual(values  = cat_cols, name = "Category") +
  scale_color_manual(values = cat_cols, name = "Category") +
  labs(
    title    = "FDA — PC1 score distribution by category",
    subtitle = "PC1 captures the dominant mode of temporal variation across sites",
    x        = NULL,
    y        = "PC1 score"
  ) +
  theme_bw(base_size = 10) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold", size = 9),
    legend.position  = "bottom",
    legend.title     = element_text(face = "bold"),
    plot.title       = element_text(face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(size = 8, hjust = 0.5,
                                    color = "grey40"),
    panel.spacing    = unit(0.7, "lines")
  )

png("FDA_PC1_scores.png",
    width = 9, height = 10, units = "in", res = 600)
print(p_pc1); dev.off()
message("FDA PC1 score chart saved.")

# ── F5. FPC1 and FPC2 eigenfunctions ─────────────────────────
# Shows WHAT temporal pattern each PC captures

fpc_df <- bind_rows(lapply(names(fda_results), function(key) {
  res     <- fda_results[[key]]
  yr_fine <- seq(min(res$yr_cols), max(res$yr_cols), length.out = 200)

  bind_rows(lapply(1:2, function(pc_n) {
    fpc_vals <- eval.fd(yr_fine, res$fpca_res$harmonics[pc_n])
    data.frame(
      Year    = yr_fine,
      fpc_val = fpc_vals,
      PC      = paste0("FPC", pc_n,
                       " (", round(res$var_exp[pc_n], 1), "%)"),
      country = res$ctry,
      index   = res$idx
    )
  }))
}))

fpc_df$country <- factor(fpc_df$country, levels = levels(dat_all$country))

p_fpc <- ggplot(fpc_df,
                aes(x = Year, y = fpc_val,
                    color = PC, linetype = PC)) +
  geom_hline(yintercept = 0, color = "grey60",
             linetype = "dashed", linewidth = 0.5) +
  geom_line(linewidth = 1.0) +
  facet_grid(country ~ index, scales = "free_y") +
  scale_color_manual(values = c("FPC1 (\\d+\\.\\d+%)" = "#1B9E77",
                                 "FPC2 (\\d+\\.\\d+%)" = "#7570B3"),
                     aesthetics = c("color")) +
  scale_x_continuous(breaks = seq(1998, 2024, by = 4)) +
  labs(
    title    = "FDA — first two functional principal components",
    subtitle = "FPC1 = dominant temporal mode  |  FPC2 = secondary mode",
    x        = "Year",
    y        = "Eigenfunction value",
    color    = NULL, linetype = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold", size = 9),
    legend.position  = "bottom",
    axis.text.x      = element_text(angle = 45, vjust = 0.5),
    plot.title       = element_text(face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(size = 8, hjust = 0.5,
                                    color = "grey40"),
    panel.spacing    = unit(0.7, "lines")
  )

png("FDA_eigenfunctions.png",
    width = 9, height = 10, units = "in", res = 600)
print(p_fpc); dev.off()
message("FDA eigenfunction chart saved.")

# ============================================================
# SECTION G — BAYESIAN GAM  (brms + Stan)
# ============================================================
# Strategy:
#   Same model structure as the frequentist GAMM in Section C,
#   but fitted in a Bayesian framework:
#
#   value ~ category + s(Year_c, by = category, k = 10)
#           + (1 | nombre)
#
#   Priors (weakly informative):
#     Intercept  ~ Normal(0.5, 0.5)   — index values in [0,1]
#     b          ~ Normal(0, 0.1)     — smooth coefficients
#     sigma      ~ Exponential(10)
#     sd(nombre) ~ Exponential(10)
#
#   Inference:
#     4 chains × 2000 iterations (1000 warmup)
#     Rhat convergence check (< 1.05)
#     Posterior predictive check
#     Credible interval comparison Reserve vs No Reserve
#
#   Run per country × index (6 models).
#   NOTE: each model takes ~5–15 min. Run overnight or
#         reduce iter to 1000 for quick checks.
# ============================================================

# ── G0. Priors ───────────────────────────────────────────────
bayes_priors <- c(
  prior(normal(0.5, 0.5), class = Intercept),
  prior(normal(0,   0.1), class = b),
  prior(exponential(10),  class = sigma),
  prior(exponential(10),  class = sd)
)

bayes_models <- list()

for (ctry in levels(dat_all$country)) {
  for (idx in c("NDVI", "EVI")) {

    key <- paste(ctry, idx, sep = " | ")
    message("\nFitting Bayesian GAM: ", key)

    sub <- dat_all %>%
      filter(country == ctry, index == idx) %>%
      mutate(nombre_f = factor(nombre))

    bayes_models[[key]] <- brm(
      formula  = bf(value ~
                      category +
                      s(Year_c, by = category, k = 10) +
                      (1 | nombre_f)),
      data     = sub,
      family   = gaussian(),
      prior    = bayes_priors,
      chains   = 4,
      iter     = 2000,
      warmup   = 1000,
      cores    = 4,          # set to number of CPU cores available
      seed     = 42,
      backend  = "cmdstanr", # or "rstan" if cmdstanr not installed
      silent   = 2,
      file     = paste0("bayes_gam_",
                        gsub("[^A-Za-z]", "_", key), ".rds")
      # file = saves/reloads the model automatically — avoids
      # refitting if you re-run the script
    )
  }
}

# ── G1. Convergence diagnostics ──────────────────────────────
conv_df <- bind_rows(lapply(names(bayes_models), function(key) {
  m     <- bayes_models[[key]]
  parts <- strsplit(key, " \\| ")[[1]]
  rh    <- rhat(m)

  data.frame(
    country     = parts[1],
    index       = parts[2],
    max_rhat    = round(max(rh, na.rm = TRUE), 4),
    n_rhat_bad  = sum(rh > 1.05, na.rm = TRUE),
    converged   = max(rh, na.rm = TRUE) < 1.05
  )
}))

cat("\n========================================================\n")
cat("  BAYESIAN GAM — CONVERGENCE (Rhat < 1.05 = good)\n")
cat("========================================================\n")
print(kable(conv_df, format = "simple"))

# ── G2. Posterior predictive check — one panel per model ─────
png("BayesGAM_pp_checks.png",
    width = 12, height = 8, units = "in", res = 300)
par(mfrow = c(2, 3))
for (key in names(bayes_models)) {
  pp_check(bayes_models[[key]],
           type  = "dens_overlay",
           ndraws = 50) +
    ggtitle(key) +
    theme_bw(base_size = 9)
}
dev.off()
message("Bayesian GAM PP checks saved.")

# ── G3. Posterior marginal predictions ───────────────────────
bayes_preds <- bind_rows(lapply(names(bayes_models), function(key) {
  m     <- bayes_models[[key]]
  parts <- strsplit(key, " \\| ")[[1]]
  ctry  <- parts[1]
  idx   <- parts[2]

  sub <- dat_all %>%
    filter(country == ctry, index == idx) %>%
    mutate(nombre_f = factor(nombre))

  yr_seq <- seq(min(sub$Year_c), max(sub$Year_c), length.out = 100)

  bind_rows(lapply(c("Reserve", "No Reserve"), function(cat_lbl) {

    nd <- data.frame(
      Year_c   = yr_seq,
      category = factor(cat_lbl, levels = levels(sub$category)),
      nombre_f = factor(levels(sub$nombre_f)[1],
                        levels = levels(sub$nombre_f))
    )

    # epred_draws: posterior of E[y] — excludes residual noise
    draws <- epred_draws(m, newdata = nd,
                         re_formula = NA,   # population-level only
                         ndraws = 500)

    draws %>%
      group_by(Year_c) %>%
      summarise(
        predicted = median(.epred),
        lo        = quantile(.epred, 0.025),
        hi        = quantile(.epred, 0.975),
        lo90      = quantile(.epred, 0.050),
        hi90      = quantile(.epred, 0.950),
        .groups   = "drop"
      ) %>%
      mutate(
        Year     = Year_c + yr_min,
        category = cat_lbl,
        country  = ctry,
        index    = idx
      )
  }))
}))

bayes_preds$country  <- factor(bayes_preds$country,
                                levels = levels(dat_all$country))
bayes_preds$category <- factor(bayes_preds$category,
                                levels = c("Reserve", "No Reserve"))

# ── G4. Bayesian GAM fit chart ───────────────────────────────
p_bayes <- ggplot() +

  # Raw annual means
  geom_line(data = annual_means,
            aes(x = Year, y = mean_val, color = category),
            linewidth = 0.35, alpha = 0.4) +

  # 95% credible interval
  geom_ribbon(data = bayes_preds,
              aes(x = Year, ymin = lo, ymax = hi,
                  fill = category),
              alpha = 0.12) +

  # 90% credible interval (inner, darker)
  geom_ribbon(data = bayes_preds,
              aes(x = Year, ymin = lo90, ymax = hi90,
                  fill = category),
              alpha = 0.18) +

  # Posterior median line
  geom_line(data = bayes_preds,
            aes(x = Year, y = predicted, color = category),
            linewidth = 1.2) +

  facet_grid(country ~ index, scales = "free_y") +
  scale_color_manual(values = cat_cols, name = "Category") +
  scale_fill_manual(values  = cat_cols, name = "Category") +
  scale_x_continuous(breaks = seq(1998, 2024, by = 4)) +
  labs(
    title    = "Bayesian GAM — posterior predictions",
    subtitle = "Median posterior  |  Inner band = 90% CrI  |  Outer band = 95% CrI",
    x        = "Year",
    y        = "Index value"
  ) +
  theme_bw(base_size = 10) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold", size = 9),
    legend.position  = "bottom",
    legend.title     = element_text(face = "bold"),
    axis.text.x      = element_text(angle = 45, vjust = 0.5),
    plot.title       = element_text(face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(size = 8, hjust = 0.5,
                                    color = "grey40"),
    panel.spacing    = unit(0.7, "lines")
  )

png("BayesGAM_posterior_fits.png",
    width = 9, height = 10, units = "in", res = 600)
print(p_bayes); dev.off()
message("Bayesian GAM posterior fit chart saved.")

# ── G5. Posterior difference: P(Reserve > No Reserve) ────────
# For each year point, compute the posterior probability that
# Reserve > No Reserve — the Bayesian analogue of a p-value.
# P > 0.95 or < 0.05 is conventionally "strong evidence".

bayes_diff <- bind_rows(lapply(names(bayes_models), function(key) {
  m     <- bayes_models[[key]]
  parts <- strsplit(key, " \\| ")[[1]]
  ctry  <- parts[1]
  idx   <- parts[2]

  sub <- dat_all %>%
    filter(country == ctry, index == idx) %>%
    mutate(nombre_f = factor(nombre))

  yr_seq <- seq(min(sub$Year_c), max(sub$Year_c), length.out = 100)

  make_nd <- function(cat_lbl) {
    data.frame(
      Year_c   = yr_seq,
      category = factor(cat_lbl, levels = levels(sub$category)),
      nombre_f = factor(levels(sub$nombre_f)[1],
                        levels = levels(sub$nombre_f))
    )
  }

  draws_res  <- epred_draws(m, newdata = make_nd("Reserve"),
                             re_formula = NA, ndraws = 500)
  draws_nres <- epred_draws(m, newdata = make_nd("No Reserve"),
                             re_formula = NA, ndraws = 500)

  # Align draws by Year_c and draw index
  diff_draws <- draws_res %>%
    rename(pred_res = .epred) %>%
    bind_cols(draws_nres %>%
                select(pred_nres = .epred)) %>%
    mutate(diff = pred_res - pred_nres)

  diff_draws %>%
    group_by(Year_c) %>%
    summarise(
      prob_reserve_gt = mean(diff > 0),
      diff_median     = median(diff),
      diff_lo         = quantile(diff, 0.025),
      diff_hi         = quantile(diff, 0.975),
      .groups         = "drop"
    ) %>%
    mutate(
      Year    = Year_c + yr_min,
      country = ctry,
      index   = idx,
      # Evidence threshold
      evidence = case_when(
        prob_reserve_gt > 0.95 ~ "Strong: Reserve > No Reserve",
        prob_reserve_gt < 0.05 ~ "Strong: No Reserve > Reserve",
        TRUE                   ~ "Uncertain"
      )
    )
}))

bayes_diff$country  <- factor(bayes_diff$country,
                               levels = levels(dat_all$country))
bayes_diff$evidence <- factor(bayes_diff$evidence,
                               levels = c("Strong: Reserve > No Reserve",
                                          "Uncertain",
                                          "Strong: No Reserve > Reserve"))

ev_cols <- c("Strong: Reserve > No Reserve" = "#1B9E77",
             "Uncertain"                     = "grey70",
             "Strong: No Reserve > Reserve"  = "#D95F02")

p_prob <- ggplot(bayes_diff,
                 aes(x = Year, y = prob_reserve_gt,
                     color = evidence)) +
  geom_hline(yintercept = c(0.05, 0.50, 0.95),
             linetype = c("dashed", "solid", "dashed"),
             color    = c("grey50", "grey30", "grey50"),
             linewidth = 0.5) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 0.8, alpha = 0.6) +
  annotate("text", x = -Inf, y = 0.97,
           label = "P = 0.95 threshold",
           hjust = -0.1, size = 2.5, color = "grey40") +
  annotate("text", x = -Inf, y = 0.03,
           label = "P = 0.05 threshold",
           hjust = -0.1, size = 2.5, color = "grey40") +
  facet_grid(country ~ index) +
  scale_color_manual(values = ev_cols, name = NULL) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = scales::percent) +
  scale_x_continuous(breaks = seq(1998, 2024, by = 4)) +
  labs(
    title    = "Bayesian GAM — P(Reserve > No Reserve) over time",
    subtitle = "Green above 95% line = strong Bayesian evidence Reserve has higher index",
    x        = "Year",
    y        = "Posterior probability Reserve > No Reserve"
  ) +
  theme_bw(base_size = 10) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold", size = 9),
    legend.position  = "bottom",
    axis.text.x      = element_text(angle = 45, vjust = 0.5),
    plot.title       = element_text(face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(size = 8, hjust = 0.5,
                                    color = "grey40"),
    panel.spacing    = unit(0.7, "lines")
  )

png("BayesGAM_posterior_probability.png",
    width = 9, height = 10, units = "in", res = 600)
print(p_prob); dev.off()
message("Bayesian posterior probability chart saved.")

# ============================================================
# SECTION H — FINAL CONSOLE SUMMARY
# ============================================================
cat("\n")
cat("========================================================\n")
cat("  ALL FILES SAVED — SECTIONS E, F, G\n")
cat("========================================================\n")
cat("  SECTION E — Segmented regression\n")
cat("    Segmented_changepoint_results.csv\n")
cat("    Segmented_changepoint_chart.png\n")
cat("    Segmented_slope_change.png\n\n")
cat("  SECTION F — Functional Data Analysis\n")
cat("    FDA_FPCA_wilcoxon_summary.csv\n")
cat("    FDA_mean_functions.png\n")
cat("    FDA_PC1_scores.png\n")
cat("    FDA_eigenfunctions.png\n\n")
cat("  SECTION G — Bayesian GAM\n")
cat("    bayes_gam_*.rds          (one per country × index)\n")
cat("    BayesGAM_pp_checks.png\n")
cat("    BayesGAM_posterior_fits.png\n")
cat("    BayesGAM_posterior_probability.png\n")
cat("========================================================\n")
cat("\nINTERPRETATION GUIDE — SECTIONS E, F, G\n")
cat("--------------------------------------------------------\n")
cat("SECTION E — Segmented regression\n")
cat("  Davies p < 0.05  → breakpoint is statistically supported\n")
cat("  bp_year          → year when trend changed direction/rate\n")
cat("  slope_b4/af      → annual rate before and after breakpoint\n\n")
cat("SECTION F — FDA\n")
cat("  FPC1             → dominant mode of variation across sites\n")
cat("  Wilcoxon p < 0.05→ Reserve and No Reserve sites differ in\n")
cat("                     their dominant temporal trajectory\n")
cat("  Eigenfunction shape → when and how sites diverge in time\n\n")
cat("SECTION G — Bayesian GAM\n")
cat("  Rhat < 1.05      → chains converged (model trustworthy)\n")
cat("  PP check         → posterior distribution covers observed data\n")
cat("  P(R > NR) > 0.95 → strong Bayesian evidence Reserve higher\n")
cat("  CrI              → credible interval (Bayesian equivalent of CI)\n")
cat("========================================================\n")
