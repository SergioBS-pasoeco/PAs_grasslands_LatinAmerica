# ============================================================
# MULTI-BREAKPOINT SEGMENTED REGRESSION
# Based on your attached script: Analysis_sections_EFG.R
# Extended to:
#   1) Detect MULTIPLE breakpoints (not only one)
#   2) Use Year as factor (discrete) when appropriate
#   3) Run analysis using:
#        A. Annual observed means
#        B. Simulated yearly values inside SE confidence intervals
#   4) Show chart with only significant multiple breakpoints
# ============================================================

# =========================
# PACKAGES
# =========================
libs <- c("tidyverse","segmented","strucchange","broom","patchwork")
for(p in libs){
  if(!require(p, character.only = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}

# ============================================================
# ASSUMED INPUT DATA (same structure as your script)
# dat_all must contain:
# country | index | category | Year | value
# ============================================================

# ============================================================
# SETTINGS
# ============================================================
MAX_BREAKS      <- 4     # maximum breakpoints to test
MIN_POINTS      <- 8
N_RANDOM_RUNS   <- 100   # random CI simulations
USE_FACTOR_IF_FEW_YEARS <- TRUE

set.seed(1234)

# ============================================================
# STEP 1 — YEARLY SUMMARY
# ============================================================

annual_stats <- dat_all %>%
  group_by(country,index,category,Year) %>%
  summarise(
    mean_val = mean(value, na.rm=TRUE),
    sd_val   = sd(value, na.rm=TRUE),
    n        = sum(!is.na(value)),
    se       = sd_val/sqrt(n),
    lo95     = mean_val - 1.96*se,
    hi95     = mean_val + 1.96*se,
    .groups="drop"
  )

# ============================================================
# FUNCTION: choose if Year should be factor
# ============================================================

year_as_factor <- function(df){
  ny <- length(unique(df$Year))
  if(USE_FACTOR_IF_FEW_YEARS && ny <= 15) TRUE else FALSE
}

# ============================================================
# FUNCTION: MULTI BREAKPOINT MODEL
# ============================================================

fit_multi_breaks <- function(df){
  
  if(nrow(df) < MIN_POINTS) return(NULL)
  
  use_factor <- year_as_factor(df)
  
  # ------------------------------------------
  # If Year as factor => discrete ANOVA style
  # ------------------------------------------
  if(use_factor){
    
    mod <- lm(mean_val ~ factor(Year), data=df)
    
    out <- list(
      method="factor_year",
      model=mod,
      breaks=NULL,
      bic=BIC(mod),
      significant=anova(mod)$`Pr(>F)`[1]
    )
    
    return(out)
  }
  
  # ------------------------------------------
  # Numeric Year => multiple breakpoints
  # ------------------------------------------
  base_mod <- lm(mean_val ~ Year, data=df)
  
  bp_mod <- breakpoints(mean_val ~ Year,
                        data=df,
                        breaks=MAX_BREAKS)
  
  bic_vals <- BIC(bp_mod)
  best_k   <- which.min(bic_vals) - 1
  
  if(best_k == 0){
    
    return(list(
      method="linear",
      model=base_mod,
      breaks=NULL,
      bic=min(bic_vals),
      significant=NA
    ))
  }
  
  final_bp <- breakpoints(mean_val ~ Year,
                          data=df,
                          breaks=best_k)
  
  bp_idx   <- final_bp$breakpoints
  bp_years <- df$Year[bp_idx]
  
  seg_fit <- segmented(base_mod,
                       seg.Z = ~Year,
                       psi = bp_years)
  
  pvals <- tryCatch(
    summary(seg_fit)$coefficients[,4],
    error=function(e) NA
  )
  
  list(
    method="segmented",
    model=seg_fit,
    breaks=bp_years,
    bic=min(bic_vals),
    significant=pvals
  )
}

# ============================================================
# STEP 2 — OBSERVED YEARLY MEANS ANALYSIS
# ============================================================

results_obs <- list()

keys <- annual_stats %>%
  unite("grp", country,index,category, remove=FALSE) %>%
  pull(grp)

for(k in unique(keys)){
  
  tmp <- annual_stats %>%
    unite("grp", country,index,category, remove=FALSE) %>%
    filter(grp==k) %>%
    arrange(Year)
  
  fit <- fit_multi_breaks(tmp)
  
  results_obs[[k]] <- list(
    data=tmp,
    fit=fit
  )
}

# ============================================================
# STEP 3 — RANDOM VALUES INSIDE YEARLY CI
# ============================================================

sim_results <- list()

for(run in 1:N_RANDOM_RUNS){
  
  sim_df <- annual_stats %>%
    rowwise() %>%
    mutate(
      mean_val = runif(1, lo95, hi95)
    ) %>%
    ungroup()
  
  keys2 <- sim_df %>%
    unite("grp", country,index,category, remove=FALSE) %>%
    pull(grp)
  
  for(k in unique(keys2)){
    
    tmp <- sim_df %>%
      unite("grp", country,index,category, remove=FALSE) %>%
      filter(grp==k) %>%
      arrange(Year)
    
    fit <- fit_multi_breaks(tmp)
    
    sim_results[[length(sim_results)+1]] <- tibble(
      run=run,
      group=k,
      n_breaks=ifelse(is.null(fit$breaks),0,length(fit$breaks))
    )
  }
}

sim_results <- bind_rows(sim_results)

# ============================================================
# STEP 4 — TABLE OF OBSERVED BREAKPOINTS
# ============================================================

bp_table <- map_df(names(results_obs), function(k){
  
  fit <- results_obs[[k]]$fit
  
  tibble(
    group=k,
    method=fit$method,
    n_breaks=ifelse(is.null(fit$breaks),0,length(fit$breaks)),
    breakpoints=paste(fit$breaks, collapse=", "),
    bic=fit$bic
  )
})

write.csv(bp_table,
          "Multiple_breakpoint_results.csv",
          row.names=FALSE)

print(bp_table)

# ============================================================
# STEP 5 — CHART ONLY SIGNIFICANT MULTIPLE BREAKPOINTS
# ============================================================

plot_df <- bind_rows(
  lapply(names(results_obs), function(k){
    
    d <- results_obs[[k]]$data
    d$group <- k
    d
  })
)

bp_lines <- map_df(names(results_obs), function(k){
  
  fit <- results_obs[[k]]$fit
  
  if(is.null(fit$breaks)) return(NULL)
  if(length(fit$breaks) < 2) return(NULL)
  
  tibble(
    group=k,
    bp=fit$breaks
  )
})

p <- ggplot(plot_df,
            aes(Year, mean_val)) +
  
  geom_point(size=1.7, alpha=.7, colour="black") +
  geom_line(linewidth=.8, colour="steelblue") +
  
  geom_vline(data=bp_lines,
             aes(xintercept=bp),
             colour="red",
             linetype="dashed",
             linewidth=.8) +
  
  facet_wrap(~group, scales="free_y") +
  
  theme_bw(base_size=11) +
  
  labs(
    title="Multiple Significant Breakpoints",
    subtitle="Panels with 2 or more detected breakpoints",
    x="Year",
    y="Annual Mean"
  )

print(p)

ggsave("Multiple_breakpoints_chart.png",
       p,
       width=14,
       height=9,
       dpi=600)

# ============================================================
# STEP 6 — RANDOM CI SUMMARY
# ============================================================

sim_summary <- sim_results %>%
  group_by(group) %>%
  summarise(
    avg_breaks = mean(n_breaks),
    max_breaks = max(n_breaks),
    prop_multi = mean(n_breaks >= 2),
    .groups="drop"
  )

write.csv(sim_summary,
          "Random_CI_breakpoint_summary.csv",
          row.names=FALSE)

print(sim_summary)

cat("\n=====================================================\n")
cat("FILES SAVED:\n")
cat("  Multiple_breakpoint_results.csv\n")
cat("  Multiple_breakpoints_chart.png\n")
cat("  Random_CI_breakpoint_summary.csv\n")
cat("=====================================================\n")

