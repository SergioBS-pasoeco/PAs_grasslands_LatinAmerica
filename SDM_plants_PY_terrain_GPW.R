### SDM plants Paraguay - with terrain, GPW, accuracy metrics, and spatial blocks
### Based on SDM_plants_PY_new.R inputs/paths, extended with terrain/GPW predictors
### and permutation-based variable importance as in SDM_plants_CO_terrain_GPW.R

library(sf)
library(terra)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(dismo)
library(geodata)
library(maxnet)

setwd("H:/SDM_project")

# ---- Inputs ---------------------------------------------------------------
plants_pa  <- terra::vect("GBIF/Plants/Plants_PYG_PA.shp")
plants_opa <- terra::vect("GBIF/Plants/Plants_PYG_outPA.shp")
dat_pa     <- as.data.frame(terra::values(plants_pa))
dat_opa    <- as.data.frame(terra::values(plants_opa))
dat_all_PY <- bind_rows(dat_pa, dat_opa)
dat_all_PY$type    <- c(rep("Reserve", nrow(plants_pa)), rep("No Reserve", nrow(plants_opa)))
dat_all_PY$country <- "Paraguay"

gadm_PY <- terra::vect("Vectors/PY/Ecoregions_study_collapse_4326.shp")

# ---- Load rasters ---------------------------------------------------------

# NDVI (250 m)
ndvi_mean_path <- "data/raw/rasters/mean_ndvi_250m_PY.tif"
ndvi_sd_path   <- "data/raw/rasters/std_dev_ndvi_250m_PY.tif"
if (!file.exists(ndvi_mean_path) || !file.exists(ndvi_sd_path)) {
  warning("250 m NDVI files not found; falling back to 500 m NDVI files.")
  ndvi_mean_path <- "data/raw/rasters/mean_ndvi_500m_PY.tif"
  ndvi_sd_path   <- "data/raw/rasters/std_dev_ndvi_500m_PY.tif"
}
m_ndvi  <- terra::rast(ndvi_mean_path)
sd_ndvi <- terra::rast(ndvi_sd_path)
names(m_ndvi)  <- "ndvi_mean"
names(sd_ndvi) <- "ndvi_sd"

# Terrain stack (Paraguay, 250 m): elevation, slope, aspect, hillshade
terrain_py <- terra::rast("data/raw/rasters/terrain_stack_paraguay_resampled_250m.tif")
names(terrain_py) <- c("elevation", "slope", "aspect", "hillshade")

# Population density – GPW (250 m)
gpw_py <- terra::rast("data/raw/rasters/Resampled_250m_GPW_PY.tif")
names(gpw_py) <- "gpw_density"

# ---- Output directories ---------------------------------------------------
out_dir  <- "H:/SDM_project/GBIF/WK_results/SDM_plants/PY_terrain_GPW_250m"
fold_dir <- file.path(out_dir, "fold_metrics")

# ---- Helpers --------------------------------------------------------------

pick_coord_cols <- function(df) {
  lon_candidates <- c("dcmlLng", "decimalLon", "decimalLongitude", "longitude", "lon")
  lat_candidates <- c("dcmlLtt", "decimalLat", "decimalLatitude", "latitude", "lat")
  lon_col <- lon_candidates[lon_candidates %in% names(df)][1]
  lat_col <- lat_candidates[lat_candidates %in% names(df)][1]
  if (is.na(lon_col) || is.na(lat_col)) stop("Coordinate columns not found")
  c(lon_col, lat_col)
}

make_stratified_folds <- function(labels, k = 5, seed = 2025) {
  set.seed(seed)
  labels  <- as.integer(labels)
  pos_idx <- which(labels == 1)
  neg_idx <- which(labels == 0)
  k_eff   <- max(2, min(k, length(pos_idx), length(neg_idx)))
  folds   <- integer(length(labels))
  if (length(pos_idx) > 0) folds[pos_idx] <- sample(rep(seq_len(k_eff), length.out = length(pos_idx)))
  if (length(neg_idx) > 0) folds[neg_idx] <- sample(rep(seq_len(k_eff), length.out = length(neg_idx)))
  list(folds = folds, k = k_eff)
}

# Simple spatial blocks via k-means on coordinates
make_spatial_blocks <- function(coords_df, k_blocks = 5, seed = 2025) {
  set.seed(seed)
  k_blocks <- max(2, min(k_blocks, nrow(coords_df)))
  jittered <- coords_df
  if (nrow(unique(coords_df)) < nrow(coords_df)) {
    jittered$lon <- jitter(jittered$lon, amount = 1e-06)
    jittered$lat <- jitter(jittered$lat, amount = 1e-06)
  }
  stats::kmeans(jittered[, c("lon", "lat")], centers = k_blocks, iter.max = 50)$cluster
}

compute_auprc <- function(labels01, scores) {
  o      <- order(scores, decreasing = TRUE)
  y      <- labels01[o]
  tp     <- cumsum(y)
  fp     <- cumsum(1 - y)
  p      <- tp / (tp + fp)
  r      <- tp / sum(y)
  if (sum(y) == 0) return(NA_real_)
  keep   <- c(TRUE, diff(r) != 0)
  sum(diff(c(0, r[keep])) * p[keep])
}

compute_tss <- function(pos_scores, neg_scores) {
  e <- tryCatch(dismo::evaluate(p = pos_scores, a = neg_scores), error = function(x) NULL)
  if (is.null(e)) return(list(tss = NA_real_, threshold = NA_real_))
  tss_vals <- e@TPR + e@TNR - 1
  if (length(tss_vals) == 0 || all(!is.finite(tss_vals))) return(list(tss = NA_real_, threshold = NA_real_))
  i <- which.max(tss_vals)
  list(tss = tss_vals[i], threshold = e@t[i])
}

# ---- Boxplot helper -------------------------------------------------------
create_metric_boxplots <- function(fold_data, output_dir) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  metrics_to_plot <- c("roc_auc", "auprc", "auprc_norm", "tss")
  metric_names    <- c("ROC-AUC", "AUPRC", "nAUPRC", "TSS")

  for (i in seq_along(metrics_to_plot)) {
    metric      <- metrics_to_plot[i]
    metric_name <- metric_names[i]
    plot_data   <- fold_data[!is.na(fold_data[[metric]]), ]
    if (nrow(plot_data) > 0) {
      p <- ggplot(plot_data, aes_string(x = "species", y = metric)) +
        geom_boxplot(fill = "lightblue", alpha = 0.7) +
        geom_jitter(width = 0.2, alpha = 0.5, size = 0.8) +
        labs(title = paste(metric_name, "Distribution Across Folds"), x = "Species", y = metric_name) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
              plot.title = element_text(hjust = 0.5, size = 14),
              legend.position = "none")
      ggsave(file.path(output_dir, paste0(metric, "_boxplot.png")), p, width = 12, height = 8, dpi = 300)
      message("Saved boxplot: ", metric_name)
    } else {
      message("No valid data for ", metric_name)
    }
  }

  # Combined plot
  plot_data_long <- fold_data %>%
    tidyr::pivot_longer(cols = c("roc_auc", "auprc", "auprc_norm", "tss"),
                        names_to = "metric", values_to = "value") %>%
    dplyr::filter(!is.na(value))

  if (nrow(plot_data_long) > 0) {
    plot_data_long$metric <- factor(plot_data_long$metric,
                                    levels = c("roc_auc", "auprc", "auprc_norm", "tss"),
                                    labels = c("ROC-AUC", "AUPRC", "nAUPRC", "TSS"))
    p_combined <- ggplot(plot_data_long, aes(x = metric, y = value)) +
      geom_boxplot(fill = "lightgreen", alpha = 0.7) +
      geom_jitter(width = 0.2, alpha = 0.5, size = 0.8) +
      labs(title = "Model Performance Metrics Across All Species and Folds", x = "Metric", y = "Value") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 14))
    ggsave(file.path(output_dir, "combined_metrics_boxplot_250.png"), p_combined, width = 10, height = 8, dpi = 300)
    message("Saved combined boxplot")
  }
}

# ---- Main function --------------------------------------------------------
run_sdm_loop <- function(gadm_sf, occurrences_df, target_species,
                         resolution = 2.5, k_folds = 5, seed = 2025, k_blocks = 5) {

  coord <- pick_coord_cols(occurrences_df)
  varimp_rows     <- list()
  processed_count <- 0
  skipped_count   <- 0
  set.seed(seed)

  dir.create(out_dir,  recursive = TRUE, showWarnings = FALSE)
  dir.create(fold_dir, recursive = TRUE, showWarnings = FALSE)

  for (i in seq_along(target_species)) {
    sp <- target_species[i]

    if (is.na(sp) || sp == "" || length(sp) == 0) {
      warning("Skipping invalid species name at index ", i)
      skipped_count <- skipped_count + 1
      next
    }

    message("Modeling [", i, "/", length(target_species), "]: ", sp)

    # 1. Filter occurrences
    sp_occ <- occurrences_df %>%
      dplyr::filter(species == sp) %>%
      dplyr::distinct(.data[[coord[1]]], .data[[coord[2]]], .keep_all = TRUE)

    if (nrow(sp_occ) > 200) sp_occ <- dplyr::slice_sample(sp_occ, n = 200)
    message("  Occurrences: ", nrow(sp_occ))

    if (nrow(sp_occ) < 10) {
      warning("Skipping ", sp, ": insufficient data (<10 points)")
      skipped_count <- skipped_count + 1
      next
    }

    # 2. Spatial filter to study area
    sp_sf <- st_as_sf(sp_occ, coords = coord, crs = 4326)
    if (inherits(gadm_sf, "SpatVector")) gadm_sf <- st_as_sf(gadm_sf)
    gadm_sf <- st_transform(gadm_sf, crs = st_crs(sp_sf))
    sp_sf   <- sp_sf[st_within(sp_sf, st_union(gadm_sf), sparse = FALSE), , drop = FALSE]
    if (nrow(sp_sf) < 10) {
      warning("Skipping ", sp, ": <10 occurrences within study area after spatial filter")
      skipped_count <- skipped_count + 1
      next
    }

    study_area <- gadm_sf

    # 3. WorldClim climatic predictors
    env        <- geodata::worldclim_global(var = "bio", res = resolution, path = tempdir())
    study_env  <- terra::project(vect(study_area), crs(env))
    env_crop   <- crop(env, study_env)
    env_masked <- mask(env_crop, study_env)

    # 4. NDVI (250 m)
    ndvi_stack   <- c(m_ndvi, sd_ndvi)
    study_ndvi   <- terra::project(vect(study_area), crs(ndvi_stack))
    ndvi_crop    <- crop(ndvi_stack, study_ndvi)
    ndvi_masked  <- mask(ndvi_crop, study_ndvi)

    # 5. Terrain (250 m, Paraguay)
    terrain_stack    <- terrain_py
    study_terrain    <- terra::project(vect(study_area), crs(terrain_stack))
    terrain_crop     <- crop(terrain_stack, study_terrain)
    terrain_masked   <- mask(terrain_crop, study_terrain)

    # 6. GPW population density (250 m)
    gpw_stack    <- gpw_py
    study_gpw    <- terra::project(vect(study_area), crs(gpw_stack))
    gpw_crop     <- crop(gpw_stack, study_gpw)
    gpw_masked   <- mask(gpw_crop, study_gpw)

    # Align all layers to the 250 m NDVI grid (bilinear for continuous variables)
    terrain_resampled <- terra::resample(terrain_masked, ndvi_masked[[1]], method = "bilinear")
    gpw_resampled     <- terra::resample(gpw_masked,     ndvi_masked[[1]], method = "bilinear")
    env_resampled     <- terra::resample(env_masked,     ndvi_masked[[1]], method = "bilinear")

    # Combine all predictors: climate + NDVI + terrain + GPW
    env_with_ndvi <- c(env_resampled, ndvi_masked, terrain_resampled, gpw_resampled)

    if (nlyr(env_with_ndvi) == 0) {
      warning("Skipping ", sp, ": no valid environmental layers")
      skipped_count <- skipped_count + 1
      next
    }

    # Remove all-NA layers
    valid_layers <- sapply(seq_len(nlyr(env_with_ndvi)), function(k) {
      !all(is.na(terra::values(env_with_ndvi[[k]])))
    })
    if (sum(valid_layers) == 0) {
      warning("Skipping ", sp, ": all environmental layers contain only NA values")
      skipped_count <- skipped_count + 1
      next
    }
    if (sum(valid_layers) < nlyr(env_with_ndvi)) {
      message("  Removing ", sum(!valid_layers), " all-NA layers for ", sp)
      env_with_ndvi <- env_with_ndvi[[which(valid_layers)]]
    }

    # 7. Occurrence coordinates
    occ_coords <- st_coordinates(sp_sf)
    occ_points <- as.data.frame(occ_coords)
    colnames(occ_points) <- c("lon", "lat")

    # 8. Background points (exclude 5 km buffer around occurrences)
    sp_proj       <- st_transform(sp_sf, 3857)
    buffer_5km    <- st_buffer(sp_proj, dist = 5000)
    buffer_union  <- st_union(buffer_5km)
    buffer_union_ll <- st_transform(buffer_union, 4326)
    buffer_vect   <- terra::makeValid(vect(buffer_union_ll))
    study_vect    <- terra::makeValid(vect(study_area))
    pd_aligned    <- tryCatch(terra::project(vect(gadm_sf), crs(buffer_vect)), error = function(e) vect(gadm_sf))
    pd_aligned    <- terra::makeValid(pd_aligned)
    pd_in_study   <- terra::intersect(pd_aligned, study_vect)
    bg_area       <- terra::erase(pd_in_study, buffer_vect)
    if (nrow(bg_area) == 0) {
      warning("Background area empty after erasing buffer; using full study area")
      bg_area <- pd_in_study
    }
    set.seed(seed)
    bg_points_pool  <- spatSample(bg_area, size = 1000, method = "random")
    bg_pool_coords  <- terra::crds(bg_points_pool)
    colnames(bg_pool_coords) <- c("lon", "lat")
    bg_pool_sf <- st_as_sf(as.data.frame(bg_pool_coords), coords = c("lon", "lat"), crs = 4326)
    occ_sf     <- st_as_sf(as.data.frame(occ_coords), coords = c("X", "Y"), crs = 4326)

    # 9. Extract environmental values
    env_vals_occ     <- terra::extract(env_with_ndvi, occ_sf)
    env_vals_bg_pool <- terra::extract(env_with_ndvi, bg_pool_sf)

    if (nrow(env_vals_occ) < 10) {
      warning("Skipping ", sp, ": insufficient valid environmental records for presences (<10)")
      skipped_count <- skipped_count + 1
      next
    }

    # Main 2000-point BG subset for model fitting
    set.seed(seed)
    bg_sub_idx_main <- seq_len(nrow(env_vals_bg_pool))
    if (length(bg_sub_idx_main) > 2000) bg_sub_idx_main <- sample(bg_sub_idx_main, size = 2000)
    env_vals_bg <- env_vals_bg_pool[bg_sub_idx_main, , drop = FALSE]
    bg_coords   <- bg_pool_coords[bg_sub_idx_main, , drop = FALSE]

    # Build unified dataset
    env_all      <- rbind(env_vals_occ, env_vals_bg)
    presence_vec <- c(rep(1L, nrow(env_vals_occ)), rep(0L, nrow(env_vals_bg)))
    coords_df    <- rbind(occ_points[, c("lon", "lat")], as.data.frame(bg_coords)[, c("lon", "lat")])
    if (ncol(env_all) > 0 && names(env_all)[1] %in% c("ID", "ID_1")) env_all <- env_all[, -1, drop = FALSE]
    keep         <- complete.cases(env_all)
    env_all      <- env_all[keep, , drop = FALSE]
    coords_df    <- coords_df[keep, , drop = FALSE]
    presence_vec <- as.integer(presence_vec[keep])
    sdm_data     <- cbind(env_all, presence = presence_vec)

    # Clean predictors: remove non-finite and zero-variance columns
    predictor_names <- setdiff(names(sdm_data), "presence")
    finite_ok <- vapply(predictor_names, function(nm) all(is.finite(sdm_data[[nm]])), logical(1))
    if (!all(finite_ok)) {
      sdm_data <- sdm_data[, c("presence", predictor_names[finite_ok]), drop = FALSE]
      predictor_names <- setdiff(names(sdm_data), "presence")
    }
    if (length(predictor_names) == 0) {
      warning("Skipping ", sp, ": no valid predictors after cleaning")
      skipped_count <- skipped_count + 1
      next
    }
    zero_var <- vapply(predictor_names, function(nm) stats::sd(sdm_data[[nm]]) == 0, logical(1))
    if (any(zero_var)) {
      sdm_data <- sdm_data[, c("presence", predictor_names[!zero_var]), drop = FALSE]
      predictor_names <- setdiff(names(sdm_data), "presence")
    }

    n_pres <- sum(sdm_data$presence == 1L)
    n_abs  <- sum(sdm_data$presence == 0L)
    if (n_pres < 10) {
      warning("Skipping ", sp, ": insufficient presences after cleaning (", n_pres, ")")
      skipped_count <- skipped_count + 1
      next
    }
    # Cap background at 10x presences
    max_abs <- min(n_abs, 10L * n_pres)
    if (n_abs > max_abs) {
      set.seed(seed)
      pres_idx  <- which(sdm_data$presence == 1L)
      abs_idx   <- which(sdm_data$presence == 0L)
      keep_abs  <- sample(abs_idx, size = max_abs)
      keep_rows <- c(pres_idx, keep_abs)
      sdm_data  <- sdm_data[keep_rows, , drop = FALSE]
      coords_df <- coords_df[keep_rows, , drop = FALSE]
    }

    # Ensure all predictors are numeric
    for (col in predictor_names) {
      if (!is.numeric(sdm_data[[col]])) sdm_data[[col]] <- as.numeric(sdm_data[[col]])
    }

    # ---- Spatial blocks (on full presence + BG pool) ----------------------
    coords_all_pool <- rbind(occ_points[, c("lon", "lat")], as.data.frame(bg_pool_coords)[, c("lon", "lat")])
    blocks_all      <- make_spatial_blocks(coords_all_pool, k_blocks = k_blocks, seed = seed)
    blocks_all      <- as.integer(as.factor(blocks_all))

    # ---- Out-of-block cross-validation (spatial blocks) -------------------
    pred_all     <- numeric(0)
    labels_all   <- integer(0)
    fold_metrics <- list()
    feature_names <- setdiff(names(sdm_data), "presence")
    pres_block_ids <- unique(blocks_all[seq_len(nrow(occ_points))])
    fold_counter  <- 0

    for (b in sort(pres_block_ids)) {
      set.seed(seed + b)

      # Fold-specific BG subsample from pool
      bg_sub_idx_fold <- seq_len(nrow(env_vals_bg_pool))
      if (length(bg_sub_idx_fold) > 2000) bg_sub_idx_fold <- sample(bg_sub_idx_fold, size = 2000)

      # Build fold dataset
      env_fold      <- rbind(env_vals_occ, env_vals_bg_pool[bg_sub_idx_fold, , drop = FALSE])
      presence_fold <- c(rep(1L, nrow(env_vals_occ)), rep(0L, length(bg_sub_idx_fold)))
      coords_fold   <- rbind(occ_points[, c("lon", "lat")], as.data.frame(bg_pool_coords)[bg_sub_idx_fold, , drop = FALSE])
      if (ncol(env_fold) > 0 && names(env_fold)[1] %in% c("ID", "ID_1")) env_fold <- env_fold[, -1, drop = FALSE]
      keep_f        <- complete.cases(env_fold)
      env_fold      <- env_fold[keep_f, , drop = FALSE]
      coords_fold   <- coords_fold[keep_f, , drop = FALSE]
      presence_fold <- as.integer(presence_fold[keep_f])
      sdm_fold      <- cbind(env_fold, presence = presence_fold)

      # Derive fold blocks from precomputed pool
      rows_all    <- c(seq_len(nrow(occ_points)), nrow(occ_points) + bg_sub_idx_fold)
      blocks_fold <- as.integer(as.factor(blocks_all[rows_all][keep_f]))
      n_f         <- nrow(sdm_fold)
      if (n_f < 20) next
      test_ids  <- which(blocks_fold == b)
      train_ids <- setdiff(seq_len(n_f), test_ids)
      if (length(test_ids) == 0 || length(train_ids) < 10) next
      if (length(unique(sdm_fold$presence[train_ids])) < 2 ||
          length(unique(sdm_fold$presence[test_ids]))  < 2) next

      mdl_b <- tryCatch({
        maxnet::maxnet(
          p    = sdm_fold$presence[train_ids],
          data = sdm_fold[train_ids, feature_names, drop = FALSE],
          f    = maxnet::maxnet.formula(p = sdm_fold$presence[train_ids],
                                        data = sdm_fold[train_ids, feature_names, drop = FALSE],
                                        classes = "default")
        )
      }, error = function(e) {
        message("  Error building model for ", sp, " fold ", b, ": ", e$message)
        NULL
      })
      if (is.null(mdl_b)) next

      X_test <- sdm_fold[test_ids, feature_names, drop = FALSE]
      if (nrow(X_test) == 0 || ncol(X_test) == 0) next
      if (!identical(colnames(X_test), feature_names)) next

      preds <- tryCatch({
        stats::predict(mdl_b, newdata = X_test, type = "cloglog")
      }, error = function(e) {
        message("  Error predicting for ", sp, " fold ", b, ": ", e$message)
        rep(NA_real_, nrow(X_test))
      })

      pred_all   <- c(pred_all,   preds)
      labels_all <- c(labels_all, sdm_fold$presence[test_ids])

      fold_counter  <- fold_counter + 1
      valid_fold    <- is.finite(preds)
      labels_fold   <- as.integer(sdm_fold$presence[test_ids][valid_fold])
      scores_fold   <- preds[valid_fold]

      if (length(scores_fold) > 0 && sum(labels_fold == 1L) > 0 && sum(labels_fold == 0L) > 0) {
        n_pres_fold  <- sum(labels_fold == 1L)
        n_bg_fold    <- sum(labels_fold == 0L)
        denom_fold   <- n_pres_fold + n_bg_fold
        baseline_fold <- if (is.finite(denom_fold) && denom_fold > 0) n_pres_fold / denom_fold else NA_real_

        roc_auc_fold <- tryCatch(
          dismo::evaluate(p = scores_fold[labels_fold == 1L], a = scores_fold[labels_fold == 0L])@auc,
          error = function(e) NA_real_
        )
        auprc_fold <- tryCatch(
          if (is.finite(baseline_fold) && n_pres_fold > 0 && n_bg_fold > 0)
            compute_auprc(as.integer(labels_fold), scores_fold) else NA_real_,
          error = function(e) NA_real_
        )
        tss_res_fold <- tryCatch(
          if (n_pres_fold > 0 && n_bg_fold > 0)
            compute_tss(scores_fold[labels_fold == 1L], scores_fold[labels_fold == 0L])
          else list(tss = NA_real_, threshold = NA_real_),
          error = function(e) list(tss = NA_real_, threshold = NA_real_)
        )
        tss_fold      <- tss_res_fold$tss
        auprc_norm_fold <- if (is.finite(auprc_fold) && is.finite(baseline_fold) && baseline_fold < 1)
          (auprc_fold - baseline_fold) / (1 - baseline_fold) else NA_real_

        fold_metrics[[fold_counter]] <- data.frame(
          species    = sp,
          fold       = fold_counter,
          roc_auc    = roc_auc_fold,
          auprc      = auprc_fold,
          auprc_norm = auprc_norm_fold,
          tss        = tss_fold,
          baseline   = baseline_fold,
          n_pres     = n_pres_fold,
          n_bg       = n_bg_fold,
          stringsAsFactors = FALSE
        )
      }
    }  # end spatial block loop

    # ---- Overall metrics on pooled out-of-block predictions ---------------
    valid_p  <- is.finite(pred_all)
    labels   <- as.integer(labels_all[valid_p])
    scores   <- pred_all[valid_p]
    n_pres_eval <- sum(labels == 1L, na.rm = TRUE)
    n_bg_eval   <- sum(labels == 0L, na.rm = TRUE)
    denom       <- n_pres_eval + n_bg_eval
    baseline    <- if (is.finite(denom) && denom > 0) n_pres_eval / denom else NA_real_

    roc_auc <- tryCatch(
      dismo::evaluate(p = scores[labels == 1L], a = scores[labels == 0L])@auc,
      error = function(e) NA_real_
    )
    auprc <- tryCatch(
      if (is.finite(baseline) && n_pres_eval > 0 && n_bg_eval > 0)
        compute_auprc(as.integer(labels), scores) else NA_real_,
      error = function(e) NA_real_
    )
    tss_res <- tryCatch(
      if (n_pres_eval > 0 && n_bg_eval > 0)
        compute_tss(scores[labels == 1L], scores[labels == 0L])
      else list(tss = NA_real_, threshold = NA_real_),
      error = function(e) list(tss = NA_real_, threshold = NA_real_)
    )
    tss     <- tss_res$tss
    tss_thr <- tss_res$threshold
    auprc_norm <- if (is.finite(auprc) && is.finite(baseline) && baseline < 1)
      (auprc - baseline) / (1 - baseline) else NA_real_

    # ---- Permutation-based variable importance (stratified k-fold CV) -----
    vars <- setdiff(names(sdm_data), "presence")
    if (length(vars) > 0) {
      folds_info  <- make_stratified_folds(sdm_data$presence, k = k_folds, seed = seed)
      folds       <- folds_info$folds
      k_eff       <- folds_info$k
      n_all       <- nrow(sdm_data)
      base_pred_k <- rep(NA_real_, n_all)
      perm_preds  <- lapply(vars, function(v) rep(NA_real_, n_all))
      names(perm_preds) <- vars

      for (f in seq_len(k_eff)) {
        test_ids_f  <- which(folds == f)
        if (length(test_ids_f) == 0) next
        train_ids_f <- setdiff(seq_len(n_all), test_ids_f)

        mdl_f <- tryCatch({
          maxnet::maxnet(
            p    = sdm_data$presence[train_ids_f],
            data = sdm_data[train_ids_f, vars, drop = FALSE],
            f    = maxnet::maxnet.formula(p    = sdm_data$presence[train_ids_f],
                                          data = sdm_data[train_ids_f, vars, drop = FALSE],
                                          classes = "default")
          )
        }, error = function(e) {
          message("  Error building CV model for ", sp, " fold ", f, ": ", e$message)
          NULL
        })
        if (is.null(mdl_f)) next

        X_test_f <- sdm_data[test_ids_f, vars, drop = FALSE]
        base_pred <- tryCatch(
          stats::predict(mdl_f, newdata = X_test_f, type = "cloglog"),
          error = function(e) rep(NA_real_, nrow(X_test_f))
        )
        base_pred_k[test_ids_f] <- base_pred

        for (v in vars) {
          Xp      <- X_test_f
          Xp[[v]] <- sample(Xp[[v]])
          perm_pred <- tryCatch(
            stats::predict(mdl_f, newdata = Xp, type = "cloglog"),
            error = function(e) rep(NA_real_, nrow(Xp))
          )
          perm_preds[[v]][test_ids_f] <- perm_pred
        }
      }

      base_auc <- tryCatch(
        dismo::evaluate(p = base_pred_k[sdm_data$presence == 1],
                        a = base_pred_k[sdm_data$presence == 0])@auc,
        error = function(e) NA_real_
      )

      if (!is.na(base_auc)) {
        drops <- lapply(vars, function(v) {
          auc_p <- tryCatch(
            dismo::evaluate(p = perm_preds[[v]][sdm_data$presence == 1],
                            a = perm_preds[[v]][sdm_data$presence == 0])@auc,
            error = function(e) NA_real_
          )
          data.frame(variable = v, auc_drop = base_auc - auc_p, stringsAsFactors = FALSE)
        })
        drops <- do.call(rbind, drops)
        drops$auc_drop[!is.finite(drops$auc_drop)] <- 0
        total              <- sum(pmax(drops$auc_drop, 0))
        drops$perc_contribution <- if (total > 0) 100 * pmax(drops$auc_drop, 0) / total else 0
        drops$species      <- sp
        drops$roc_auc      <- roc_auc
        drops$baseline     <- baseline
        drops$auprc        <- auprc
        drops$AUPRC_norm   <- auprc_norm
        drops$tss          <- tss
        drops$tss_thr      <- tss_thr
        drops <- drops[, c("species", "variable", "auc_drop", "perc_contribution",
                           "roc_auc", "baseline", "auprc", "AUPRC_norm", "tss", "tss_thr")]
        varimp_rows[[length(varimp_rows) + 1]] <- drops
        processed_count <- processed_count + 1
      } else {
        message("  CV baseline AUC could not be computed for ", sp)
        skipped_count <- skipped_count + 1
      }
    } else {
      message("  No predictors found for ", sp)
    }

    # ---- Save per-species fold metrics ------------------------------------
    if (length(fold_metrics) > 0) {
      fold_metrics_df   <- do.call(rbind, fold_metrics)
      fold_metrics_path <- file.path(fold_dir, paste0(gsub("[ /]", "_", sp), "_fold_metrics.csv"))
      dir.create(dirname(fold_metrics_path), showWarnings = FALSE, recursive = TRUE)
      utils::write.csv(fold_metrics_df, fold_metrics_path, row.names = FALSE)
      message("  Saved fold metrics for ", sp, " (", nrow(fold_metrics_df), " folds)")
    }

    # ---- Final model and raster prediction --------------------------------
    model <- tryCatch({
      maxnet::maxnet(
        p    = sdm_data$presence,
        data = sdm_data[, vars, drop = FALSE],
        f    = maxnet::maxnet.formula(p = sdm_data$presence, data = sdm_data[, vars, drop = FALSE], classes = "default")
      )
    }, error = function(e) {
      message("  Error building main model for ", sp, ": ", e$message)
      NULL
    })

    if (!is.null(model)) {
      env_for_pred         <- env_with_ndvi
      names(env_for_pred)  <- vars
      prediction <- tryCatch(
        terra::predict(env_for_pred, model, type = "cloglog", na.rm = TRUE),
        error = function(e) {
          message("  Error predicting raster for ", sp, ": ", e$message)
          NULL
        }
      )

      if (!is.null(prediction)) {
        out_path <- file.path(out_dir, paste0(gsub("[ /]", "_", sp), "_SDM.tif"))
        dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
        writeRaster(prediction, out_path, overwrite = TRUE)
        message("  Saved: ", out_path)

        # Two-panel plot: continuous + TSS-binary
        bg_points_main <- vect(as.data.frame(bg_coords), geom = c("lon", "lat"), crs = "EPSG:4326")
        par(mfrow = c(1, 2))
        plot(prediction, main = paste0(i, "_", gsub(" ", "_", sp), " - Continuous"))
        plot(vect(study_area), add = TRUE, border = "white", lwd = 1.2)
        plot(bg_points_main, add = TRUE, pch = 20, cex = 0.4, col = "white")
        plot(vect(occ_sf), add = TRUE, pch = 19, cex = 0.6, col = "red")

        if (is.finite(tss_thr)) {
          binary_map <- prediction >= tss_thr
          plot(binary_map, main = paste0(i, "_", gsub(" ", "_", sp), " - Binary (TSS)"))
          plot(vect(study_area), add = TRUE, border = "white", lwd = 1.2)
          plot(bg_points_main, add = TRUE, pch = 20, cex = 0.4, col = "white")
          plot(vect(occ_sf), add = TRUE, pch = 19, cex = 0.6, col = "red")
        } else {
          plot(prediction, main = paste0(i, "_", gsub(" ", "_", sp), " - Binary (TSS unavailable)"))
          plot(vect(study_area), add = TRUE, border = "white", lwd = 1.2)
          plot(bg_points_main, add = TRUE, pch = 20, cex = 0.4, col = "white")
          plot(vect(occ_sf), add = TRUE, pch = 19, cex = 0.6, col = "red")
        }
      }
    }
  }  # end species loop

  message("Loop completed: ", processed_count, " species processed, ", skipped_count, " skipped")

  # ---- Combine fold metrics from saved files --------------------------------
  result <- list()
  result$varimp <- if (length(varimp_rows) > 0) do.call(rbind, varimp_rows) else NULL

  if (dir.exists(fold_dir)) {
    fold_files <- list.files(fold_dir, pattern = "_fold_metrics.csv$", full.names = TRUE)
    if (length(fold_files) > 0) {
      fold_metrics_all <- lapply(fold_files, function(f) {
        tryCatch(read.csv(f, stringsAsFactors = FALSE),
                 error = function(e) { message("Error reading ", f, ": ", e$message); NULL })
      })
      fold_metrics_all <- fold_metrics_all[!sapply(fold_metrics_all, is.null)]
      result$fold_metrics <- if (length(fold_metrics_all) > 0) do.call(rbind, fold_metrics_all) else NULL
      message("Combined fold metrics from ", length(fold_files), " files")
    } else {
      result$fold_metrics <- NULL
    }
  } else {
    result$fold_metrics <- NULL
  }

  result
}

# ---- Run ------------------------------------------------------------------
target_species_test <- unique(dat_all_PY$species)
length(target_species_test)
#5678
message("Total species to process: ", length(target_species_test))

results_PY <- run_sdm_loop(
  gadm_sf         = gadm_PY,
  occurrences_df  = dat_all_PY,
  target_species  = target_species_test[1:5],   # adjust range as needed
  resolution      = 2.5,
  k_folds         = 5,
  seed            = 2025,
  k_blocks        = 5
)

# ---- Post-processing ------------------------------------------------------
varimp_all_PY   <- results_PY$varimp
fold_metrics_PY <- results_PY$fold_metrics

message("Total species processed: ", ifelse(is.null(varimp_all_PY), 0, length(unique(varimp_all_PY$species))))

# Variable importance + metrics wide table
if (!is.null(varimp_all_PY)) {
  wide_PY <- tidyr::pivot_wider(varimp_all_PY,
                                id_cols     = species,
                                names_from  = variable,
                                values_from = perc_contribution,
                                values_fill = 0)
  metrics_df <- varimp_all_PY %>%
    dplyr::group_by(species) %>%
    dplyr::summarise(
      roc_auc    = dplyr::first(roc_auc),
      baseline   = dplyr::first(baseline),
      auprc      = dplyr::first(auprc),
      AUPRC_norm = dplyr::first(AUPRC_norm),
      tss        = dplyr::first(tss),
      tss_thr    = dplyr::first(tss_thr)
    )
  wide_PY <- dplyr::left_join(wide_PY, metrics_df, by = "species")

  # Friendly column name aliases
  rename_map <- c(ndvi_mean   = "NDVI_m",
                  ndvi_sd     = "NDVI_sd",
                  gpw_density = "GPW_density",
                  elevation   = "Elevation",
                  slope       = "Slope",
                  aspect      = "Aspect",
                  hillshade   = "Hillshade")
  for (old in names(rename_map)) {
    if (old %in% names(wide_PY)) names(wide_PY)[names(wide_PY) == old] <- rename_map[[old]]
  }

  out_csv <- file.path(out_dir, "PY_plants_varimp_with_metrics_250.csv")
  dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(wide_PY, out_csv, row.names = FALSE)
  message("Saved variable importance + metrics: ", out_csv)
}

# Fold metrics: save combined CSV and produce boxplots
if (!is.null(fold_metrics_PY)) {
  message("Fold metrics available for ", length(unique(fold_metrics_PY$species)), " species")

  fold_metrics_out <- file.path(out_dir, "combined_fold_metrics_250.csv")
  utils::write.csv(fold_metrics_PY, fold_metrics_out, row.names = FALSE)
  message("Saved combined fold metrics to: ", fold_metrics_out)

  boxplot_dir <- file.path(out_dir, "boxplots")
  create_metric_boxplots(fold_metrics_PY, boxplot_dir)
} else {
  message("No fold metrics available for boxplots")
}
