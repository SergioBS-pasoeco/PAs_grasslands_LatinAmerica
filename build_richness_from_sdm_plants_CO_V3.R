
suppressPackageStartupMessages({
  library(terra)
})

# Directory with SDM prediction rasters (GeoTIFF)
PRED_DIR <- "F:/SDM_project/results/SDM_plants/CO_terrain_GPW_250m"

# Output directory for derived products
OUT_DIR <- file.path(PRED_DIR, "derived")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Optional CLI override: Rscript R/build_richness_from_sdm.R "path/to/preds"
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1 && nchar(args[[1]]) > 0) {
  PRED_DIR <- args[[1]]
  OUT_DIR  <- file.path(PRED_DIR, "derived")
  dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
}

message("Reading predictions from: ", PRED_DIR)

tif_files <- list.files(PRED_DIR, pattern = "\\.tif$", full.names = TRUE)
if (length(tif_files) == 0) {
  stop("No .tif files found in ", PRED_DIR)
}

# ---------------------------------------------------------------------------
# Read fold-level metrics and derive per-species mean TSS threshold
# ---------------------------------------------------------------------------
metrics_csv <- file.path(PRED_DIR, "combined_fold_metrics_250.csv")
if (!file.exists(metrics_csv)) {
  stop("Metrics file not found: ", metrics_csv)
}

metrics <- utils::read.csv(metrics_csv, stringsAsFactors = FALSE)

if (!"tss" %in% names(metrics)) {
  stop("Column 'tss' not found in ", metrics_csv,
       ". Available columns: ", paste(names(metrics), collapse = ", "))
}

if (!"species" %in% names(metrics)) {
  stop("Column 'species' not found in ", metrics_csv,
       ". Available columns: ", paste(names(metrics), collapse = ", "))
}

# Compute mean TSS per species
mean_tss <- aggregate(tss ~ species, data = metrics, FUN = mean, na.rm = TRUE)
names(mean_tss)[names(mean_tss) == "tss"] <- "threshold"

message("Mean TSS computed for ", nrow(mean_tss), " species.")

# ---------------------------------------------------------------------------
# Build the thresholds table pre-populated from mean TSS
# ---------------------------------------------------------------------------

# Derive species names from filenames:
# e.g. "Zanthoxylum_fagara_SDM.tif" -> "Zanthoxylum fagara"
sp_names_raw <- tools::file_path_sans_ext(basename(tif_files))
sp_names     <- gsub("_", " ",                        # underscores -> spaces
                  trimws(
                    gsub("(?i)\\s*SDM\\s*$", "",      # remove trailing "SDM"
                      gsub("_SDM$", "", sp_names_raw, ignore.case = TRUE),
                    perl = TRUE)
                  )
                )

thresholds <- data.frame(
  species   = sp_names,
  threshold = NA_real_,
  stringsAsFactors = FALSE
)

# Match species names and assign mean TSS values
matched <- match(thresholds$species, mean_tss$species)
thresholds$threshold <- mean_tss$threshold[matched]

# Warn about any species with no TSS entry
missing_sp <- thresholds$species[is.na(thresholds$threshold)]
if (length(missing_sp) > 0) {
  warning("No TSS entry found for ", length(missing_sp), " species; ",
          "they will be skipped during binarization:\n",
          paste(missing_sp, collapse = "\n"))
}

message("Thresholds assigned for ",
        sum(!is.na(thresholds$threshold)), " of ",
        nrow(thresholds), " species.")


thresholds$threshold[which(is.na(thresholds$threshold))] = mean(thresholds$threshold, na.rm = TRUE)


# ---------------------------------------------------------------------------
# Build richness raster by binarizing each .tif with its TSS threshold
# ---------------------------------------------------------------------------

template <- terra::rast(tif_files[1])
richness  <- template * 0

for (f in tif_files) {
  sp_raw <- tools::file_path_sans_ext(basename(f))
  sp     <- gsub("_", " ", trimws(gsub("_SDM$", "", sp_raw, ignore.case = TRUE)))
  thr <- thresholds$threshold[thresholds$species == sp]

  # Skip species with no valid threshold
  if (length(thr) == 0 || is.na(thr)) {
    message("Skipping ", sp, " (no TSS threshold available).")
    next
  }

  r <- terra::rast(f)
  if (!terra::compareGeom(r, template, stopOnError = FALSE)) {
    r <- terra::resample(r, template, method = "bilinear")
  }

  bin      <- terra::ifel(r >= thr, 1, 0)
  bin      <- terra::ifel(is.na(bin), 0, bin)
  richness <- richness + bin

  message(sprintf("  %-50s  thr = %.4f", sp, thr))
}

# ---------------------------------------------------------------------------
# Save outputs
# ---------------------------------------------------------------------------

thr_csv <- file.path(OUT_DIR, "species_thresholds_V2.csv")
utils::write.csv(thresholds, thr_csv, row.names = FALSE)
message("Saved thresholds table: ", thr_csv)

rich_tif <- file.path(OUT_DIR, "species_richness_plants_CO_V2.tif")
terra::writeRaster(richness, rich_tif, overwrite = TRUE)
message("Saved richness raster: ", rich_tif)

# Plot with a nice palette
pal <- tryCatch({
  if (requireNamespace("viridis", quietly = TRUE)) viridis::magma(100) else NULL
}, error = function(e) NULL)
if (is.null(pal)) pal <- terrain.colors(100)

plot(richness, col = pal, main = "Species richness (binary sum)")
