#!/usr/bin/env Rscript
# Extract the 19 WorldClim v2.1 bioclimatic variables (2.5 arc-minute) at the
# coordinates of each population in the sample sheet (R packages terra and
# geodata; requires internet access on first use).
#   Rscript scripts/teosinte/03_extract_climate.R [metadata_dir]
# Output: <metadata_dir>/climate_table_hilo.tsv (population, lat, lon, bio1..bio19)
suppressMessages({library(terra); library(geodata)})
args <- commandArgs(trailingOnly = TRUE)
metadir <- if (length(args) >= 1) args[1] else "metadata"

ss <- read.delim(file.path(metadir, "sample_sheet_hilo_mexicana.tsv"), stringsAsFactors = FALSE)
pops <- unique(ss[, c("population", "lat", "lon")])
cat("populations:", nrow(pops), "\n")

cache <- file.path(metadir, "worldclim")
dir.create(cache, showWarnings = FALSE, recursive = TRUE)
bio <- worldclim_global(var = "bio", res = 2.5, path = cache)
bn <- as.integer(sub(".*bio_?", "", names(bio)))
names(bio) <- paste0("bio", bn)
bio <- bio[[order(bn)]]

pts  <- vect(data.frame(lon = pops$lon, lat = pops$lat), geom = c("lon", "lat"), crs = "EPSG:4326")
vals <- terra::extract(bio, pts)
clim <- cbind(pops, vals[, -1])
out <- file.path(metadir, "climate_table_hilo.tsv")
write.table(clim, out, sep = "\t", quote = FALSE, row.names = FALSE)
cat("written", out, "\n")
cat(sprintf("bio1 (mean annual temperature): %.1f to %.1f C\n", min(clim$bio1), max(clim$bio1)))
cat(sprintf("bio12 (annual precipitation): %.0f to %.0f mm\n", min(clim$bio12), max(clim$bio12)))
