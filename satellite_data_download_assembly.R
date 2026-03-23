# ============================================================================
# Satellite Data Download and Assembly
# Lake Erie Central Basin Ecosystem Metabolism Study
#
# Downloads and processes daily satellite-derived water quality data from the
# NOAA GLERL ERDDAP server (https://apps.glerl.noaa.gov/erddap/):
#   - Chromophoric dissolved organic matter (cDOM, m-1)
#   - Chlorophyll-a (Chla, ug L-1)
#   - Diffuse attenuation coefficient (Kd, m-1)
#   - Suspended minerals (SM, mg L-1)
#
# MODIS data used for 2002-2017; VIIRS data used for 2018-present.
# For each site, values are extracted as the mean of the target pixel and its
# eight surrounding pixels (nearest-neighbour mean) to reduce the impact of
# cloud cover and processing artifacts. Temporal gaps are filled by linear
# interpolation where gap length does not exceed the product-specific limit.
#
# Outputs (saved as .rds):
#   lake_erie_cdom_satellite.rds     (gap fill: <= 3 days)
#   lake_erie_chla_satellite.rds     (gap fill: <= 5 days)
#   lake_erie_kd_satellite.rds       (gap fill: <= 10 days)
#   lake_erie_sm_satellite.rds       (gap fill: <= 3 days)
# ============================================================================

library(ncdf4)       # Reading NetCDF files
library(terra)       # Spatial raster operations
library(sf)          # Spatial vector operations
library(raster)      # Legacy raster (for extract/adjacent functions)
library(dplyr)
library(lubridate)
library(tidyr)
library(RCurl)       # File download

mean_na <- function(x) mean(x, na.rm = TRUE)

# ============================================================================
# PART 1: Downloading raw NetCDF files from NOAA GLERL ERDDAP
# ============================================================================
# URL lists are downloaded from the ERDDAP interface as CSV files.
# The function converts the image (.png) URLs to science-quality NetCDF URLs.
# Run once for each product, pointing to the appropriate URL list CSV and
# destination directory.

curl_noaa <- function(x, dest_dir) {
  x <- gsub("_IMG", "_SQ",  x, fixed = TRUE)
  x <- gsub(".png",  ".nc", x, fixed = TRUE)
  x <- gsub(" ",     "",    x, fixed = TRUE)
  fname <- basename(x)
  download.file(x, paste0(dest_dir, fname))
}

# --- Suspended minerals (example; repeat for cDOM, Chla, Kd) ---------------
url_list <- read.csv(
  "/Users/jdh/lakeerie/data/LE_SM_VIIRS_IMG_0a08_57d9_4e0f_urls.csv",
  col.names = "url"
)
for (i in url_list$url) {
  curl_noaa(i, dest_dir = "/Users/jdh/lakeerie/data/lake_erie_sm/")
}

# ============================================================================
# PART 2: Shared Helper Functions
# ============================================================================

# Read a single NetCDF water-quality file and return a correctly oriented
# SpatRaster.  The function handles the Band1 variable structure used by
# NOAA GLERL CPA products, replaces fill values, and applies a sanity check
# on the data range.
read_netcdf_wq <- function(file_path) {
  nc          <- nc_open(file_path)
  data_3d     <- ncvar_get(nc, "Band1")
  lat         <- ncvar_get(nc, "lat")
  lon         <- ncvar_get(nc, "lon")
  fill_value  <- ncatt_get(nc, "Band1", "_FillValue")$value
  color_min   <- ncatt_get(nc, "Band1", "colorBarMinimum")
  color_max   <- ncatt_get(nc, "Band1", "colorBarMaximum")
  nc_close(nc)

  # Extract first time slice if array is 3-D [time, lat, lon]
  wq <- if (length(dim(data_3d)) == 3) data_3d[1, , ] else data_3d

  # Replace fill values with NA and apply reasonable data bounds
  wq[wq == fill_value] <- NA
  if (!is.null(color_min$value) && !is.null(color_max$value)) {
    reasonable_max <- color_max$value * 5   # allow up to 5× colorBarMaximum
    wq[wq < 0 | wq > reasonable_max] <- NA
  }

  # Transpose and create SpatRaster; flip vertically to correct N-S orientation
  wq <- t(wq)
  r  <- rast(wq,
             extent = c(min(lon), max(lon), min(lat), max(lat)),
             crs    = "EPSG:4326")
  r  <- flip(r, direction = "vertical")
  names(r) <- "wq"
  return(r)
}

# Extract the acquisition date from a NetCDF file's global attributes or
# the time variable (seconds since 1970-01-01).
extract_time_coverage_start <- function(file_path) {
  nc         <- nc_open(file_path)
  time_start <- ncatt_get(nc, 0, "time_coverage_start")

  if (!time_start$hasatt) {
    time_var   <- ncvar_get(nc, "time")
    time_units <- ncatt_get(nc, "time", "units")$value
    if (grepl("seconds since 1970", time_units)) {
      t0 <- as.POSIXct(time_var[1], origin = "1970-01-01", tz = "UTC")
      time_start_formatted <- format(t0, "%Y-%m-%dT%H:%M:%SZ")
    } else {
      time_start_formatted <- "Unknown"
    }
  } else {
    time_start_formatted <- time_start$value
  }
  nc_close(nc)
  return(time_start_formatted)
}

# Extract Band1 raster values at each CHRP site using a nearest-neighbour
# mean (target pixel + 8 adjacent pixels).
extract_nn_mean <- function(r_raster, site_df, lat_col = "Latitude",
                            lon_col = "Longitude", out_col = "nn_mean") {
  site_df[[out_col]] <- NA
  for (j in seq_len(nrow(site_df))) {
    cells <- cellFromXY(r_raster, site_df[j, c(lon_col, lat_col)])
    adj   <- adjacent(r_raster, cells[1], 8, include = TRUE)
    site_df[[out_col]][j] <- mean_na(extract(r_raster, adj[, 2]))
  }
  return(site_df)
}

# Linear interpolation of a site-level time series across gaps of
# <= max_gap_days.  Returns the original data frame with an additional
# column named "<value_col>_interp".
interpolate_with_constraint <- function(data, value_col,
                                        date_col = "date",
                                        max_gap_days = 10) {
  data[[date_col]] <- as.Date(data[[date_col]])
  result_list      <- list()

  for (site_name in unique(data$site)) {
    site_data <- data %>% filter(site == site_name)
    merged    <- merge(
      data.frame(date      = site_data[[date_col]],
                 site      = site_data$site,
                 Latitude  = unique(site_data$Latitude)[1],
                 Longitude = unique(site_data$Longitude)[1]),
      site_data,
      by = c("site", date_col), all.x = TRUE
    )

    non_na <- which(!is.na(merged[[value_col]]))
    merged[[paste0(value_col, "_interp")]] <- merged[[value_col]]

    if (length(non_na) > 1) {
      gaps <- diff(merged[[date_col]][non_na])
      for (i in seq_len(length(non_na) - 1)) {
        gap_days <- as.numeric(gaps[i])
        if (gap_days <= max_gap_days && gap_days > 1) {
          idx <- non_na[i]:non_na[i + 1]
          merged[[paste0(value_col, "_interp")]][idx] <- approx(
            x    = c(1, length(idx)),
            y    = merged[[value_col]][c(non_na[i], non_na[i + 1])],
            xout = seq_along(idx)
          )$y
        }
      }
    }
    result_list[[site_name]] <- merged
  }
  bind_rows(result_list)
}

# Load CHRP mooring coordinates (derived from Compiled_metabolism.csv)
chrp_sites <- readRDS("/Users/jdh/lakeerie/data/chrp_sites.rds")

# ============================================================================
# PART 3: Process cDOM (coloured dissolved organic matter absorption, m-1)
# Instrument: MODIS (2002-2017), VIIRS (2018-present)
# Temporal interpolation: <= 3 days
# ============================================================================
gisfiles_cdom <- list.files("/Users/jdh/lakeerie/data/lake_erie_cdom",
                             pattern = "*.nc")
cdom_raw <- data.frame()

for (i in gisfiles_cdom) {
  file_path <- paste0("/Users/jdh/lakeerie/data/lake_erie_cdom/", i)
  r         <- raster(file_path, varname = "Band1")
  proj4string(r) <- CRS("+init=EPSG:4326")

  chrp_sites$cdom    <- extract(r, chrp_sites[, c("Longitude", "Latitude")])
  chrp_sites         <- extract_nn_mean(r, chrp_sites, out_col = "nn_mean")
  chrp_sites$gis_dates <- extract_time_coverage_start(file_path)
  cdom_raw           <- bind_rows(cdom_raw, chrp_sites)
  chrp_sites$nn_mean <- NA
  chrp_sites$cdom    <- NA
}

cdom_raw$dtp  <- as.POSIXct(cdom_raw$gis_dates, format = "%Y-%m-%dT%H:%M:%SZ")
cdom_raw$date <- as.Date(cdom_raw$dtp)

cdom_interp <- interpolate_with_constraint(cdom_raw, "nn_mean", max_gap_days = 3)
saveRDS(cdom_interp, "/Users/jdh/lakeerie/data/lake_erie_cdom_satellite.rds")

# ============================================================================
# PART 4: Process Chlorophyll-a (Chla, ug L-1)
# Instrument: MODIS (2002-2017), VIIRS (2018-present)
# Temporal interpolation: <= 5 days
# ============================================================================
gisfiles_chla <- list.files("/Users/jdh/lakeerie/data/lake_erie_chla",
                              pattern = "*.nc")
chla_raw <- data.frame()

for (i in gisfiles_chla) {
  file_path    <- paste0("/Users/jdh/lakeerie/data/lake_erie_chla/", i)
  chlor_raster <- read_netcdf_wq(file_path)
  gis_date     <- extract_time_coverage_start(file_path)

  # Extract pixel value and 8-neighbour mean at each CHRP site
  chrp_sites$chlor_a <- extract(
    raster(chlor_raster), chrp_sites[, c("Longitude", "Latitude")]
  )
  chrp_sites <- extract_nn_mean(raster(chlor_raster), chrp_sites,
                                out_col = "chlor_a_nn")

  results      <- bind_cols(chrp_sites,
                            data.frame(gis_dates = gis_date))
  chla_raw     <- bind_rows(chla_raw, results)
  chrp_sites$chlor_a    <- NA
  chrp_sites$chlor_a_nn <- NA
}

chla_raw$dtp  <- as.POSIXct(chla_raw$gis_dates, format = "%Y-%m-%dT%H:%M:%SZ")
chla_raw$date <- as.Date(chla_raw$dtp)

chla_interp <- interpolate_with_constraint(chla_raw, "chlor_a_nn",
                                            max_gap_days = 5)
chla_interp$month <- month(chla_interp$date)
chla_interp$year  <- year(chla_interp$date)
saveRDS(chla_interp, "/Users/jdh/lakeerie/data/lake_erie_chla_satellite.rds")

# ============================================================================
# PART 5: Process Kd (diffuse attenuation coefficient of PAR, m-1)
# Instrument: MODIS (2002-2017), VIIRS (2018-present)
# Temporal interpolation: <= 10 days
# ============================================================================
gisfiles_kd <- list.files("/Users/jdh/lakeerie/data/lake_erie_kd",
                           pattern = "*.nc")
kd_raw <- data.frame()

for (i in gisfiles_kd) {
  file_path <- paste0("/Users/jdh/lakeerie/data/lake_erie_kd/", i)
  r         <- raster(file_path, varname = "Band1")
  proj4string(r) <- CRS("+init=EPSG:4326")

  chrp_sites$kd        <- extract(r, chrp_sites[, c("Longitude", "Latitude")])
  chrp_sites           <- extract_nn_mean(r, chrp_sites, out_col = "kd_nn_mean")
  chrp_sites$gis_dates <- extract_time_coverage_start(file_path)
  kd_raw               <- bind_rows(kd_raw, chrp_sites)
  chrp_sites$kd_nn_mean <- NA
  chrp_sites$kd         <- NA
}

kd_raw$dtp  <- as.POSIXct(kd_raw$gis_dates, format = "%Y-%m-%dT%H:%M:%SZ")
kd_raw$date <- as.Date(kd_raw$dtp)

kd_interp <- interpolate_with_constraint(kd_raw, "kd_nn_mean",
                                          max_gap_days = 10)
saveRDS(kd_interp, "/Users/jdh/lakeerie/data/lake_erie_kd_satellite.rds")

# ============================================================================
# PART 6: Process Suspended Minerals (SM, mg L-1)
# Instrument: MODIS (2002-2017), VIIRS (2018-present)
# Temporal interpolation: <= 3 days
# ============================================================================
gisfiles_sm <- list.files("/Users/jdh/lakeerie/data/lake_erie_sm",
                           pattern = "*.nc")
sm_raw <- data.frame()

for (i in gisfiles_sm) {
  file_path <- paste0("/Users/jdh/lakeerie/data/lake_erie_sm/", i)
  r         <- raster(file_path, varname = "Band1")
  proj4string(r) <- CRS("+init=EPSG:4326")

  chrp_sites$sm        <- extract(r, chrp_sites[, c("Longitude", "Latitude")])
  chrp_sites           <- extract_nn_mean(r, chrp_sites, out_col = "sm_nn_mean")
  chrp_sites$gis_dates <- extract_time_coverage_start(file_path)
  sm_raw               <- bind_rows(sm_raw, chrp_sites)
  chrp_sites$sm_nn_mean <- NA
  chrp_sites$sm         <- NA
}

sm_raw$dtp  <- as.POSIXct(sm_raw$gis_dates, format = "%Y-%m-%dT%H:%M:%SZ")
sm_raw$date <- as.Date(sm_raw$dtp)

sm_interp <- interpolate_with_constraint(sm_raw, "sm_nn_mean",
                                          max_gap_days = 3)
saveRDS(sm_interp, "/Users/jdh/lakeerie/data/lake_erie_sm_satellite.rds")
