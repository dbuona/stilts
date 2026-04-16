library(prism)
library(terra)
prism_set_dl_dir("~/prism_data")


years <- c(2020,2021)

get_prism_monthlys(
  type = "tmean",
  years = years,
  keepZip = FALSE
)


dirs <- prism_archive_ls()

# Keep monthly tmean datasets
#dirs <- dirs[grepl("tmean", dirs)]
dirs <- dirs[grepl("tmean", dirs)]
# Build full path to the .bil files

bil_files <- file.path(
  prism_get_dl_dir(),
  dirs,
  paste0(dirs, ".bil")
)

# Sanity check (important)
bil_files <- bil_files[file.exists(bil_files)]





# Load rasters
tmean <- rast(bil_files)
#ppt <- rast(bil_files)

# Add time index
dates <- as.Date(
  paste0(
    sub(".*_(\\d{4})(\\d{2})$", "\\1", dirs),
    "-",
    sub(".*_(\\d{4})(\\d{2})$", "\\2", dirs),
    "-15"
  )
)
time(tmean) <- dates
#time(ppt) <- dates

months <- as.integer(format(time(tmean), "%m"))
years  <- format(time(tmean), "%Y")

tmean_march_sept <- tmean[[months >= 2 & months <= 3]]

tmean_annual <- tapp(
  tmean_march_sept,
  years[months >= 2 & months <= 3],
  mean,
  na.rm = TRUE
)

df <- dplyr::select(dat,decimalLatitude,decimalLongitude, year) ## your coordinates

pts <- vect(df, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")
pts <- project(pts, crs(tmean_annual))

vals <- extract(tmean_annual, pts)


names(tmean_annual) <- gsub("^X", "", names(tmean_annual))
years_available <- as.integer(names(tmean_annual))

df_out <- df |>
  mutate(
    col_id = match(year, years_available),
    tmean = vals[cbind(seq_len(nrow(vals)), col_id)]
  ) |>
  dplyr::select(decimalLongitude, decimalLatitude, year, tmean)


