###make R polygon
rm(list=ls()) #removes everything in R environment
options(stringsAsFactors = FALSE) #stops R from automatically treating string as factors
options(mc.cores = parallel::detectCores()) # for Bayesian models, paralleled chain sampling
graphics.off() # removes anything in graphic environemnt
#options(device = "RStudioGD")
library(dplyr)
if(FALSE){
library(tigris)
library(rnaturalearth)
library(sf)
library(dplyr)

us_states <- states(cb = TRUE) %>%
  filter(NAME %in% c(
    "New York",
    "Maine",
    "Vermont",
    "Maryland",
    "Massachusetts",
    "Minnesota",
    "Indiana",
    "Michigan",
    "Wisconsin",
    "Virginia",
    "Pennsylvania",
    "North Carolina",
    "Rhode Island",
    "Missouri",
    "Ohio",
    "Tennessee",
    "West Virginia",
    "New Hampshire",
    "New Jersey",
    "Connecticut",
    "Kentucky",
    "Kansas"
  ))

canada <- ne_states(
  country = "Canada",
  returnclass = "sf"
) %>%
  filter(name_en %in% c(
    "Ontario",
    "Quebec"
  ))


library(sf)
library(dplyr)

# US states
us_states2 <- us_states %>%
  transmute(
    region_name = NAME,
    geometry = geometry
  )

# Canadian provinces
canada2 <- canada %>%
  transmute(
    region_name = name_en,
    geometry = geometry
  )

# Check current CRS
st_crs(us_states2)
st_crs(canada2)

# Reproject both to WGS84
us_states2 <- st_transform(us_states2, 4326)
canada2    <- st_transform(canada2, 4326)

# Now combine
regions <- rbind(us_states2, canada2)




regions$region_name <- ifelse(
  is.na(regions$NAME),
  regions$name_en,
  regions$NAME
)

st_write(
  regions,
  "~/Documents/git/stilts/data/spring_regions.geojson",
  delete_dsn = TRUE
)
}

###now make the data flile by state/province
setwd("~/Documents/git/stilts/data/greenup/")
if(FALSE){
library(terra)

files <- list.files(
  pattern = "\\.nc$",
  full.names = TRUE
)

length(files)
library(terra)

files <- list.files(
  pattern = "\\.nc$",
  full.names = TRUE
)

length(files)


x <- rast(files[1])

names(x)
library(ncdf4)

nc <- nc_open(files[1])

names(nc$var)
stats <- read.csv("MCD12Q2-061-Statistics.csv")

names(stats)

head(stats)
list.files(pattern = "csv")

length(unique(stats$aid))

unique(stats$aid)

table(substr(stats$Date, 1, 4))


library(dplyr)
library(stringr)

dat <- stats %>%
  mutate(
    cycle = str_extract(File.Name, "Greenup_[01]")
  )

table(dat$cycle)

library(dplyr)
library(lubridate)

dat2 <- dat %>%
  filter(cycle == "Greenup_0") %>%   # primary spring greenup
  mutate(
    greenup_date = as.Date("1970-01-01") + Mean,
    spring_year = year(greenup_date),
    doy = yday(greenup_date)
  )

range(dat2$doy, na.rm = TRUE)

summary(dat2$doy)
library(sf)

regions <- st_read("..//spring_regions.geojson")
regions
names(regions)

lookup <- data.frame(
  aid = sort(unique(dat$aid)),
  region = regions$region_name
)

lookup
regions$region_name
dat2 %>%
  filter(aid == "aid0001") %>%
  select(greenup_date, doy) %>%
  head()

dat2 %>%
  filter(aid == "aid0024") %>%
  select(greenup_date, doy) %>%
  head()
dat3 <- left_join(dat2, lookup, by = "aid")

print(dat3 %>%
  group_by(region) %>%
  summarise(
    mean_doy = mean(doy, na.rm = TRUE)
  ) %>%
  arrange(mean_doy),n=24)


dat %>%
  group_by(cycle) %>%
  summarise(
    mean_date = mean(Mean)
  )

as.Date("1970-01-01") + 15639
as.Date("1970-01-01") + 15791
stats[1, c("Date", "Mean")]
as.Date("1970-01-01") + 11513.03
head(stats$File.Name, 20)


dat %>%
  group_by(cycle) %>%
  summarise(
    min_mean = min(Mean),
    max_mean = max(Mean),
    mean_mean = mean(Mean)
  )

dat %>%
  filter(aid == "aid0016") %>%  # New York
  select(File.Name, Date, Mean) %>%
  head(20)
library(ncdf4)

nc <- nc_open(files[1])

print(nc)

library(ncdf4)

nc <- nc_open(files[1])

print(nc)

nc$dim
tmp <- dat %>%
  filter(aid == "aid0016",
         cycle == "Greenup_0") %>%
  mutate(
    greenup_date = as.Date("1970-01-01") + Mean,
    doy = lubridate::yday(greenup_date)
  )

tmp %>%
  select(Date, greenup_date, doy)


dat2 <- dat %>%
  filter(cycle == "Greenup_0") %>%
  mutate(
    analysis_year = as.integer(substr(Date, 1, 4)),
    greenup_date = as.Date("1970-01-01") + Mean,
    doy = lubridate::yday(greenup_date)
  )

range(dat2$doy, na.rm = TRUE)

dat3 <- left_join(dat2, lookup, by = "aid")

state_shift <- dat3 %>%
  group_by(region) %>%
  summarise(
    slope = coef(lm(doy ~ analysis_year))[2],
    total_shift_days =
      slope * (max(analysis_year) - min(analysis_year)),
    p_value =
      summary(lm(doy ~ analysis_year))$coefficients[2,4],
    mean_doy = mean(doy),
    .groups = "drop"
  ) %>%
  arrange(total_shift_days)

library(ggplot2)

ggplot(
  filter(dat3, region == "New York"),
  aes(analysis_year, doy)
) +
  geom_point() +
  geom_smooth(method = "lm")

dat3 %>%
  filter(region == "New York") %>%
  count(analysis_year)
dat3 %>%
  filter(region == "New York") %>%
  select(analysis_year, cycle, doy) %>%
  arrange(analysis_year)

dat3%>%filter(region == "New York") %>%
  count(analysis_year)

``
dat %>%
  filter(aid == "aid0016") %>%
  distinct(File.Name) %>%
  head(20)
``
table(grepl("Greenup_0_", dat$File.Name))
table(grepl("Greenup_1_", dat$File.Name))


dat %>%
  filter(aid == "aid0016") %>%
  mutate(
    greenup_date = as.Date("1970-01-01") + Mean
  ) %>%
  select(File.Name, greenup_date)



library(dplyr)
library(stringr)

dat2 <- stats %>%
  mutate(
    metric = case_when(
      str_detect(File.Name, "MidGreenup_0") ~ "MidGreenup_0",
      str_detect(File.Name, "MidGreenup_1") ~ "MidGreenup_1",
      str_detect(File.Name, "Greenup_0") ~ "Greenup_0",
      str_detect(File.Name, "Greenup_1") ~ "Greenup_1",
      TRUE ~ NA_character_
    )
  )

table(dat2$metric)

greenup <- dat2 %>%
  filter(metric == "Greenup_0") %>%
  mutate(
    analysis_year = as.integer(substr(Date, 1, 4)),
    greenup_date = as.Date("1970-01-01") + Mean,
    doy = lubridate::yday(greenup_date)
  ) %>%
  left_join(
    lookup,
    by = "aid"
  )

greenup %>%
  filter(region == "New York") %>%
  count(analysis_year)

library(ggplot2)

ggplot(
  filter(greenup, region == "New York"),
  aes(analysis_year, doy)
) +
  geom_point() +
  geom_smooth(method = "lm")

state_shift <- greenup %>%
  group_by(region) %>%
  summarise(
    slope = coef(
      lm(doy ~ analysis_year)
    )[2],
    
    total_shift_days =
      slope *
      (max(analysis_year) -
         min(analysis_year)),
    
    p_value =
      summary(
        lm(doy ~ analysis_year)
      )$coefficients[2,4],
    
    mean_doy =
      mean(doy),
    
    .groups = "drop"
  ) %>%
  arrange(total_shift_days)

state_shift <- greenup %>%
  group_by(region) %>%
  summarise(
    slope = coef(
      lm(doy ~ analysis_year)
    )[2],
    
    total_shift_days =
      slope *
      (max(analysis_year) -
         min(analysis_year)),
    
    p_value =
      summary(
        lm(doy ~ analysis_year)
      )$coefficients[2,4],
    
    mean_doy =
      mean(doy),
    
    .groups = "drop"
  ) %>%
  arrange(total_shift_days)
}
# ============================================================
# SPRING GREENUP SHIFT ANALYSIS (MODIS MCD12Q2 via AppEEARS)
#
# Goal:
# Calculate the average shift in spring onset (Greenup_0)
# from 2001-2024 for:
# - 22 US states
# - Ontario
# - Quebec
#
# Input files:
# 1. MCD12Q2-061-Statistics.csv
# 2. spring_regions.geojson
#
# Output:
# Spring_Shift_2001_2024.csv
# ============================================================

# Load packages
library(dplyr)
library(stringr)
library(lubridate)
library(sf)
library(ggplot2)

# ------------------------------------------------------------
# Read AppEEARS statistics file
# ------------------------------------------------------------

stats <- read.csv("MCD12Q2-061-Statistics.csv")

# ------------------------------------------------------------
# Read the GeoJSON used in the AppEEARS request
# ------------------------------------------------------------

regions <- st_read("..//spring_regions.geojson")

# ------------------------------------------------------------
# Create aid-to-region lookup table
#
# AppEEARS assigns:
# aid0001, aid0002, ...
# in the same order polygons were uploaded
# ------------------------------------------------------------

lookup <- data.frame(
  aid = sort(unique(stats$aid)),
  region = regions$region_name
)

# ------------------------------------------------------------
# Identify phenology metric
#
# The statistics file contains:
# Greenup_0
# Greenup_1
# MidGreenup_0
# MidGreenup_1
#
# Greenup_0 = primary spring greenup
# (this is the metric we want)
# ------------------------------------------------------------

dat <- stats %>%
  mutate(
    metric = case_when(
      str_detect(File.Name, "MidGreenup_0") ~ "MidGreenup_0",
      str_detect(File.Name, "MidGreenup_1") ~ "MidGreenup_1",
      str_detect(File.Name, "Greenup_0") ~ "Greenup_0",
      str_detect(File.Name, "Greenup_1") ~ "Greenup_1",
      TRUE ~ NA_character_
    )
  )

# Optional check
table(dat$metric)

# ------------------------------------------------------------
# Keep only spring greenup (Greenup_0)
#
# Mean values are stored as:
# days since 1970-01-01
#
# Convert to actual dates and then DOY
# ------------------------------------------------------------

greenup <- dat %>%
  filter(metric == "Greenup_0") %>%
  mutate(
    analysis_year = as.integer(substr(Date, 1, 4)),
    greenup_date = as.Date("1970-01-01") + Mean,
    doy = yday(greenup_date)
  ) %>%
  left_join(
    lookup,
    by = "aid"
  )

# ------------------------------------------------------------
# Sanity checks
# ------------------------------------------------------------

# Should be:
# 24 regions x 24 years = 576 rows
nrow(greenup)

# One observation per year for New York
greenup %>%
  filter(region == "New York") %>%
  count(analysis_year)

# Range of spring DOY values
range(greenup$doy, na.rm = TRUE)

# Average spring onset by region
greenup %>%
  group_by(region) %>%
  summarise(mean_doy = mean(doy)) %>%
  arrange(mean_doy)

# ------------------------------------------------------------
# Plot New York time series
# ------------------------------------------------------------

ggplot(
  filter(greenup, region == "New York"),
  aes(analysis_year, doy)
) +
  geom_point() +
  geom_smooth(
    method = "lm",
    se = TRUE
  ) +
  labs(
    x = "Year",
    y = "Spring Greenup (DOY)",
    title = "New York Spring Greenup (2001-2024)"
  )

# ------------------------------------------------------------
# Calculate spring onset trends
#
# Negative values = earlier spring
# Positive values = later spring
# ------------------------------------------------------------
# ------------------------------------------------------------
# Calculate spring onset summary statistics and trends
#
# Negative values = earlier spring
# Positive values = later spring
# ------------------------------------------------------------
state_shift <- greenup %>%
  group_by(region) %>%
  summarise(
    
    mean_doy =
      mean(doy, na.rm = TRUE),
    
    sd_doy =
      sd(doy, na.rm = TRUE),
    
    se_doy =
      sd(doy, na.rm = TRUE) /
      sqrt(sum(!is.na(doy))),
    
    slope =
      coef(
        lm(doy ~ analysis_year)
      )[2],
    
    slope_se =
      summary(
        lm(doy ~ analysis_year)
      )$coefficients[2, 2],
    
    total_shift_days =
      slope *
      (
        max(analysis_year) -
          min(analysis_year)
      ),
    
    total_shift_se =
      slope_se *
      (
        max(analysis_year) -
          min(analysis_year)
      ),
    
    .groups = "drop"
  ) %>%
  mutate(
    total_shift_lower =
      total_shift_days -
      1.96 * total_shift_se,
    
    total_shift_upper =
      total_shift_days +
      1.96 * total_shift_se
  ) %>%
  arrange(total_shift_days)

# ------------------------------------------------------------
# Save results
# ------------------------------------------------------------

write.csv(
  state_shift,
  "..//..//Input/Spring_Shift_2001_2024.csv",
  row.names = FALSE
)

