############################################################
# Download GHCN Hourly Weather Data
# Mid-Atlantic / Northeastern US
############################################################
 
rm(list = ls())

library(worldmet)
library(data.table)
library(dplyr)
library(lubridate)

options(stringsAsFactors = FALSE)

############################################################
# Station metadata
############################################################

stations <- import_ghcn_stations()

setDT(stations)

# Study region
study_states <- c(
  "ME","NH","VT","NY","MA","CT","RI",
  "NJ","PA","DE","MD","VA","KY"
)

# Keep only hourly-capable station networks
candidate_stations <- stations[
  country == "US" &
    state %in% study_states &
    network %in% c("I", "W")
]

cat(
  "Candidate stations:",
  nrow(candidate_stations),
  "\n"
)

############################################################
# Function to download one station
############################################################

download_station <- function(station,
                             years = 1960:2025){
  
  station_data <- rbindlist(
    
    lapply(years, function(y){
      
      cat(
        "Station:",
        station,
        "Year:",
        y,
        "\n"
      )
      
      tryCatch({
        
        x <- import_ghcn_hourly(
          station = station,
          year = y
        )
        
        # Skip missing years
        if(is.null(x) || nrow(x) == 0){
          return(NULL)
        }
        
        setDT(x)
        
        # Fix occasional type mismatches
        if("hour" %in% names(x)){
          x[, hour := as.character(hour)]
        }
        
        x
        
      }, error = function(e){
        
        message(
          "Failed: ",
          station,
          " ",
          y
        )
        
        NULL
        
      })
      
    }),
    
    fill = TRUE
  )
  
  station_data
}

############################################################
# Download all stations
############################################################

all_weather <- rbindlist(
  
  lapply(candidate_stations$id, function(station){
    
    cat(
      "\n=================================\n"
    )
    
    cat(
      "Downloading station:",
      station,
      "\n"
    )
    
    download_station(
      station = station,
      years = 1940:2025
    )
    
  }),
  
  fill = TRUE
)

############################################################
# Keep only variables likely needed
############################################################

weather_keep <- all_weather[
  ,
  .(
    station_id,
    station_name,
    date,
    
    # Temperature / humidity
    air_temp,
    dew_point,
    rh,
    
    # Wind
    ws,
    wd,
    
    # Pressure
    atmos_pres,
    
    # Visibility
    visibility,
    
    # Precipitation
    precip,
    precip_03,
    precip_06,
    precip_12,
    precip_24
  )
]

############################################################
# Save files
############################################################

saveRDS(
  weather_keep,
  file = "northeast_hourly_weather_1960_2025.rds",
  compress = "xz"
)

saveRDS(
  candidate_stations,
  file = "northeast_station_metadata.rds"
)

fwrite(
  candidate_stations,
  file = "northeast_station_metadata.csv"
)

############################################################
# Quick quality check
############################################################

weather_keep[
  ,
  .(
    first_date = min(date, na.rm = TRUE),
    last_date = max(date, na.rm = TRUE),
    n_hours = .N
  ),
  by = station_id
]