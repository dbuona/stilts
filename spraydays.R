####calculate changes in number of "spray days" due to climate
###house keeping
rm(list=ls()) #removes everything in R environment
options(stringsAsFactors = FALSE) #stops R from automatically treating string as factors
options(mc.cores = parallel::detectCores()) # for Bayesian models, paralleled chain sampling
graphics.off() # removes anything in graphic environemnt
dev.off()
library(worldmet)
library(data.table)
library(dplyr)

stations <- getMeta()

head(stations)
names(stations)

library(data.table)

# Convert to data.table
stations <- as.data.table(stations)

# Look at date formats
head(stations[, .(station, ctry, begin, end, code)])

# Keep only USA stations
us_stations <- stations[
  ctry == "US"
]

nrow(us_stations)

# Candidate stations with long records
candidate_stations <- us_stations[
  begin <= as.Date("1960-01-01") &
    end >= as.Date("2020-01-01")
]

sort(table(candidate_stations$st))
candidate_stations<-filter(candidate_stations,!st %in% c("AK","HI"))

if(FALSE){
test <- importNOAA(
  code = candidate_stations$code[1],
  year = 1945:2025
)

summary(test$air_temp)
}

# Function to test hourly temperature availability
check_station <- function(code){
  
  x <- tryCatch(
    importNOAA(
      code = code,
      year = 1980:2020
    ),
    error = function(e) NULL
  )
  
  # Skip stations that fail
  if(is.null(x)){
    return(NULL)
  }
  
  data.table(
    code = code,
    start = min(x$date, na.rm = TRUE),
    end = max(x$date, na.rm = TRUE),
    hours = nrow(x),
    temp_hours = sum(!is.na(x$air_temp)),
    temp_percent = 100 * sum(!is.na(x$air_temp)) / nrow(x)
  )
}

temperature_screen <- rbindlist(
  lapply(
    candidate_stations$code,
    check_station
  ),
  fill = TRUE
)

candidate_stations<-filter(candidate_stations,code!="999999-13752")

###no precip
check_precip <- function(code){
  
  x <- tryCatch(
    importNOAA(
      code = code,
      year = 1980:2020
    ),
    error = function(e) NULL
  )
  
  if(is.null(x))
    return(NULL)
  
  precip_vars <- grep(
    "precip",
    names(x),
    value = TRUE
  )
  
  if(length(precip_vars) == 0){
    return(
      data.table(
        code = code,
        precip_hours = 0
      )
    )
  }
  
  data.table(
    code = code,
    precip_hours =
      sum(
        rowSums(
          !is.na(x[, ..precip_vars])
        ) > 0
      ),
    total_hours = nrow(x)
  )
}
precip_screen <- rbindlist(
  lapply(
    temperature_screen$code,
    check_precip
  ),
  fill = TRUE
)

precip_screen[
  ,
  precip_percent :=
    precip_hours / total_hours * 100
]



fwrite(
  temperature_screen,
  "hourly_temperature_screen.csv"
)


precip_screen <- rbindlist(
  lapply(
    temperature_screen$code,
    check_precip
  ),
  fill = TRUE
)

precip_screen[
  ,
  precip_percent :=
    precip_hours / total_hours * 100
]



x <- importNOAA(
  code = "724040-13721",
  year = 2020
)

names(x)
##############################maybe start here
###extract from stations
library(worldmet)
library(data.table)

# Stations you decided to keep
stations_use<-candidate_stations$usaf

#<- c(
#  "724010-13740",
#  "724037-93728",
#  "724040-13721",
#  "724050-13743",
#  "724060-93721",
#  "745940-13705",
#  "745980-13702"
#)

out_dir <- "weather_station_data"

dir.create(out_dir, showWarnings = FALSE)

library(worldmet)
library(data.table)

all_weather <- rbindlist(
  lapply(stations_use, function(code) {
    
    cat("Downloading:", code, "\n")
    
    x <- tryCatch(
      importNOAA(
        code = code,
        year = 1945:2025
      ),
      error = function(e) {
        message("Failed: ", code)
        return(NULL)
      }
    )
    
    if (is.null(x)) return(NULL)
    
    setDT(x)
    
    x[, station_code := code]
    
    # Add missing precip columns if needed
    if (!"precip_6" %in% names(x)) {
      x[, precip_6 := NA_real_]
    }
    
    if (!"pwc" %in% names(x)) {
      x[, pwc := NA_character_]
    }
    
    # Binary precipitation indicator from precip_6 when available
    x[, precip_any := fifelse(
      !is.na(precip_6) & precip_6 > 0,
      1L,
      0L,
      na = NA_integer_
    )]
    
    x[, .(
      station_code,
      date,
      temp_c = air_temp,
      precip_6,
      precip_any,
      pwc,
      RH,
      ws,
      wd
    )]
  }),
  fill = TRUE
)

fwrite(
  all_weather,
  "hourly_weather_1945_2025.csv"
)

all_weather<-read.csv("~/md_va_hourly_weather_1945_2025.csv")

all_weather[
  ,
  .(
    temp_pct = mean(!is.na(temp_c))*100,
    precip_pct = mean(!is.na(precip_6))*100,
    pwc_pct = mean(!is.na(pwc))*100
  ),
  by = station_code
]


library(data.table)

setDT(all_weather)

all_weather_subset <- all_weather[
  # Month-day filter
  format(date, "%m-%d") >= "04-01" &
    format(date, "%m-%d") <= "10-31" &
    
    # Hour filter
    as.integer(format(date, "%H")) >= 7 &
    as.integer(format(date, "%H")) <= 17
]


range(all_weather_subset$date)

table(format(all_weather_subset$date, "%H"))

all_weather_subset[
  ,
  spray_ok := (
    temp_c >= 10 &
      temp_c <= 29 &
      RH >= 40 &
      ws >= 1.3 &       # ~3 mph
      ws <= 4.5         # ~10 mph
  )
]


all_weather_subset[
  ,
  period := fcase(
    year(date) < 1980, "Pre-1980",
    year(date) < 2000, "1980-1999",
    year(date) <= 2025, "2000-2025"
  )
]


all_weather_subset[
  ,
  .(
    total_hours = .N,
    spray_hours = sum(spray_ok, na.rm = TRUE),
    spray_pct = mean(spray_ok, na.rm = TRUE) * 100
  ),
  by = .(period)
]

library(data.table)
library(ggplot2)

# Month number and label
all_weather_subset[, month := as.integer(format(date, "%m"))]

# Sprayable hours by station-year-month
monthly <- all_weather_subset[
  ,
  .(
    spray_hours = sum(spray_ok, na.rm = TRUE)
  ),
  by = .(
    station_code,
    period,
    year = as.integer(format(date, "%Y")),
    month
  )
]

# Average across all station-years within each period
monthly_avg <- monthly[
  ,
  .(
    avg_spray_hours = mean(spray_hours, na.rm = TRUE)
  ),
  by = .(period, month)
]

monthly_avg[
  ,
  month_name := factor(
    month.abb[month],
    levels = month.abb[4:10]
  )
]

ggplot(
  monthly_avg,
  aes(
    x = month_name,
    y = avg_spray_hours,
    color = period,
    group = period
  )
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  labs(
    x = "Month",
    y = "Average Sprayable Hours",
    color = "Period"
  ) +
  theme_minimal()

library(data.table)
library(ggplot2)

# Monthly sprayable hours by station-year
monthly <- all_weather_subset[
  ,
  .(
    spray_hours = sum(spray_ok, na.rm = TRUE)
  ),
  by = .(
    station_code,
    period,
    year = as.integer(format(date, "%Y")),
    month = as.integer(format(date, "%m"))
  )
]

# Mean and SE across station-years
monthly_summary <- monthly[
  ,
  .(
    mean_hours = mean(spray_hours, na.rm = TRUE),
    se = sd(spray_hours, na.rm = TRUE) / sqrt(.N),
    n = .N
  ),
  by = .(period, month)
]

monthly_summary[
  ,
  month_name := factor(
    month.abb[month],
    levels = month.abb[4:10]
  )
]

ggplot(
  monthly_summary,
  aes(
    x = month_name,
    y = mean_hours,
    fill = period
  )
) +
  geom_col(
    position = position_dodge(width = 0.9)
  ) +
  geom_errorbar(
    aes(
      ymin = mean_hours - se,
      ymax = mean_hours + se
    ),
    position = position_dodge(width = 0.9),
    width = 0.2
  ) +
  labs(
    x = NULL,
    y = "Average Sprayable Hours per Month",
    fill = "Period"
  ) +
  theme_minimal(base_size = 14)
