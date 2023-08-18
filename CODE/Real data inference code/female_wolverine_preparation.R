##### WOLVERINE DATA PREPARATION#####

# Loading libraries
library(raster)
library(tidyverse)
library(lubridate)

# Setting working directory
setwd("/import/ecoc9/data-jeltsch/arceguillen/female_wolverines") 

# Loading spatial covariates
landscape_lakes = raster("F6F11F12_lake edge distance.tif")
landscape_rivers = raster("F6F11F12_streams rivers distance.tif")
landscape_log_lakes = raster("F6F11F12_lake edge distance.tif")
landscape_log_rivers = raster("F6F11F12_streams rivers distance.tif")
landscape_terrain = raster("F6F11F12_TRI.tif")

# Loading wolverine tracks from Movebank
movebank_download <- read.csv("Wolverine - Glass - Arctic Alaska_part-gps.csv")

df1 <- movebank_download %>%
  filter(gps.dop < 4 | is.na(gps.dop)) %>%
  mutate(timestamp = ymd_hms(timestamp),
         season = as.factor(if_else(yday(timestamp) < 130, "Spring", "Summer"))) %>%
  filter(season == "Spring", individual.local.identifier %in% c("F6", "F11", "F12")) %>%
  select(x=location.long, y=location.lat, time_gmt = timestamp, id=individual.local.identifier)

# Processing spatial covariates
names(landscape_log_lakes) = "log_lakes"

names(landscape_log_rivers) = "log_rivers"

names(landscape_terrain) = "terrain"

names(landscape_lakes) = "lakes"

names(landscape_rivers) = "rivers"

# Logarithmizing covariates
values(landscape_log_lakes) = log(values(landscape_log_lakes))

values(landscape_log_rivers) = log(values(landscape_log_rivers))

# Removing infinite values
values(landscape_log_lakes)[is.infinite(values(landscape_log_lakes))] = NA

values(landscape_log_rivers)[is.infinite(values(landscape_log_rivers))] = NA

# Normalizing landscape values
values(landscape_log_lakes) = (values(landscape_log_lakes)-mean(values(landscape_log_lakes), na.rm = T))/ sd(values(landscape_log_lakes), na.rm = T)

values(landscape_log_rivers) = (values(landscape_log_rivers)-mean(values(landscape_log_rivers), na.rm = T))/ sd(values(landscape_log_rivers), na.rm = T)

values(landscape_terrain) = (values(landscape_terrain)-mean(values(landscape_terrain), na.rm = T))/ sd(values(landscape_terrain), na.rm =T)

values(landscape_lakes) = (values(landscape_lakes)-mean(values(landscape_lakes), na.rm = T))/ sd(values(landscape_lakes), na.rm = T)

values(landscape_rivers) = (values(landscape_rivers)-mean(values(landscape_rivers), na.rm = T))/ sd(values(landscape_rivers), na.rm = T)

# Squared terrain ruggedness index
landscape_terrain_2 = landscape_terrain
names(landscape_terrain_2) = "terrain_2"
values(landscape_terrain_2) = values(landscape_terrain_2)^2

# Transforming layers to SpatialPixelsDataFrame objects
gcov <- list(
  log_lakes = as(landscape_log_lakes$log_lakes, "SpatialPixelsDataFrame"),
  log_rivers = as(landscape_log_rivers$log_rivers, "SpatialPixelsDataFrame"),
  terrain = as(landscape_terrain$terrain, "SpatialPixelsDataFrame"),
  terrain_2 = as(landscape_terrain_2$terrain_2, "SpatialPixelsDataFrame"),
  lakes = as(landscape_lakes$lakes, "SpatialPixelsDataFrame"),
  rivers = as(landscape_rivers$rivers, "SpatialPixelsDataFrame")
)

log_lakes <- gcov$log_lakes
log_rivers <- gcov$log_rivers
terrain = gcov$terrain
terrain_2 = gcov$terrain_2
lakes = gcov$lakes
rivers = gcov$rivers

# Saving spatial covariates and female tracks as RData object
save(list = c("df1", "landscape_log_lakes" , "landscape_log_rivers" , "landscape_terrain" , "landscape_terrain_2", "landscape_lakes", "landscape_rivers", "log_lakes", "log_rivers", "terrain", "terrain_2", "lakes", "rivers"), file = "wolverine_f_data.RData")

