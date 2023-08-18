####### ANIMAL TRACKS SIMULATION #######

# set working directory
setwd("/import/ecoc9j/data-jeltsch/arceguillen/")

# Loading libraries
library(raster)
library(CircStats)
library(geoR)
library(tictoc)
library(foreach)
library(doParallel)
library(circular)

# Parallel kernels specification
no_cores <- 25

cl <- makeCluster(no_cores, type = "FORK")

registerDoParallel(cl)

# source function to simulate data
source("fun_complete.R")

## Defining Landscapes ##

# Dimension of space
ext <- 100 # area of 2*ext x 2*ext

resol <- 0.2 # grid resolution

# Coordinate reference system
proj.utm <- CRS("+proj=utm +zone=33 +datum=WGS84")

refrast <- raster(
  xmn = -ext - resol / 2, ymn = -ext - resol / 2,
  xmx = ext + resol / 2, ymx = ext + resol / 2,
  resolution = resol, crs = proj.utm
)

cells <- ncell(refrast) # number of cells

values(refrast) <- 1

refrast

plot(refrast)

coord_refrast <- rasterToPoints(refrast) # ccordinates of raster


# Defining centralizing tendency
# Setting starting point of the track
cen <- c(0, 0)

cenmap <- distanceFromPoints(refrast, cen)

names(cenmap) <- "cen"

plot(cenmap)


# Generating continuous covariate "x1" as a Gaussian random field
nx <- ny <- sqrt(cells)

set.seed(999)

temp <- raster(grf(cells,
  grid = "reg", nx = nx, ny = nx, cov.model = "exponential",
  xlims = c(-ext - resol / 2, ext + resol / 2),
  ylims = c(-ext - resol / 2, ext + resol / 2),
  cov.pars = c(1, 20)
))

plot(temp)

landscape_1 <- refrast

values(landscape_1) <- values(temp)

names(landscape_1) <- "x1"


# Generating discrete covariate "x2" as a Gaussian random field
set.seed(123)

temp <- raster(grf(cells,
  grid = "reg", nx = nx, ny = nx, cov.model = "exponential",
  xlims = c(-ext - resol / 2, ext + resol / 2),
  ylims = c(-ext - resol / 2, ext + resol / 2),
  cov.pars = c(1, 10)
))

quant <- quantile(values(temp), probs = 0.60)

values(temp)[values(temp) > quant] <- 1

values(temp)[values(temp) <= quant] <- 0

plot(temp)

landscape_2 <- refrast

values(landscape_2) <- values(temp)

names(landscape_2) <- "x2"

# Generating GAUSSIAN FIELD scenarios and tracks
for (range in c(30, 40, 50)) {
  for (sigma_2 in c(1, 4, 9)) {
    # Suitable seeds to avoid high correlated landscapes
    if (range == 30) {
      set.seed(911)
    }

    if (range == 40) {
      set.seed(1909)
    }

    if (range == 50) {
      set.seed(919)
    }

    temp <- raster(grf(cells,
      grid = "reg", nx = nx, ny = nx, cov.model = "matern",
      xlims = c(-ext - resol / 2, ext + resol / 2),
      ylims = c(-ext - resol / 2, ext + resol / 2),
      cov.pars = c(sigma_2, range)
    ))

    plot(temp)

    landscape_3 <- refrast

    values(landscape_3) <- values(temp)

    names(landscape_3) <- "x3"

    ## Simulation specifications ##

    # Time series length
    Ts <- 1020

    # Making fake times - arbitrary
    time.begin <- as.POSIXct("2017-11-25 12:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "Europe/Berlin")

    dates <- time.begin + seq(0, (Ts - 1) * 360, 360)

    parameter <- list()

    parameter$shape <- 4 # Gamma movement distribution

    parameter$rate <- 2 # Gamma movement distribution

    parameter$beta <- c(1.5, 1, -0.04, 1) # Selection coefficients

    parameter$kappa <- 1 # Von Mises movement distribution

    ##### Simulation #####

    # stack the landscapes that should be used for the simulation
    landscape <- stack(landscape_1, landscape_2, cenmap, landscape_3)

    my_sim <- function(i) {
      result <- simdata_generic(
        parameter = parameter, center = cen, Ts = Ts, seed = i,
        vzero = 1e-12, landscape = landscape
      )

      result <- result[21:Ts, ] # first 20 steps for initialisation

      result$t <- dates[21:Ts]

      return(result)
    }

    my_data <- foreach(i = 1:25) %dopar% my_sim(i)

    setwd("/import/ecoc9/data-jeltsch/arceguillen/")

    save(list = c("my_data", "landscape"), file = paste("paper_tracks", range, "_", sigma_2, ".RData", sep = ""))
  }
}
