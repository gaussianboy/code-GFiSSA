##### simulation code based on code from Ulrike Schlägel #####

# Called by script "paper_simulation.R"

# This function simulates movement based on iSSF movement densities

# Parameter list
#** movement kernel (Gamma distribution): shape, scale
#** movement kernel (Von Mises distribution): concentration
#** weighting function: Resource selection function

# Further arguments:
#** center : where individual starts (later: centralising tendency)
#** T: number of time points
#** seed: seed for set.seed
#** vzero: small number replacing step length=0 (to avoid numerical problems
#   due to discretisation, i.e.  dividing by 0 - in original script 0.5 was used)
#** landscape: landscape map (raster object with values describing the presence or
#   absence of the resource the individual is attracted to

# Modified function from MoveHMM package (in order to calculate the turning angles)
turnAngle <- function(x, y, z, LLangle) {
  # NA angle if zero step length
  if (all(x == y) | all(y == z)) {
    return(NA)
  }

  if (LLangle) {
    angle <- (bearing(x, y) - bearing(y, z)) / 180 * pi
  } else {
    v <- c(y[1] - x[1], y[2] - x[2])
    w <- matrix(c(z[, 1] - y[1], z[, 2] - y[2]), ncol = 2, byrow = FALSE)
    angle <- atan2(w[, 2], w[, 1]) - atan2(v[2], v[1])
  }

  for (i in 1:length(angle)) {
    while (angle[i] <= -pi) {
      angle[i] <- angle[i] + 2 * pi
    }
    while (angle[i] >= pi) {
      angle[i] <- angle[i] - 2 * pi
    }
  }

  return(angle)
}


simdata_generic <- function(parameter, center, Ts, seed, vzero, landscape) {
  loc <- matrix(NA, ncol = 2, nrow = Ts) # location of the animal

  colnames(loc) <- c("x", "y")

  loc <- as.data.frame(loc)

  # Note: need to convert to data frame so that plotting the points in raster works
  steps <- numeric(Ts) # save selected step length

  cells <- ncell(landscape) # number of cells in the raster

  # Set seed
  set.seed(seed)

  # Initial 2 locations are random
  ext_hr <- 20

  loc[1, ] <- c(
    sample((center[1] - ext_hr):(center[1] + ext_hr), 1),
    sample((center[2] - ext_hr):(center[2] + ext_hr), 1)
  )

  loc[2, ] <- c(
    sample((center[1] - ext_hr):(center[1] + ext_hr), 1),
    sample((center[2] - ext_hr):(center[2] + ext_hr), 1)
  )

  # Coordinates from raster
  xy <- xyFromCell(landscape, 1:length(landscape$x1))

  # Calculations SLs and TAs at each time point
  for (ts in 3:Ts) {
    # SL based on the last observed location
    sl1 <- distanceFromPoints(landscape, loc[ts - 1, ])

    # TA based on the last observed location
    ta <- turnAngle(as.numeric(loc[ts - 2, ]), as.numeric(loc[ts - 1, ]), xy, LLangle = FALSE)

    # Setting value of current location to a small number because otherwise gamma is NA (at zero)
    sl1[cellFromXY(landscape, loc[ts - 1, ])] <- vzero

    # Calculating whole movement kernel
    kern1 <- dgamma(values(sl1), shape = parameter$shape, rate = parameter$rate) *
      dvonmises(circular(ta), mu = circular(0), kappa = parameter$kappa) /
      (values(sl1))

    # Calculating the resource selection function
    wfun <- exp(as.matrix(values(landscape)) %*% parameter$beta)

    # Calculating step probability
    step.prob <- (kern1 * wfun) / sum(kern1 * wfun)

    # Sampling a new step
    temp <- sample(1:cells, 1, prob = step.prob)
    loc[ts, ] <- xyFromCell(landscape, temp)
    steps[ts - 1] <- sl1[cellFromXY(landscape, loc[ts, ])]
    #
  }
  return(cbind(loc, steps))
}
