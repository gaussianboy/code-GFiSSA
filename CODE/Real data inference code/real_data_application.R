##### WOLVERINE DATA APPLICATION #####

# Parallel libraries
library(foreach)
library(parallel)

cl <- makeCluster(5)

# Modified function from the MoveHMM package (in order to calculate the turning angles)
turnAngle <- function(x, y, z, LLangle) {
  if (is.empty(x)) {
    return(rep(NA, nrow(z)))
  }

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

# Specifying working directory
setwd("/import/ecoc9/data-jeltsch/arceguillen/female_wolverines")

load("wolverine_f_data.RData")

# Loading packages
library(dplyr)
library(tidyverse)
library(lubridate)
library(raster)
library(rgdal)
library(sf)
library(rgeos)
library(sp)
library(INLA)
library(deldir)
library(maptools)
library(spatstat)
library(inlabru)
library(circular)
library(purrr)

inla.setOption(mkl = TRUE)

### DATA PREPARATION ###
my_df <- df1

# Transforming the time column
my_df$time_gmt <- as.POSIXlt(my_df$time_gmt, format = "%Y-%m-%d %H:%M:%OS")

my_df <- my_df %>% dplyr::filter(!(is.na(time_gmt)))

# Sorting by data by dates and animal ID
my_df <- my_df %>%
  dplyr::group_by(id) %>%
  dplyr::arrange(time_gmt, .by_group = TRUE)

# Transforming data frame to an sp object
coordinates(my_df) <- ~ lon + lat

# Setting coordinate reference system
crs(my_df) <- "+init=epsg:4326 +proj=longlat +ellps=WGS84"

# Defining the coordinate reference system
coor_rs <- landscape_terrain@crs

my_df <- spTransform(my_df, CRSobj = coor_rs)

# Transforming back to a data frame object
my_df <- as.data.frame(my_df)

# Defining list of animal identifiers
id_list <- as.list(c("F11", "F12", "F6"))

# Calculating the step lengths
filter_step_lengths <- lapply(id_list, function(u) {
  xyset <- my_df %>% filter(id == u)
  xyset <- cbind(xyset$lon, xyset$lat)
  step_lengths <- c(NA, sqrt(rowSums((xyset[-1, ] - xyset[-nrow(xyset), ])^2)))
  step_lengths <- ifelse(step_lengths == 0, 1e-6, step_lengths)
})

all_filtered <- c(filter_step_lengths[[1]], filter_step_lengths[[2]], filter_step_lengths[[3]])

my_df$my_sl <- all_filtered


# Removing outliers from the data frame
sl99 <- quantile(my_df$my_sl, 0.99, na.rm = T)

my_df <- my_df %>% filter(is.na(my_sl) | my_sl <= sl99)

## DATA PREPARATION WITH THE AMT PACKAGE ##
library(amt)

my_df <- make_track(my_df, .x = lon, .y = lat, .t = time_gmt, id = id)

my_df_nested <- my_df %>% nest(data = -"id")

issa_df <- my_df_nested %>%
  mutate(steps = map(data, function(x) {
    x %>%
      track_resample(rate = minutes(40), tolerance = minutes(3)) %>%
      steps_by_burst()
  }))

burst_df <- my_df_nested %>%
  mutate(steps = map(data, function(x) {
    x %>%
      track_resample(rate = minutes(40), tolerance = minutes(3))
  })) %>%
  select(id, steps) %>%
  unnest(cols = steps)

issa_df <- issa_df %>%
  select(id, steps) %>%
  unnest(cols = steps)

issa_df <- issa_df %>% random_steps(n_control = 200)

issa_df <- issa_df %>% dplyr::arrange(id)

issa_df <- issa_df %>%
  extract_covariates(landscape_lakes) %>%
  extract_covariates(landscape_rivers) %>%
  extract_covariates(landscape_terrain) %>%
  extract_covariates(landscape_terrain_2) %>%
  extract_covariates(landscape_log_lakes) %>%
  extract_covariates(landscape_log_rivers)

issa_df$id <- as.factor(issa_df$id)

issa_df$ANIMAL_ID1 <- issa_df$id

issa_df$case <- as.numeric(issa_df$case_)

# Using kilometers instead of meters
issa_df$sl_ <- issa_df$sl_ / 1000

issa_df$log_sl_ <- log(issa_df$sl_)

# Defining formula for iSSA model
formula.random <- case ~ -1 +
  lakes +
  rivers +
  sl_ +
  cos(ta_) +
  f(step_id_, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T))) +
  f(ANIMAL_ID1, rivers,
    model = "iid",
    hyper = list(theta = list(initial = log(1), fixed = F, prior = "pc.prec"))
  )

# Fitting the iSSA model
m0 <- inla(formula.random, family = "Poisson", data = issa_df, verbose = TRUE)

## DATA PREPARATION FOR THE GF-ISSA METHODS FOR MANY ANIMALS ##
spde_df <- as.data.frame(burst_df)

names(spde_df) <- c("id", "lon", "lat", "time", "burst")

spde_df <- do.call("rbind", lapply(id_list, function(u) {
  spde_df %>%
    dplyr::filter(id == u) %>%
    dplyr::group_by(burst) %>%
    mutate(
      select = ifelse(row_number() %in% c(1), FALSE, TRUE),
      na_ta = ifelse(row_number() %in% c(2), TRUE, FALSE)
    ) %>%
    dplyr::ungroup()
}))

# Defining data frame in order to calculate the radius of the disks of availability
lonlat_set <- data.frame(lon = spde_df$lon, lat = spde_df$lat, id = spde_df$id)

# Calculating step lengths
step_lengths <- lapply(id_list, function(u) {
  xyset <- lonlat_set %>% filter(id == u)
  xyset <- cbind(xyset$lon, xyset$lat)
  step_lengths <- sqrt(rowSums((xyset[-1, ] - xyset[-nrow(xyset), ])^2))
  step_lengths <- ifelse(step_lengths == 0, 1e-6, step_lengths)
})

radius_vector <- c(step_lengths[[1]], step_lengths[[2]], step_lengths[[3]])

# Maximum step length to be used as radius of the disks of availability
max_sl <- max(sapply(id_list, function(u) {
  xyset <- lonlat_set %>% filter(id == u)
  xyset <- cbind(xyset$lon, xyset$lat)
  step_lengths <- sqrt(rowSums((xyset[-1, ] - xyset[-nrow(xyset), ])^2))
  step_lengths <- ifelse(step_lengths == 0, 1e-6, step_lengths)
  filter_df <- spde_df %>% filter(id == u)
  step_lengths_filtered <- step_lengths[which(filter_df$select) - 1]
  max_sl <- max(step_lengths_filtered)
  return(max_sl)
}))

radius_vector[1:length(radius_vector)] <- max_sl

# Defining observed locations
observed <- lapply(id_list, function(u) {
  id_df <- spde_df %>% dplyr::filter(id == u)
  lonlat <- cbind(id_df$lon, id_df$lat)[2:nrow(id_df), ]
  return(lonlat)
})

# Defining list of consecutive time points for different animals
time_list <- lapply(id_list, function(u) {
  id_df <- spde_df %>% dplyr::filter(id == u)
  result <- as.list(2:nrow(id_df))
  return(result)
})

# Observed locations as an sp object
observed_sp <- do.call("rbind", mapply(function(u, v, s, o) {
  id_df <- spde_df %>% dplyr::filter(id == v)

  SpatialPointsDataFrame(
    coords = o,
    data = data.frame(
      case = 1,
      time = 1:nrow(o),
      sl = s[1:length(s)],
      log_sl = log(s[1:length(s)]),
      ta = sapply(u, function(t) {
        result <- turnAngle(cbind(id_df$lon[t - 2], id_df$lat[t - 2]),
          cbind(id_df$lon[t - 1], id_df$lat[t - 1]),
          cbind(id_df$lon[t], id_df$lat[t]),
          LLangle = FALSE
        )
        return(result)
      }),
      cos_ta = cos(sapply(u, function(t) {
        result <- turnAngle(cbind(id_df$lon[t - 2], id_df$lat[t - 2]),
          cbind(id_df$lon[t - 1], id_df$lat[t - 1]),
          cbind(id_df$lon[t], id_df$lat[t]),
          LLangle = FALSE
        )
        return(result)
      })),
      id = v
    ),
    proj4string = coor_rs
  )
}, o = observed, s = step_lengths, u = time_list, v = id_list, SIMPLIFY = FALSE))

# Defining samplers (domains of availability)
samplers <- st_buffer(st_as_sf(observed_sp), radius_vector)

samplers_sp <- as_Spatial(samplers)

# Defining boundary
boundary <- st_union(samplers)

boundary_sp <- as_Spatial(boundary)

# Defining fine mesh (for the approximation of the Gaussian Field)
mesh_inner <- inla.mesh.2d(
  boundary = boundary_sp,
  max.edge = 600,
  cutoff = 600,
  crs = coor_rs
)

boundary_outer <- inla.nonconvex.hull(mesh_inner$loc, convex = max_sl * 1.10)

mesh <- inla.mesh.2d(
  boundary = list(boundary_sp, boundary_outer),
  max.edge = c(600, 800),
  crs = coor_rs
)

# Defining integration points and weights
compute_ips <- function() {
  mapply(function(u, v) {
    integration_mesh <- mesh
    id_samplers_sp <- samplers_sp[samplers_sp$id == v, ]

    int_points <- ipoints(
      samplers = id_samplers_sp[unlist(u) - 1, ],
      domain = integration_mesh,
      group = "time"
    )

    calc_info <- function(t) {
      id_df <- spde_df %>% dplyr::filter(id == v)
      time_coords <- int_points[int_points$time == t - 1, ]@coords[, 1:2]

      x_t_1 <- cbind(id_df$lon, id_df$lat)[t - 1, ]
      x_t_2 <- cbind(id_df$lon, id_df$lat)[t - 2, ]

      my_sl <- pointDistance(
        x_t_1,
        time_coords,
        lonlat = FALSE
      )

      data.frame(
        sl = my_sl,
        id = v,
        log_sl = log(my_sl),
        ta = turnAngle(x_t_2, as.numeric(x_t_1), time_coords,
          LLangle = FALSE
        ),
        cos_ta = cos(turnAngle(x_t_2, as.numeric(x_t_1), time_coords,
          LLangle = FALSE
        ))
      )
    }

    result <-
      cbind(
        int_points,
        do.call("rbind", lapply(u, calc_info))
      )
    return(result)
  }, u = time_list, v = id_list, SIMPLIFY = FALSE)
}

ips <- compute_ips()

# Removing first element of the bursts (steps which were not consecutive)
ips_many <- do.call("rbind", mapply(function(i, v) {
  id_ips <- i
  id_spde_df <- spde_df %>% dplyr::filter(id == v)
  result <- id_ips[id_ips$time %in% (which(id_spde_df$select) - 1), ]
  return(result)
}, i = ips, v = id_list, SIMPLIFY = FALSE))

observed_many <- do.call("rbind", lapply(id_list, function(v) {
  id_observed <- observed_sp[observed_sp$id == v, ]
  id_spde_df <- spde_df %>% dplyr::filter(id == v)
  id_observed$ta[id_observed$time %in% (which(id_spde_df$na_ta) - 1)] <- NA
  result <- id_observed[id_observed$time %in% (which(id_spde_df$select) - 1), ]
}))

observed_many$newtimes <- 1:nrow(observed_many)

probeobs <- as.data.frame(ips_many) %>% mutate(newtimes = group_indices(., id, time))

# Adapting consecutive time points for many animals
ips_many$newtimes <- probeobs$newtimes

# Defining the priors for the GF hyperparameters
gf_model <- inla.spde2.pcmatern(
  mesh = mesh,
  alpha = 2,
  prior.range = c(500, 0.05),
  prior.sigma = c(2, 0.05)
)

# Working with kilometers instead of meters
ips_many$sl <- ips_many$sl / 1000

ips_many$log_sl <- log(ips_many$sl)

observed_many$sl <- observed_many$sl / 1000

observed_many$log_sl <- log(observed_many$sl)

# Defining model components for the GF-iSSA
comp <- coordinates ~
  -1 +
  rivers(rivers, model = "linear", mean.linear = -1, prec.linear = 0.1) +
  rivers_slope(id, model = "iid", weights = rivers, hyper = list(theta = list(initial = log(1), fixed = F, prior = "pc.prec", param = c(2, 0.05)))) +
  lakes(lakes, model = "linear", mean.linear = -1, prec.linear = 0.1) +
  sl(sl, model = "linear") +
  sl_radial_term(-log(sl), model = "const") +
  cos_ta(cos_ta, model = "linear") +
  mySmooth(coordinates, model = gf_model) +
  random_intercept(newtimes, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T)))

# Defining the likelihood structure
lik <- try(like(coordinates ~ .,
  family = "cp",
  data = observed_many,
  ips = ips_many,
  E = 1e6 # Convert the integration units to m^2 instead of km^2
), silent = TRUE)

# Fitting the model
fit <- try(bru(comp,
  lik,
  options = list(
    num.threads = 5,
    verbose = TRUE,
    control.predictor = list(compute = FALSE),
    control.compute = list(
      hyperpar = TRUE,
      return.marginals = FALSE,
      mlik = FALSE
    ),
    control.inla = list(int.strategy = "eb"),
    safe = TRUE
  )
), silent = TRUE)

# Defining model components for the NHPP model
comp2 <- coordinates  ~
  -1 +
  rivers(rivers, model = "linear", mean.linear = -1, prec.linear = 0.1) +
  rivers_slope(id, model = "iid", weights = rivers, hyper = list(theta = list(initial = log(1), fixed = F, prior = "pc.prec", param = c(2, 0.05)))) +
  lakes(lakes, model = "linear", mean.linear = -1, prec.linear = 0.1) +
  sl(sl, model = "linear") +
  sl_radial_term(-log(sl), model = "const") +
  cos_ta(cos_ta, model = "linear") +
  random_intercept(newtimes, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T)))


# Fitting the NHPP model
fit2 <- try(bru(comp2,
  lik,
  options = list(
    num.threads = 5,
    verbose = TRUE,
    control.predictor = list(compute = FALSE),
    control.compute = list(
      hyperpar = TRUE,
      return.marginals = FALSE,
      mlik = FALSE
    ),
    control.inla = list(int.strategy = "eb"),
    safe = TRUE
  )
), silent = TRUE)

# Defining model components for the NHPP_terrain model
comp3 <- coordinates ~
  -1 +
  terrain(terrain, model = "linear") +
  terrain_2(terrain_2, model = "linear") +
  rivers(rivers, model = "linear", mean.linear = -1, prec.linear = 0.1) +
  rivers_slope(id, model = "iid", weights = rivers, hyper = list(theta = list(initial = log(1), fixed = F, prior = "pc.prec", param = c(2, 0.05)))) +
  lakes(lakes, model = "linear", mean.linear = -1, prec.linear = 0.1) +
  sl(sl, model = "linear") +
  sl_radial_term(-log(sl), model = "const") +
  cos_ta(cos_ta, model = "linear") +
  random_intercept(newtimes, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = T)))

# Fitting the NHPP_terrain model
fit3 <- try(bru(comp3,
  lik,
  options = list(
    num.threads = 5,
    verbose = TRUE,
    control.predictor = list(compute = FALSE),
    control.compute = list(
      hyperpar = TRUE,
      return.marginals = FALSE,
      mlik = FALSE
    ),
    control.inla = list(int.strategy = "eb"),
    safe = TRUE
  )
), silent = TRUE)

# Table GF-iSSA
wolverine_table_f1 <- data.frame(
  Parameter = c(rownames(fit$summary.fixed), rownames(fit$summary.hyperpar)[1]),
  Estimate = c(fit$summary.fixed[, 1], fit$summary.hyperpar$mean[1]),
  lower = c(fit$summary.fixed[, 3], fit$summary.hyperpar[1, 3]),
  upper = c(fit$summary.fixed[, 5], fit$summary.hyperpar[1, 5])
)

# Table NHPP
wolverine_table_f2 <- data.frame(
  Estimate = c(fit2$summary.fixed[, 1], fit2$summary.hyperpar$mean[1]),
  lower = c(fit2$summary.fixed[, 3], fit2$summary.hyperpar[1, 3]),
  upper = c(fit2$summary.fixed[, 5], fit2$summary.hyperpar[1, 5])
)

# Table NHPP_terrain
wolverine_table_f3 <- data.frame(
  Parameter = c(rownames(fit3$summary.fixed), rownames(fit3$summary.hyperpar)[1]),
  Estimate = c(fit3$summary.fixed[, 1], fit3$summary.hyperpar$mean[1]),
  lower = c(fit3$summary.fixed[, 3], fit3$summary.hyperpar[1, 3]),
  upper = c(fit3$summary.fixed[, 5], fit3$summary.hyperpar[1, 5])
)

wolverine_table <- cbind(wolverine_table_f1, wolverine_table_f2)

## PLOTS ##
library(wesanderson)
library(ggthemes)

# Defining palette
colsc <- function(...) {
  scale_fill_gradientn(
    colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu")),
    limits = range(..., na.rm = TRUE)
  )
}

# Defining lattice object
df_pixels <- pixels(mesh, mask = boundary_sp, nx = 200, ny = 200)

# Linear predictor estimation
gf <- predict(fit, df_pixels, ~ list(
  fixed = rivers + lakes,
  field = mySmooth,
  complete = rivers + lakes + mySmooth
), include = c("rivers", "lakes", "mySmooth"))

csc <- colsc(gf$complete$mean, gf$fixed$mean)

csc2 <- colsc(gf$field$mean)

my_theme <-
  theme(
    plot.title = element_text(color = "black", size = 35, face = "bold.italic", hjust = 0.5),
    axis.text.x = element_text(color = "grey20", size = 25, angle = 45, hjust = .5, vjust = .5, face = "plain"),
    axis.text.y = element_text(color = "grey20", size = 25, angle = 0, hjust = 1, vjust = 0, face = "plain"),
    axis.title.x = element_text(color = "grey20", size = 35, angle = 0, hjust = .5, vjust = 0, face = "bold.italic"),
    axis.title.y = element_text(color = "grey20", size = 35, angle = 90, hjust = .5, vjust = .5, face = "bold.italic"),
    legend.text = element_text(color = "grey20", size = 35), legend.direction = "vertical", legend.position = "right",
    legend.title = element_text(color = "grey20", size = 35, face = "bold.italic"), legend.key.size = unit(1, "cm"),
    strip.text = element_text(size = 35)
  )

gf_field <- ggplot(boundary) +
  gg(gf$field, aes(fill = mean), size = 5) +
  geom_sf(linewidth = 2, color = "black", alpha = 0.0001) +
  csc2 +
  labs(fill = "Mean") +
  theme_minimal() +
  xlab("Eastings") +
  ylab("Northings") +
  ggtitle("GF: Mean estimate") +
  my_theme

gf_complete <- ggplot(boundary) +
  gg(gf$complete, aes(fill = mean), size = 5) +
  geom_sf(linewidth = 2, color = "black", alpha = 0.0001) +
  csc +
  labs(fill = "Mean") +
  gg(observed_many, alpha = 0.2) +
  theme_minimal() +
  xlab("Eastings") +
  ylab("Northings") +
  ggtitle("Spatial covariates and GF") +
  my_theme

gf_fixed <- ggplot(boundary) +
  gg(gf$fixed, aes(fill = mean), size = 5) +
  geom_sf(linewidth = 2, color = "black", alpha = 0.0001) +
  csc +
  labs(fill = "Mean") +
  gg(observed_many, alpha = 0.2) +
  theme_minimal() +
  xlab("Eastings") +
  ylab("Northings") +
  ggtitle("Spatial covariates") +
  my_theme

# Saving plots
ggsave("gf_field_finer_wolverine.png", gf_field, width = 15, height = 10)

ggsave("gf_complete_finer_wolverine.png", gf_complete, width = 15, height = 10)

ggsave("gf_fixed_finer_wolverine.png", gf_fixed, width = 15, height = 10)

library(cowplot)

gf_all <- plot_grid(gf_fixed,
  gf_field,
  gf_complete,
  ncol = 3, nrow = 1
)

ggsave("gf_all_finer_wolverine.png", gf_all, width = 33.5, height = 10)

# Saving summary tables
saveRDS(wolverine_table, "wolverine_finer_table.rds")
saveRDS(wolverine_table_f3, "wolverine_finer_table_terrain.rds")

# Saving all objects as an RData object
save.image(file = "wolverine_finer_GFiSSA.RData")

parallel::stopCluster(cl)
