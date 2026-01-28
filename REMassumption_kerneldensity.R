#     R code to replicate the simulations developed on 
#    "Simulations illustrate the usefulness of camera trapping to estimate 
#     wolf densities using the random encounter model"
#
#    Assessing REM assumption when applied to wolves

# please note that wolf locations are sensible data, so we can't publish this info.
# However, the most relevant aspect on this section is how "random", "attracted"
# and "avoided" point patterns were simulated, and how they were compared against
# wolf locations

# 1. Study area definition (example) ----
x_min <- 0
x_max <- 10000
y_min <- 0
y_max <- 10000

lp <- 50 # pixel size (m)
nc <- 64 # number of cameras
n_x <- ceiling((x_max - x_min) / lp)
n_y <- ceiling((y_max - y_min) / lp)

# 2. Camera placements ----
# 2.1 Regular-grid cameras ----
sqrt_nc <- sqrt(nc)
regridCTs_camera_locations <- data.frame(
  X = rep(round(seq(x_min + 50, x_max - 50, length.out = sqrt_nc), 4), times = sqrt_nc),
  Y = rep(round(seq(y_min + 50, y_max - 50, length.out = sqrt_nc), 4), each = sqrt_nc)
)

# 2.2 Random cameras ----
randomCTs_camera_locations <- data.frame(
  X = runif(nc, min = x_min + 50, max = x_max - 50),
  Y = runif(nc, min = y_min + 50, max = y_max - 50)
)

# 3. Point patterns ----
# 3.1 Random points ----
n_wolf_locs <- 252 # Let's work with a seven-days-long period in which 252 wolf locations were recorded.
random_points <- data.frame(
  X = runif(n_wolf_locs, x_min, x_max),
  Y = runif(n_wolf_locs, y_min, y_max)
)
kde_random <- kde2d(random_points$X, random_points$Y, n = c(n_x, n_y), lims = c(x_min, x_max, y_min, y_max), h = c(1000, 1000))

random_density_grid <- interp.surface(kde_random, regridCTs_camera_locations) # kernel density at grid CTs
random_density_rand <- interp.surface(kde_random, randomCTs_camera_locations) # kernel density at random CTs

# 3.2 Attracted points ----
# CTs random
near_camera_points_rand <- data.frame(
  X = rep(randomCTs_camera_locations$X, length.out = n_wolf_locs) + 
    rnorm(n_wolf_locs, mean = 0, sd = 100),
  Y = rep(randomCTs_camera_locations$Y, length.out = n_wolf_locs) + 
    rnorm(n_wolf_locs, mean = 0, sd = 100)
)
kde_attracted <- kde2d(near_camera_points_rand$X, near_camera_points_rand$Y, n = c(n_x, n_y), lims = c(x_min, x_max, y_min, y_max), h = c(1000, 1000))
attracted_density_rand <- interp.surface(kde_attracted, randomCTs_camera_locations) # kernel density at random CTs

# CTs grid
near_camera_points_grid <- data.frame(
  X = rep(regridCTs_camera_locations$X, length.out = n_wolf_locs) + 
    rnorm(n_wolf_locs, mean = 0, sd = 100),
  Y = rep(regridCTs_camera_locations$Y, length.out = n_wolf_locs) + 
    rnorm(n_wolf_locs, mean = 0, sd = 100)
)

kde_attracted <- kde2d(near_camera_points_grid$X, near_camera_points_grid$Y, n = c(n_x, n_y), lims = c(x_min, x_max, y_min, y_max), h = c(1000, 1000))
attracted_density_grid <- interp.surface(kde_attracted, regridCTs_camera_locations) # kernel density at grid CTs

# 3.3 Avoided points ----

# CTs random
candidates <- data.frame(
  X = runif(n_wolf_locs*10, x_min, x_max),
  Y = runif(n_wolf_locs*10, y_min, y_max)
) # candidate points

distances_to_cameras <- apply(candidates, 1, function(p) {
  min(spDistsN1(as.matrix(randomCTs_camera_locations), p, longlat = FALSE))
}) # distance to the CTs
weights <- exp(0.001 * distances_to_cameras)
indices <- sample(seq_len(n_wolf_locs*10), size = n_wolf_locs, replace = TRUE, prob = weights)
far_points <- candidates[indices, ]

kde_far_rand <- kde2d(far_points$X, far_points$Y, n = c(n_x, n_y), lims = c(x_min, x_max, y_min, y_max), h = c(1000, 1000))
avoided_density_rand <- interp.surface(kde_far_rand, randomCTs_camera_locations) # kernel density at random CTs

# CTs grid
candidates <- data.frame(
  X = runif(n_wolf_locs*10, x_min, x_max),
  Y = runif(n_wolf_locs*10, y_min, y_max)
) # candidate points

distances_to_cameras <- apply(candidates, 1, function(p) {
  min(spDistsN1(as.matrix(regridCTs_camera_locations), p, longlat = FALSE))
}) # distance to the CTs
weights <- exp(0.001 * distances_to_cameras)
indices <- sample(seq_len(n_wolf_locs*10), size = n_wolf_locs, replace = TRUE, prob = weights)
far_points <- candidates[indices, ]

kde_far_grid <- kde2d(far_points$X, far_points$Y, n = c(n_x, n_y), lims = c(x_min, x_max, y_min, y_max), h = c(1000, 1000))
avoided_density_grid <- interp.surface(kde_far_grid, regridCTs_camera_locations) # kernel density at grid CTs

# 4. Storing results ----

# CTs random
result_df_rand <- data.frame(
  period = 1,  # dummy
  npoints = n_wolf_locs,
  lado_pixel = lp,
  n_cts = nc,
  design = "random",
  mean_random_density = mean(random_density_rand, na.rm = TRUE),
  mean_attracted_density = mean(attracted_density_rand, na.rm = TRUE),
  mean_avoided_density = mean(avoided_density_rand, na.rm = TRUE)
)
results_list <- list()
results_list[[length(results_list) + 1]] <- result_df_rand

# CTs grid
result_df_grid <- data.frame(
  period = 1,  # dummy
  npoints = n_wolf_locs,
  lado_pixel = lp,
  n_cts = nc,
  design = "grid",
  mean_random_density = mean(random_density_grid, na.rm = TRUE),
  mean_attracted_density = mean(attracted_density_grid, na.rm = TRUE),
  mean_avoided_density = mean(avoided_density_grid, na.rm = TRUE)
)
results_list[[length(results_list) + 1]] <- result_df_grid
results_df <- do.call(rbind, results_list) # CTs grids & random 
