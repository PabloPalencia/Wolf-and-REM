
#     R code to replicate the simulations developed on 
#    "Simulations illustrate the usefulness of camera trapping to estimate 
#     wolf densities using the random encounter model"
#
#     Simulated population sampled with camera traps

# simulation
library(sp)
library(dplyr)
library(sf)
library(CircStats)

# analysis
library(Distance)
library(stringr)
library(trappingmotion)

# functions (not available in R packages)
source("required_functions.R")

# 1. Defining simulations parameters ----
## 1.1 Survey design ----
areas_km2 <- 800 # study area size
puntos <- 196 # camera placements
dayssimulated <- 90 # sampling period (days)
sector_radius <- 20 # camera detection zone - radius
fov_angle <- pi / 4 # camera detection zone - angle
sigma <- 6 # imperfect detection - sigma parameter normal distribution
buffer <- 50  # buffer from the study area border (m) to place cameras

## 1.2 Wolf population ----
n_groups <- 16 # number of groups
groupSize <- 10 # individuals per group
HR_simulated <- 50 #km2
groupCohesion <- 0.6 # proportion of the time travelling together
groupSpacing <- 50 # sd for the normal distribution (m)
floaters_perc <- 0.25 # proportion of the time travelling together
stepsDailyPerAnimal <- 12*3 # (1loc/20min) - Spanish wolves-
HR_simulated <- 50 # home range simulated (km2)
homeRangeRadius <- round(sqrt(HR_simulated/pi), 2)*1000 # home range radius (m)
step_sec <- 1200  #seconds between consecutive locations (i.e 1 loc/20 min)

# Parameters of the step length and turning angle distribution. These parameters
# were obtained from analysing the data from 13 GPS collared wolves using 
# moveHMM:fitHMM
stepPar_breedingpair <- matrix(c(
  0.008503947, 0.006226183,
  0.2049300,   0.2034753,
  0.8704988,   0.4539068
), nrow=2, byrow=FALSE)

anglePar_breedingpair <- matrix(c(
  -3.1047632,  0.5366011,
  0.23830675, 0.09991957,
  -0.02550855, 1.29108773
), nrow=2, byrow=FALSE)

stepPar_helper <- matrix(c(
  0.006074573, 0.003998430,
  0.04127050,  0.03698423,
  0.6591254,   0.4305531
), nrow=2, byrow=FALSE)

anglePar_helper <- matrix(c(
  -3.1074527,  0.6336442,
  -2.4653042,  0.4268319,
  -0.1014530,  0.9437155
), nrow=2, byrow=FALSE)

stepPar_floater <- matrix(c(
  0.008326256, 0.005653181,
  0.1234742,   0.1223342,
  0.6700017,   0.4129365
), nrow=2, byrow=FALSE)

anglePar_floater <- matrix(c(
  -2.9360288,  0.5825553,
  -3.070981,   0.308978,
  0.05339519, 0.94337257
), nrow=2, byrow=FALSE)

# 2. Simulating the population ----

## 2.1 Simulating the wolves' trajectories ----
total_population <- n_groups * groupSize # group-living individuals 
floater <- round(floaters_perc * total_population) # Number of floaters
nbAnimals <- groupSize * n_groups + floater # total number of individuals
obsPerAnimal <- dayssimulated * stepsDailyPerAnimal # locations per individual

# simulating breeding pairs
simulated_data_breedingpair <- simDataHR_group(nbAnimals = n_groups*2, 
                                               nbStates = 3, 
                                               stepDist = "gamma", 
                                               angleDist = "vm", 
                                               stepPar = stepPar_breedingpair,
                                               anglePar = anglePar_breedingpair,
                                               obsPerAnimal = obsPerAnimal, 
                                               groupSize = 2, 
                                               groupSpacing = groupSpacing,
                                               groupCohesion = groupCohesion,
                                               springEffect = TRUE, 
                                               k_spring = 1, 
                                               homeRangeRadius = homeRangeRadius)
simulated_data_helper <- simulated_data_breedingpair[0, ]

# simulating helpers
if ((groupSize - 2) > 0) {
              simulated_data_helper <- simDataHR_group(
                nbAnimals = (groupSize - 2) * n_groups, 
                nbStates = 3, 
                stepDist = "gamma", 
                angleDist = "vm", 
                stepPar = stepPar_helper,
                anglePar = anglePar_helper,
                obsPerAnimal = obsPerAnimal, 
                groupSize = (groupSize - 2), 
                groupSpacing = groupSpacing,
                groupCohesion = groupCohesion,
                springEffect = TRUE, 
                k_spring = 1, 
                homeRangeRadius = homeRangeRadius)
              
              simulated_data_helper$ID <- vapply(simulated_data_helper$ID, function(id) {
                match <- regexec("(_A)([0-9]+)", id)
                parts <- regmatches(id, match)[[1]]
                if (length(parts) == 3) {
                  new_num <- as.numeric(parts[3]) + 2
                  sub("(_A)[0-9]+", paste0(parts[2], new_num), id)
                } else {
                  id
                }
              }, character(1))
              
            } else {
              message("No helpers. Group size 2.")
            }
            
simulated_data_floater <- simulated_data_breedingpair[0, ] # in case floaters are not simulated
            
# simulating floaters
if ((floater) > 0) {
              simulated_data_floater <- simDataHR_group(nbAnimals = floater, 
                                                        nbStates = 3, 
                                                        stepDist = "gamma", 
                                                        angleDist = "vm", 
                                                        stepPar = stepPar_floater,
                                                        anglePar = anglePar_floater,
                                                        obsPerAnimal = obsPerAnimal, 
                                                        groupSize = 1, 
                                                        groupSpacing = groupSpacing,
                                                        groupCohesion = groupCohesion,
                                                        springEffect = TRUE, 
                                                        k_spring = 1, 
                                                        homeRangeRadius = homeRangeRadius)
              
              summary(simulated_data_floater$step)
              unique_ids_floater <- unique(simulated_data_floater$ID)
              new_ids_floater <- paste0("G", n_groups + seq_along(unique_ids_floater), "_A1")
              simulated_data_floater$ID <- new_ids_floater[match(simulated_data_floater$ID, unique_ids_floater)]
              
            } else {
              message("Floaters not simulated")
            }
            
simulated_data_breedingpair_and_helper <- rbind(simulated_data_breedingpair, simulated_data_helper)

# Study area borders
studyAreas <- lapply(areas_km2, function(area) {
  half_side_m <- round(sqrt(area) * 1000 / 2)
  c(-half_side_m, half_side_m)
})
area <- studyAreas[[seq_along(studyAreas)]]
xmin <- area[1]
xmax <- area[2]
ymin <- area[1]
ymax <- area[2]
            
# Parameters
n_groups_sim <- n_groups  

# Generating home ranges centers
ncol <- ceiling(sqrt(n_groups_sim))
nrow <- ceiling(n_groups_sim / ncol)
x_range <- xmax - xmin
y_range <- ymax - ymin
quadrant_width <- x_range / ncol
quadrant_height <- y_range / nrow
            
group_centers <- list()
            
for (i in 1:n_groups_sim) {
      col_index <- (i - 1) %% ncol
      row_index <- floor((i - 1) / ncol)
              
      quad_xmin <- xmin + col_index * quadrant_width
      quad_xmax <- quad_xmin + quadrant_width
      quad_ymin <- ymin + row_index * quadrant_height
      quad_ymax <- quad_ymin + quadrant_height
              
  if ((quad_xmax - quad_xmin) < 2 * homeRangeRadius || (quad_ymax - quad_ymin) < 2 * homeRangeRadius) {
                center_x <- (quad_xmin + quad_xmax) / 2
                center_y <- (quad_ymin + quad_ymax) / 2
              } else {
                center_x <- runif(1, quad_xmin + homeRangeRadius, quad_xmax - homeRangeRadius)
                center_y <- runif(1, quad_ymin + homeRangeRadius, quad_ymax - homeRangeRadius)
              }
              
              group_centers[[i]] <- c(center_x, center_y)
            }

group_centers_df <- as.data.frame(do.call(rbind, group_centers))
colnames(group_centers_df) <- c("x", "y")
group_centers_df$group <- paste0("G", 1:n_groups_sim)
delta_xy <- group_centers_df[, c("group", "x", "y")]
simulated_data_breedingpair_and_helper$group <- sub("(_.*)$", "", simulated_data_breedingpair_and_helper$ID)

simulated_data_breedingpair_and_helper <- simulated_data_breedingpair_and_helper %>%
  left_join(delta_xy, by = "group", suffix = c("", "_delta")) %>%
  mutate(
    x = ifelse(!is.na(x_delta), x + x_delta, x),
    y = ifelse(!is.na(y_delta), y + y_delta, y)
  ) %>%
  dplyr::select(-ends_with("_delta"))

simulated_data_floater$group <- sub("(_.*)$", "", simulated_data_floater$ID)
unique_groups <- unique(simulated_data_floater$group)

for(g in unique_groups) {
  delta_x <- runif(1, xmin+(homeRangeRadius), xmax-(homeRangeRadius))
  delta_y <- runif(1, ymin+(homeRangeRadius), ymax-(homeRangeRadius))
  
  idx <- simulated_data_floater$group == g
  simulated_data_floater$x[idx] <- simulated_data_floater$x[idx] + delta_x
  simulated_data_floater$y[idx] <- simulated_data_floater$y[idx] + delta_y
}

simulated_data_breedingpair_and_helper <- simulated_data_breedingpair_and_helper[, !(names(simulated_data_breedingpair_and_helper) == "group")]
simulated_data_floater <- simulated_data_floater[, !(names(simulated_data_floater) == "group")]

# joining all the individuals locations
simulated_data <- rbind(simulated_data_breedingpair_and_helper, 
                        simulated_data_floater)

## 2.2.Simulating the camera trap placements ----
# regular grids
mallas <- lapply(puntos, crear_malla, 
                 xmin = xmin, xmax = xmax, 
                 ymin = ymin, ymax = ymax,
                 buffer = buffer)
malla_completa <- do.call(rbind, mallas)
names(malla_completa)[names(malla_completa) == "puntos"] <- "n_cts"
names(malla_completa)[names(malla_completa) == "id"] <- "cameraId"
malla_completa <- malla_completa[, !names(malla_completa) %in% "densidad"]
malla_completa$type <- "grid"

# random cameras
cameras_random <- generar_cameras_aleatorios(196, xmin, xmax, ymin, ymax, buffer = buffer, min_distancia = 100)
names(cameras_random)[names(cameras_random) == "n"] <- "n_cts"
names(cameras_random)[names(cameras_random) == "id"] <- "cameraId"
cameras_random <- cameras_random[, !names(cameras_random) %in% "buffer"]
cameras_random <- cameras_random[, !names(cameras_random) %in% "min_distancia"]
cameras_random$type <- "random"

# Joining random and grid cameras
cameras_total <- rbind(cameras_random, malla_completa)
cameras_total$orientation <- runif(nrow(cameras_total), 0, 2 * pi)  # random heading -radians-
cameras_total$radius <- sector_radius
cameras_total$fov <- fov_angle

# List to store camera detection zone
camera_polygons <- vector("list", nrow(cameras_total))

for (i in 1:nrow(cameras_total)) {
  center <- c(cameras_total$x[i], cameras_total$y[i])
  start_angle <- cameras_total$orientation[i] - cameras_total$fov[i] / 2
  end_angle <- cameras_total$orientation[i] + cameras_total$fov[i] / 2
  angles <- seq(start_angle, end_angle, length.out = 50)
  
  sector_coords <- cbind(center[1] + cameras_total$radius[i] * cos(angles),
                         center[2] + cameras_total$radius[i] * sin(angles))
  poly_coords <- rbind(center, sector_coords, center)
  
  camera_polygons[[i]] <- poly_coords
  
}

# Step-ID
simulated_data <- simulated_data %>%
  group_by(ID) %>%
  mutate(stepID = row_number()) %>%
  ungroup()

# Creating a dataframe to store the results
results <- data.frame(cameraId = character(),
                      animalID = character(),
                      cam_x = numeric(),
                      cam_y = numeric(),
                      animal_x = numeric(),
                      animal_y = numeric(),
                      stepID = numeric(),
                      step = numeric(),
                      angle = numeric(),
                      stringsAsFactors = FALSE)

# Creating trajectories
n_points <- nrow(simulated_data)
segments_list <- vector("list", n_points - 1)

for (j in 1:(n_points - 1)) {
  coords <- matrix(c(simulated_data$x[j], simulated_data$y[j],
                     simulated_data$x[j + 1], simulated_data$y[j + 1]),
                   ncol = 2, byrow = TRUE)
  segments_list[[j]] <- st_linestring(coords)
}

# Create sf object
segments_sf <- st_sf(geometry = st_sfc(segments_list))
st_crs(segments_sf) <- 32630  # EPSG:32630 (UTM 30N)

segments_sf$ID       <- simulated_data$ID[1:(n_points - 1)]
segments_sf$stepID   <- simulated_data$stepID[1:(n_points - 1)]
segments_sf$step     <- simulated_data$step[1:(n_points - 1)]
segments_sf$angle    <- simulated_data$angle[1:(n_points - 1)]
segments_sf$animal_x <- simulated_data$x[1:(n_points - 1)]
segments_sf$animal_y <- simulated_data$y[1:(n_points - 1)]

results_sf <- data.frame() # to store the results

# Encounters at the camera level
for (i in 1:nrow(cameras_total)) {
  cam_id <- cameras_total$cameraId[i]
  
  poly_coords <- camera_polygons[[i]]
  if (!all(poly_coords[1, ] == poly_coords[nrow(poly_coords), ])) {
    poly_coords <- rbind(poly_coords, poly_coords[1, ])
  }
  cam_poly <- st_polygon(list(poly_coords))
  cam_poly_sf <- st_sf(geometry = st_sfc(cam_poly), crs = st_crs(segments_sf))
  
  intersects <- st_intersects(segments_sf, cam_poly_sf, sparse = FALSE)[, 1]
  idx <- which(intersects)
  
  if (length(idx) > 0) {
    for (j in idx) {
      seg <- segments_sf[j, ]
      inter_geom <- st_intersection(seg$geometry, cam_poly_sf$geometry)
      
      if (length(inter_geom) > 0) {
        midpoint <- st_line_sample(inter_geom, sample = 0.5)  
        mid_coords <- st_coordinates(midpoint)
        inter_x <- mid_coords[1, "X"]
        inter_y <- mid_coords[1, "Y"]
      } else {
        inter_x <- NA
        inter_y <- NA
      }
      
      temp <- data.frame(
        cameraId    = cam_id,
        animalID    = seg$ID,
        cam_x       = cameras_total$x[i],
        cam_y       = cameras_total$y[i],
        
        animal_x    = inter_x,   
        animal_y    = inter_y,  
        stepID      = seg$stepID,
        step        = seg$step,
        angle       = seg$angle,
        
        stringsAsFactors = FALSE
      )
      
      results_sf <- rbind(results_sf, temp)
    }
  }
}

results <- results_sf

# Creating sequences
results <- results %>%
  arrange(animalID, stepID, cameraId) %>%
  mutate(
    seqID = cumsum(
      c(TRUE, 
        diff(stepID) != 1 | 
          cameraId[-1] != cameraId[-n()] |
          animalID[-1] != animalID[-n()]
      )
    )
  )

results <- results %>%
  mutate(
    type = str_extract(cameraId, "(?<=P)[a-z]+"),  # "grid" vs "random"
    n_cts = str_extract(cameraId, "(?<=_)\\d+")   # IDs
  )
results$n_cts <- as.numeric(results$n_cts)

results_1month <- subset(results, stepID <=30*stepsDailyPerAnimal)
results_1month$oper <- 30
results_2month <- subset(results, stepID <=2*30*stepsDailyPerAnimal)
results_2month$oper <- 60
results_3month <- subset(results, stepID <=3*30*stepsDailyPerAnimal)
results_3month$oper <- 90

super_studyarea_km2 <- round(((ymax-ymin)*(xmax-xmin))/1000000, 0)
true_density <- nbAnimals/super_studyarea_km2

## Imperfect detection
results_impdet <- results

# Euclidean distance between animals and CTs "distancetoCT"
results_impdet$distancetoCT <- sqrt((results_impdet$animal_x - results_impdet$cam_x)^2 +
                                      (results_impdet$animal_y - results_impdet$cam_y)^2)

# Probability of detection decreased as distance increased
#sigma <- 6 # 4.44 m -EDR - wolf Palencia et al., 2022 - Transboundary

# Estimating probability of detection
results_impdet$det_prob <- exp(- (results_impdet$distancetoCT^2) / (2 * sigma^2))
results_impdet$detection <- rbinom(n = nrow(results_impdet), size = 1, prob = results_impdet$det_prob)

# Selecting encounters detected
results_impdet_detected <- subset(results_impdet, detection == 1)

# Three sampling efforts (1, 2 and 3 months)
results_1month <- subset(results_impdet_detected, stepID <=30*stepsDailyPerAnimal)
results_1month$oper <- 30
results_2month <- subset(results_impdet_detected, stepID <=2*30*stepsDailyPerAnimal)
results_2month$oper <- 60
results_3month <- subset(results_impdet_detected, stepID <=3*30*stepsDailyPerAnimal)
results_3month$oper <- 90

### 1 month ----
# unique combinations type(random vs grid) & number of cameras (n_cts)
combinaciones_unicas <- results_1month %>%
  distinct(type, n_cts)

resultados <- list()

w_rad <- 10 # truncation distance for distance sampling models
for(i in 1:nrow(combinaciones_unicas)) {
  current_type <- combinaciones_unicas$type[i]
  current_n_cts <- combinaciones_unicas$n_cts[i]
  
  results_filtrado <- results_1month %>%
    filter(type == current_type, n_cts == current_n_cts)
  
  cameras_total_filtrado <- cameras_total %>%
    filter(type == current_type, n_cts == current_n_cts)
  
  # Encounters per camera 
  conteo <- results_filtrado %>%
    add_count(cameraId, name = "seq")
  
  results_filtrado2 <- cameras_total_filtrado %>%
    left_join(results_filtrado %>%
                group_by(cameraId) %>%
                summarise(seq = n()), by = "cameraId")
  results_filtrado2$seq[is.na(results_filtrado2$seq)] <- 0 # camera without encounters
  results_filtrado2$oper <- 30
  
  # Speed
    results_filtrado_detdistance <- results_filtrado %>%
    group_by(seqID) %>%
    summarise(
      speed     = mean(step, na.rm = TRUE),  
      across(-step, ~ first(.)),            
      .groups = "drop"
    )
  
  
  df_speed_detimp <- results_filtrado_detdistance 
  df_speed <- subset(df_speed_detimp, speed > 0 & speed < 1000)  
  
  if(nrow(df_speed) >= 20) {
    
    df_speed$speed <- df_speed$speed/step_sec #transform speed from m/20min to m/s
    
    df_speed <- subset(df_speed, speed > 0.02) 
    harmonic_mean_detimp <- length(df_speed$speed)/sum(1/df_speed$speed)
    harmonic_SE_detimp <- harmonic_mean_detimp^2 * sqrt(var(1/df_speed$speed)/length(df_speed$speed))
    DR_detimp <- harmonic_mean_detimp*60*60*12/1000
    DR_SE_detimp <- harmonic_SE_detimp*60*60*12/1000
    
    # radius
    results_filtrado_detdistance <- results_filtrado %>%
      group_by(seqID) %>%
      slice_head(n = 1) %>%
      ungroup()
    
    best_modRad <- ds(results_filtrado_detdistance$distancetoCT, transect = "point", key="hn", adjustment = NULL, truncation=w_rad)
    
    # Estimating effective detection radius and (SE)
    EfecRad <- EDRtransform(best_modRad)

    # Average values
    param <- list(DR = DR_detimp,
                  r = EfecRad$EDR/1000,
                  theta = fov_angle)
    
    # Standard error values
    paramse <- list(DR = DR_SE_detimp,
                    r = EfecRad$se.EDR/1000,
                    theta = 0)
    # density
    REMdensity <- bootTRD(results_filtrado2$seq, results_filtrado2$oper, param, paramse)
    
    # Store results
    resultados[[paste(current_type, current_n_cts, sep = "_")]] <- list(
      true_dens = true_density,
      dens = REMdensity[1],
      SE_dens = REMdensity[2],
      dens_type = "ImPerfDet", 
      n_cts = current_n_cts,
      type = current_type,
      sur_per_days = unique(results_filtrado2$oper),
      dayrange = DR_detimp,
      SE_dayrange = DR_SE_detimp,
      radius = EfecRad$EDR/1000,
      SE_radius = EfecRad$se.EDR/1000,
      angle = fov_angle,
      n_seq = nrow(df_speed),
      studyarea_km2 = super_studyarea_km2,
      n_cts_seq_ma0 = sum(results_filtrado2$seq > 0),
      porcen_max_seq = (max(results_filtrado2$seq)/sum(results_filtrado2$seq)) * 100,
      scenario = paste0("groupSize", groupSize, "_n_groups", n_groups, "_groupCohesion", groupCohesion, "_floater", floater) # ANADIR LOBOS SOLITARIOS
    )
    
  } else {
    # When no encounters (or few)
    resultados[[paste(current_type, current_n_cts, sep = "_")]] <- list(
      true_dens = true_density,
      dens = NA,
      SE_dens = NA,
      dens_type = "ImPerfDet",
      n_cts = current_n_cts,
      type = current_type,
      sur_per_days = unique(results_filtrado2$oper),
      dayrange = NA,
      SE_dayrange = NA,
      radius = NA,
      SE_radius = NA,
      angle = fov_angle,
      n_seq = nrow(df_speed),
      studyarea_km2 = super_studyarea_km2,
      n_cts_seq_ma0 = sum(results_filtrado2$seq > 0),
      porcen_max_seq = (max(results_filtrado2$seq)/sum(results_filtrado2$seq)) * 100,
      scenario = paste0("groupSize", groupSize, "_n_groups", n_groups, "_groupCohesion", groupCohesion, "_floater", floater)
    )
  }
  
}

# joining with previous results
resultados_df <- bind_rows(resultados, .id = "combinacion")

### 2 months ----
combinaciones_unicas <- results_2month %>%
  distinct(type, n_cts)

resultados <- list()

for(i in 1:nrow(combinaciones_unicas)) {
  current_type <- combinaciones_unicas$type[i]
  current_n_cts <- combinaciones_unicas$n_cts[i]
  
  results_filtrado <- results_2month %>%
    filter(type == current_type, n_cts == current_n_cts)
  
  cameras_total_filtrado <- cameras_total %>%
    filter(type == current_type, n_cts == current_n_cts)
  
  conteo <- results_filtrado %>%
    add_count(cameraId, name = "seq")
  
  results_filtrado2 <- cameras_total_filtrado %>%
    left_join(results_filtrado %>%
                group_by(cameraId) %>%
                summarise(seq = n()), by = "cameraId")
  results_filtrado2$seq[is.na(results_filtrado2$seq)] <- 0 # camera without detections
  results_filtrado2$oper <- 60
  
  # 3.0 speed
  results_filtrado_detdistance <- results_filtrado %>%
    group_by(seqID) %>%
    summarise(
      speed     = mean(step, na.rm = TRUE), 
      across(-step, ~ first(.)),            
      .groups = "drop"
    )
  
  
  df_speed_detimp <- results_filtrado_detdistance 
  
  df_speed <- subset(df_speed_detimp, speed > 0 & speed < 1000) 
  
  if(nrow(df_speed) >= 20) {
    
    df_speed$speed <- df_speed$speed/step_sec 
    df_speed <- subset(df_speed, speed > 0.02) 
    
    harmonic_mean_detimp <- length(df_speed$speed)/sum(1/df_speed$speed)
    harmonic_SE_detimp <- harmonic_mean_detimp^2 * sqrt(var(1/df_speed$speed)/length(df_speed$speed))
    DR_detimp <- harmonic_mean_detimp*60*60*12/1000
    DR_SE_detimp <- harmonic_SE_detimp*60*60*12/1000
    
    # radius
    results_filtrado_detdistance <- results_filtrado %>%
      group_by(seqID) %>%
      slice_head(n = 1) %>%
      ungroup()
    
    best_modRad <- ds(results_filtrado_detdistance$distancetoCT, transect = "point", key="hn", adjustment = NULL, truncation=w_rad)
    
    # Estimating effective detection radius and (SE)
    EfecRad <- EDRtransform(best_modRad)
    
    # Average values
    param <- list(DR = DR_detimp,
                  r = EfecRad$EDR/1000,
                  theta = fov_angle)
    
    # Standard error values
    paramse <- list(DR = DR_SE_detimp,
                    r = EfecRad$se.EDR/1000,
                    theta = 0)
    # density
    REMdensity <- bootTRD(results_filtrado2$seq, results_filtrado2$oper, param, paramse)
    
    resultados[[paste(current_type, current_n_cts, sep = "_")]] <- list(
      true_dens = true_density,
      dens = REMdensity[1],
      SE_dens = REMdensity[2],
      dens_type = "ImPerfDet", 
      n_cts = current_n_cts,
      type = current_type,
      sur_per_days = unique(results_filtrado2$oper),
      dayrange = DR_detimp,
      SE_dayrange = DR_SE_detimp,
      radius = EfecRad$EDR/1000,
      SE_radius = EfecRad$se.EDR/1000,
      angle = fov_angle,
      n_seq = nrow(df_speed),
      studyarea_km2 = super_studyarea_km2,
      n_cts_seq_ma0 = sum(results_filtrado2$seq > 0),
      porcen_max_seq = (max(results_filtrado2$seq)/sum(results_filtrado2$seq)) * 100,
      scenario = paste0("groupSize", groupSize, "_n_groups", n_groups, "_groupCohesion", groupCohesion, "_floater", floater) 
    )
  } else {
    resultados[[paste(current_type, current_n_cts, sep = "_")]] <- list(
      true_dens = true_density,
      dens = NA,
      SE_dens = NA,
      dens_type = "ImPerfDet",
      n_cts = current_n_cts,
      type = current_type,
      sur_per_days = unique(results_filtrado2$oper),
      dayrange = NA,
      SE_dayrange = NA,
      radius = NA,
      SE_radius = NA,
      angle = fov_angle,
      n_seq = nrow(df_speed),
      studyarea_km2 = super_studyarea_km2,
      n_cts_seq_ma0 = sum(results_filtrado2$seq > 0),
      porcen_max_seq = (max(results_filtrado2$seq)/sum(results_filtrado2$seq)) * 100,
      scenario = paste0("groupSize", groupSize, "_n_groups", n_groups, "_groupCohesion", groupCohesion, "_floater", floater)
    )
  }
  
}

resultados_df <- bind_rows(resultados_df, bind_rows(resultados, .id = "combinacion"))


### 3 months ----
combinaciones_unicas <- results_3month %>%
  distinct(type, n_cts)

resultados <- list()

for(i in 1:nrow(combinaciones_unicas)) {
  current_type <- combinaciones_unicas$type[i]
  current_n_cts <- combinaciones_unicas$n_cts[i]
  
  results_filtrado <- results_3month %>%
    filter(type == current_type, n_cts == current_n_cts)
  
  cameras_total_filtrado <- cameras_total %>%
    filter(type == current_type, n_cts == current_n_cts)
  
  conteo <- results_filtrado %>%
    add_count(cameraId, name = "seq")
  
  results_filtrado2 <- cameras_total_filtrado %>%
    left_join(results_filtrado %>%
                group_by(cameraId) %>%
                summarise(seq = n()), by = "cameraId")
  results_filtrado2$seq[is.na(results_filtrado2$seq)] <- 0
  results_filtrado2$oper <- 90
  
  # 3.0 speed
  results_filtrado_detdistance <- results_filtrado %>%
    group_by(seqID) %>%
    summarise(
      speed     = mean(step, na.rm = TRUE),  
      across(-step, ~ first(.)),            
      .groups = "drop"
    )
  
  df_speed_detimp <- results_filtrado_detdistance 
  
  df_speed <- subset(df_speed_detimp, speed > 0 & speed < 1000) 
  
  if(nrow(df_speed) >= 20) {
    
    df_speed$speed <- df_speed$speed/step_sec 
    
    df_speed <- subset(df_speed, speed > 0.02) 
    
    harmonic_mean_detimp <- length(df_speed$speed)/sum(1/df_speed$speed)
    harmonic_SE_detimp <- harmonic_mean_detimp^2 * sqrt(var(1/df_speed$speed)/length(df_speed$speed))
    DR_detimp <- harmonic_mean_detimp*60*60*12/1000
    DR_SE_detimp <- harmonic_SE_detimp*60*60*12/1000
    
    # radius
    results_filtrado_detdistance <- results_filtrado %>%
      group_by(seqID) %>%
      slice_head(n = 1) %>%
      ungroup()
    
    best_modRad <- ds(results_filtrado_detdistance$distancetoCT, transect = "point", key="hn", adjustment = NULL, truncation=w_rad)
    
    # Estimating effective detection radius and (SE)
    EfecRad <- EDRtransform(best_modRad)
    
    # Average values
    param <- list(DR = DR_detimp,
                  r = EfecRad$EDR/1000,
                  theta = fov_angle)
    
    # Standard error values
    paramse <- list(DR = DR_SE_detimp,
                    r = EfecRad$se.EDR/1000,
                    theta = 0)
    # density
    REMdensity <- bootTRD(results_filtrado2$seq, results_filtrado2$oper, param, paramse)
    
    # Almacenar resultados
    resultados[[paste(current_type, current_n_cts, sep = "_")]] <- list(
      true_dens = true_density,
      dens = REMdensity[1],
      SE_dens = REMdensity[2],
      dens_type = "ImPerfDet", 
      n_cts = current_n_cts,
      type = current_type,
      sur_per_days = unique(results_filtrado2$oper),
      dayrange = DR_detimp,
      SE_dayrange = DR_SE_detimp,
      radius = EfecRad$EDR/1000,
      SE_radius = EfecRad$se.EDR/1000,
      angle = fov_angle,
      n_seq = nrow(df_speed),
      studyarea_km2 = super_studyarea_km2,
      n_cts_seq_ma0 = sum(results_filtrado2$seq > 0),
      porcen_max_seq = (max(results_filtrado2$seq)/sum(results_filtrado2$seq)) * 100,
      scenario = paste0("groupSize", groupSize, "_n_groups", n_groups, "_groupCohesion", groupCohesion, "_floater", floater) 
    )
  } else {
    resultados[[paste(current_type, current_n_cts, sep = "_")]] <- list(
      true_dens = true_density,
      dens = NA,
      SE_dens = NA,
      dens_type = "ImPerfDet",
      n_cts = current_n_cts,
      type = current_type,
      sur_per_days = unique(results_filtrado2$oper),
      dayrange = NA,
      SE_dayrange = NA,
      radius = NA,
      SE_radius = NA,
      angle = fov_angle,
      n_seq = nrow(df_speed),
      studyarea_km2 = super_studyarea_km2,
      n_cts_seq_ma0 = sum(results_filtrado2$seq > 0),
      porcen_max_seq = (max(results_filtrado2$seq)/sum(results_filtrado2$seq)) * 100,
      scenario = paste0("groupSize", groupSize, "_n_groups", n_groups, "_groupCohesion", groupCohesion, "_floater", floater, "_ncts", current_n_cts, "_cts", current_type)
    )
  }
}

# joining results
resultados_df <- bind_rows(resultados_df, bind_rows(resultados, .id = "combinacion"))

resultados_df$groupSize <- groupSize
resultados_df$floater <- floater
resultados_df$groupCohesion <- groupCohesion
resultados_df$n_groups <- n_groups
resultados_df$HRsimulated <- HR_simulated

# Performance measurements: relative bias, coefficient of variation & 95%CI coverage 
resultados_df <- cbind(resultados_df, lnorm_confint(resultados_df$dens, resultados_df$SE_dens))
resultados_df$CV <- resultados_df$SE_dens/resultados_df$dens
resultados_df <- resultados_df %>%
  mutate(
    bias = dens - true_dens,
    relative_bias = (bias / true_dens) * 100,
    coverage = ifelse(
      (true_dens >= lcl) & 
        (true_dens <= ucl), 1, 0)
  )
