
# Loading some functions developed by M Rowcliffe but not available in R packages
#source("https://raw.githubusercontent.com/MarcusRowcliffe/camtools/master/camtools.R")

# Modified to work with DR instead of  speed and activity
TRD <- function(P, T, param, strata=NULL, areas=NULL){
  if(length(P)!=length(T)) stop("P and T have unequal lengths")
  if(!("g" %in% names(param))) param <- c(param,g=1)
  if(!("p" %in% names(param))) param <- c(param,p=1)
  if(is.null(strata))
  {	res <- pi * param$g * sum(P) / (sum(T) * param$DR * param$r * (2+param$theta))
  names(res) <- NULL
  } else
  {	if(is.null(areas)) stop("areas are missing")
    if(length(strata)!=length(P)) stop("strata vector is a different length to P/T")
    if(sum(names(areas) %in% levels(strata)) != length(names(areas)) |
       sum(levels(strata) %in% names(areas)) != length(levels(strata)))
      stop("strata levels do not match areas names")
    nstrata <- length(areas)
    areas <- unlist(areas[order(names(areas))])
    locdens <- sapply(1:nstrata, STRD, P, T, param, strata)
    res <- sum(locdens * areas) / sum(areas)
  }
  res
} 
TRDsample <- function(i, P, T, param, strata=NULL, areas=NULL){
  if(is.null(strata)){
    x <- sample(1:length(T), replace=TRUE)
  } else{
    nstrata <- length(areas)
    x <- NULL
    for (i in 1:nstrata) x <- c(x, sample(which(strata==levels(strata)[i]), replace=TRUE))
    strata <- sort(strata)
  }
  TRD(P[x], T[x], param, strata, areas)
}
bootTRD <- function(P, T, param, paramSE, strata=NULL, areas=NULL, its=1000){
  BSdens <- sapply(1:its, TRDsample, P, T, param, strata, areas)
  BSse <- sd(BSdens)
  Dens <- TRD(P,T,param,strata,areas)
  prms <- length(param)
  addn <- rep(0,prms)
  addn[which(names(param)=="theta")] <- 2
  Es <- c(Dens,addn+unlist(param))
  SEs <- c(BSse,unlist(paramSE))
  SE <- Dens * sqrt(sum((SEs/Es)^2))
  cbind(Density=Dens, SE=SE)
  
}

lnorm_confint <- function(estimate, se, percent=95){
  if(length(estimate) != length(se)) 
    stop("estimate and se must have the same number of values")
  z <- qt((1 - percent/100) / 2, Inf, lower.tail = FALSE)
  w <- exp(z * sqrt(log(1 + (se/estimate)^2)))
  data.frame(lcl=estimate/w, ucl=estimate*w)
}


# Loading an version of moveHMM::simData to simulates groups accounting for 
# homeranges and group cohesion
simDataHR_group <- function(nbAnimals = 1, nbStates = 2, 
                            stepDist = c("gamma", "weibull", "lnorm", "exp"), 
                            angleDist = c("vm", "wrpcauchy", "none"), 
                            stepPar = NULL, anglePar = NULL, beta = NULL, 
                            covs = NULL, nbCovs = 0, zeroInflation = FALSE, 
                            obsPerAnimal = 500,   
                            springEffect = FALSE, 
                            k_spring = 1, 
                            homeRangeRadius = 100, 
                            states = FALSE,
                            groupSize = 1,
                            groupSpacing = 10,
                            groupCohesion = 1
) {
  
  if (nbAnimals < 1) stop("nbAnimals debe ser al menos 1.")
  if (nbStates < 1) stop("nbStates debe ser al menos 1.")
  if (is.null(stepPar)) stop("'stepPar' necesita especificarse.")
  if (groupSize < 1) stop("groupSize debe ser al menos 1.")
  if (groupCohesion < 0 || groupCohesion > 1) 
    stop("groupCohesion range should be [0,1].")
  
  stepDist  <- match.arg(stepDist)
  angleDist <- match.arg(angleDist)
  nbGroups  <- ceiling(nbAnimals / groupSize)
  
  p <- moveHMM:::parDef(stepDist, angleDist, nbStates, TRUE, zeroInflation)
  if (length(stepPar) != p$parSize[1] * nbStates ||
      length(anglePar) != p$parSize[2] * nbStates) {
    stop("Wrong number of parameters")
  }
  
  if (length(obsPerAnimal) == 1) {
    allNbObs <- rep(obsPerAnimal, nbGroups)
  } else if (length(obsPerAnimal) == nbGroups) {
    allNbObs <- obsPerAnimal
  } else if (length(obsPerAnimal) == 2) {
    allNbObs <- sample(obsPerAnimal[1]:obsPerAnimal[2], nbGroups, replace = TRUE)
  } else {
    stop("obsPerAnimal length should be 1, 2 or nbGroups")
  }
  
  data <- data.frame(ID = character(), step = numeric(),
                     angle = numeric(), x = numeric(),
                     y = numeric(), stringsAsFactors = FALSE)
  
  for (grp in seq_len(nbGroups)) {
    nbObs <- allNbObs[grp]
    Z <- sample(1:nbStates, nbObs, replace = TRUE)
    
    # Leader
    X_leader <- matrix(0, nrow = nbObs, ncol = 2)
    phi <- 0
    s <- numeric(nbObs)
    a <- numeric(nbObs)
    
    for (k in seq_len(nbObs-1)) {
      stepArgs <- c(1, stepPar[, Z[k]])
      angleArgs <- c(1, anglePar[, Z[k]])
      
      if (stepDist == "gamma") {
        shape <- stepArgs[2]^2 / stepArgs[3]^2
        scale <- stepArgs[3]^2 / stepArgs[2]
        stepArgs <- c(1, shape, 1/scale)
      }
      
      s[k] <- switch(stepDist,
                     gamma = rgamma(1, shape = stepArgs[2], scale = stepArgs[3]),
                     weibull = rweibull(1, shape = stepArgs[2], scale = stepArgs[3]),
                     lnorm = rlnorm(1, meanlog = stepArgs[2], sdlog = stepArgs[3]),
                     exp = rexp(1, rate = stepArgs[2]))
      
      if (angleDist != "none" && s[k] > 0) {
        a[k] <- do.call(paste0("r", angleDist),
                        list(n = 1, mean = angleArgs[2], k = angleArgs[3]))
        phi <- phi + a[k]
      }
      
      disp <- c(s[k] * cos(phi), s[k] * sin(phi))
      
      if (springEffect) {
        predicted_pos <- X_leader[k, ] + disp
        predicted_distance <- sqrt(sum(predicted_pos^2))
        if (predicted_distance > homeRangeRadius) {
          unit_vec <- predicted_pos / predicted_distance
          force_magnitude <- k_spring * (predicted_distance - homeRangeRadius)
          force <- -force_magnitude * unit_vec
        } else {
          force <- c(0, 0)
        }
        X_leader[k + 1, ] <- X_leader[k, ] + disp + force
      } else {
        X_leader[k + 1, ] <- X_leader[k, ] + disp
      }
    }
    a[1] <- NA
    
    data <- rbind(data, data.frame(
      ID = paste0("G", grp, "_A1"),
      step = s, angle = a,
      x = X_leader[, 1], y = X_leader[, 2],
      stringsAsFactors = FALSE
    ))
    
    # "followers"
    if (groupSize > 1) {
      # group cohesion
      cohesion_steps <- runif(nbObs) < groupCohesion
      
      for (i in 2:groupSize) {
        X_i <- matrix(0, nrow = nbObs, ncol = 2)
        phi_i <- 0
        s_i <- numeric(nbObs)
        a_i <- numeric(nbObs)
        
        # starting point (same as the leader)
        X_i[1, ] <- X_leader[1, ]
        
        for (k in seq_len(nbObs-1)) {
          if (cohesion_steps[k]) {
            # follow the leader - but not exactly the same trajectory-
            new_pos <- c(
              X_leader[k+1, 1] + rnorm(1, 0, groupSpacing),
              X_leader[k+1, 2] + rnorm(1, 0, groupSpacing)
            )
            
            dx <- new_pos[1] - X_i[k, 1]
            dy <- new_pos[2] - X_i[k, 2]
            
            # Real steps and angles
            s_i[k] <- sqrt(dx^2 + dy^2)
            angle_abs <- atan2(dy, dx)  
            
            if (k == 1) {
              a_i[k] <- 0  
            } else {
              turn_angle <- angle_abs - phi_i
              a_i[k] <- atan2(sin(turn_angle), cos(turn_angle))
            }
            
            X_i[k+1, ] <- new_pos
            phi_i <- angle_abs
            
          } else {
            # now following the leader
            z <- sample(1:nbStates, 1)
            stepArgs <- c(1, stepPar[, z])
            angleArgs <- c(1, anglePar[, z])
            
            if (stepDist == "gamma") {
              shape <- stepArgs[2]^2 / stepArgs[3]^2
              scale <- stepArgs[3]^2 / stepArgs[2]
              stepArgs <- c(1, shape, 1/scale)
            }
            
            s_i[k] <- switch(stepDist,
                             gamma = rgamma(1, shape = stepArgs[2], scale = stepArgs[3]),
                             weibull = rweibull(1, shape = stepArgs[2], scale = stepArgs[3]),
                             lnorm = rlnorm(1, meanlog = stepArgs[2], sdlog = stepArgs[3]),
                             exp = rexp(1, rate = stepArgs[2]))
            
            if (angleDist != "none" && s_i[k] > 0) {
              a_i[k] <- do.call(paste0("r", angleDist),
                                list(n = 1, mean = angleArgs[2], k = angleArgs[3]))
              phi_i <- phi_i + a_i[k]
            } else {
              a_i[k] <- 0
            }
            
            disp_i <- c(s_i[k] * cos(phi_i), s_i[k] * sin(phi_i))
            
            if (springEffect) {
              predicted_pos <- X_i[k, ] + disp_i
              predicted_distance <- sqrt(sum(predicted_pos^2))
              if (predicted_distance > homeRangeRadius) {
                unit_vec <- predicted_pos / predicted_distance
                force_magnitude <- k_spring * (predicted_distance - homeRangeRadius)
                force <- -force_magnitude * unit_vec
              } else {
                force <- c(0, 0)
              }
              X_i[k + 1, ] <- X_i[k, ] + disp_i + force
            } else {
              X_i[k + 1, ] <- X_i[k, ] + disp_i
            }
          }
        }
        a_i[1] <- NA
        
        data <- rbind(data, data.frame(
          ID = paste0("G", grp, "_A", i),
          step = s_i, angle = a_i,
          x = X_i[, 1], y = X_i[, 2],
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  return(moveHMM:::moveData(data))
}

# regular grids
crear_malla <- function(n, xmin, xmax, ymin, ymax, buffer = 0) {
  if (buffer < 0) stop("Buffer should be positive")
  if (buffer >= (xmax - xmin)/2 || buffer >= (ymax - ymin)/2) {
    stop("The buffer is too large for the defined space")
  }
  
  xmin_eff <- xmin + buffer
  xmax_eff <- xmax - buffer
  ymin_eff <- ymin + buffer
  ymax_eff <- ymax - buffer
  
  n_lado <- sqrt(n)
  
  if (n_lado != round(n_lado)) {
    stop("The number of points must be a perfect square for a regular grid")
  }
  
  x_seq <- seq(from = xmin_eff, to = xmax_eff, length.out = n_lado)
  y_seq <- seq(from = ymin_eff, to = ymax_eff, length.out = n_lado)
  
  malla <- expand.grid(x = x_seq, y = y_seq)
  
  malla$puntos <- n
  malla$densidad <- n / ((xmax_eff - xmin_eff) * (ymax_eff - ymin_eff))
  malla$id <- paste0("Pgrid", 1:n, "_", n)
  
  return(malla)
}

# random cameras
generar_cameras_aleatorios <- function(n, xmin, xmax, ymin, ymax, buffer = 0, 
                                       min_distancia = 0, max_intentos = 100) {
  if (buffer < 0) stop("Buffer should be positive")
  if (min_distancia < 0) stop("Minimum distance should be positive")
  if (buffer >= (xmax - xmin)/2 || buffer >= (ymax - ymin)/2) {
    stop("The buffer is too large for the defined space")
  }
  
  xmin_eff <- xmin + buffer
  xmax_eff <- xmax - buffer
  ymin_eff <- ymin + buffer
  ymax_eff <- ymax - buffer
  ancho_eff <- xmax_eff - xmin_eff
  alto_eff <- ymax_eff - ymin_eff
  
  cameras <- data.frame(x = numeric(0), y = numeric(0))
  intentos <- 0
  
  while (nrow(cameras) < n && intentos < max_intentos) {
    x_candidato <- runif(1, xmin_eff, xmax_eff)
    y_candidato <- runif(1, ymin_eff, ymax_eff)
    
    if (min_distancia > 0 && nrow(cameras) > 0) {
      distancias <- sqrt((cameras$x - x_candidato)^2 + (cameras$y - y_candidato)^2)
      if (min(distancias) < min_distancia) {
        intentos <- intentos + 1
        next
      }
    }
    
    cameras <- rbind(cameras, data.frame(x = x_candidato, y = y_candidato))
    intentos <- 0
  }
  
  cameras$n <- n
  cameras$buffer <- buffer
  cameras$min_distancia <- min_distancia
  cameras$id <- paste0("Prandom", 1:nrow(cameras), "_", n)
  
  if (nrow(cameras) < n) {
    warning("Only simulated ", nrow(cameras), " cameras from ", n, " cameras. ",
            "Increase max_intentos o reduce min_distancia.")
  }
  
  return(cameras)
}

# estimate effective detection zone - developed by Distance Sampling team - CREEM
EDRtransform <- function(dsobject, alpha=0.05) {
  if(class(dsobject) != "dsmodel") stop("First argument must be a dsmodel object")
  if(!dsobject$ddf$meta.data$point) stop("EDR can only be computed for point transect data")
  summary.ds.model <- summary(dsobject)
  p_a <- summary.ds.model$ds$average.p
  se.p_a <- summary.ds.model$ds$average.p.se
  cv.p_a <- se.p_a / p_a
  w <- summary.ds.model$ds$width
  edr <- sqrt(p_a * w^2)
  se.edr <- cv.p_a/2 * edr
  degfree <- summary.ds.model$ds$n - length(summary.ds.model$ddf$par)
  t.crit <- qt(1 - alpha/2, degfree)
  se.log.edr <- sqrt(log(1 + (cv.p_a/2)^2))
  c.mult <- exp(t.crit * se.log.edr)
  ci.edr <- c(edr / c.mult, edr * c.mult)
  return(list(EDR=edr, se.EDR=se.edr, ci.EDR=ci.edr))
}

