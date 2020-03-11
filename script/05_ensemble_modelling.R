#...................................................
#...................................................
# Model species distribution using bioclimatic variables
# 
# Kaue de Sousa
# Inland Norway University
# 
# ...................................................
# ...................................................
# Packages ####
library("data.table")
library("rgdal")
library("raster")
library("dismo")
library("rgeos")
library("gbm")
library("glmnet")
library("BiodiversityR")
library("PresenceAbsence")

#...................................................
#...................................................
# Data ####

# bioclimatic variables
bio <- list.files("data/bioclim",
                  pattern = ".tif$",
                  full.names = TRUE)

bio <- stack(bio)

# define projection and extension
myproj <- proj4string(bio)
myext  <- extent(bio) 
myres  <- res(bio)

# gcm models
gcm <- list.dirs("data/gcm")[-1]
gcm <- do.call("rbind", strsplit(gcm, "/"))[,3]


# passport data
df <- fread("data/passport_data.csv")

sp <- sort(unique(df$acronym))

#...................................................
#...................................................
# Run ensemble modelling ####
# the BiodiversityR saves outputs on the current working directory
# get the parent wd to return here if needed
parentwd <- getwd()
  
# check for species already processed 
filepattern <- paste(c(gcm, "current"), collapse = "|")
nfiles <- length(c(gcm, "current")) * 2

for (i in seq_along(sp)) {
  
  pres <- paste0("processing/enm/", sp[i], "/ensembles/presence/")
  pres <- grepl(filepattern, list.files(pres))
  
  suit <- paste0("processing/enm/", sp[i], "/ensembles/suitability/")
  suit <- grepl(filepattern, list.files(suit))
  
  # if folder doesnt contain these files then it should be removed
  if (sum(suit) != nfiles | sum(pres) != nfiles) {
    
    unlink(paste0("processing/enm/", sp[i]), recursive = TRUE)
    
  }
  
}

done <- list.dirs("processing/enm")[-1]
done <- strsplit(done, "/")
done <- suppressWarnings(do.call("rbind", done)[,3])
done <- unique(done)


# filter and run the model for the remnant 
sp <- sp[!sp %in% done] 

for (i in seq_along(sp)) {
  
  message("\n Ensemble modelling for ", sp[i], "\n Time: ", date(), "\n")
  
  # create a dir for the species and work in that dir
  output <- paste0("processing/enm/", sp[i], "/")
  dir.create(output, showWarnings = FALSE, recursive = TRUE)
  
  setwd(output)
  
  # sampling data
  coord <- df[df$acronym == sp[i], ]
  
  coord <- coord[, c("lon","lat")]
  
  coord <- as.data.frame(coord)
  
  # create a convexHull to limit the model to the 
  # areas where the species is actually present
  # calculate largest distance
  largedist <- pointDistance(coord, longlat = FALSE)
  largedist <- max(largedist, na.rm = TRUE)

  # make a convex hull 
  hull <- convHull(coord, lonlat = TRUE)
  
  # extent convex hull
  hull <- gBuffer(hull@polygons, 
                  width = 0.1 * largedist)
  
  crs(hull) <- myproj
  
  # convert into a raster
  r <- raster(hull, res = myres, ext = myext)

  hull <- rasterize(hull, r, field = 1,
                    background = NA)
  
  hull[hull == 0] <- NA
  
  # crop bioclim layers to fit this extention
  bio_i <- mask(bio, hull)
  
  bio_i <- stack(bio_i)

  # create background points over the area
  nbg <- nrow(coord) * 3
  if (nbg > 1000) nbg <- 1000
  set.seed(123)
  bg <- randomPoints(bio_i, n = nbg, p = coord)
  bg <- as.data.frame(bg)
  names(bg) <- c("lon", "lat")
  
  # Run ensemble modelling
  # step 1: model calibration
  # here the function tests for the best algorithms
  # since the algorithms were previously selected,
  # a 3-fold cross-validation is performed to make sure that all 
  # pass the output.weights threshold
  message("\n Step 1: Calibrating algorithms \n", "Time: ", date(), "\n")
  
  enm_step1 <- ensemble.calibrate.weights(x = bio_i, 
                                          p = coord, 
                                          a = bg,
                                          k = 3,
                                          layer.drops = NULL,
                                          SINK = TRUE, 
                                          species.name = sp[[i]],
                                          GBM = 1,
                                          GAM = 1,
                                          GLM = 1,
                                          SVM = 1,
                                          BIOCLIM = 1,
                                          MAHAL = 0,
                                          RPART = 0, 
                                          GBMSTEP = 0,
                                          MAXENT = 0, 
                                          NNET = 0, 
                                          RF = 0, 
                                          EARTH = 0,
                                          GLMSTEP = 0, 
                                          GAMSTEP = 0, 
                                          MGCV = 0, 
                                          MGCVFIX = 0, 
                                          CF = 0, 
                                          FDA = 0,
                                          SVME = 0,
                                          DOMAIN = 0, 
                                          ENSEMBLE.tune = TRUE, 
                                          PROBIT = TRUE,
                                          # see Liu et al (2013) doi:10.1111/jbi.12058
                                          threshold.method = "MaxSens+Spec", 
                                          threshold.sensitivity = 0.9,
                                          threshold.PresenceAbsence = TRUE, 
                                          ENSEMBLE.min = 0.7,
                                          Yweights = "BIOMOD")
  
  # step 2: create models that will be used for the raster predictions
  # models with output.weights <0.05 are excluded
  output_weights <- enm_step1$output.weights
  output_weights[output_weights < 0.05] <- 0
  
  message("Step 2: Model species distribution with selected ENM algorithms \n")
  
  enm_step2 <- ensemble.calibrate.models(x = bio_i, 
                                         p = coord, 
                                         a = bg,
                                         k = 10, 
                                         layer.drops = NULL,
                                         SINK = TRUE, 
                                         species.name = sp[[i]],
                                         models.keep = TRUE,
                                         input.weights = output_weights,
                                         # see Liu et al (2013) doi:10.1111/jbi.12058
                                         threshold.method = "MaxSens+Spec", 
                                         threshold.sensitivity = 0.9,
                                         threshold.PresenceAbsence = TRUE, 
                                         ENSEMBLE.tune = FALSE, 
                                         ENSEMBLE.min = 0.7,
                                         PROBIT = TRUE,
                                         Yweights = "BIOMOD", 
                                         models.save = FALSE)
  
  
  # save AUCs
  auc <- data.frame(auc = enm_step2$AUC.testing)
  write.csv(auc, file = "outputs/auc_testing.csv")
  
  message("Step 3.1: Generate map of current distribution \n")
  #step3: use previously calibrated models to construct consensus layers
  ensemble_current <- ensemble.raster(xn = bio_i,
                                      models.list = enm_step2$models,
                                      input.weights = output_weights,
                                      thresholds = enm_step2$models$thresholds,
                                      SINK = TRUE,
                                      RASTER.species.name = sp[[i]], 
                                      RASTER.stack.name = "current")
  
  # project suitability under climate change
  for (k in seq_along(gcm)){
    message("Step 3.2: Predict future distribution, GCM ", toupper(gcm[[k]]), "\n")

    #load GCM layers
    gcmfiles <- list.files(paste0(parentwd, "/data/gcm/", gcm[[k]], "/"),
                           pattern = ".tif$",
                           full.names = TRUE)

    gcm_model <- stack(gcmfiles)

    crs(gcm_model) <- myproj

    gcm_model <- mask(gcm_model, hull)

    gcm_model <- stack(gcm_model)

    ensemble_gcm_model <- ensemble.raster(xn = gcm_model,
                                          models.list = enm_step2$models,
                                          input.weights = output_weights,
                                          thresholds = enm_step2$models$thresholds,
                                          SINK = TRUE,
                                          RASTER.species.name = sp[[i]],
                                          RASTER.stack.name = gcm[[k]])

  }
  
  # remove working files created in the third step
  unlink("models", recursive = TRUE)
  unlink("ensembles/count", recursive = TRUE)
  file.remove(list.files(pattern = "working", full.names = TRUE))
  
  
  # return to parent wd
  setwd(parentwd)
  
  
  
}

message("Done at ", Sys.time())
