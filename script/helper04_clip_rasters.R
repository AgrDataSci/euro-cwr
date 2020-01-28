# clip bioclim rasters to fit in 
# the target research area
# source data is available in an external hard drive or 
# http://worldclim.org/version1
library("raster")
library("rgdal")

source("script/helper02_functions.R")

gadm <- readOGR("data/gadm/europe", "europe")

ext <- floor(extent(bbox(gadm)))

r <- raster(ext, res = 0.04166666)

gadm <- rasterize(gadm, r, field = 1)

# ...........................................
# ...........................................
# worldclim 1.0

# read rasters from an external hard drive
r <- stack(list.files("data/bioclim", 
                      pattern = ".bil",
                      full.names = TRUE))
# crop it
r <- crop(r, ext)

# mask it
r <- mask(r, gadm)

# check the first layer
plot(r[[1]])

# rename rasters
names(r) <- gsub("bio","bio_",names(r))
names(r) <- nr

# save it on the project folder
output <- "data/bioclim/"

dir.create(output,
           showWarnings = FALSE,
           recursive = TRUE)

r <- stack(r)

names_r <- paste0(output, names(r))


file.remove(list.files("data/bioclim/",
                       pattern = ".bil|.hdr", 
                       full.names = TRUE))


writeRaster(r, 
            filename = names_r, 
            format = "GTiff",
            bylayer = TRUE,
            overwrite = TRUE)



# ..................................
# ..................................
# future scenarios
gcm <- list.dirs("data/gcm")[-1]

for (i in seq_along(gcm)) {
  print(gcm[[i]])
  # read rasters from an external hard drive
  r <- stack(list.files(gcm[[i]], 
                        pattern = ".tif",
                        full.names = TRUE))
  
  # crop it
  r <- crop(r, ext)
  
  # mask it
  r <- mask(r, gadm)
  
  model <- substr(names(r), 1, 2)[1]
  
  bio <- gsub(paste0(model,"45bi7"),"", names(r))
  
  bio[nchar(bio) > 2] <- gsub("01","1", bio[nchar(bio) > 2])
  
  bio <- paste0("bio_", bio)
  
  names(r) <- bio
  
  names_r <- paste0(gcm[[i]], "/", names(r))
  
  file.remove(list.files(gcm[[i]], full.names = TRUE))
  
  writeRaster(r, 
              filename = names_r, 
              format = "GTiff",
              bylayer = TRUE,
              overwrite = TRUE)
  

}





# ..................................
# ..................................
# evapotranspiration


r <- raster("E:/rasters/et0_yr/et0_yr.tif")

r <- crop(r, ext)

# reduce resolution
r <- aggregate(r, fact = 5)

r <- mask(r, gadm)


plot(r)

writeRaster(r,
            filename = paste0(output, "eto.tif"),
            format = "GTiff",
            overwrite = TRUE)