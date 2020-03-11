# Process layers from species distribution models

# ...............................................
# ...............................................
# Packages
library("raster")
library("sp")

# ...............................................
# ...............................................
# Parameters 
# define projection
myproj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# GCM models
gcm <- c("cn", "mp", "no")

# ...............................................
# ...............................................
# get species names from processed models
sp <- list.dirs("processing/enm")[-1]
sp <- strsplit(sp, "/")
sp <- suppressWarnings(do.call("rbind", sp)[,3])
sp <- unique(sp)

n <- max(seq_along(sp))

# each species layer has its own bbox based on the max hull for the 
# presence points used take a raster that includes all the regional
# area used here to create a new layer with the same bbox for all species
# use one of the bioclim layers and set all values as zero
eur <- raster("data/bioclim/bio_01.tif")
eur[eur[] != 0] <- 0
  
# ...............................................
# ...............................................
# Run over current presence 
sp_r <- list()

pb <- txtProgressBar(min = 0, max = n, initial = 0)
for(i in seq_along(sp)) {
  
  path_i <- paste0("processing/enm/", 
                   sp[i], 
                   "/ensembles/presence/")
  
  r_i <- stack(list.files(path_i,
                          pattern = "current.grd",
                          full.names = TRUE))

  # presence is defined as 1,
  # set all values different tha 1 as NA
  r_i[r_i[] != 1 ] <- NA

  crs(r_i) <- myproj
  
  # reconstruct layer using the regional layer as baseline
  r_i <- mosaic(eur, r_i, fun = sum)
  # set the regional layer as a mask
  r_i <- mask(r_i, eur)

  sp_r[[sp[i]]] <- r_i
  
  setTxtProgressBar(pb, i)
  
}

# put all layers together as a raster stack object
sp_r <- stack(sp_r)

# sum all layers from the stack
x <- calc(sp_r, fun = sum)

plot(x)


# ...............................................
# ...............................................
# Run over future presence 


