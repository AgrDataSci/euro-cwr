

library("raster")
library("sp")

myproj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

sp <- list.dirs("processing/enm")[-1]
sp <- strsplit(sp, "/")
sp <- do.call("rbind", sp)[,3]
sp <- unique(sp)

eur <- raster("data/bioclim/bio_01.tif")

eur[eur[] != 0] <- 0
  
sp_r <- list()

for(i in seq_along(sp)) {
  
  path_i <- paste0("processing/enm/", 
                   sp[i], 
                   "/ensembles/presence/")
  
  r_i <- stack(list.files(path_i,
                          pattern = "current.grd",
                          full.names = TRUE))

  r_i[r_i[] != 1 ] <- NA

  crs(r_i) <- myproj
  
  r_i <- mosaic(eur, r_i, fun = sum)
  
  r_i <- mask(r_i, eur)

  sp_r[[sp[i]]] <- r_i
  
}


sp_r <- stack(sp_r)


x <- calc(sp_r, fun = sum)

plot(x)
