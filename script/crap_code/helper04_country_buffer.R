library("rgeos")
library("sp")
library("rgdal")
library("raster")


# Country borders from GADM
country <- c("NO","SE","FI","DK","RU")
gadm <- lapply(country, function(x) {
  getData("GADM", 
          country = x,
          path = "data/gadm/",
          level = 0)
})

names(gadm) <- country

# create a clipping polygon for RU
ext <- as(extent(5, 42, 59, 71),
          "SpatialPolygons")
proj <- proj4string(gadm[[1]])

proj4string(ext) <- CRS(proj)

# clip the map
rus <- gadm[["RU"]]
rus <- gIntersection(rus, ext, byid = TRUE)

gadm[["RU"]] <- rus

rm(rus)

# Create a buffer around country borders
B <- list()

for(i in seq_along(country)) {
  print(country[i])
  
  X <- gadm[[i]]
  
  X <- SpatialPolygons(X@polygons, proj4string = CRS(proj) )
  
  b <- gBuffer(X, width = 0.25)
  
  e <- gDifference(b, X, byid = TRUE)
  
  B[[i]] <- e
  
}

# Write shapefiles
output <- "data/buffer"
dir.create(output,
           showWarnings = FALSE,
           recursive = TRUE)

for(i in seq_along(country)) {
  print(country[i])
  
  e <- B[[i]]
  
  ids <- sapply(slot(e, "polygons"), function(x) slot(x, "ID"))
  df <- data.frame(rep(0, length(ids)), row.names=ids)
  
  # Create a spatial polygon data frame (includes shp attributes)
  spdf <- SpatialPolygonsDataFrame(e, df)
  
  writeOGR(spdf, 
           dsn = output,
           layer = country[i],
           driver = "ESRI Shapefile",
           overwrite_layer = TRUE)
  
}


