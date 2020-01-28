

library(rgdal)
library(raster)
library(rgeos)
library(maptools)


fenno <- readOGR("data/gadm", "fennoscandia")

proj <- proj4string(fenno)

ext <- floor(extent(fenno))

r <- raster(ext, res = 0.05)

fenno <- rasterize(fenno, r, field = 1)


inverse <- readOGR("data/buffer", "inverse")

inverse <- rasterize(inverse, r, field = 1)


m <- mask(fenno, inverse, inverse = TRUE)

p <- as.data.frame(m, xy = TRUE)

p <- p[!is.na(p[,3]),]

p <- p[,-3]

p <- SpatialPoints(p, proj4string = CRS(proj))

df <- as.data.frame(p@coords)

pdf <- SpatialPointsDataFrame(p, df)


writeOGR(pdf, 
         dsn = "data/buffer",
         layer = "coast_points",
         driver = "ESRI Shapefile",
         overwrite_layer = TRUE)





