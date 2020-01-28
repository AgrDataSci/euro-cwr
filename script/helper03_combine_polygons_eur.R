library("sf")
library("rgdal")
library("rgeos")
library("raster")
library("maptools")
library("tidyverse")
library("magrittr")


# read file with contries iso codes
iso <- read_csv("data/country_iso.csv")

iso %<>% 
  filter(region == "Europe") %>% 
  select(alpha3) %>% 
  t()



# .............................
# .............................
# read world polygon and keep european countries
# except russia, which will be added later
list.files("data/gadm/other")

eur <- readOGR(dsn = "data/gadm/world", layer = "world_borders_adm0")

eur <- eur[eur$ADM0_A3 %in% iso, ]

# pick focal columns 
eur@data <- eur@data[,c("ADMIN","ADM0_A3", "GEOUNIT")]

# set a new bbox to remove countries areas outside Europe
e <- as(raster::extent(-25, 60, 34, 71), "SpatialPolygons")
proj4string(e) <- proj4string(eur)

# intersect
eur <- raster::intersect(eur, e)

plot(eur)


# # ...............................
# # ...............................
# # dissolve russia layer
# rus <- read_sf(dsn = "data/gadm/other", layer = "gadm36_RUS_2")
# 
# rus %>% 
#   st_set_geometry(NULL) %>% 
#   glimpse()
# 
# rus$area <- st_area(rus)
# 
# rus %>% 
#   summarise(area = sum(area)) ->
#   rus_eur
# 
# rus <- as(rus_eur, "Spatial")
# 
# rus_df <- eur@data[1,]
# 
# rus_df[1,] <- c("Russia","RUS","Russia")
# 
# rus@data <- rus_df
# 
# plot(rus)
# 
# # combine the datasets into a single polygon
# gadm <- list(eur,
#              rus)
# 
# gadm <- do.call("bind", gadm)
# 
# plot(gadm)

writeOGR(eur, 
         dsn = "data/gadm/europe",
         layer = "europe_adm0",
         driver = "ESRI Shapefile",
         overwrite_layer = TRUE)


# ................................
# ................................
# dissolve eur polygon 

# now dissolve
eur@data$eur <- 1

eur <- aggregate(eur, c("eur"))

plot(eur)

writeOGR(eur, 
         dsn = "data/gadm/europe",
         layer = "europe",
         driver = "ESRI Shapefile",
         overwrite_layer = TRUE)

# ......................................
# ......................................
# create buffer area
buff <- SpatialPolygons(eur@polygons, 
                        proj4string = CRS(proj4string(eur)) )

buff <- gBuffer(buff, width = 0.12)

# buffer into polygons dataframe
ids <- sapply(slot(buff, "polygons"), function(x) slot(x, "ID"))
df <- data.frame(region="Europe", 
                 adm = "buffer_010arcmin", 
                 row.names=ids)

# Create a spatial polygon data frame (includes shp attributes)
spdf <- SpatialPolygonsDataFrame(buff, df)

plot(spdf)

writeOGR(spdf, 
         dsn = "data/gadm/europe",
         layer = "europe_buffer_010arcmin",
         driver = "ESRI Shapefile",
         overwrite_layer = TRUE)

# .........................................
# ........................................
# create points around coast

buff <- SpatialPolygons(eur@polygons, 
                        proj4string = CRS(proj4string(eur)))

buff <- gBuffer(buff, width = -0.15)

e <- raster::intersect(buff, eur)

# set eur as raster
ext <- floor(extent(eur))
proj <- proj4string(buff)

r <- raster(ext, res = 0.05)

eur <- rasterize(eur, r, field = 1)

# and the inverse buffer
e <- rasterize(e, r, field = 1)

# now mask eur using the inverse raster
m <- mask(eur, e, inverse = TRUE)

p <- as.data.frame(m, xy = TRUE)

p <- p[!is.na(p[,3]),]

# remove russia east
p <- p[p[[1]] <= xmax(ext)-1, ]

p <- p[,-3]

p <- SpatialPoints(p, proj4string = CRS(proj))

df <- as.data.frame(p@coords)

pdf <- SpatialPointsDataFrame(p, df)

plot(pdf)

writeOGR(pdf, 
         dsn = "data/gadm/europe",
         layer = "europe_coast_points",
         driver = "ESRI Shapefile",
         overwrite_layer = TRUE)



