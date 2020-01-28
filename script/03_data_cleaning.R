# Clean occurrence data for distribution analysis 
# updated by K. de Sousa 
# Inland Norway University
#..........................................
#..........................................
# Packages ####
library("tidyverse")
library("magrittr")
library("maptools")
library("raster")
library("rgdal")
library("dismo")
library("alphahull")
library("rgeos")
library("sf")


# extra function 
source("script/helper02_functions.R")

# ........................................
# ........................................
# Data ####

# gadm adm0 europe
list.files("data/gadm/europe", pattern = ".shp$")

gadm <- readOGR(dsn = "data/gadm/europe", layer = "europe")

# europe 0.10 arc-min buffer area 
b <- readOGR("data/gadm/europe", "europe_buffer_010arcmin")

# europe coast points
coast <- readOGR("data/gadm/europe", "europe_coast_points")

# europe as a whole
eur <- readOGR("data/gadm/europe", "europe")

# list of target countries
country <- "data/country_iso.csv"

country %<>% 
  read_csv() %>% 
  filter(region == "Europe") %>% 
  dplyr::select(alpha2) %>% 
  t()


# define projection
proj <- proj4string(gadm)

#........................................
#........................................
# Read occurrence data ####

# Read genesys data
gen <- "data/raw/genesys_occurrences.csv"
gen %<>% 
  read_csv()

names(gen) <- c("accession_number","acquisition_date","lat","lon","collection_site",
                "country_origin","id","in_svalbard","institute","genus",
                "species","subtaxa","institute_name","country_name","country")

gen %<>% 
  mutate(scientific_name = paste(genus, species),
         source = "genesys") %>% 
  dplyr::select(scientific_name, lon, lat, country, source)


# Read GBIF data
gbif <- "data/raw/gbif_occurrences.csv"
gbif %<>% 
  read_csv()

# rename the dataframe
names(gbif) <- c("name","scientific_name","country","lon",
               "lat","key","record","publishing_org",
               "year","gbif_id","taxa","acronym")

# remove entries older than 1950
keep <- gbif$year >= 1950 & !is.na(gbif$year)

gbif %<>% 
  filter(keep) 

gbif %<>% 
  mutate(source = "gbif") %>% 
  dplyr::select(scientific_name, lon, lat, country, source, acronym)


# Check occurrences in both datasets
both <- "data/raw/in_both_databases.csv"
both %<>% 
  read_csv() %>% 
  mutate(scientific_name = paste(genus, species)) %>% 
  dplyr::select(scientific_name, acronym)


gen <- right_join(gen, both, by = "scientific_name")


df <- bind_rows(gen, gbif)

#........................................
#........................................
# Data cleaning 1 ####

# no NAs in lon
keep <- !is.na(df$lon)

# and no NAs in lat
keep <- !is.na(df$lat) & keep

# apply filter
df %<>% 
  filter(keep)

# remove duplicated coordinates within species
df %<>% 
  group_by(acronym, source) %>% 
  distinct(lon,lat, .keep_all = TRUE) %>% 
  ungroup()

# remove coordinates without decimals
keep <- unlist(lapply(df$lon, .decimalplaces)) != 0

keep <- unlist(lapply(df$lat, .decimalplaces)) != 0 & keep

df %<>% 
  filter(keep)

#........................................
#........................................
# Spatial cleaning 1 ####

# remove points outside country borders
df %>% 
  dplyr::select(lon,lat) %>% 
  SpatialPoints(., proj4string = CRS(proj)) ->
  coord

keep <- over(coord, b)

keep <- !is.na(keep[[1]])

df %<>% 
  filter(keep)

rm(coord)

#........................................
#........................................
# Spatial cleaning 2 ####

# revise points in the sea and correct those placed within 10 arc-min
# from the coastal border
df %>%
  dplyr::select(lon, lat) %>%
  SpatialPoints(., proj4string = CRS(proj)) ->
  coord


# keep sea area in a buffer of 10 arcmin
sea <- gDifference(b, eur, byid = TRUE)
ids <- sapply(slot(sea, "polygons"), function(x) slot(x, "ID"))
ids <- data.frame(sea = "sea", row.names = ids)
sea <- SpatialPolygonsDataFrame(sea, ids)

# identify points in the sea
sea <- over(coord, sea)
sea <- !is.na(sea[[1]])

# subset coordinates in the sea
sea_coord <- df[sea, c("lon","lat","acronym")]

sea_coord <- split(sea_coord, rownames(sea_coord))

# transform dataframes into spatial points then set the original projection
sea_buffer <- lapply(sea_coord, function(x) {
  # df into spatial points object
  x <- SpatialPoints(x[, c("lon","lat")], proj4string = CRS(proj))
  # set the projection
  x <- gBuffer(x, byid = TRUE, width = 0.5)
  # return the result
  x
})

# run over coordinates and find the nearest point in the land
nearest <- list()
pb <- txtProgressBar(min = 0, max = length(sea_buffer), initial = 1) 
for(i in seq_along(sea_buffer)) {
  
  x <- over(coast, sea_buffer[[i]])
  x <- coast@coords[!is.na(x),]

  if (nrow(x) > 0) {
    x <- SpatialPoints(x, proj4string = CRS(proj))
    
    # the original point
    y <- sea_coord[[i]][,c("lon","lat")]
    y <- SpatialPoints(y, proj4string = CRS(proj))
    
    # find the nearest point
    d <- pointDistance(x, y, lonlat = TRUE)
    d <- which.min(d)
    d <- x@coords[d,]
    names(d) <- c("lon","lat")
    
    nearest[[i]] <- d 
    
  } else {
    nearest[[i]] <- data.frame(lon = -999, lat = -999)
  }
  
  setTxtProgressBar(pb,i)
}

rm(x,y,d,i)

nearest <- do.call("rbind", nearest)


# replace sea points by their nearest points in the land
df[sea,"lon"] <- nearest[,1]
df[sea,"lat"] <- nearest[,2]

# remove possible -999 values
keep <- df$lon != -999

df <- df[keep, ]

# ....................................
# ....................................
# remove points within the same grid cell

# place genesys data in a different dataset
gen <- df[df$source == "genesys", ]

df <- df[df$source == "gbif", ]

bio <- raster("data/bioclim/eto.tif")

acronym <- unique(df$acronym)

df_coords <- data.frame()

pb <- txtProgressBar(min = 0, max = length(acronym), initial = 1) 
for (i in seq_along(acronym)) {
  # sampling data
  df_i <- df[df$acronym == acronym[i], ]
  
  coord <- df_i[, c("lon","lat")]
  
  largedist <- coord %>%
    raster::pointDistance(., longlat = FALSE) %>%
    max(., na.rm = TRUE)
  
  # make a convex hull and remove duplicated coordinates
  # in the same grid-cell
  hull <- convHull(coord, lonlat = TRUE)
  # extent convex hull
  ext_hull <- gBuffer(hull@polygons, width = 0.1 * largedist)
  crs(ext_hull) <- proj4string(eur)
  
  # define raster
  r <- raster(ext_hull)
  # set the resolution of the cells to 
  # 30 arc-sec
  res(r) <- res(bio)
  
  coord %<>%
    as.matrix() %>% 
    as.data.frame() %>%
    dismo::gridSample(., r, n=1)
  
  
  coord$acronym <- acronym[i]
  
  coord$scientific_name <- df_i$scientific_name[1]
  
  coord$source <- df_i$source[1]
  
  df_coords <- rbind(df_coords, coord)
  
  setTxtProgressBar(pb,i)
}


# put genesys data back to the main dataset
gen %<>% 
  dplyr::select(-country)

df_coords %<>% 
  bind_rows(df_coords, gen) %>% 
  as_tibble()

# ......................................
# ......................................
# write outputs

# occurrence in both datasets
df_coords %>% 
  group_by(acronym) %>% 
  distinct(source, .keep_all = TRUE) %>% 
  count() %>% 
  mutate(keep = n == 2) %>% 
  filter(keep) -> 
  keep


keep <- df_coords$acronym %in% keep$acronym

df_coords %<>% 
  filter(keep)


# count number of points per species and add to species names
df_coords %>% 
  group_by(acronym) %>%  
  count(acronym) %>% 
  mutate(keep = n > 29) %>% 
  filter(keep) ->
  keep 

keep <- df_coords$acronym %in% keep$acronym

df_coords %<>% 
  filter(keep)



df_coords %>% 
  group_by(acronym, source) %>% 
  count() -> x


# write file with cleaned points
write.csv(df_coords, "data/passport_data.csv", row.names = FALSE)


df_coords %>%
  dplyr::select(lon,lat) %>%
  SpatialPoints(., proj4string = CRS(proj)) ->
  coord

plot(eur, col = "lightgrey")
plot(coord, col = "blue", cex = 1,
     bg = "Steelblue1", pch = 21, add = TRUE)


# # Country borders from GADM
# country <- c("NO","SE","FI","DK")
# gadm <- lapply(country, function(x) {
#   getData("GADM", 
#           country = x,
#           path = "data/gadm/",
#           level = 0)
# })
# #gadm <- do.call("bind", gadm)
# 
# # add Russia
# rus <- getData("GADM", 
#                country = "RU",
#                path = "data/gadm/",
#                level = 0)
# 
# # create a clipping polygon
# ext <- as(extent(5, 42, 59, 71),
#           "SpatialPolygons")
# 
# proj <- proj4string(gadm[[1]])
# 
# proj4string(ext) <- CRS(proj)
# 
# # clip the map
# rus <- gIntersection(rus, ext, byid = TRUE)
# 
# # add dataframe again
# ids <- sapply(slot(rus, "polygons"), function(x) slot(x, "ID"))
# d <- data.frame(GID_0 = rep("RU", length(ids)), 
#                 NAME_0 = rep("Russia", length(ids)),
#                 row.names = ids)
# 
# # Create a spatial polygon data frame (includes shp attributes)
# rus <- SpatialPolygonsDataFrame(rus, d)
# 
# gadm[["RU"]] <- rus
# 
# #gadm <- do.call("bind", gadm)
# 
# rm(rus, d, ids)
# 
# country <- c(country, "RU")