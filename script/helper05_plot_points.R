library("tidyverse")
library("magrittr")
library("rworldmap")
library("sp")
source("script/helper02_functions.R")

list.files("data/raw")

df <- read_csv("data/raw/gbif_occurrences.csv")

names(df) <- c("name","scientific_name","country","lon",
               "lat","key","record","publishing_org",
               "year","gbif_id","taxa","acronym")

sp <- unique(df$acronym)

scname <- unique(df$taxa)

scname <- gsub("<d7>","x", scname)


# get coarse resolution world from rworldmap
layer <- getMap() 

# simple cleaning steps 

# no NAs in lon
keep <- !is.na(df$lon)

# and no NAs in lat
keep <- !is.na(df$lat) & keep

# apply filter
df %<>% 
  filter(keep) %>% 
  dplyr::select(-name, -publishing_org)

# remove duplicated coordinates within species
df %<>% 
  group_by(acronym) %>% 
  distinct(lon,lat, .keep_all = TRUE) %>% 
  ungroup()

# remove coordinates without decimals
keep <- unlist(lapply(df$lon, decimalplaces)) != 0

keep <- unlist(lapply(df$lat, decimalplaces)) != 0 & keep

df %<>% 
  filter(keep)


# write plots
output <- "output/original_points/"

dir.create(output, showWarnings = FALSE, recursive = TRUE)

for(i in seq_along(sp)) {
  cat(i, "...")
  df %>% 
    filter(acronym == sp[i]) %>% 
    dplyr::select(lon,lat) %>%
    SpatialPoints(., proj4string = CRS(proj4string(layer))) ->
    coord
  

  png(filename = paste0(output, sp[i],".png"),
      width = 30,
      height = 25, 
      units = "cm",
      res = 100)
  plot(layer, col = "lightgrey", main = scname[i])
  plot(coord, col = "blue", cex = 1,
       bg = "Steelblue1", pch = 21, add = TRUE)
  dev.off()
}



