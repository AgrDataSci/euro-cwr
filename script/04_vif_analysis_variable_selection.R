#...................................................
# Select bioclimatic variables for modelling
# using VIF analysis

#...................................................
#...................................................
# Packages ####
library("tidyverse")
library("magrittr")
library("raster")
library("dismo")
library("BiodiversityR")
library("car")

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

# species acronyms
list.files("data")

# passport data
df <- "data/passport_data.csv"

df %<>%
  read_csv()

# .......................................
# .......................................
# Set background points ####
xy <-
  df %>%
  dplyr::select(lon, lat) %>%
  distinct(lon, lat, .keep_all = TRUE) %>%
  as.data.frame()

set.seed(123)
bg <- randomPoints(bio[[1]], 10000)

#...................................................
#...................................................
# Variable selection with VIF ####

vif <- ensemble.VIF(
  x = bio,
  a = xy,
  an = bg,
  VIF.max = 10,
  keep = NULL,
  layer.drops = "bio_07",
  factors = NULL,
  dummy.vars = NULL
)

# save outputs
output <- "processing/vif/"

dir.create(output,
           recursive = TRUE,
           showWarnings = FALSE)

save(vif, file = paste0(output, "vif_results.rda"))

# remove files not selected by vif analysis
out <- vif$var.drops

file.remove(paste0("data/bioclim/", out, ".tif"))


# ..................................
# ..................................
# future scenarios
gcm <- list.dirs("data/gcm")[-1]

for (i in seq_along(gcm)) {
  print(gcm[[i]])
  
  file.remove(paste0(gcm[[i]], "/", out, ".tif"))
  
}
