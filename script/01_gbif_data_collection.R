# Retrieve GBIF occurrence data for distribution analysis 
# updated by K. de Sousa 
# Inland Norway University
# Adapted from M. van Zonneveld
# World Vegetable Center
#..........................................
#..........................................
# Packages ####
library("rgbif")
library("tidyverse")
library("magrittr")

#..........................................
#..........................................
# Data ####

df <- "data/species_names.csv"
df %<>% 
  read_csv()

taxa <- df$taxa

# ....................................
# ....................................
# Data collection ####

# use that taxonconceptkey to get occurrence data for spp and its synonyms
# fields that will be retrieved by rgbif::occ_search
fields <- c("name","scientificName",
           "countryCode","decimalLongitude",
           "decimalLatitude","key",
           "basisOfRecord","publishingOrg",
           "year","gbifID")

gbif <- NULL

# run over species names
for (i in seq_along(taxa)){
  cat(i, "out of", length(taxa), taxa[i],"\n")
  
  # check for all scientific names that include the species names  
  name <- name_suggest(q = taxa[i])
  
  # remove subspecies and varieties
  name <- name[name$rank != "SUBSPECIES", ]
  
  g <- NULL
  
  if (nrow(name) == 0) next 
  
  cat("fetching GBIF data...\n")
  
  for (j in seq_len(nrow(name))) {  
    
    o <- occ_search(taxonKey = name[[j,1]],
                    hasCoordinate = TRUE, 
                    fields = fields)
    
    if (length(o$data) == 0) next
    
    o$data$taxa <- taxa[i]
    
    o$data$acronym <- df$acronym[i]
    
    # sometimes it dont come with all desired field 
    # create an empty data.frame and match wich fields
    # are provided by the sourch
    p <- as.data.frame(matrix(NA, nrow=nrow(o$data), ncol=12))
    names(p) <- c(fields,"taxa","acronym")
    index <- match(names(o$data), names(p))
    
    p[,index] <- o$data
    
    g <- rbind(g, p)
  
  } 
    
  gbif <- rbind(gbif, g)

}

# keep these occurrences as raw data
write.csv(gbif, "data/raw/gbif_occurrences.csv", 
          row.names = FALSE)

