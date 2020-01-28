# Get Genesys occurrence data for distribution analysis 
# updated by K. de Sousa 
# Inland Norway University
#..........................................
#..........................................
# Packages ####
library("genesysr")
library("data.table")

#..........................................
#..........................................
# Data ####

df <- fread("data/species_names.csv")

taxa <- df[ ,c("genus","species", "acronym")]

#..........................................
#..........................................
# Set Genesys query ####

setup_sandbox()

user_login()

taxa <- split(taxa, 1:nrow(taxa))

gen_df <- list()

for (i in seq_along(taxa)) {
  
  cat(i, " ",  as.vector(unlist(taxa[[i]])), "\n")
  
  genus <- as.character(as.matrix(taxa[[i]][,"genus"]))
  sp <- as.character(as.matrix(taxa[[i]][,"species"]))
  
  tax <- list(genus = genus, species = sp)
  
  call <- get_accessions(filters = list(taxonomy = tax), 
                         at.least = 1000)
  
  gen_df[[i]] <- call
  
}

#..........................................
#..........................................
# Combine data from genesys query ####

vars <- c("accessionNumber","acquisitionDate","geo.latitude","geo.longitude",
          "coll.collSite","countryOfOrigin.code2","id","inSvalbard","institute.acronym",
          "taxonomy.genus","taxonomy.species","taxonomy.subtaxa","institute.fullName",
          "institute.country.name","institute.country.code2")

gen_sub <- data.frame(matrix(NA, 
                             ncol = length(vars),
                             nrow = 1,
                             dimnames = list(1, vars)))


for (i in seq_along(gen_df)) {
  
  if (dim(gen_df[[i]])[1] == 0) { 
    next()
  }
  
  # take the data
  x <- gen_df[[i]]
  
  # select the variable available in the data
  in_x <- vars %in% names(x)
  
  # if any missing variable, then add it as NAs
  if (any(!in_x)) {
    
    miss <- vars[!in_x]
    
    miss <- data.frame(matrix(NA, 
                              ncol = length(miss),
                              nrow = nrow(x),
                              dimnames = list(1:nrow(x), miss)))
    
    x <- cbind(x, miss)
    
  }
  
  x <- x[, vars]
  
  # bind with the main data
  gen_sub <- rbind(gen_sub, x)
  
  
}

# remove first line 
gen_sub <- gen_sub[-1, ]

in_gen <- unique(with(gen_sub, paste(taxonomy.genus, taxonomy.species)))

target <- unique(with(bind_rows(taxa), paste(genus, species)))

taxa <- bind_rows(taxa)

in_both <- taxa[in_gen %in% target,]

write.csv(in_both, "data/raw/in_both_databases.csv", row.names = FALSE)


keep <- with(gen_sub, paste(taxonomy.genus, taxonomy.species)) %in% in_both

gen_sub <- gen_sub[keep, ]

write.csv(gen_sub, file = "data/raw/genesys_occurrences.csv", row.names = FALSE)

save(gen_df, file = "data/raw/genesys_occurrences.rdm")



# plot(gen_sub$geo.longitude, gen_sub$geo.latitude)
# 
# lonlat <- gen_sub[,c("geo.longitude", "geo.latitude")]
# 
# names(lonlat) <- c("lon","lat")
# 
# library("rnaturalearth")
# library("rnaturalearthdata")
# 
# world <- ne_countries(scale = "medium", returnclass = "sf")
# class(world)
# 
# map1 <- 
#   ggplot(data = world) +
#   geom_sf() +
#   geom_point(data = lonlat, aes(x = lon, y = lat), size = 1, 
#              shape = 23, fill = "darkred") + 
#   labs(x="",y="")
# map1
# ggsave("map.png", map1, dpi = 500, width = 20, height = 15, units = "cm")
# 
