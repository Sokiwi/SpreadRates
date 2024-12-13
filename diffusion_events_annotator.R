library(pastclim)
# it is assumed that the entire Beyer2020 dataset has been downloaded
# using something like download_dataset(dataset = "Beyer2020", bio_variables = get_vars_for_dataset(dataset = "Beyer2020"))
# pastclim has functions for extraction data on net primary productivity (npp), 
# biome4, altitude, rugosity
library(sf)  # for reading the shape files
# library(rgeos)  # for checking if a point is within a polygon # no longer maintained
library(sp)  #  point.in.polygon()
library(geosphere)  # distHaversine()

# This produces a file called diffusion_events_annotated.txt from diffusion_events.txt
# Given two homeland coordinates and dates for each, get the following
# information about diffusion event:
# - distance as the crow flies from A to B, |AB|
# - speed (distance/time difference)
# - estimate of error on speed, err_speed
# - world area	
# - subsistence	
# - bearing
# - main biome4 traversed from A to B with X stations between, where X = |AB| / 10
# - percentage of cases where the main biome is involved ("perc_b4")
# - modern biomes (called biome98 because of Olson & Dinerstein (1998)), 
#   mainly to verify that NA for bio variables in the pastclim data
#   can be interpreted as water bodies, but also as a general sanity check
# - percentage of cases where the main biome98 is involved ("perc_b98")
# - average altitude across stations from A to B and the same X stations between as for biomes
# - average rugosity across stations from A to B and the same X stations between as for biomes
# - average net primary productivity (npp)

de <- read.table(file="diffusion_events.txt", header=TRUE)

# reduce it taking out cases where the mother is younger or equally old as
# the daughter or where age BP is undefined
time_diff_tmp <- de$BP_mo - de$BP_da
da_not_older <- which(time_diff_tmp <= 0)
age_undefined <- which(is.na(time_diff_tmp))
remove <- c(da_not_older, age_undefined)
if (length(remove) > 0) {
  dea <- de[-remove,]  # stands for diffusion events annotated, reduced by 498
}
# save(dea, file="dea.RData")

# get errors on homeland estimates for mothers and daughters
# takes around 2 mins
g <- read.table(file="glottolog_classification_strings.txt", header=TRUE, sep="\t", quote="")
error <- function(area) .4182 * sqrt(area) + 77.075
err_home_mo <- c()
err_home_da <- c()
area_mo <- c()
area_da <- c()
for (i in 1:nrow(dea)) {
  mo <- dea$mo[i]
  w_g_mo <- grep(paste0("^", mo), g$Classification)
  lons_mo <- g$lon[w_g_mo]
  lats_mo <- g$lat[w_g_mo]
  na_lons_mo <- which(is.na(lons_mo))
  na_lats_mo <- which(is.na(lats_mo))
  if (length(na_lons_mo) > 0 | length(na_lats_mo) > 0) {
    nas <- unique(c(na_lons_mo,na_lats_mo))
    lons_mo <- lons_mo[-nas]
    lats_mo <- lats_mo[-nas]
  }
  coors_mo <- data.frame(lons_mo,lats_mo)
  area_mo <- areaPolygon(coors_mo)/1000000
  err_home_mo[i] <- round(error(area_mo))
  da <- dea$da[i]
  w_g_da <- grep(paste0("^", da), g$Classification)
  lons_da <- g$lon[w_g_da]
  lats_da <- g$lat[w_g_da]
  na_lons_da <- which(is.na(lons_da))
  na_lats_da <- which(is.na(lats_da))
  if (length(na_lons_da) > 0 | length(na_lats_da) > 0) {
    nas_da <- unique(c(na_lons_da,na_lats_da))
    lons_da <- lons_da[-nas]
    lats_da <- lats_da[-nas]
  }
  coors_da <- data.frame(lons_da,lats_da)
  area_da <- areaPolygon(coors_da)/1000000
  err_home_da[i] <- round(error(area_da))
}
# save(err_home_mo, file="err_home_mo.RData")
# save(err_home_da, file="err_home_da.RData")

# sum_err_home <- err_home_mo + err_home_da

# for each diffusion event, make a list of stations between A and B
# geosphere also has distHaversine for the GCD
make_line <- function(lat1, lon1, lat2, lon2) {
  d <- distHaversine(c(lon1, lat1), c(lon2, lat2))/1000
  N <- round(d/10)
  points = gcIntermediate(c(lon1, lat1), c(lon2, lat1), 
                          n=N, addStartEnd=T)
  return(points)
}

lines <- list()
for (i in 1:nrow(dea)) {
  lines[[i]] <- make_line(dea$lat_mo[i], dea$lon_mo[i], dea$lat_da[i], dea$lon_da[i])
}

# distance as the crow flies from A to B, |AB|
distance <- c()
for (i in 1:nrow(dea)) {
  distance[i] <- distHaversine(c(dea$lon_mo[i], dea$lat_mo[i]), c(dea$lon_da[i], dea$lat_da[i]))/1000
}
# save(distance, file="distance.RData")

# speed
speed <- round(distance/(dea$BP_mo - dea$BP_da), 3)
# save(speed, file="speed.RData")

# read prepared data on families, areas, and subsistence from Hammarström
fas <- read.table(file="families_areas_subsistence.txt", header=TRUE, 
                  sep="\t", quote="", stringsAsFactors=FALSE, na.strings="", comment.char="")

# world area, subsistence
area <- c()
subsistence <- c()
for (i in 1:nrow(dea)) {
  curr_fam <- strsplit(dea$mo[i], ",")[[1]][1]
  w_fam <- which(fas$family==curr_fam)
    area[i] <- fas$area[w_fam]
    subsistence[i] <- fas$subsistence[w_fam]
  if ( area[i] == "Papunesia" & curr_fam == "Austronesian" ) {
    area[i] <- "Austronesian"
  }
  if ( area[i] == "Papunesia" & curr_fam != "Austronesian" ) {
    area[i] <- "Papuan"
  }
}
# save(area, file="area.RData")
# save(subsistence, file="subsistence.RData")

# bearing
bearing <- c()
for (i in 1:nrow(dea)) {
  if ( distance[i] > 0 ) {
    b <- geosphere::bearing(c(dea$lon_mo[i], dea$lat_mo[i]), c(dea$lon_da[i], dea$lat_da[i]))
    if ( b >= 45 & b < 135 ) { bearing[i] <- "EW" }
    if ( b < -45 & b >= -135 ) { bearing[i] <- "EW" }
    if ( b >= -45 & b < 45 ) { bearing[i] <- "NS" }
    if ( b < -135 & b >= -180 ) { bearing[i] <- "NS" }
    if ( b >= 135 & b < 180 ) { bearing[i] <- "NS" }
  } else {
    bearing[i] <- NA
  }
}
# save(bearing, file="bearing.RData")

# general function for getting data from Beyer et al. (2020)
# (see  https://doi.org/10.1038/s41597-020-0552-1)
pcd <- function(name, lon, lat, bp, biovar) {  # stands for pastclim data
  # options for biovar include npp, biome, altitude, rugosity, can be a vector of several
  locations <- data.frame(
    name = name,
    longitude = lon,
    latitude = lat,
    time_bp = -1 * bp
  )
  
  ls <- location_slice(
    x = locations, bio_variables = biovar,
    dataset = "Beyer2020", nn_interpol = TRUE
  )
}

# pcd() was modified to use the location_series() function since 
# location_slice() no longer seems to work for biomes
pcd2 <- function(taxonname, lon, lat, bp, biovar) {  # pcd stands for pastclim data
  # options for biovar include npp, biome, altitude, rugosity, can be a vector of several
  locations <- data.frame(
    name = taxonname,
    longitude = lon,
    latitude = lat,
    time_bp = -1 * bp
  )

  # rather than doing a slice a whole series of all periods is done (which 
  # works, unlike slices) and then the relevant period is extracted
  lseries <- location_series(
    x = locations, bio_variables = biovar,
    dataset = "Beyer2020", nn_interpol = TRUE
  )
  bp_period <- round(bp/1000)*-1000
  w_bp_period <- which(lseries$time_bp==bp_period)
  return(lseries[w_bp_period,])
}

# main reconstructed biome traversed from A to B with X stations between, 
# where X = |AB| / 10 percentage of cases where the main biome is involved ("percent")
# and percentage of subparts of trip pertaining to this main biome
# could take some 24 hours to run
biome4 <- c()
perc_b4 <- c()
for (i in 1:nrow(dea)) {
 if (i %% 10 == 0) {
   cat("finding biome4 no.", i, "out of", nrow(dea), "\n")
 }
 pcd_out <- pcd2(dea$da[i], lines[[i]][,1], lines[[i]][,2], round((dea$BP_mo[i] + dea$BP_da[i])/2), "biome")
 w_na <- which(is.na(pcd_out$biome))
 if( length(w_na) > 0 ) {
   pcd_out$biome <- "Sea"
 }
 biome_table <- sort(table(as.vector(pcd_out$biome)), decreasing=TRUE)
 mb <- names(biome_table[1])  # stands for main biome
 perc <- round(100*unname(biome_table[1]/sum(biome_table)),1)
 if (length(mb) > 0) {
   biome4[i] <- mb
 } else {
   biome4[i] <- NA
 }
 if (!is.na(perc)) {
   perc_b4[i] <- perc
 } else {
   perc_b4[i] <- NA
 }
}
save(perc_b4, file="perc_b4.RData")
save(biome4, file="biome4.RData")

# function given two geographical coordinates find the proportion of 
# the type of modern biome that most of the connecting line goes through
# NA output means that there is no movement
# Olson & Dinerstein (1998) biomes for WWF
shape <- read_sf(dsn = "C:/Wichmann/ASJP/Homelands/Homeland2025/official_teow/official", 
                 layer = "wwf_terr_ecos")
biome98finder <- function(i) {
  points <- lines[[i]]
  
  # function for finding if a given point is within one of the polygons
  is_in <- function(lon, lat, shape, s=0) {
    if ( s == 0 ) {
      for (i in 1:length(shape[[22]])) {
        m <- as.matrix(shape[[22]][i][[1]])
        x <- m[,1]
        y <- m[,2]
        if ( point.in.polygon(lon, lat, x, y) > 0 ) {
          return(i)
        }
      }
    }
    # s is the index in the shape file to look at first
    if ( s != 0 ) {
      m <- as.matrix(shape[[22]][s][[1]])
      x <- m[,1]
      y <- m[,2]
      if ( point.in.polygon(lon, lat, x, y) > 0 ) {
        return(s)
      } else {
        for (i in 1:length(shape[[22]])) {
          m <- as.matrix(shape[[22]][i][[1]])
          x <- m[,1]
          y <- m[,2]
          if ( point.in.polygon(lon, lat, x, y) > 0 ) {
            return(i)
          }
        }
      }
    }
    if ( i == length(shape[[22]]) ) {
      return(0)
    }
  }
  b <- c()
  for (p in 1:length(points[,1])) {
    if ( p==1 ) {
      pgon <- is_in(points[p,1], points[p,2], shape, s=0)
      spgon <- pgon
    } else {
      pgon <- is_in(points[p,1], points[p,2], shape, s=spgon)
      spgon <- pgon
    }
    if ( pgon == 0 ) {
      biome <- 100
    } else {
      biome <- shape[[6]][pgon]
    }
    b[p] <- biome
  }
  biomes_n <- as.vector(table(b))
  max_biomes <- which(biomes_n==max(biomes_n))[1]
  percent <- round(100*max(biomes_n)/sum(biomes_n),1)
  most_freq_biome <- as.numeric(names(table(b))[max_biomes])
  return(list(most_freq_biome, percent))
}

# running the following could take 1-2 days
# output is optionally saved
biome98 <- c()
perc_b98 <- c()
for (i in 1:nrow(dea)) {
  if (i %% 100 == 0) {
    cat("finding biome98 no.", i, "out of", nrow(dea), "\n")
  }
  biome98out <- biome98finder(i)
  biome98[i] <- biome98out[[1]]
  perc_b98[i] <- biome98out[[2]]
}
save(biome98, file="biome98.RData")
save(perc_b98, file="perc_b98.RData")

# mean altitude, rugosity, and net primary productivity across stations 
# from A to B and the same X stations between as for biomes
altitude <- c()
rugosity <- c()
npp <- c()
for (i in 1:nrow(dea)) {
  pcd_out <- pcd(dea$da[i], lines[[i]][,1], lines[[i]][,2], 
                 round((dea$BP_mo[i] + dea$BP_da[i])/2), c("altitude", "rugosity", "npp"))
  altitude[i] <- round(mean(pcd_out$altitude),2)
  rugosity[i] <- round(mean(pcd_out$rugosity),2)
  npp[i] <- round(mean(pcd_out$npp),2)
}
save(altitude, file="altitude.RData")
save(rugosity, file="rugosity.RData")
save(npp, file="npp.RData")

biome4_names <- biome4
save(biome4_names, file="biome4_names.RData")
biome4 <- c()
biome4_key <- get_biome_classes("Example")
for (i in 1:length(biome4_names)) {
  if (!is.na(biome4_names[i])) {
    w_biome4_name <- which(biome4_key$category==biome4_names[i])
    if (length(w_biome4_name) > 0) {
      biome4[i] <- biome4_key$id[w_biome4_name]
    } else {
      biome4[i] <- NA
    }
  } else {
    biome4[i] <- NA
  }
}
save(biome4, file="biome4.RData")

biome98_names <- c()
biome98_key <- read.table(file="key biomes98.txt", header=TRUE, sep="\t")
for (i in 1:length(biome98)) {
  if (!is.na(biome98[i])) {
    w_biome98 <- which(biome98_key$id==as.numeric(biome98[i]))
    if (length(w_biome98) > 0) {
      biome98_names[i] <- biome98_key$category[w_biome98]
    } else {
      biome98_names[i] <- NA
    }
  } else {
    biome98_names[i] <- NA
  }
}
save(biome98_names, file="biome98_names.RData")

# in case only some elements were updated as others were saved from 
# previous sessions the saved ones can be loaded as needed
load("dea.RData")
load("err_home_mo.RData")
load("err_home_da.RData")
load("distance.RData")
load("speed.RData")
load("area.RData")
load("subsistence.RData")
load("bearing.RData")
load("biome4.RData")
load("biome4_names.RData")
load("perc_b4.RData")
load("biome98.RData")
load("biome98_names.RData")
load("perc_b98.RData")
load("altitude.RData")
load("rugosity.RData")
load("npp.RData")

dea_full <- data.frame(dea, err_home_mo, err_home_da, distance, speed, area, subsistence, bearing, 
                       biome4, biome4_names, perc_b4, biome98, biome98_names, 
                       perc_b98, altitude, rugosity, npp)

write.table(dea_full, file="diffusion_events_annotated.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")

