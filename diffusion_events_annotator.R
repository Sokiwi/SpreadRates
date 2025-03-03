library(pastclim)
# it is assumed that the entire Beyer2020 dataset has been downloaded
# using something like download_dataset(dataset = "Beyer2020", bio_variables = get_vars_for_dataset(dataset = "Beyer2020"))
# pastclim has functions for extraction data on net primary productivity (npp), 
# biome4, altitude, rugosity
library(sf)  # for reading the shape files
# library(rgeos)  # for checking if a point is within a polygon # no longer maintained
library(sp)  #  point.in.polygon()
library(geosphere)  # distHaversine()
library(swaRm)  # chull_area()
library(tibble)

# This produces a file called diffusion_events_annotated.txt from diffusion_events.txt
# Given two homeland coordinates and dates for each, add the following
# information about diffusion event:
# - distance as the crow flies from A to B, |AB|
# - estimate of error on distance
# - estimate of error on time differences
# - velocity (distance/time difference)
# - estimate of error on velocity, err_vel
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
# in the following two functions the intercepts and slopes come
# from the script error_bar_homelands

# formula for estimating error from area of language group for BT method
get.error.bt <- function(area) {
  # error.bt = slope.bt * area + intercept.bt
  error.bt <- 0.0002481033 * area + 108.947
  return(error.bt)
}

# formula for estimating error from area of language group for MD method
get.error.md <- function(area) {
  # error.md = slope.md * area + intercept.md
  error.md <- 0.0002619432 * area + 111.4923
  return(error.md)
}

# given a classification string access the pruned ASJP database,
# extract coordinates for the languages, compute the area, and then the 
# expected error

# extracting coordinates requires the database
db <- read.table(file="listss20_pruned_updated.tab", sep="\t", header=TRUE, 
                 quote = "", na.strings="", comment.char="")

# function for calculating the area of a convex hull (envelope) of a set of 
# geographical coordinates (in square kilometers)
area <- function(x, y) {  # x and y are vectors of longitudes and latitudes
  chull_area(x, y, geo=TRUE)/1000000
}

# function for computing the area and then the error
get.exp.err <- function(clstring, method) {
  w_clstring <- grep(clstring, db$hh, fixed=TRUE)
  if (length(w_clstring)==0) {
    return(NA)
  }
  lons <- db$lon[w_clstring]
  lats <- db$lat[w_clstring]
  a <- area(lons, lats)
  if(method=="MD") {
    return(get.error.md(a))
  }
  if(method=="BT") {
    return(get.error.bt(a))
  }
}

# given two homelands get the error on them using get.exp.err
# (which takes the method as a variable),
# and subsequently an error on the distance;
# mo and da are classification strings;
# the function takes as variable a row in the diffusion_events.txt file
# or a data frame based on this file, such as dea
# (dea is assumed to be present in memory)

err.bar.dist <- function(i) {
  mo <- dea$mo[i]
  method.mo <- dea$method_mo[i]
  da <- dea$da[i]
  method.da <- dea$method_da[i]
  lon.mo <- dea$lon_mo[i]
  lat.mo <- dea$lat_mo[i]
  lon.da <- dea$lon_da[i]
  lat.da <- dea$lat_da[i]
  exp.err.mo <- get.exp.err(mo, method.mo)
  exp.err.da <- get.exp.err(da, method.da)
  dist.mo.da <- distHaversine(c(lon.mo,lat.mo),c(lon.da,lat.da))/1000
  if (is.na(exp.err.mo) | is.na(exp.err.da)) {
    dist.min <- NA
    dist.max <- NA
    print("NA case found")  # not expected to happen
    return(c(dist.mo.da,dist.min,dist.max))
  }
  dist.max <- exp.err.mo + exp.err.da + dist.mo.da
  if (dist.mo.da >= exp.err.mo + exp.err.da) {  # non-overlapping
    dist.min <- dist.mo.da - (exp.err.mo + exp.err.da)
  } else {  # overlapping or one inside the other
    dist.min <- 0
  }
  return(c(dist.mo.da,dist.min,dist.max))
}

# run all diffusion events, filling the below vectors for distances and error bars
distance <- c()
dist_min <- c()
dist_max <- c()

for (i in 1:nrow(dea)) {
  if(i %% 100 == 0) {
    cat("doing", i, "out of", nrow(dea), "\n")
  }
  res <- err.bar.dist(i)
  distance[i] <- res[1]
  dist_min[i] <- res[2]
  dist_max[i] <- res[3]
}

# save(distance, file="distance.RData")
# save(dist_min, file="dist_min.RData")
# save(dist_max, file="dist_max.RData")

# get error bars on dates for mothers and daughters using the plus/minus 29%
# error on calibration in Holman et al.
err.bar.date <- function(i) {
  BP.mo <- dea$BP_mo[i]
  BP.da <- dea$BP_da[i]
  BP.mo.from <- round(BP.mo + BP.mo * .29, 0)
  BP.mo.to <- round(BP.mo - BP.mo * .29, 0)
  if (!is.na(BP.mo.to) & BP.mo.to < 0) {BP.mo.to <- 0}
  BP.da.from <- round(BP.da + BP.da * .29, 0)
  BP.da.to <- round(BP.da - BP.da * .29, 0)
  if (!is.na(BP.da.to) & BP.da.to < 0) {BP.da.to <- 0}
  return(c(BP.mo.from, BP.mo.to, BP.da.from, BP.da.to))
}

BP_mo_from <- c()
BP_mo_to <- c()
BP_da_from <- c()
BP_da_to <- c()

# run all diffusion events and fill the above vectors for error bars on dates
for (i in 1:nrow(dea)) {
  res <- err.bar.date(i)
  BP_mo_from[i] <- res[1]
  BP_mo_to[i] <- res[2]
  BP_da_from[i] <- res[3]
  BP_da_to[i] <- res[4]
}

# save(BP_mo_from, file="BP_mo_from.RData")
# save(BP_mo_to, file="BP_mo_to.RData")
# save(BP_da_from, file="BP_da_from.RData")
# save(BP_da_to, file="BP_da_to.RData")

# compute speed and error bars on speed;
# this uses propagation of uncertainties for division
# where errV^2 = errD^2 + errT^2,
# respectively error on velocity, distance, and time
err.bar.speed <- function(i) {
  D <- distance[i]
  # when areas of uncertainty are overlapping the minimal distance is zero
  # in this situation the proportion of possible error to distance can be large 
  # or infinite, but in practice the measurement on distance is probably
  # going to be quite accurate, so complications in the calculation of
  # propagation of uncertainty are avoided by simply setting the error to zero
  if (dist_min[i]==0) {
    errD <- 0
  } else {
    abs.errD <- (dist_max[i] - dist_min[i])/2
    errD <- abs.errD/D
  }
  if (is.na(dea$BP_mo[i]) | is.na(dea$BP_da[i])) {
    return(c(NA, NA, NA))
  }
  T <- dea$BP_mo[i] - dea$BP_da[i]
  if (T <= 0) {
    return(c(NA, NA, NA))
  }
  errT <- .29
  V <- D/T
  errV <- V * sqrt(errD^2 + errT^2)
  V.from <- V - errV
  if (V.from < 0) {
    V.from <- 0
  }
  V.to <- V + errV
  return(c(V, V.from, V.to, errV))
}

speed <- c()
speed_from <- c()
speed_to <- c()
speed_err <- c()

# run all diffusion events and fill the above vectors for error bars on speed
for (i in 1:nrow(dea)) {
  res <- err.bar.speed(i)
  speed[i] <- res[1]
  speed_from[i] <- res[2]
  speed_to[i] <- res[3]
  speed_err[i] <- res[4]
}

# save(speed, file="speed.RData")
# save(speed_from, file="speed_from.RData")
# save(speed_to, file="speed_to.RData")

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
# save(perc_b4, file="perc_b4.RData")
# save(biome4, file="biome4.RData")

# function given two geographical coordinates find the proportion of 
# the type of modern biome that most of the connecting line goes through
# NA output means that there is no movement
# Olson & Dinerstein (1998) biomes for WWF
shape <- read_sf(dsn = "E:/Wichmann/ASJP/Homelands/Homeland2025/official_teow/official", 
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
# save(biome98, file="biome98.RData")
# save(perc_b98, file="perc_b98.RData")

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
# save(altitude, file="altitude.RData")
# save(rugosity, file="rugosity.RData")
# save(npp, file="npp.RData")

biome4_names <- biome4
# save(biome4_names, file="biome4_names.RData")
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
# save(biome4, file="biome4.RData")

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
# save(biome98_names, file="biome98_names.RData")

# in case only some elements were updated and others were saved from 
# previous sessions the saved ones can be loaded as needed
load("dea.RData")
load("distance.RData")
load("dist_min.RData")
load("dist_max.RData")
load("speed.RData")
load("speed_from.RData")
load("speed_to.RData")
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

dea_full <- data.frame(dea, distance, dist_min, dist_max, 
                       speed, speed_from, speed_to, area, subsistence, bearing, 
                       biome4, biome4_names, perc_b4, biome98, biome98_names, 
                       perc_b98, altitude, rugosity, npp)

write.table(dea_full, file="diffusion_events_annotated.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")

