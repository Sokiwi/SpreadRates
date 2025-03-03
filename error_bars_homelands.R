library(swaRm)  # chull_area()

# This script defines linear formulas for computing expected errors (distance 
# from true homeland) from the area covered by a language family for the two 
# methods BayesTraits with fixed rates (BT) and Minimum Distance (MD). The
# estimates are based on simulated data from Wichmann and Rama (2021).

# Download supplementary materials from Wichmann and Rama (2021) from
# https://figshare.com/articles/dataset/Testing_methods_of_linguistic_homeland_detection_SI/13308926?file=25642538
# relevant files are in SI-03 language diffusion simulation output;
# file names indicate true homelands and the last line are coordinates
# for language locations at the end of simulations
# these should be read and areas and errors inferred

# Replace my path with your own location of the SI-03 files
# path to simulated locations for families
my.path.sim <- "C:/Wichmann/Library/GeneralLibrary/My publications/Wichmann and Rama 2021 SI-01-13/SI-03 language diffusion simulation output/"
# path to results for baseline methods, including md
my.path.md <- "C:/Wichmann/Library/GeneralLibrary/My publications/Wichmann and Rama 2021 SI-01-13/SI-04 baseline methods software and output/"
# path to results for BayesTraits
my.path.bt <- "C:/Wichmann/Library/GeneralLibrary/My publications/Wichmann and Rama 2021 SI-01-13/SI-05 bayestraits fixed rates/"

# read output files for MD and BT
baselines.read <- read.table(file=paste0(my.path.md, "baselines.txt"), sep="\t", header=TRUE, strip.white = TRUE)
bt.read <- read.table(file=paste0(my.path.bt, "BTF_results.txt"), sep="\t", header=TRUE, strip.white = TRUE)

# get prefix for the 1000 simulation files
prefixes <- baselines.read$file_prefix

# function for calculating the area of a convex hull (envelope) of a set of 
# geographical coordinates (in square kilometers)
area <- function(x, y) {  # x and y are vectors of longitudes and latitudes
  chull_area(x, y, geo=TRUE)/1000000
}

# function for reading a pair of simulation output files given a file name prefix
# and extracting a set of either longitudes or latitudes
# example of prefix: F_20_wanej_-0.3659_121.948
read.coors <- function(prefix) {
  lon.file <- paste0(my.path.sim, prefix, ".longs")
  lat.file <- paste0(my.path.sim, prefix, ".lats")
  lon.read <- read.table(lon.file, header=TRUE, row.names = 1, sep="\t")
  lat.read <- read.table(lat.file, header=TRUE, row.names = 1, sep="\t")
  lons <- as.vector(unlist(lon.read[nrow(lon.read),]))
  lats <- as.vector(unlist(lat.read[nrow(lat.read),]))
  return(list(lons,lats))
}

# function for getting the absolute md error given a certain file prefix
# the file gives errors in km
md.abs.err <- function(prefix) {
  w_prefix <- grep(prefix, baselines.read$file_prefix)
  md.err <- baselines.read$error_mindist[w_prefix]
  return(md.err)
}

# function for getting the absolute BT error given a certain file prefix
# the file gives errors in km
bt.abs.err <- function(prefix) {
  w_prefix <- grep(prefix, bt.read$fileprefix)
  bt.err <- bt.read$error[w_prefix]
  return(bt.err)
}

# loop through files for MD and BT and get the errors and areas
areas <- c()
md.errors <- c()
bt.errors <- c()
for (i in 1:length(prefixes)) {
  if (i %% 100 == 0) {
    cat("doing", i, "out of 1000\n")
  }
  md.errors[i] <- md.abs.err(prefixes[i])
  bt.errors[i] <- bt.abs.err(prefixes[i])
  coors <- read.coors(prefixes[i])
  areas[i] <- area(coors[[1]], coors[[2]])
}

summary(lm(md.errors~areas))
# Residual standard error: 154.6 on 998 degrees of freedom
# Multiple R-squared:  0.222,	Adjusted R-squared:  0.2212 
# F-statistic: 284.8 on 1 and 998 DF,  p-value: < 2.2e-16

intercept.md <- as.vector(lm(md.errors~areas)$coefficients[1])  # 111.4923
slope.md <- as.vector(lm(md.errors~areas)$coefficients[2])  # 0.0002619432

summary(lm(bt.errors~areas))
# Residual standard error: 147.9 on 998 degrees of freedom
# Multiple R-squared:  0.2185,	Adjusted R-squared:  0.2177 
# F-statistic:   279 on 1 and 998 DF,  p-value: < 2.2e-16

intercept.bt <- as.vector(lm(bt.errors~areas)$coefficients[1])  # 108.947
slope.bt <- as.vector(lm(bt.errors~areas)$coefficients[2])  # 0.0002481033
