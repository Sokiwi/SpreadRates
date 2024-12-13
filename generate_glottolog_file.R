## This script generates a file from Glottolog which can subsequently
## used by the script md_fams_subgroups.R to compute
## homelands using the minimal distance (md) method of Wichmann
## and Rama (2021) [doi: 10.1098/rstb.2020.0202]
## Here the file is just used to get classification strings
library(stringi)

# read languoid.csv from https://glottolog.org/meta/downloads
# put it in the same folder as the present script
g <- read.csv(file="languoid.csv", encoding="UTF-8")

# generate a file with the headers "ISO639-3", "Classification", "lat", and "lon" (i.e. coordinates)

## functions used in generating the file

# function for creating a pseudo-code when ISO-code is missing
NOCODEstring <- function(s) {
	s1 <- gsub(" ", "", s)
	s2 <- stri_trans_general(s1, "Latin-ASCII")
	s3 <- paste("NOCODE_", s2, sep="")
	return(s3)
}

# function for formatting a classification string
# changing it to ascii and taking out the last comma
class_string_format <- function(s) {
	s1 <- gsub(" ", "", s)
	s2 <- stri_trans_general(s1, "Latin-ASCII")
	s3 <- strsplit(s2, ",")[[1]]
	s4 <- paste(s3, collapse=",")
	return(s4)
}

# function for composing a classification string
class_string <- function(ind) {
	ancestors <- ""
	parent_id <- "master_Wu"
	while ( parent_id != "" ) {
		parent_id <- g$parent_id[ind]
		if ( parent_id != "" ) { 
			w <- which(g$id==parent_id)
			ancestors <- paste(g$name[w], ",", ancestors, sep="")
			ind <- w
		}
	}
	if ( ancestors == "" ) {
		ancestors <- g$name[ind]
	}
	return(ancestors)
}

# routine for composing a classification string for each language
# including ISO-code languages considered dialects in Glottolog
# this is used for ASJP purposes, but not for homelands, so 
# put behind comment characters here
# cat("ISO639-3\tClassification\n", file="glottolog.tab")
# for (i in 1:nrow(g)) {
# 	if ( g$level[i]=="language" | nchar(g$iso639P3code[i])==3 ) {
# 		iso <- g$iso639P3code[i]
# 		if ( iso == "" ) {
# 			iso <- NOCODEstring(g$name[i])
# 		}
# 		clstr <- class_string(i)
# 		clstrf <- class_string_format(clstr)
# 		cat(iso, "\t", clstrf, "\n", sep="", file="glottolog.tab", append=TRUE)
# 		cat(iso, "\t", clstrf, "\n")
# 	}
# }


# routine for composing a classification string for each language
# NOT including ISO-code languages considered dialects in Glottolog
# AND also including coordinates,
# which can be used to get homelands for each
# family and subgroup
cat("ISO\tGlottologID\tClassification\tlat\tlon\n", file="glottolog_classification_strings.txt")
for (i in 1:nrow(g)) {
	if ( g$level[i]=="language" ) {
		iso <- g$iso639P3code[i]
		if ( iso == "" ) {
			iso <- NOCODEstring(g$name[i])
		}
		glotid <- g$id[i]
		lat <- g$latitude[i]
		lon <- g$longitude[i]
		clstr <- class_string(i)
		clstrf <- class_string_format(clstr)
		cat(iso, "\t", glotid, "\t", clstrf, "\t", lat, "\t", lon, "\n", sep="", file="glottolog_classification_strings.txt", append=TRUE)
		# cat(iso, "\t", clstrf, "\n")
	}
}

