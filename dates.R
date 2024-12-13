# Prior to the run a selection was made of lists that have at least
# 26 words and describe a language attested after 1700.
# All languages that are creoles, pidgins, and constructed
# according to the WALS classification were removed.
# This filtering was done using the interactive R program at 
# https://github.com/Sokiwi/InteractiveASJP02

# The Glottolog classification in listss20_pruned.txt is changed to reflect 
# version 5.1
# The script requires listss20_pruned_updated.txt, preamble.txt, asjp62c.exe
# A large part of the script, involving the function combine() is not used,
# and is commented out
# This was used in earlier studies to characterize diffusion events
# in relation to modern biomes
# To run the dates type just read in the script, type run_all() and press return

# wrapper for all dates and homelands
run_all <- function() {
	clean_up_before()
	initiate_files()
	fams <- families("listss20_pruned_updated.txt")
	for (i in 1:length(fams)) {
		cat("\n\ndoing", fams[i], "which is", i, "out of", length(fams), "\n")
		family <- fams[i]
		dates(family)
	}
	# combine(walking_distance=TRUE)
	clean_up_after()
}

# clean away old files and start afresh
clean_up_before <- function() {
	# function for removing files if they exist
	delete <- function(to_remove) {
		if ( length(which(dir()==to_remove)) > 0 ) {
			invisible(file.remove(to_remove))
		}
	}
	all_to_delete <- c("dates.txt", "family_dates.txt", 
	  "dates.bat", "data_family.txt")
	for (f in 1:length(all_to_delete)) {
		delete(all_to_delete[f])
	}
}

# initiate dates file
initiate_files <- function() {
	# initiate dates file
	cat(c("group", "pairs", "BP"), sep="\t", file="dates.txt")
	cat("\n", file="dates.txt", append=TRUE)
}

# get a list of Glottolog families excluding singletons and some unwanted ones
families <- function(listss="listss20_pruned_updated.txt") {
	x <- readLines(listss)
	w_meta_1 <- grep("@", x)
	w_meta_2 <- grep("\\{", x)
	w_meta <- intersect(w_meta_1, w_meta_2)
	get_fam <- function(line) {
		f1 <- strsplit(line, "@")[[1]][2]
		f2 <- strsplit(f1, ",")[[1]][1]
		f3 <- strsplit(f2, "\\}")[[1]][1]
		return(f3)
	}
	fams_raw1 <- as.vector(unlist(lapply(x[w_meta], get_fam)))
	fam_table <- table(fams_raw1)
	fam_table_numbers <- as.vector(fam_table)
	fam_table_names <- names(fam_table)
	w_singletons <- which(fam_table_numbers==1)
	singletons <- fam_table_names[w_singletons]
	fams_raw2 <- sort(unique(fams_raw1))
	unwanted <- unique(c("", "ArtificialLanguage", "NA", "MixedLanguage", 
	  "Pidgin", "SpeechRegister", "Spurious", singletons))
	fams_raw3 <- setdiff(fams_raw2, unwanted)
	w_uncl <- grep("Unclassifiable", fams_raw3)
	if ( length(w_uncl) > 0 ) {
		fams_raw3 <- fams_raw3[-w_uncl]
	}
	return(fams_raw3)
}

clean_up_after <- function() {
	# function for removing files if they exist
	delete <- function(to_remove) {
		if ( length(which(dir()==to_remove)) > 0 ) {
			invisible(file.remove(to_remove))
		}
	}
	all_to_delete <- c("family_dates.txt", "dates.bat")
	for (f in 1:length(all_to_delete)) {
		delete(all_to_delete[f])
	}
}

# run dates for families and subgroups and collect them in dates.txt
dates <- function(family, listss="listss20_pruned_updated.txt", shape) {
	x <- readLines(listss)
	w_start <- grep("\\{", x)[1]
	x <- x[-c(1:(w_start-1))]
	# search for lists pertaining to a family and extract the set
	# this should output a file and a number of levels
	# then the function group_set goes through the file and makes
	# appropriate subset, runs the Holman softwares

	# vector of beginnings of each list pertaining to the family
	w_s <- grep(paste("@", family, sep=""), x)

	# vector of all ends of lists, including the last of the file
	first <- grep("\\{", x)
	last <- first[-1] - 1
	very_last <- grep("^     ", x) - 1
	last <- c(last, very_last)

	# function for finding the closest number in the list called last
	# after a given number in the list called first
	find_last <- function(start, last) {
		last[which(last > start)[1]]
	}

	# now look for the family members and generate the input file
	first_fam_1 <- grep(paste("@", family, "}", sep=""), x)
	first_fam_2 <- grep(paste("@", family, ",", sep=""), x)
	first_fam <- sort(union(first_fam_1, first_fam_2))
	pa <- readLines("preamble.txt")
	cat(pa, file="data_family.txt", sep="\n")
	last_fam <- as.vector(unlist(lapply(first_fam, find_last, last)))
	for (i in 1:length(first_fam)) {
		cat(x[first_fam[i]:last_fam[i]], file="data_family.txt", 
		  sep="\n", append=TRUE)
	}
	cat("     \n", file="data_family.txt", append=TRUE)

	# change the format so that the Glottolog classification
	# appears in the Ethnologue slot
	x <- readLines("data_family.txt")
	w_meta <- grep("\\|", x)
	# function for deleting the Ethnologue classification
	# so the Glottolog classification appears in its place
	for (i in 1:length(w_meta)) {
		x[w_meta[i]] <- gsub("\\|.*@", "|", x[w_meta[i]])
	}
	cat(x, file="data_family.txt", sep="\n")

	# now run Holman's dating software on the top (family) level
	cat(paste("asjp62c <", "data_family.txt", "> family_dates.txt\n"), 
	  file="dates.bat")
	# system("dates.bat", show.output.on.console=FALSE)
	system("dates.bat")
	
	# function for reading the family_dates.txt file and
	# putting the results in dates.txt
	# Then repeat for as many levels as can be identified
	get_date_and_group_info <- function(level) {
		d <- readLines("family_dates.txt")
		if ( length(d) > 4 ) {
			cat("  level", level, "\n")
			# make one called w_group, of which there are 
			# as many as w_info
			w_info <- grep("PAIRS", d)
			w_group <- which(d == "") + 1
			w_group <- w_group[-(length(w_group))]
			classification_strings <- c()
			for (k in 1:length(w_info)) {
				info1 <- gsub(" ", "", d[w_info[k]])
				pairs <- strsplit(info1, "PAIRS")[[1]][1]
				info2 <- strsplit(info1, "%")[[1]][2]
				years <- strsplit(info2, "BP")[[1]][1]
				group1 <- strsplit(d[w_group[k]], "[[:digit:]]  ")[[1]][2]
				group2 <- strsplit(group1, ",")[[1]][1:level]
				group3 <- paste(group2, collapse=",")
				group <- strsplit(group3, "\\}")[[1]][1]
				cat(group, "\t", pairs, "\t", years, "\n", file="dates.txt", 
				  sep="", append=TRUE)
				classification_strings[k] <- group
			}
			return(classification_strings)
		}
	}
	get_date_and_group_info(1)

	# Get the number of levels
	# and get dates

	# collect all G classifications in a vector
	get_g <- function(meta) {
		g1 <- strsplit(meta, "\\|")[[1]][2]
		g2 <- strsplit(g1, "\\}")[[1]][1]
		return(g2)
	}
	g <- as.vector(unlist(lapply(x[w_meta], get_g)))
	# code for getting the number of commas in the g classification
	count_commas <- function(x) {
		length(strsplit(x, ",")[[1]])
	}
	commas <- as.vector(unlist(lapply(g, count_commas)))
	levs <- max(commas)

	# code for going through lower levels and getting dates
	if ( levs > 1 & length(unique(g)) > 1 ) {
		for (L in 2:levs) {
			if ( L < 10 ) {
				newfirst1 <- strsplit(x[1],"")[[1]]
				newfirst1[23] <- " "
				newfirst1[24] <- L
				x[1] <- paste(newfirst1, collapse="")
				cat(x, file="data_family.txt", sep="\n")
			}
			if ( L >= 10 ) {
				firstdigit <- (L - (L %% 10))/10
				seconddigit <- L - (L - (L %% 10))
				newfirst1 <- strsplit(x[1],"")[[1]]
				newfirst1[23] <- firstdigit
				newfirst1[24] <- seconddigit
				x[1] <- paste(newfirst1, collapse="")
				cat(x, file="data_family.txt", sep="\n")
			}
			cat(paste("asjp62c <", "data_family.txt", "> family_dates.txt\n"), 
			  file="dates.bat")
			# system("dates.bat", show.output.on.console=FALSE)
			system("dates.bat")
			cs <- get_date_and_group_info(L)  # classification strings
		}
	}
}

combine <- function(walking_distance=TRUE) {
	# combine the outputs for dates and homelands such that for all homelands
	# for which there are also dates a row is produced
	d <- read.table(file="dates.txt", header=TRUE, sep="\t", quote="", 
	  stringsAsFactors=FALSE, na.strings="", comment.char="")
	h <- read.table(file="homelands_glottolog_groups.txt", header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE, na.strings="", comment.char="")
	# all _ and spaces are converted to nothing by hand in both tables
	remove_us <- function(s) gsub("_", "", s)
	h[,"group"] <- as.vector(unlist(lapply(h[,"group"],trimws)))
	h[,"group"] <- as.vector(unlist(lapply(h[,"group"],remove_us)))
	d[,"group"] <- as.vector(unlist(lapply(d[,"group"],trimws)))
	d[,"group"] <- as.vector(unlist(lapply(d[,"group"],remove_us)))
	cat("group\tpairs_date\tBP\tlat\tlon\tN_homeland\n", 
	  file="dates_homelands.txt")
	for (i in 1:length(h[,1])) {
		w_h <- which(d$group==h$group[i])
		if ( length(w_h) > 0 ) {
			comb <- c(d$group[w_h], d$pairs[w_h], d$BP[w_h], 
			  h$lat[i], h$lon[i], h$N[i])
			cat(comb, sep="\t", file="dates_homelands.txt", append=TRUE)
			cat("\n", file="dates_homelands.txt", append=TRUE)
		}
	}
	cat("Done combining info for", i, "dates and homelands\n\n")

 	# read the file with combined date and homeland info from the preceding
	# and make a list of homeland events: cases where a mother-daughter
	# relationship can be identified; measure distance using both GCD and,
	# if walking_distance==TRUE (the default), 
	# the best walking distance measure from the Physica A paper;
	# that can be excluded because it is time consuming
	cat("Now identifying migration events and calculating distances\n")
	if (walking_distance==TRUE) {
		source("DD-.25.R")
	}
	source("distance_function.R")
	f_a <- read.table(file="families_areas_subsistence.txt", header=TRUE, 
	  sep="\t", quote="", stringsAsFactors=FALSE, na.strings="", comment.char="")
	x <- read.table(file="dates_homelands.txt", header=TRUE, sep="\t", quote="", 
	  stringsAsFactors=FALSE, na.strings="", comment.char="")
	# get rid of unclassifiables among the groups
	w_u <- grep("Unclassified", x$group)
	if ( length(w_u) > 1 ) {
		x <- x[-w_u,]
	}
	cat("M\tpairs_date_M\tBP_M\tlat_M\tlon_M\tN_homeland_M\tD\tpairs_date_D\tBP_D\tlat_D\tlon_D\tN_homeland_D\tGCD\tDD\trate_GCD\trate_DD\tarea\tsubsistence\tbearing\tbiome\tpercent\n", file="migration_events.txt")

	# function given two geographical coordinates find the proportion of 
	# the type of biome that most of the connecting line goes through
	# NA output means that there is no movement
	library(sf)  # for reading the shape file
	source("distance_function.R")  # distance() for computing distances
	library(rgeos)  # for checking if a point is within a polygon
	library(sp)  # for point.in.polygon function
	shape <- read_sf(dsn = "C:/Wichmann/ASJP/Homelands/Homeland2022/official_teow/official", 
	  layer = "wwf_terr_ecos")

	biomes <- function(lat1, lon1, lat2, lon2, shape) {
	library(geosphere)  # for function gcIntermediate
		make_line <- function(lat1, lon1, lat2, lon2) {
			d <- distance(lat1, lat2, lon1, lon2)
			N <- round(d/10)
			points = gcIntermediate(c(lon1, lat1), c(lon2, lat1), 
			  n=N, addStartEnd=T)
			return(points)
		}
		points <- make_line(lat1, lon1, lat2, lon2)

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

	for (i in 1:length(x[,1])) {
		if ( i %% 10 == 0 ) {
			cat("doing", i, "out of", length(x[,1]), "\n")
		}
		D_cand <- x[i,1]
		L <- length(strsplit(D_cand, ",")[[1]])
		if ( L > 1 ) {
			M_search1 <- strsplit(D_cand, ",")[[1]][-L]
			M_search2 <- paste(M_search1, collapse=",")
			w_M_search2 <- which(x[,1]==M_search2)
			if ( length(w_M_search2) > 1 ) {
				cat("too many mothers, decide what to do\n")
			}
			if ( length(w_M_search2) == 1 ) {
				M_info <- x[w_M_search2,]
				# now distances are computed (rounded off to digits)
				# lats and longs are spelled out here for more clarity
				M_lat <- M_info$lat
				M_lon <- M_info$lon
				D_lat <- x$lat[i]
				D_lon <- x$lon[i]
				GCD <- round(distance(M_lat, D_lat, M_lon, D_lon))
				# The DD always assumes a movement so is set to 0
				# when it actually is zero
				if ( walking_distance==TRUE ) {
					DD <- suppressWarnings(DDquarter(M_lat,M_lon,D_lat,D_lon))
				} else {
					DD <- NA
				}
				if ( GCD == 0 ) {
					DD <- 0
				}
				years <- M_info$BP_M - x$BP_D[i]
				curr_fam1 <- M_info$group
				curr_fam2 <- strsplit(curr_fam1, ",")[[1]][1]
				w_curr_fam2 <- which(f_a[,1]==curr_fam2)
				if ( length(w_curr_fam2) > 0 ) {
					area <- f_a[w_curr_fam2,2]
					subsistence <- f_a[w_curr_fam2,3]
				} else {
					area <- "unidentified"
					subsistence <- "unknown"
				}
				if ( area == "Papunesia" & curr_fam2 == "Austronesian" ) {
					area <- "Austronesian"
				}
				if ( area == "Papunesia" & curr_fam2 != "Austronesian" ) {
					area <- "Papuan"
				}
				years <- as.numeric(M_info[3]) - as.numeric(x[i,3])
				# add the rate information and areas to the file
				rate_GCD <- round(GCD/years, 3)
				if ( walking_distance==TRUE ) {
					rate_DD <- round(DD/years, 3)
				} else {
					rate_DD <- NA
				}
				if ( GCD > 0 ) {
					b <- argosfilter::bearing(M_lat, D_lat, M_lon, D_lon)
					if ( !is.na(b) ) {
						if ( b >= 45 & b < 135 ) { bear <- "EW" }
						if ( b < -45 & b >= -135 ) { bear <- "EW" }
						if ( b >= -45 & b < 45 ) { bear <- "NS" }
						if ( b < -135 & b >= -180 ) { bear <- "NS" }
						if ( b >= 135 & b < 180 ) { bear <- "NS" }
					} else {
						bear <- NA
					}
					biomes_out <- biomes(M_lat, M_lon, D_lat, D_lon, shape)
					biome <- biomes_out[[1]]
					percent <- biomes_out[[2]]
				} else {
					biome <- NA
					percent <- NA
					bear <- NA
				}
				cat(paste(c(M_info, x[i,], GCD, DD, rate_GCD, rate_DD, 
				  area, subsistence, bear, biome, percent), 
				  collapse="\t"), file="migration_events.txt", append=TRUE)
				cat("\n", file="migration_events.txt", append=TRUE)
			}
		}
	}
}
