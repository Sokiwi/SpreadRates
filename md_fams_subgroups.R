library(fields)

## routine for finding identical classification strings and getting homelands
hd <- read.table(file="glottolog_classification_strings.txt", header=TRUE, sep="\t", 
  na.strings="", quote="", comment.char="")  # hd for homeland data

# function used for geographical distances
distance <- function (lat1, lat2, lon1, lon2) 
{
    if (lat1 == lat2 & lon1 == lon2) 
        distance <- 0
    else {
        rlat1 = radian(lat1)
        rlat2 = radian(lat2)
        rlon = radian(lon2 - lon1)
        distance <- 60 * (180/pi) * acos(sin(rlat1) * sin(rlat2) + 
            cos(rlat1) * cos(rlat2) * cos(rlon))
        distance <- distance * 1852/1000
    }
    distance
}

radian <- function (degree) 
{
    radian <- degree * (pi/180)
    radian
}

# change ( to ° and ) to § since they are reserved strings in R
from_parens <- function(sa) {
	sa1 <- gsub("\\(", "°", sa)
	sa2 <- gsub("\\)", "§", sa1)
	return(sa2)
}

# change ° and § back again later
to_parens <- function(sb) {
	sb1 <- gsub("°", "(", sb)
	sb2 <- gsub("§", ")", sb1)
	return(sb2)
}

hd$Classification <- as.vector(unlist(sapply(hd$Classification, from_parens)))

# find and eliminate unwanted languages
unwanted <- c("Artificial_Language", "Bookkeeping", "Pidgin",
  "Speech Register", "Sign Language", "Unattested", "Unclassifiable")
w_unwanted <- which(hd$Classification %in% unwanted)
if ( length(w_unwanted) > 0 ) {
	hd <- hd[-w_unwanted,]
}

w_nas <- which(is.na(hd$lat))
if ( length(w_nas) > 0 ) {
	hd <- hd[-w_nas,]
}

w_nas_string <- which(hd$lat=="NA")
if ( length(w_nas_string) > 0 ) {
	hd <- hd[-w_nas_string,]
}

w_pidgin <- grep("Pidgin", hd$Classification)
if ( length(w_pidgin) > 0 ) {
	hd <- hd[-w_pidgin,]
}

w_Unattested <- grep("^Unattested", hd$Classification)
if ( length(w_Unattested) > 0 ) {
	hd <- hd[-w_Unattested,]
}

w_Artificial_Language <- grep("Artificial_Language", hd$Classification)
if ( length(w_Artificial_Language) > 0 ) {
	hd <- hd[-w_Artificial_Language,]
}

w_Bookkeeping <- grep("Bookkeeping", hd$Classification)
if ( length(w_Bookkeeping) > 0 ) {
	hd <- hd[-w_Bookkeeping,]
}

w_Mixed_Language <- grep("Mixed_Language", hd$Classification)
if ( length(w_Mixed_Language) > 0 ) {
	hd <- hd[-w_Mixed_Language,]
}

w_Sign_Language <- grep("Sign_Language", hd$Classification)
if ( length(w_Sign_Language) > 0 ) {
	hd <- hd[-w_Sign_Language,]
}

w_Speech_Register <- grep("Speech_Register", hd$Classification)
if ( length(w_Speech_Register) > 0 ) {
	hd <- hd[-w_Speech_Register,]
}

# for each unique classification string generate all the subgroups that it contains
# and make a classification string for each
# put this set of 'folded-out' classification strings in a vector;
# that vector, when made to only contain unique strings will contain all
# the possible family and subgroup classification strings in Glottolog
groups <- unique(hd$Classification)
all_c_strings <- c()

# function for generating subgroups
unfold_subgroups <- function(s) {
	subgroups <- c()
	s1 <- strsplit(s, ",")[[1]]
	L <- length(s1)
	if ( L==1 ) {
		subgroups <- c(subgroups, s1)
	} else {
		subgroups <- c(subgroups, s1[1])
		for (i in 2:L) {
			s2 <- paste(s1[1:i], collapse=",")
			subgroups <- c(subgroups, s2)
		}
	}
	return(subgroups)
}

for (i in 1:length(groups)) {
	unfold_subgroups_out <- unfold_subgroups(groups[i])
	all_c_strings <- c(all_c_strings, unfold_subgroups_out)
	all_c_strings <- unique(all_c_strings)
}
all_c_strings <- sort(all_c_strings)


# now search for each of all the possible families and subgroups
# in the Glottolog language set, make a geographical distance
# matrix for each (using the distance() function for Great
# Circle Distance), and get the homeland as the location with
# the smallest average distance to all the other languages

cat("group\tlat\tlon\tN\n", file="md_homelands_glottolog_groups.txt")


for (i in 1:length(all_c_strings)) {
	cat("doing ", i, " out of ", length(all_c_strings), ": ", all_c_strings[i], "\n", sep="")
	# if the classification string has no commas, it is a family, and so
	# compute all distances
	if ( length(grep(",", all_c_strings[i]))==0 ) {
		lat_lon <- hd[grep(paste("^", all_c_strings[i], sep=""), hd$Classification), c(4,5,2)]
		L <- nrow(lat_lon)
		if ( L==1 ) {
			home_lat  <- lat_lon[1,1]
			home_lon <- lat_lon[1,2]
		} else if ( L==2) {
			lat_lon <- lat_lon[order(lat_lon$GlottologID),]
			home_lat  <- lat_lon[1,1]
			home_lon <- lat_lon[1,2]
		} else {
			lat_lon <- lat_lon[order(lat_lon$GlottologID),]
			m <- lat_lon[,c(2,1)]
			m <- as.matrix(m)
			m <- apply(m, c(1,2), as.numeric)
			distances <- rdist.earth(m, m, miles=FALSE)
			means <- apply(distances, 2, mean)
			w_min <- which(means==min(means))[1]
			home_lat <- lat_lon[w_min,1]
			home_lon <- lat_lon[w_min,2]
		}
		cat(to_parens(all_c_strings[i]), "\t", home_lat, "\t", home_lon, "\t", L, "\n", file="md_homelands_glottolog_groups.txt", sep="", append=TRUE)
	# if the classification string has commas, it is a subgroup, and so use distances
	# already computed
	} else {
		lat_lon_sub <- hd[grep(paste("^", all_c_strings[i], sep=""), hd$Classification), c(4,5,2)]
		L <- nrow(lat_lon_sub)
		if ( L==1 ) {
			home_lat  <- lat_lon_sub[1,1]
			home_lon <- lat_lon_sub[1,2]
		} else if ( L==2) {
			lat_lon_sub <- lat_lon_sub[order(lat_lon$GlottologID),]
			home_lat  <- lat_lon_sub[1,1]
			home_lon <- lat_lon_sub[1,2]
		} else {
			# subset the family distance matrix
			w_sub <- match(row.names(lat_lon_sub), row.names(distances))
			distances_sub <- distances[w_sub, w_sub]
			means <- apply(distances_sub, 2, mean)
			w_min <- which(means==min(means))[1]
			name <- names(w_min)
			home_lat <- lat_lon[name,1]
			home_lon <- lat_lon[name,2]
		}
		cat(to_parens(all_c_strings[i]), "\t", home_lat, "\t", home_lon, "\t", L, "\n", file="md_homelands_glottolog_groups.txt", sep="", append=TRUE)
	}
}
