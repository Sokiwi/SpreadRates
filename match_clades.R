# initiate file for holding the matched data on dates and homelands
cat("clade\tpairs\tBP\tlat\tlon\tNtaxa\tmethod\n", file="dates_homelands.txt")

# data for dates, homelands, and the key to match them are read in
d <- read.table(file="dates.txt", sep="\t", header=TRUE, quote="", na.strings = "", comment.char = "")
h <- read.table(file="homelands_BT.txt", sep="\t", header=TRUE, quote="", na.strings = "", comment.char = "")
k <- read.table(file="clstrings_ids.txt", sep="\t", header=TRUE, quote="", na.strings = "", comment.char = "")  # k stands for key
md <- read.table(file="md_homelands_glottolog_groups.txt", sep="\t", header=TRUE, quote="", na.strings = "", comment.char = "")

# for the homeland data a column is added with the number of taxa in each clade,
# which will help to narrow down searches for matches
ntaxa <- function(s) {
  length(strsplit(s, "_")[[1]])
}
Ntaxa <- unname(sapply(h$glotto, ntaxa))
h <- data.frame(h, Ntaxa)

# routine for looking up a classification string in 
# md_homelands_glottolog_groups.txt in cases where results for the BT
# method could not be matched
find_md_homeland <- function(clstring) {
  w_s <- which(md$group==clstring)
  if (length(w_s) > 0) {
    md.info <- c(md$lat[w_s], md$lon[w_s], md$N[w_s])
  } else {
    md.info <- c(NA, NA, NA)
  }
  return(md.info)
}

# now go from data on dates via the key to data on homelands, matching up
# clades and outputting results to a file called dates_homelands.txt
for (i in 1:nrow(d)) {
  found <- FALSE
  clstring <- d$group[i]
  w_clstring <- which(k$clades==clstring)
  if (length(w_clstring) > 0) {
    ids <- k$ids[w_clstring]
    ids_v <- strsplit(ids, "_")[[1]]
    L <- length(ids_v)
    w_L <- which(h$Ntaxa==L)
    if (length(w_L) > 0) {
      hc <- h[w_L,]
      hc_groups <- lapply(hc$glotto, function(x) strsplit(x, "_")[[1]])
      for (j in 1:length(hc_groups)) {
        if (setequal(hc_groups[[j]], ids_v)==TRUE) {
          found <- TRUE
          relevant.info <- c(d$group[i], d$pairs[i], d$BP[i], hc$lat[j], hc$lon[j], hc$Ntaxa[j])
          cat(relevant.info[1], "\t", file="dates_homelands.txt", sep="", append=TRUE)
          cat(relevant.info[2], "\t", file="dates_homelands.txt", sep="", append=TRUE)
          cat(relevant.info[3], "\t", file="dates_homelands.txt", sep="", append=TRUE)
          cat(relevant.info[4], "\t", file="dates_homelands.txt", sep="", append=TRUE)
          cat(relevant.info[5], "\t", file="dates_homelands.txt", sep="", append=TRUE)
          cat(relevant.info[6], "\t", file="dates_homelands.txt", sep="", append=TRUE)
          cat("BT", file="dates_homelands.txt", append=TRUE)
          cat("\n", file="dates_homelands.txt", append=TRUE)
          break
        }
      }
      if (found==FALSE) {
        md.info <- find_md_homeland(clstring)
        relevant.info <- c(d$group[i], d$pairs[i], d$BP[i], md.info[1], md.info[2], md.info[3])
        cat(relevant.info[1], "\t", file="dates_homelands.txt", sep="", append=TRUE)
        cat(relevant.info[2], "\t", file="dates_homelands.txt", sep="", append=TRUE)
        cat(relevant.info[3], "\t", file="dates_homelands.txt", sep="", append=TRUE)
        cat(relevant.info[4], "\t", file="dates_homelands.txt", sep="", append=TRUE)
        cat(relevant.info[5], "\t", file="dates_homelands.txt", sep="", append=TRUE)
        cat(relevant.info[6], "\t", file="dates_homelands.txt", sep="", append=TRUE)
        cat("MD", file="dates_homelands.txt", append=TRUE)
        cat("\n", file="dates_homelands.txt", append=TRUE)
      }
    } else {
      md.info <- find_md_homeland(clstring)
      relevant.info <- c(d$group[i], d$pairs[i], d$BP[i], md.info[1], md.info[2], md.info[3])
      cat(relevant.info[1], "\t", file="dates_homelands.txt", sep="", append=TRUE)
      cat(relevant.info[2], "\t", file="dates_homelands.txt", sep="", append=TRUE)
      cat(relevant.info[3], "\t", file="dates_homelands.txt", sep="", append=TRUE)
      cat(relevant.info[4], "\t", file="dates_homelands.txt", sep="", append=TRUE)
      cat(relevant.info[5], "\t", file="dates_homelands.txt", sep="", append=TRUE)
      cat(relevant.info[6], "\t", file="dates_homelands.txt", sep="", append=TRUE)
      cat("MD", file="dates_homelands.txt", append=TRUE)
      cat("\n", file="dates_homelands.txt", append=TRUE)
    }
  } else {
    md.info <- find_md_homeland(clstring)
    relevant.info <- c(d$group[i], d$pairs[i], d$BP[i], md.info[1], md.info[2], md.info[3])
    cat(relevant.info[1], "\t", file="dates_homelands.txt", sep="", append=TRUE)
    cat(relevant.info[2], "\t", file="dates_homelands.txt", sep="", append=TRUE)
    cat(relevant.info[3], "\t", file="dates_homelands.txt", sep="", append=TRUE)
    cat(relevant.info[4], "\t", file="dates_homelands.txt", sep="", append=TRUE)
    cat(relevant.info[5], "\t", file="dates_homelands.txt", sep="", append=TRUE)
    cat(relevant.info[6], "\t", file="dates_homelands.txt", sep="", append=TRUE)
    cat("MD", file="dates_homelands.txt", append=TRUE)
    cat("\n", file="dates_homelands.txt", append=TRUE)
  }
}
