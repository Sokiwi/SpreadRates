# read information about dates and homelands
dh <- read.table(file="dates_homelands.txt", sep="\t", header=TRUE, quote="")
cat("mo\tpairs_mo\tBP_mo\tlat_mo\tlon_mo\tNtaxa_mo\tmethod_mo\tda\tpairs_da\tBP_da\tlat_da\tlon_da\tNtaxa_da\tmethod_da\n", 
    file="diffusion_events.txt")

# function for extracting the mother in a classification string
extract_mo <- function(clade) {
  s1 <- strsplit(clade, ",")[[1]]
  L <- length(s1)
  s2 <- s1[-L]
  s3 <- paste(s2, collapse=",")
  return(s3)
}

for (i in 1:nrow(dh)) {
  clade <- dh$clade[i]
  if(length(grep(",", clade)) > 0) {
    mo <- extract_mo(clade)
    w_mo <- which(dh$clade==mo)
    if (length(w_mo) > 0) {
      cat(dh[w_mo,1], "\t", dh[w_mo,2], "\t", dh[w_mo,3], "\t", 
          dh[w_mo,4], "\t", dh[w_mo,5], "\t", dh[w_mo,6], "\t", 
          dh[w_mo,7], "\t", dh[i,1], "\t", dh[i,2], "\t", 
          dh[i,3], "\t", dh[i,4], "\t", dh[i,5], "\t", 
          dh[i,6], "\t", dh[i,7], "\n", file="diffusion_events.txt", 
          sep="", append=TRUE)
    }
  }
}
