cat("", file="clstrings_ids.txt")

clades <- c()
ids <- c()
count <- 0
clstrings <- read.table(file="glottolog_classification_strings.txt", header=TRUE, sep="\t", quote="")
ustrings <- unique(clstrings$Classification)
for (i in 1:length(ustrings)) {
  L <- length(strsplit(ustrings[i], ",")[[1]])
  for (j in L:1) {
    count <- count + 1
    s <- paste(strsplit(ustrings[i], ",")[[1]][1:j], collapse=",")
    clades[count] <- s
  }
}
uclades <- unique(clades)
for (i in 1:length(uclades)) {
  s <- paste("^", uclades[i], sep="")
  w_c <- grep(s, clstrings$Classification)
  ids[i] <- paste(clstrings$GlottologID[w_c], collapse="_")
}
clstrings_ids <- data.frame(uclades, ids)
names(clstrings_ids) <- c("clades", "ids")
clstrings_ids <- clstrings_ids[order(clstrings_ids$clades),]
write.table(clstrings_ids, file="clstrings_ids.txt", sep="\t", quote=FALSE, row.names=FALSE)
