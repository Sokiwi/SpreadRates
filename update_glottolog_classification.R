g <- read.table(file="glottolog_classification_strings.txt", header=TRUE, sep="\t", quote="")
asjp <- readLines("listss20_pruned.txt")
cat("", file="classification_changes_log.txt")

get_iso_code <- function(s) {
  s1 <- strsplit(s, "")[[1]][40:42]
  s2 <- paste(s1, collapse="")
  return(s2)
}

replace_classification <- function(s, new.clas) {
  s1 <- strsplit(s, "@")[[1]][2]
  s2 <- strsplit(s1, "}")[[1]][1]  # this is the old classification
  part1 <- strsplit(s, "@")[[1]][1]  # this is the text preceding the classification
  part2 <- "@"
  part3 <- "}"
  new.s <- paste(part1, part2, new.clas, part3, sep="")
  if (s2 != new.clas) {
    cat("\n", s2, "=>", "\n", new.clas, "\n", file="classification_changes_log.txt", append=TRUE)
  }
  return(new.s)
}

metalines1 <- grep("@", asjp)
metalines2 <- metalines1 + 1
for (i in 1:length(metalines2)) {
  iso <- get_iso_code(asjp[metalines2[i]])
  w_iso <- which(g$ISO==iso)
  if(length(w_iso) > 0) {
    new.clas <- g$Classification[w_iso]
  } else {
    new.clas <- ""
  }
  asjp[metalines1[i]] <- replace_classification(asjp[metalines1[i]], new.clas)
}

writeLines(asjp, con = "listss20_pruned_updated.txt")
