library(TreeTools)  # as.Newick
library(ape)  # read.nexus, which.edge, node.depth.edgelength
library(caper)  # clade.members.list
library(TreeTools)  # NewickTree
library(phytools)  # getDescendants
library(sp)  # SpatialPoints
library(spatstat)  # density.ppp
library(spatstat.geom)  # ppp
library(assertthat)  # assert_that
library(dplyr)  # inner_join
library(purrr)  # map_dfr
library(igraph)  # graph_from_data_frame

# parse the ...AncStates.txt files of the BT output

# get a list of family names corresponding to their
# indices in the BT output folder (BT_outfiles)
getfams <- function() {
  glotto.info <- read.csv("languoid.csv", header=TRUE)
  outfiles <- dir("BT_outfiles")
  fams <- c()
  for (i in 1:length(outfiles)) {
    glottocode <- strsplit(outfiles[i], "_")[[1]][1]
    w_id <- which(glotto.info$id==glottocode)
    fams[i] <- glotto.info$name[w_id]
  }
  return(fams)
}

#  makes vectors out of a line in the tree info beginning of the BT outfile
tab_delim_text_to_vec <- function(textline) {
  vec <- strsplit(textline, "\\\t")[[1]]
  vec <- vec[3:length(vec)]
}

get_languages <- function(all.text) {
  textline <- all.text[1]
  lgs <- tab_delim_text_to_vec(textline)
  return(lgs)
}

get_phylo <- function(filename) {
  s1 <- strsplit(filename, "\\/")[[1]][2]
  famcode <- strsplit(s1, "_")[[1]][1]
  nexfile <- paste0("nexfiles/", famcode, ".nex")
  phylo <- read.nexus(file=nexfile)
}

sanity.check <- function(lgs,phylo) {
  ok <- identical(sort(lgs),sort(phylo$tip.label))
  if (ok) {
    out <- "BT output files and the tree used here have the same languages"
  }
  if (!ok) {
    out <- "WARNING: BT output files and the tree used here do not have the same languages"
  }
  return(out)
}  

get_newick <- function(phylo) {
  phylo2 <- phylo
  phylo2$edge.length <- NULL
  newick <- NewickTree(phylo2)  # function from TreeTools
  return(newick)
}

# match clade identifiers in phylo and BT output
# add a column with names of taxa
match_nodes <- function(phylo, all.text, filename) {
  BT.tree.info <- all.text[grep("^Node", all.text)]
  BT.clades <- lapply(BT.tree.info, tab_delim_text_to_vec)
  BT.nodes <- unlist(lapply(BT.tree.info, function(x) strsplit(x, "\t")[[1]][1]))
  # change - to _ in node names, for later compatibility
  BT.nodes <- unlist(lapply(BT.nodes, function(x) gsub("-", "_", x)))
  Nl <- length(BT.clades[[1]])  # no. of languages
  Nc <- length(BT.clades)  # no. of clades, including individual languages
  Ne <- nrow(phylo$edge)  # no. of edges
  if (Nc - 1 == Ne) {
    cat("Number of clades (", Nc, ") and edges (", Ne, ") as expected\n", sep="")
  } else {
    cat("Number of clades (", Nc, ") and edges (", Ne, ") unexpected\n", sep="")
  }
  # make a data frame to hold matches between clade identifiers
  clades.tmp <- as.data.frame(matrix(NA, nrow=Nc, ncol=6))
  names(clades.tmp) <- c("glotto", "name", "phylo", "BT", "lat", "lon")
  clades.tmp$phylo <- sort(unique(as.vector(phylo$edge)))
  # get clade members numbered as in the phylo object
  cml <- list()
  for (i in 1:nrow(clades.tmp)) {
    cml[[i]] <- getDescendants(phylo, i)
  }
  # get clade members in phylo object named
  cmn <- list()  # clade members named
  for (i in 1:length(cml)) {
    taxa.incl.na <- phylo$tip.label[cml[[i]]]
    taxa.excl.na <- as.vector(na.omit(taxa.incl.na))
    cmn[[i]] <- taxa.excl.na
  }
  # match the clade members in cmn with clade members in
  # the BT clades, thereby matching BT and phylo clades
  for (i in 1:length(cmn)) {
    for (j in 1:length(BT.clades)) {
      if (setequal(cmn[[i]], BT.clades[[j]])) {
        clades.tmp$BT[i] <- BT.nodes[j]
        clades.tmp$glotto[i] <- paste(BT.clades[[j]], collapse="_")
      }
    }
  }
  # add names of languages and the proto-family
  # but for subgroups only ""
  glotto.info <- read.csv("languoid.csv", header=TRUE)
  for (i in 1:nrow(clades.tmp)) {
    if(length(grep("_", clades.tmp$glotto[i]))==0) {
      w_id <- which(glotto.info$id==clades.tmp$glotto[i])
      clades.tmp$name[i] <- glotto.info$name[w_id]
    } else if (clades.tmp$BT[i]=="Node_00000") {
      s1 <- strsplit(filename, "/")[[1]][2]
      glottocode <- strsplit(s1, "_")[[1]][1]
      w_fam <- grep(glottocode, glotto.info$id)
      fam.name <- glotto.info$name[w_fam]
      clades.tmp$name[i] <- fam.name 
    } else {
      clades.tmp$name[i] <- ""
    }
  }    
  return(clades.tmp)
}

# get the best supported single locations for each node and
# add that to the clades.tmp data frame
get_clades_dots <- function(clades.tmp, all.text) {
  BT.node.text <- all.text[-grep("^Node", all.text)]
  writeLines(BT.node.text, con="tmp.txt")
  BT.ni <- read.table(file="tmp.txt", header=TRUE, sep="\t") # BT node info
  BT.ni <- BT.ni[,-ncol(BT.ni)]
  names(BT.ni) <- unlist(lapply(names(BT.ni), function(x) gsub("\\.", "_", x)))
  file.remove("tmp.txt")
  w_max_lh <- which(BT.ni$Lh==max(BT.ni$Lh))[1]
  for (i in 1:nrow(clades.tmp)) {
    column.name.lat <- paste0(clades.tmp$BT[i], "___Lat")
    column.name.lon <- paste0(clades.tmp$BT[i], "___Long")
    klat <- which(names(BT.ni)==column.name.lat)
    klon <- which(names(BT.ni)==column.name.lon)
    clades.tmp$lat[i] <- BT.ni[w_max_lh,klat]
    clades.tmp$lon[i] <- BT.ni[w_max_lh,klon]
  }
  return(clades.tmp)
}

# the following loops over BT output files collecting homeland data
# this is then put in a file called homelands_BT.txt
cat("glotto\tname\tphylo\tBT\tlat\tlon\n", file="homelands_BT.txt")
for (i in 1:length(dir("BT_outfiles"))) {
  if ("fams" %in% ls()==FALSE) fams <- getfams()  # only needed once
  filename <- paste0("BT_outfiles/", dir("BT_outfiles")[i])
  all.text <- readLines(filename)
  lgs <- get_languages(all.text)
  phylo <- get_phylo(filename)
  sanity <- sanity.check(lgs,phylo)
  newick <- get_newick(phylo)
  # output of the following only used to feed the next 
  clades.tmp <- match_nodes(phylo, all.text, filename)
  clades <- get_clades_dots(clades.tmp, all.text); rm(clades.tmp)
  write.table(clades, file="homelands_BT.txt", sep="\t", quote=FALSE, 
              row.names=FALSE, col.names=FALSE, append=TRUE)
}

