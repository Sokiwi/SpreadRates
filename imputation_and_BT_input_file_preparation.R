## LIBRARIES
# libraries needed
library(ape)  # read.tree(), write.tree(), drop.tip(), write.nexus(), di2multi()
# see https://ladal.edu.au/phylo
library(phangorn)  # nnls.tree
library(stringr)  # str_split, str_sub
library(phytools)  # force.ultrametric
# library(ade4)  # mantel.rtest

## LOCALLY SAVED SCRIPTS
# from the glottoTrees library relevant functions have been copied into 
# glottoTrees_functions.R, which obviates the need for an installation
# from GitHub and which somehow works better than using the installed version
# the relevant function is keep_as_tip(), but it calls other functions also
# included in glottoTrees_functions.R
source("glottoTrees_functions.R")

## SOFTWARE
# Download BayesTraits 4.0.1 from
# http://www.evolution.reading.ac.uk/BayesTraitsV4.0.1/BayesTraitsV4.0.1.html

## DATA PREPARATION PRIOR TO THIS SCRIPT
# this study uses version 20 of ASJP
# unlike the dating procedure it uses the whole database
# Using the interactive ASJP software at https://github.com/Sokiwi/InteractiveASJP02
# produce a tab-delimited version of the database; here supplied for convenience
# as listss20_formatted.tab

## DOWNLOADING FILES
# this study uses version 5.1 of Glottolog
# download the file tree_glottolog_newick.txt from https://glottolog.org/meta/downloads
# by right-clicking and downloading; if opened in a browser and then copied
# to a file the encoding may be upset
# download asjp-v20.zip from https://zenodo.org/record/7079637/files/lexibank/asjp-v20.zip
# unzip asjp-v20.zip and retain only the file "languages.csv" 
# download languages_and_dialects_geo.csv from https://glottolog.org/meta/downloads
# download languoid.csv from https://glottolog.org/meta/downloads
# the above download may not work, and will then have to be done manually

## READ FILES, PREPARE DATA
# Prepare Glottolog trees, reading tree_glottolog_newick.txt
trees <- suppressWarnings(readLines("tree_glottolog_newick.txt"))
# two apostrophes in names of taxa symbolize one apostrophe; these are
# converted to "§" here
for (i in 1:length(trees)) {
  trees[i] <- gsub("\'\'\'", "'§", trees[i])
}
for (i in 1:length(trees)) {
  trees[i] <- gsub("\'\'", "§", trees[i])
}

# Prepare index to Glottolog tree data
fam.names.trees <- c()
for (i in 1:length(trees)) {
  fam.names.trees[i] <- read.tree(text=trees[i])$node.label[1]
}
# Extract the family names from fam.names.trees
efn <- function(x) {
  a1 <- strsplit(x, "\\'")[[1]][2]
  a2 <- strsplit(a1, " \\[")[[1]][1]
}
tree.index <- unlist(lapply(fam.names.trees, efn))
# Extract the family glottocodes from fam.names.trees
efg <- function(x) {
  a1 <- strsplit(x, "\\[")[[1]][2]
  a2 <- strsplit(a1, "\\]")[[1]][1]
}
tree.index.gcodes <- unlist(lapply(fam.names.trees, efg))

# Prepare ASJP metadata
# revise asjp_meta such that a dialect or family glottocode match to an ASJP
# doculect is converted to a language-level languoid
# takes a couple of minutes
unify.level <- function(asjp_meta) {
  # having gotten the ASJP languages.csv file with metadata now read it
  asjp_meta <- read.csv("languages.csv")
  languoid <- read.csv("languoid.csv")
  for (i in 1:nrow(asjp_meta)) {
    gcode_old <- asjp_meta$Glottocode[i]
    if (gcode_old != "") {
      w_gcode <- grep(gcode_old, languoid$id)
      if(languoid$level[w_gcode]=="dialect") {
        parent1 <- languoid$parent_id[w_gcode]
        w_parent1 <- which(languoid$id==parent1)
        if(languoid$level[w_parent1]=="language") {
          asjp_meta$Glottocode[i] <- languoid$id[w_parent1]
        } else if(languoid$level[w_parent1]=="dialect") {
          parent2 <- languoid$parent_id[w_parent1]
          w_parent2 <- which(languoid$id==parent2)
          if(languoid$level[w_parent2]=="language") {
            asjp_meta$Glottocode[i] <- languoid$id[w_parent2]
          } else if(languoid$level[w_parent2]=="dialect") {
            parent3 <- languoid$parent_id[w_parent2]
            w_parent3 <- which(languoid$id==parent3)
            if(languoid$level[w_parent3]=="language") {
              asjp_meta$Glottocode[i] <- languoid$id[w_parent3]
            } else if (languoid$level[w_parent3]=="dialect") {
              parent4 <- languoid$parent_id[w_parent3]
              w_parent4 <- which(languoid$id==parent4)
              if(languoid$level[w_parent4]=="language") {
                asjp_meta$Glottocode[i] <- languoid$id[w_parent4]
              } else if (languoid$level[w_parent3]=="dialect") {
                cat("at", i, "there is more than 4 levels of dialect")
              }
            }
          }
        }
      } else if(languoid$level[w_gcode]=="family") {
        if(gcode_old=="alba1267") {
          # Standard Albanian fixed by hand
          asjp_meta$Glottocode[i] <- "tosk1239"
        } else {
          # other cases fixed arbitrarily as the first encountered
          # daughter language
          asjp_meta$Glottocode[i] <-
            languoid$id[which(languoid$parent_id==gcode_old
                              & languoid$level=="language")][1]
        }
      }
    }
  }
  return(asjp_meta)
}
asjp_meta <- unify.level(asjp_meta)
# save(asjp_meta, file="asjp_meta.RData")
# load("asjp_meta.RData")

# add a column to asjp_meta with counts of attested words
asjp <- read.table(file="listss20_formatted.tab", header=TRUE, sep="\t",
                   quote="", na.strings="", comment.char="")
forty <- c("I", "you", "we", "one", "two", "person", "fish", "dog", "louse",
           "tree", "leaf", "skin", "blood", "bone", "horn", "ear",
           "eye", "nose", "tooth", "tongue", "knee", "hand", "breast",
           "liver", "drink", "see", "hear", "die", "come", "sun", "star",
           "water", "stone", "fire", "path", "mountain", "night", "full",
           "new", "name")

attestations_asjp <- 40 - apply(asjp[,forty], 1, function(x) length(which(x=="XXX")))
m1 <- match(asjp_meta$Name, asjp$names)
attestations <- attestations_asjp[m1]
asjp_meta <- data.frame(asjp_meta, attestations)

#  prune Glottolog trees so as to only contain "language"-level taxa
languoid <- read.csv("languoid.csv")
languoid <- languoid[languoid$level=="language",]
# the following excludes isolates, which have the empty string as family_id
# and non-genuine families (artificial languages etc.)
languoid <- languoid[!languoid$family_id%in%c("", "arti1236", "book1242", 
                                              "mixe1287", "pidg1258", "spee1234", 
                                              "unat1236", "uncl1493", "sign1238"),]
# furthermore families with less than 3 members are excluded
glotfams <- names(table(languoid$family_id)[table(languoid$family_id) >= 3])
m <- match(glotfams, tree.index.gcodes)
trees.select <- trees[m]
trees.select.gcodes <- tree.index.gcodes[m]
trees.select.names <- tree.index[m]
languoid <- languoid[which(languoid$family_id%in%trees.select.gcodes),]
trees.select.pruned <- c()

# make patristic matrices from the Glottolog trees
pat.mat <- function() {
  patristic.matrices <- list()
  L <- length(trees.select)
  for (i in 1:L) {
    if (i%%10==0) {
      cat("doing", i, "out of", L, "\n")
    }
    phylo <- read.tree(text=trees.select[i])
    efg <- function(x) {
      a1 <- strsplit(x, "\\[")[[1]][2]
      a2 <- strsplit(a1, "\\]")[[1]][1]
    }
    phylo$tip.label <- unlist(lapply(phylo$tip.label, efg))
    phylo$node.label <- unlist(lapply(phylo$node.label, efg))
    # the keep_as_tip function from https://github.com/erichround/glottoTrees 
    # mysteriously didn't initially work. Seems best to make it part of
    # this script
    pruned.phylo <- keep_as_tip(phylo, languoid$id[which(languoid$family_id==trees.select.gcodes[i])])
    # make non-ultrametric Glottolog trees ultrametric
    if (is.ultrametric(phylo)==FALSE) {
      pruned.phylo <- force.ultrametric(pruned.phylo, method="extend")
    }
    pruned.newick <- write.tree(pruned.phylo)
    trees.select.pruned[i] <- pruned.newick
    patristic.matrices[[i]] <- cophenetic(pruned.phylo)
  }
  names(patristic.matrices) <- trees.select.gcodes
  return(patristic.matrices)
}
patristic.matrices <- pat.mat()  # takes less than a minute
# save(patristic.matrices, file="patristic_matrices.RData")
# load("patristic_matrices.RData")

# create LDN matrices for the best attested ASJP doculects per Glottolog
# language
# first reduce ASJP to doculects with at least 21 attested words
# and then to the best attested glottocode representative
m1 <- match(asjp_meta$Name[which(asjp_meta$attestations > 20)], asjp$names)
asjp <- asjp[m1,]

selected.doculects <- rep(NA, length(languoid$id))
for (i in 1:length(languoid$id)) {
  w_g_lg <- which(asjp_meta$Glottocode==languoid$id[i])
  if (length(w_g_lg) > 0) {
    o <- order(asjp_meta$attestations[w_g_lg])
    w_max <- match(max(o), o)
    candidate <- asjp_meta$Name[w_g_lg[w_max]]
    if (candidate %in% asjp$names) {
      selected.doculects[i] <- candidate
    }
  }
}
s.d <- selected.doculects[which(!is.na(selected.doculects))]
m2 <- match(s.d, asjp$names)
asjp <- asjp[m2,]
languoid <- data.frame(languoid, selected.doculects)

fd <- list()  # stands for families and doculects
u <- unique(languoid$family_id)
for (i in 1:length(u)) {
  w_u <- which(languoid$family_id==u[i])
  doculects <- languoid$selected.doculects[w_u]
  w_na <- which(is.na(doculects))
  if (length(w_na) > 0) {
    doculects <- doculects[-w_na]
  }
  fd[[i]] <- doculects
}
names(fd) <- u

prepare_data <- function(data_all) {
  data <- data_all[,c("names",forty)]
  # replace NAs with XXX
  XXX <- function(x) {
    if ( is.na(x) ) {
      x <- "XXX"
    }
    return(x)
  }
  data[,2:41] <- apply(data[,2:41], c(1,2), XXX)
  # cat("getting rid of the loanword symbol...\n")
  # get rid of the loanword symbol
  eliminate.loans <- function(x) {
    where.loan <- grep("%", x)
    if(length(where.loan) > 0) {
      for (i in 1:length(where.loan)) {
        x[where.loan[i]] <- gsub("%", "", x[where.loan[i]])
      }
    }
    return(x)
  }
  data[,2:41] <- apply(data[,2:41], c(1,2), eliminate.loans)
  # trim whitespace
  # cat("trimming whitespace...\n")
  trim <- function(x) {
    whites <- which(strsplit(x, "")[[1]]==" ")
    if( length(whites) > 0 ) {
      x <- paste(strsplit(x, " ")[[1]][1:(length(whites)+1)],collapse="")
    } else {
      x <- x
    }
    return(x)
  }
  data <- apply(data, c(1,2), trim)
  # get rid of more than two synonyms
  # cat("getting rid of more than two synonyms...\n")
  synred <- function(x) {
    if( length(strsplit(x,",")[[1]]) > 2 ) {
      x <- paste(strsplit(x,",")[[1]][1:2],collapse=",")
    } else {
      x <- x
    }
    return(x)
  }
  data <- apply(data, c(1,2), synred)
  # get a vector of unique symbols, including strings
  # to be treated as one symbol
  # first a function for splitting words, sw
  sw <- function(x) {
    ws <- strsplit(x,"")[[1]]
    # QUOTES
    wq <- grep("\"",ws)
    if ( length(wq) > 0 ) {
      for (i in 1:length(wq)) {
        ws[wq[i]-1] <- paste(ws[wq[i]-1],ws[wq[i]],sep="")
      }
      ws <- ws[-wq]
    }
    # STARS
    wst <- grep("\\*",ws)
    if ( length(wst) > 0 ) {
      for (i in 1:length(wst)) {
        ws[wst[i]-1] <- paste(ws[wst[i]-1],ws[wst[i]],sep="")
      }
      ws <- ws[-wst]
    }
    # TILDE
    wt <- grep("\\~",ws)
    if ( length(wt) > 0 ) {
      for (i in 1:length(wt)) {
        ws[wt[i]-2] <- paste(ws[wt[i]-2],ws[wt[i]-1],sep="")
      }
      ws <- ws[-c(wt,wt-1)]
    }
    # DOLLAR
    wd <- grep("\\$",ws)
    if ( length(wd) > 0 ) {
      for (i in 1:length(wd)) {
        ws[wd[i]-3] <- paste(ws[wd[i]-3],ws[wd[i]-2],ws[wd[i]-1],sep="")
      }
      ws <- ws[-c(wd,wd-1,wd-2)]
    }
    return(unique(ws))
  }
  l <- lapply(data[,2:41],sw)
  uni <- unique(unlist(l))
  rm(l)
  # the following is a trick to get the comma in position 44
  # this will make the value of the comma stay as a comma after being
  # transformed to a number in the list of unique symbols and then
  # transformed to unicode
  L_uni <- length(uni)
  wcomma <- which(uni==",")
  if ( length(wcomma) > 0 ) {
    if ( L_uni > 44 ) {
      uni <- uni[-wcomma]
      uni <- c(uni[1:43],",",uni[44:length(uni)])
    } else if ( L_uni <= 44 ) {
      uni <- uni[-wcomma]
      uni <- c(uni[1:43],",")
    }
  }
  # the following is a trick to get the X in position 88
  # this will make the value of the X stay as an X after being
  # transformed to a number in the list of unique symbols and then
  # transformed to unicode
  # in order not to disturb the position of the comma
  # an occurrence of X before 44 is replaced by the dummy _
  L_uni <- length(uni)
  wX <- which(uni=="X")
  if ( length(wX) > 0 ) {
    if ( L_uni > 87 ) {
      uni[wX] <- "_"
      uni <- c(uni[1:87],"X",uni[88:length(uni)])
    } else if ( L_uni <= 87 ) {
      uni[wX] <- "_"
      uni <- c(uni[1:87],"X")
    }
  }
  assign("uni", uni, envir=.GlobalEnv)
  # translate sound-encoding sequences to numbers
  # then to utf8 and then to words
  # cat("transforming sound representations to single utf8 symbols...\n")
  transf <- function(x) {
    # first split the string into individual symbols
    ws <- strsplit(x,"")[[1]]  # stands for word split
    # now start joining the ones that should be treated as
    # single symbols
    # find quotes and join them with their preceding symbol
    wq <- grep("\"",ws)  # stands for where quote
    if ( length(wq) > 0 ) {
      for (i in 1:length(wq)) {
        ws[wq[i]-1] <- paste(ws[wq[i]-1],ws[wq[i]],sep="")
      }
      ws <- ws[-wq]
    }
    # find stars and join them with their preceding symbol
    wst <- grep("\\*",ws)  # stands for where star
    if ( length(wst) > 0 ) {
      for (i in 1:length(wst)) {
        ws[wst[i]-1] <- paste(ws[wst[i]-1],ws[wst[i]],sep="")
      }
      ws <- ws[-wst]
    }
    # find tildes and join them with their preceding symbol
    wt <- grep("\\~",ws)  # stands for where tilde
    if ( length(wt) > 0 ) {
      for (i in 1:length(wt)) {
        ws[wt[i]-2] <- paste(ws[wt[i]-2],ws[wt[i]-1],sep="")
      }
      ws <- ws[-c(wt,wt-1)]
    }
    # find dollar signs and join them with their two preceding symbols
    wd <- grep("\\$",ws)  # stands for where dollar
    if ( length(wd) > 0 ) {
      for (i in 1:length(wd)) {
        ws[wd[i]-3] <- paste(ws[wd[i]-3],ws[wd[i]-2],ws[wd[i]-1],sep="")
      }
      ws <- ws[-c(wd,wd-1,wd-2)]
    }
    # for each symbol or symbol combination in the word
    # look up its place in the vector of unique symbols;
    # its place there is a number which is subsequently
    # transformed into one of the more than 1 million unicode symbols
    for (i in 1:length(ws)) {
      ws[i] <- which(uni==ws[i])
    }
    # the last value, that of out, will not be printed but if the
    # output is assigned to a variable this value will pass on
    out <- intToUtf8(ws)
  }
  data[,2:41] <- apply(data[,2:41], c(1,2), transf)
  return(data)
}
data <- prepare_data(asjp) # takes a couple of minutes
# save(data, file="data.RData")
# load("data.RData")

listwise_LDN <- function(A, B, data) {
  # reduce the matrix to the two relevant languages
  whereA <- which(data[,1]==A)
  whereB <- which(data[,1]==B)
  data.red <- data[c(whereA, whereB),2:41]
  # do LDN for individual word pairs in data.red
  LDN <- function(x,y) {
    LDN <- adist(x,y)/max(nchar(x),nchar(y))
    return(LDN)
  }
  # reduce further when something is missing for the purpose
  # of LDN between words referring to the same concept
  missing <- union(which(data.red[1,]=="XXX"),which(data.red[2,]=="XXX"))
  if(length(missing) > 0 & length(missing) < 40) {
    data.red1 <- data.red[,-missing]
  } else if ( length(missing)==0 ) {
    data.red1 <- data.red
  } else {
    cat(A, "and", B, "are non-overlapping,\n")
    cat("so the LDN could not be calculated\n")
    return()
  }
  # do average LDN for all pairs of words referring to the same concept,
  # take into account synonyms
  LDNs <- rep(0,length(data.red1[1,]))
  for (i in 1:length(data.red1[1,])) {
    L1 <- length(grep(",",data.red1[1,i]))
    L2 <- length(grep(",",data.red1[2,i]))
    if ( L1==0 & L2==0 ) {
      LDNs[i] <- LDN(data.red1[1,i],data.red1[2,i])
    } else if ( L1 > 0 & L2==0 ) {
      w1 <- strsplit(data.red1[1,i],",")[[1]][1]
      w2 <- strsplit(data.red1[1,i],",")[[1]][2]
      w3 <- data.red1[2,i]
      LDNs[i] <- mean(c(LDN(w1,w3),LDN(w2,w3)))
    } else if ( L1==0 & L2 > 0 ) {
      w1 <- data.red1[1,i]
      w2 <- strsplit(data.red1[2,i],",")[[1]][1]
      w3 <- strsplit(data.red1[2,i],",")[[1]][2]
      LDNs[i] <- mean(c(LDN(w1,w2),LDN(w1,w3)))
    } else {
      w1 <- strsplit(data.red1[1,i],",")[[1]][1]
      w2 <- strsplit(data.red1[1,i],",")[[1]][2]
      w3 <- strsplit(data.red1[2,i],",")[[1]][1]
      w4 <- strsplit(data.red1[2,i],",")[[1]][2]
      LDNs[i] <- mean(c(LDN(w1,w3),LDN(w1,w4),LDN(w2,w3),LDN(w2,w4)))
    }
  }
  LDN.lg <- mean(LDNs)
  return(LDN.lg)
}

# now run LDN matrices for the individual families;
# those that are represented by 3 or
# more doculects in ASJP will be output as NA
ldn.mat <- function() {
  ldn.matrices <- list()
  for (i in 1:length(fd)) {
    cat("doing", names(fd)[i], "which is", i, "out of", length(fd), "\n")
    L <- length(fd[[i]])
    if (L >= 3) {
      m <- matrix(0, nrow=L, ncol=L, dimnames=list(fd[[i]],fd[[i]]))
      for (j in 1:(L-1)) {
        for (k in (j+1):L) {
          m[k,j] <- listwise_LDN(fd[[i]][j],fd[[i]][k],data)
        }
      }
      m <- m + t(m)
      ldn.matrices[[i]] <- m
    } else {
      ldn.matrices[[i]] <- NA
    }
  }
  names(ldn.matrices) <- names(fd)
  return(ldn.matrices)
}
ldn.matrices <- ldn.mat()
# save(ldn.matrices, file="ldn_matrices.RData")
# load("ldn_matrices.RData")

# function for aligning a patristic Glottolog matrix (pm) and an LDN matrix (lm)
# - when doing these comparisons remember that they sit in different
# places in the lists patristic.matrices and ldn.matrices, but can be
# called by their Glottolog family ID's, which are the same,
# the object glotfams can be used; this is ordered like patristic.matrices
align.matrices <- function(pm, lm) {
  lm.asjp.names <- colnames(lm)
  m <- match(lm.asjp.names, asjp_meta$Name)
  colnames(lm) <- rownames(lm) <- asjp_meta$Glottocode[m]
  lm.new <- apply(pm, c(1,2), function(x) x <- NA)
  diag(lm.new) <- 0
  lm.new[colnames(lm),colnames(lm)] <- lm[colnames(lm),colnames(lm)]
  return(lm.new)
}

# align  patristic Glottolog matries and LDN matrices
ext.ldn.matrices <- list()
for (i in 1:length(glotfams)) {
  if ( length(ldn.matrices[[glotfams[i]]]) > 1 ) {  # excludes NA matrices
    ext.ldn.matrices[[i]] <- align.matrices(patristic.matrices[[glotfams[i]]], ldn.matrices[[glotfams[i]]])
  } else {
    ext.ldn.matrices[[i]] <- NA
  }
}
names(ext.ldn.matrices) <- glotfams

## exclude families with no more than three levels in the classification,
## i.e. with unique patristic distances being 3 or less
## this will exlude all families that have at most one patristic distance
## with associated LDN values;
## the 86 families thus excluded have 3-11 languages
#excluded <- c()
#N <- c()
#for (i in 1:length(patristic.matrices)) {
#  if (length(unique(as.vector(patristic.matrices[[i]]))) <= 4) {
#    excluded <- c(excluded, i)
#    N <- c(N, nrow(patristic.matrices[[i]]))
#  }
#}

#elm.red <- ext.ldn.matrices[-excluded]
#pam.red <- patristic.matrices[-excluded]
#glotfams.red <- glotfams[-excluded]

# the following just done because subsequent code refers to 
# objects thus named; above code may eventually be
# used to exclude some families in which case the
# three lines below should go out
elm.red <- ext.ldn.matrices
pam.red <- patristic.matrices
glotfams.red <- glotfams

# do imputations of missing LDNs based on a function fitted to the entire
# PAM-LDN correlation for a family
imputation <- function() {
  imput <- list()
  for (i in 1:length(glotfams.red)) {
    cat("doing", i, "\n")
    elm <- elm.red[[i]]
    # imputation not done if there aren't missing points
    if (length(which(is.na(elm))) > 1) {
      pam <- pam.red[[i]]
      elm.data <- elm[lower.tri(elm)]
      pam.data <- pam[lower.tri(pam)]
      u.pam <- sort(unique(pam.data))
      
      if (length(u.pam)==1) {
        imput[[i]] <- elm
        L <- nrow(elm)
        for (p in 1:(L-1)) {
          for (q in (p+1):L) {
            if (is.na(elm[p,q])) {
              imput[[i]][q,p] <- imput[[i]][p,q] <- mean(elm.data, na.rm=TRUE)
            }
          }
        }
      }
      if (length(u.pam)==2) {
        model.fam <- lm(elm.data ~ poly(pam.data, 1))
        predicted <- predict(model.fam,newdata=data.frame(pam.data=u.pam))
        imput[[i]] <- elm
        L <- nrow(elm)
        for (p in 1:(L-1)) {
          for (q in (p+1):L) {
            if (is.na(elm[p,q])) {
              imput[[i]][q,p] <- imput[[i]][p,q] <- predicted[which(u.pam==pam[p,q])]
            }
          }
        }
      }
      if (length(u.pam) == 3) {
        model.fam <- lm(elm.data ~ poly(pam.data, 2))
        predicted <- predict(model.fam,newdata=data.frame(pam.data=u.pam))
        imput[[i]] <- elm
        L <- nrow(elm)
        for (p in 1:(L-1)) {
          for (q in (p+1):L) {
            if (is.na(elm[p,q])) {
              imput[[i]][q,p] <- imput[[i]][p,q] <- predicted[which(u.pam==pam[p,q])]
            }
          }
        }
      }
      if (length(u.pam) > 3) {
        model.fam <- lm(elm.data ~ poly(pam.data, 3))
        predicted <- predict(model.fam,newdata=data.frame(pam.data=u.pam))
        imput[[i]] <- elm
        L <- nrow(elm)
        for (p in 1:(L-1)) {
          for (q in (p+1):L) {
            if (is.na(elm[p,q])) {
              imput[[i]][q,p] <- imput[[i]][p,q] <- predicted[which(u.pam==pam[p,q])]
            }
          }
        }
      }
    } else {
      imput[[i]] <- elm
    }
  }
  return(imput)
}
# imput <- imputation()
# save(imput, file="imput.RData")
load("imput.RData")

# for the purpose of an illustrative graph something like the following 
# can be run (handmade Chocoan example)
predicted.x <- c(0,2,4,6,8)
predicted.y <- c(0,0.2018,0.4018,0.6216,0.7307)
observed.x <- c(0,4,6,6,6,8,6,6,6,8,2,4,8,4,8,8) 
observed.y <- c(0,0.6164,0.3115,0.3392,0.4434,0.7221,0.6424,0.6311,0.6089,0.7827,0.2018,0.3138,0.7124,0.2752,0.7297,0.7067)
plot(predicted.x, predicted.y, main="", xlim=c(0,8), ylim=c(0,0.8), type="l", lwd=2, xlab="", ylab="", col="red")
par(new=TRUE)
plot(observed.x, observed.y, main="Observed and predicted LDN values", pch=16, xlim=c(0,8), ylim=c(0,0.8), xlab="patristic distance", ylab="Levenshtein distance normalized (LDN)", col="pink", cex=1.5)

# at this point relevant objects are:
# glotfams.red, listing relevant Glottolog families
# imput, containing imputed distance matrices
# tree.index.gcodes, listing family glottocodes in the order of the trees object
# trees, containing newick trees from Glottolog
# now add branch lengths to Glottolog trees from imput

# match glotfams.red to get a corresponding vector of family names
fams <- tree.index[match(glotfams.red, tree.index.gcodes)]

## PIPELINE FOR ADDING ASJP LDN BRANCH LENGTHS TO GLOTTOLOG TREES
# The pipeline contains the following:
#  find.index(fam)
#  -- looks up a tree index given the name of a Glottolog family
# simplify(tr)
#  -- takes out internal node names and branch lengths
#  fix.non.br(tr)
#  -- removes non-branching nodes in a pruned tree
#  pruning(fam, matrix.file)
#  -- outputs a list of a phylo object with the pruned tree and a distance matrix
#  bl(fam)
#  -- calling function producing a distance matrix and a pruned Glottolog tree
#     with ASJP LDN branch lengths using the above functions

# function for looking up a tree index given the name of
# a Glottolog family
find.index <- function(fam) {which(tree.index==fam)}
# for instance: find.index("Austroasiatic") # 329

# function for taking branch lengths (digits)
# and internal node labels out of a newick tree
# takes a newick string as input
simplify <- function(tr) {
  new1 <- gsub("\\:[[:digit:]]*", "", tr)
  phy1 <- read.tree(text=new1)
  phy1$node.label <- c()
  new2 <- write.tree(phy1)
  return(new2)
}

# function for detecting and fixing cases of an internal
# non-branching node left "hanging" after tree pruning
# as done with Round's keep_as_tip() function
fix.non.br <- function(tr) {
  library(stringr)
  # tr is a string representing a newick tree
  # pattern ()
  baddo1 <- "\\([[:alnum:]]{8}\\)"
  while ( length(grep(baddo1, tr))==1 ) {
    to.fix <- regmatches(tr,regexpr(baddo1,tr))
    before <- strsplit(tr, baddo1)[[1]][1]
    after <- unlist(str_split(tr,baddo1,n=2))[2]
    fixed <- str_sub(to.fix,2,-2)
    tr <- paste0(before, fixed, after)
  }
  # pattern ((...))
  s1 <- strsplit(tr,"")[[1]]
  to.remove <- c()
  left.friend <- sort(which(s1=="("), decreasing=TRUE)  # positions of (
  right.indices <- which(s1==")")  # positions of )
  right.friend <- c()  # to hold ) pair members
  # pair members
  for (i in 1:length(left.friend)) {
    right.friend[i] <- right.indices[which(right.indices > left.friend[i])[1]]
    right.indices <- right.indices[-which(right.indices > left.friend[i])[1]]
  }
  # find cases where the right friends of two consecutive left guys are
  # also consecutive and store the indices of the outermost friend pair
  # for later removal
  if (length(left.friend)>1) {
    for (i in length(left.friend):2) {
      if ( left.friend[i-1]-left.friend[i] == 1 & right.friend[i]-right.friend[i-1] == 1 ) {
        to.remove <- c(left.friend[i], right.friend[i], to.remove)
      }
    }
    if ( length(to.remove) > 0 ) {
      s1 <- s1[-to.remove]
    }
  }
  tr <- paste(s1, collapse="")
  return(tr)
}

# given a Glottolog tree and an LDN matrix, prune the tree to
# only contain the languages in the matrix
# this will get rid of dialects
# Needs a family name, e.g. Sino-Tibetan, and a distance matrix
pruning <- function(fam, m) {
  index.number <- find.index(fam)
  t.s <- read.tree(text = trees[index.number])  # stands for tree string
	# turn labels into bare glottocodes
	get.gc <- function(s) {  # stands for get glottocode
		a1 <- strsplit(s, "\\[")[[1]][2]
		a2 <- strsplit(a1, "\\]")[[1]][1]
	}
	t.s$node.label <- unlist(lapply(t.s$node.label, get.gc))
	t.s$tip.label <- unlist(lapply(t.s$tip.label, get.gc))
	Round.tree.phylo <- keep_as_tip(t.s, rownames(m))  # Erics Round's function
	Round.tree.newick <- write.tree(Round.tree.phylo)
  simplified.tree.newick <- simplify(Round.tree.newick)
  corrected.tree.newick <- fix.non.br(simplified.tree.newick)
	pruned.tree.phylo <- read.tree(text=corrected.tree.newick)
  return(list(pruned.tree.phylo,m))
}

# calling function: given a language family name (fam)
# selects a distance matrix, prunes the corresponding Glottolog tree,
# and adds ASJP LDN branch lengths to the Glottolog tree,
# outputting a phylo object
bl <- function(fam) {  # stands for branch lengths
  library(ape)  # drop.tip function
  index <- which(fams==fam)
  m <- imput[[index]]
  pruning.out <- pruning(fam, m)
  pruned.tree.phylo <- pruning.out[[1]]
  pruned.matrix <- pruning.out[[2]]
  bltree.phylo <- nnls.tree(pruned.matrix, pruned.tree.phylo, method="unrooted")
  return(bltree.phylo)
}

bl.prune <- function(bltree.phylo, languages.in.data) {
  # drop tips that are not in the data
  languages.in.tree <- bltree.phylo$tip.label
  languages.not.in.data <- setdiff(languages.in.tree, languages.in.data)
  bltree.phylo.pruned <- drop.tip(bltree.phylo, languages.not.in.data)
  return(bltree.phylo.pruned)
}

## PREPARING FOR RUNNING BAYESTRAITS
# prepare .nex and .dat files
# relevant objects: glotfams.red, imput, tree.index.gcodes, trees, fams
# takes a couple of minutes
geo.data <- read.csv("languages_and_dialects_geo.csv", header=TRUE)
# make empty folders to hold the files
dir.create("datfiles")
dir.create("nexfiles")
for (i in 1:length(fams)) {
  cat("preparing BT input files for", fams[i], "\n")
  liim <- rownames(imput[[i]]) # stands for languages in imputed matrix
  # languages.in.data are in this case normally all languages, so excluding
  # dialects, and also languages without coordinates in Glottolog are
  # excluded
  if (length(liim) > 0) {
    phylo.full <- bl(fams[i])
    nexfile <- paste0("nexfiles/", glotfams.red[i], ".nex")
    datfile <- paste0("datfiles/", glotfams.red[i], ".dat")
    cat("", file=datfile)
    include <- c()
    for(j in 1:length(liim)) {
      w_code <- which(geo.data$glottocode==liim[j])
      if (length(w_code) > 0) {
        lon <- geo.data$longitude[w_code]
        lat <- geo.data$latitude[w_code]
        if (is.na(lon)==FALSE & is.na(lat)==FALSE) {
          cat(liim[j], "\t", lon, "\t", lat, "\n", sep="", file=datfile, append=TRUE)
          include <- c(include, liim[j])
        }
      }
    }
    phylo <- bl.prune(phylo.full, include)
    phylo.no.zero.branches <- di2multi(phylo)
    write.nexus(phylo.no.zero.branches, file=nexfile)
  }
}

# check for which families BayesTraits input files have become available
# and the number of taxa in each
check_input_files <- function() {
  d <- dir("nexfiles")
  fs <- c()  # vector to hold family names
  Ns <- c()  # vector to hold number of languages
  for (i in 1:length(d)) {
    f <- paste(strsplit(d[i],"")[[1]][1:8], collapse="")  # family name in Glottolog
    fs[i] <- f
    r <- readLines(paste0("nexfiles/",d[i]))
    a1 <- r[5]
    a2 <- strsplit(a1, "\tDIMENSIONS NTAX = ")[[1]][2]
    a3 <- strsplit(a2, ";")[[1]][1]
    Ns[i] <- as.numeric(a3)
  }
  fams_N <- data.frame(fs,Ns)
  names(fams_N) <- c("fam", "N")
  cat("\nBT input files now available for", length(d), "families,\n")
  cat("where the smallest, with", min(fams_N$N), "taxa, are:\n\n")
  print(fams_N$fam[which(fams_N$N==min(fams_N$N))])
  cat("\n")
  return(fams_N)
}
# optionally run the above function; typing cf will list families and numbers of taxa
# cf <- check_input_files()

# after BayesTraits has run, the below can optionally be used to check 
# for which families BayesTraits output files have become available,
# the number of taxa in each, families not done and the number of taxa in those
check_output_files <- function() {
  d <- dir("BT_outfiles")
  get_name <- function(z) strsplit(z, "_AncStates.txt")[[1]][1]
  outfams <- unlist(lapply(d, get_name))
  data_input_fams <- check_input_files()
  w_out <- match(outfams, data_input_fams$fam)
  missing <- data_input_fams[-w_out,]
  missing
}
# check_output_files()
