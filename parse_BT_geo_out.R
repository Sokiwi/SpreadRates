library(TreeTools)  # as.Newick, NewickTree
library(ape)  # read.nexus, which.edge, node.depth.edgelength
library(caper)  # clade.members.list
library(phytools)  # getDescendants
library(sp)  # SpatialPoints
library(spatstat)  # density.ppp
library(spatstat.geom)  # ppp
library(plotly)  # fpr preparing maps
library(maps)  # map.where
library(assertthat)  # assert_that
library(dplyr)  # inner_join
library(purrr)  # map_dfr
library(igraph)  # graph_from_data_frame
library(ggplot2)  # various plotting functions
library(viridis)  # colors for maps

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

# make a data frame guiding the drawing of edges on a physical map
# get the embedding of nodes, add a value to an edge
# corresponding to the log depth of the most recent of the two
# nodes defining an edge; this corresponds to relative time
# values are scaled 1-12, corresponding to heat colors, cf.
# barplot(rep(1,12), col = heat.colors(12))
# the color will correspond to the age of the "most recent end" of
# the branch

get_drawing_guides <- function(phylo) {
  # first a guide that will help assigning colors to points
  # corresponding to their depths
  # to enable log transform zero nodes (if any) will be assigned a small depth
  nd <- node.depth.edgelength(force.ultrametric(phylo, method="extend"))  # node depths
  nda <- nd  # node depths adjusted
  w_zero <- which(nda==0)
  if (length(w_zero) > 0) {
    # the small depth is the square of the minimum of the other depths,
    # making certain that the root is the oldest, but maybe not out of line
    nda[w_zero] <- min(nda[-w_zero])^2
  }
  ndl <- log(nda)
  nd.scaled <- round(1 + 11*((ndl - min(ndl)) / (max(ndl) - min(ndl))))
  drawing.guide.nodes <- cbind(1:length(nd.scaled),nd.scaled)
  colnames(drawing.guide.nodes) <- c("node","node.age.scaled")
  
  drawing.guide.edges <- as.data.frame(phylo$edge)
  names(drawing.guide.edges) <- c("node1", "node2")
  edge.age.scaled <- c()
  for (i in 1:nrow(drawing.guide.edges)) {
    node1.age <- drawing.guide.nodes[drawing.guide.edges[i,1],2]
    node2.age <- drawing.guide.nodes[drawing.guide.edges[i,2],2]
    max.age <- max(node1.age, node2.age)  # max age since root
    edge.age.scaled[i] <- max.age
  }  
  drawing.guide.edges <- cbind(drawing.guide.edges,edge.age.scaled)
  return(list(drawing.guide.edges,drawing.guide.nodes))
}

get_coords <- function(clades, all.text) {
  BT.node.text <- all.text[-grep("^Node", all.text)]
  writeLines(BT.node.text, con="tmp.txt")
  BT.ni <- read.table(file="tmp.txt", header=TRUE, sep="\t") # BT node info
  BT.ni <- BT.ni[,-ncol(BT.ni)]
  BT.ni <- BT.ni[501:1000,]
  names(BT.ni) <- unlist(lapply(names(BT.ni), function(x) gsub("\\.", "_", x)))
  file.remove("tmp.txt")
  lats <- list()
  lons <- list()
  for (i in 1:nrow(clades)) {
    column.name.lat <- paste0(clades$BT[i], "___Lat")
    column.name.lon <- paste0(clades$BT[i], "___Long")
    klat <- which(names(BT.ni)==column.name.lat)
    klon <- which(names(BT.ni)==column.name.lon)
    ulat <- unique(BT.ni[,klat])
    ulon <- unique(BT.ni[,klon])
    # tips have only one coordinate, which is extracted here
    if (length(ulat)==1 & length(ulon)==1) {
      lats[[i]] <- ulat
      lons[[i]] <- ulon
    } else {
      # find 95% HPD interval, put the points in swarms.lat and swarms.lon
      x <- as.vector(BT.ni[,klon])
      y <- as.vector(BT.ni[,klat])
      # p <- SpatialPoints(coords = matrix(c(x, y, ncol = 2, nrow=length(x))))
      p <- SpatialPoints(coords = as.matrix(cbind(x,y)))
      w=owin(c(min(x),max(x)),c(min(y),max(y)))
      pp = ppp(x,y, window=w)
      dp <- density.ppp(pp, sigma = 0.1, at="points")
      # how many to exclude? 25 in case of 500 in the posterior
      biteoff <- nrow(BT.ni) - .95 * nrow(BT.ni)  
      dp.sorted <- sort(dp)
      exclude <- dp.sorted[1:biteoff]
      lons[[i]] <- x[-match(exclude,dp)]  # 95% HPD longitudes
      lats[[i]] <- y[-match(exclude,dp)]  # 95% HPD latitudes
    }
  }
  coordinates <- (list(lats, lons))
  return(coordinates)
}

# the timeframe is the number of time steps from the roots to include,
# with the default set to 12, which is the maximum
map.homelands <- function(clades, drawing.guide.edges, 
         lats, lons, filename, format="pdf", timeframe=12) {
  glotto.info <- read.csv("languoid.csv", header=TRUE)
  s1 <- strsplit(filename, "/")[[1]][2]
  glottocode <- strsplit(s1, "_")[[1]][1]
  w_fam <- grep(glottocode, glotto.info$id)
  fam.name <- glotto.info$name[w_fam]
  
  w_family.homeland <- which(clades$BT=="Node_00000")
  homelat <- clades$lat[w_family.homeland]
  homelon <- clades$lon[w_family.homeland]

  edges <- drawing.guide.edges
  edges <- edges[which(drawing.guide.edges$edge.age.scaled <= timeframe),]
  weight <- rep(1, nrow(edges))
  edges <- cbind(edges, weight)
  names(edges) <- c("from", "to", "timesteps", "weight")
  edges <- edges[,c("from","to","weight","timesteps")]
  if (timeframe < 12) {
    edges <- edges[which(edges$timesteps <= timeframe),]
    relevant.ids <- unique(union(edges$from, edges$to))
  }
  
  nodes <- data.frame(clades$phylo, clades$lon, clades$lat, clades$name)
  names(nodes) <- c('id', 'lon', 'lat', 'name')
  if (timeframe < 12) {
    nodes <- nodes[which(nodes$id %in% relevant.ids),]
    relevant.ids <- unique(union(edges$from, edges$to))
  }
  
  edges <- edges %>% mutate(timesteps = as.factor(timesteps))
  g <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
  edges_for_plot <- edges %>%
    inner_join(nodes %>% select(id, lon, lat), by = c('from' = 'id')) %>%
    rename(x = lon, y = lat) %>%
    inner_join(nodes %>% select(id, lon, lat), by = c('to' = 'id')) %>%
    rename(xend = lon, yend = lat)
  
  assert_that(nrow(edges_for_plot) == nrow(edges))
  
  nodes$weight = degree(g)
  
  maptheme <- theme(panel.grid = element_blank()) +
    theme(axis.text = element_blank()) +
    theme(axis.ticks = element_blank()) +
    theme(axis.title = element_blank()) +
    theme(legend.position = "bottom") +
    theme(panel.grid = element_blank()) +
    theme(panel.background = element_rect(fill = "#CECECE")) +
    # theme(panel.background = element_rect(fill = "#596673")) +
    theme(plot.margin = unit(c(0, 0, 0.5, 0), 'cm'))
  
  country_shapes <- geom_polygon(aes(x = long, y = lat, group = group),
                                 data = map_data('world'),
                                 fill = "#596673", color = "#515151",
  #                              fill = "#CECECE", color = "#515151",
                                 size = 0.15)
  if (timeframe==12) {
    minlat <- round(min(unlist(lapply(lats, min))) - .2, 1)
    maxlat <- round(max(unlist(lapply(lats, max))) + .2, 1)
    minlon <- round(min(unlist(lapply(lons, min))) - .2, 1)
    maxlon <- round(max(unlist(lapply(lons, max))) + .2, 1)
    mapcoords <- coord_fixed(xlim = c(minlon, maxlon), ylim = c(minlat, maxlat))
  }
  if (timeframe < 12) {
    minlat <- round(min(unlist(lapply(nodes$lat, min))) - .2, 1)
    maxlat <- round(max(unlist(lapply(nodes$lat, max))) + .2, 1)
    minlon <- round(min(unlist(lapply(nodes$lon, min))) - .2, 1)
    maxlon <- round(max(unlist(lapply(nodes$lon, max))) + .2, 1)
    mapcoords <- coord_fixed(xlim = c(minlon, maxlon), ylim = c(minlat, maxlat))
  }
  
  p <- ggplot(nodes) + country_shapes +
    geom_curve(aes(x = x, y = y, xend = xend, yend = yend,     # draw edges as arcs
                   # color = category, size = weight),
                   color = timesteps, size = .5),
               # data = edges_for_plot, curvature = 0.33,
               data = edges_for_plot, curvature = 0.1,
               alpha = 0.5) +
    # scale_color_brewer(palette = "RdYlBu")
    scale_color_viridis(discrete = TRUE, option = "A") +
    scale_size_continuous(guide = FALSE, range = c(0.25, 2)) + # scale for edge widths
    #  geom_point(aes(x = lon, y = lat, size = weight),           # draw nodes
    #             shape = 21, fill = 'white',
    #             color = 'black', stroke = 0.5) +
    geom_point(aes(x = homelon, y = homelat, size = weight),           # draw nodes
               shape = 21, fill = 'white',
               color = 'black', stroke = 0.5) +
    scale_size_continuous(guide = FALSE, range = c(1, 6)) +    # scale for node size
    # geom_text(aes(x = lon, y = lat, label = name),             # draw text labels
    #          hjust = 0, nudge_x = 1, nudge_y = 4,
    #          size = 3, color = "white", fontface = "bold") +
    # geom_point(data=lons[[1]],lats[[1]])
    ggtitle(fam.name) +
    mapcoords + maptheme
    if (format=="png") {
      if (timeframe==12) {
        name.plotfile <- paste0("mapfiles/", fam.name, ".png")
        png(name.plotfile)
      }
      if (timeframe<12) {
        name.plotfile <- paste0("mapfiles/", fam.name, "_timeframe", timeframe, ".png")
        png(name.plotfile)
      }
    }  
  if (format=="pdf") {
    if (timeframe==12) {
      name.plotfile <- paste0("mapfiles/", fam.name, ".pdf")
      pdf(name.plotfile)
    }
    if (timeframe < 12) {
      name.plotfile <- paste0("mapfiles/", fam.name, "_timeframe", timeframe, ".pdf")
      pdf(name.plotfile)
    }
  }  
  print(p)
    dev.off()
}

# the following loops over BT output files making pictures for individual
# families; options in map.homelands to change timeframe and format
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
  drawing.guides <- get_drawing_guides(phylo)  # output accessed in next
  drawing.guide.edges <- drawing.guides[[1]]
  drawing.guide.nodes <- drawing.guides[[2]]; rm(drawing.guides)
  coords <- get_coords(clades, all.text)  # output accessed in next
  lats <- coords[[1]]
  lons <- coords[[2]]; rm(coords)
  map.homelands(clades, drawing.guide.edges, lats, lons, filename, format="png", timeframe=12)
}

# Here the process of making maps for BT output stops and the 
# process of comparisons with md output starts
# (using files copied from C:\Wichmann\Current\Archaeology and language conference\md\
# but not completely replicable without some other data preparation files from
# that folder)

# make maps similar to the one in the Phylogeography chapter for the
# Oxford Handbook of Archaeology and Linguistics
# (non-Pacific and Pacific), displaying results for md and BT for families
# in one map
# This uses code and files prepared for the OHAL chapter, but
# copied and modified to the present working directory

# using some of the above code make a data frame with just families and 
# single-coordinate homelands, BT.family.homelands
# takes 2-3 mins, so saved as an R object
get.BT.family.homelands <- function() {
  BTgroup <- c(); BTlat <- c(); BTlon <- c()
  for (i in 1:length(dir("BT_outfiles"))) {
    if ("fams" %in% ls()==FALSE) fams <- getfams()  # only needed once
    filename <- paste0("BT_outfiles/", dir("BT_outfiles")[i])
    all.text <- readLines(filename)
    phylo <- get_phylo(filename)
    clades.tmp <- match_nodes(phylo, all.text, filename)
    clades <- get_clades_dots(clades.tmp, all.text); rm(clades.tmp)
    w_node0 <- which(clades$BT=="Node_00000")
    BTgroup[i] <- clades$name[w_node0]
    BTlat[i] <- clades$lat[w_node0]
    BTlon[i] <- clades$lon[w_node0]
  }
  BTgroup <- unlist(lapply(BTgroup, function(x) gsub(" ", "_", x)))
  BT.family.homelands <- data.frame(BTgroup, BTlat, BTlon)
  return(BT.family.homelands)
}
# can be run as BT.family.homelands <- get.BT.family.homelands()
# save(BT.family.homelands, file="BT.family.homelands.RData")
load("BT.family.homelands.RData")

# read output of md_fams_subgroups_2.R with homeland data
all_homelands <- read.table(file="homelands_Glottolog_groups.txt", sep="\t", header=TRUE, quote="")
# eliminate data for subgroups
w_subgroups <- grep(",", all_homelands$group)
f <- all_homelands[-w_subgroups,]
abbreviation_guide <- read.table(file="family_abbreviations.txt", sep="\t", header=TRUE, quote="")
mfam <- match(f$group, abbreviation_guide$group)
abb <- abbreviation_guide$abbreviation[mfam]
area <- abbreviation_guide$area[mfam]
f <- cbind(f, abb, area)
f$group <- unlist(lapply(f$group, trimws))

# add column to hold names similar to those in the BT output
f$BTgroup <- f$group
# remove white space in the BTgroup column, but otherwise use this as the
# model for the names
f$BTgroup <- unlist(lapply(f$BTgroup, trimws))
# for matching fams to BT names use underscore in fams instead of spaces
fams <- getfams()
fams2 <- unlist(lapply(fams, function(x) gsub(" ", "_", x)))
# now the names of md and BT output have the same format

# join md output (f) and BT output (BT.family.homelands)
library(plyr)  # inner_join function
library(geosphere)  # greatCircle function
# compute distances between the family homelands for the two methods
names(BT.family.homelands)[1] <- "BTgroup"
f2 <- inner_join(f, BT.family.homelands, "BTgroup")
m1 <- as.matrix(data.frame(f2$lon, f2$lat))
m2 <- as.matrix(data.frame(f2$BTlon, f2$BTlat))
dist.km <- round(distHaversine(m1, m2)/1000)
f2 <- cbind(f2, dist.km)

# map homelands differing by more than mindist km (defaults to 100 km)
library(maps)

# X <- x
# Y <- y

library(maptools)
library(rgeos)

# define Greater Mesoamerican families in case it is made use of in
# the simple_map() function below 
GreatMes <- c("Mayan", "Misumalpan", "Mixe-Zoque", "Otomanguean", "Uto-Aztecan")
GreatMes2 <- c("Chitimacha", "Coahuilteco", "Comecrudan", "Cotoname", 
  "Cuitlatec", "Huavean", "Jicaquean", "Lencan", "Tarascan", "Tequistlatecan", 
  "Totonacan", "Xincan", "Atakapa")

# options for target_area: "ot" (other) and "NG"
# and for choosing Greater Mesoamerica ("GM")
simple_map <- function(mindist=300, target_area="ot", GM=FALSE) {
  if (GM==TRUE) {
    target_area <- "ot"
    mindist <- -1
  }
  f2$group <- unlist(lapply(f2$group, trimws))
  wa <- which(f2$area==target_area)
  f2 <- f2[wa,]
  gm <- which(f2$group %in% GreatMes)
  f2 <- f2[gm,]
  w_mindist <- which(f2$dist.km > mindist)
  f2 <- f2[w_mindist,]
  # for the Greater Mesoamerica families for which there are no BT results
  # it is posited, as a trick, that they are the same as for md; then
  # that can be added as if there were two sets of coordinates to plot,
  # and the red dots for md are plotted last so they will overwrite the 
  # blue dots for BT
  if (GM==TRUE) {
    fadd <- data.frame(f$lat, f$lon, rep(0, nrow(f)))
    names(fadd) <- c("BTlat","BTlon","dist.km")
    f3 <- cbind(f, fadd)
    f3 <- f3[match(GreatMes2, f3$group),]
    f2 <- rbind(f2, f3)
  }
  x_md <- f2$lat
  y_md <- f2$lon
  x_BT <- f2$BTlat
  y_BT <- f2$BTlon
  lb <- f2$abb
  md <- cbind(as.numeric(x_md), as.numeric(y_md))
  BT <- cbind(as.numeric(x_BT), as.numeric(y_BT))
  all <- rbind(md, BT)
  # get max min long and lat and add a frame of 10% around the points
  lat_range <- range(all[,1])
  lng_range <- range(all[,2])
  lat_extend <- 0.1 * diff(lat_range)
  lng_extend <- 0.1 * diff(lng_range)
  lat_range <- c(lat_range - lat_extend, lat_range + lat_extend)
  lng_range <- c(lng_range - lng_extend, lng_range + lng_extend)
  data(wrld_simpl)
  plot(wrld_simpl, col = 'grey', border = NA, axes = T, xlim=c(min(lng_range), 
    max(lng_range)), ylim=c(min(lat_range), max(lat_range)), 
    main="Comparison of homelands", 
    xlab="Blue: BayesTraits  Red: minimal distance")
  for (i in 1:length(md[,1])) {
    points(BT[i,2],BT[i,1],col="blue",cex=1.2,pch=i)
    points(md[i,2],md[i,1],col="red",cex=1.2,pch=i)
  }
  text(x = y_BT, y = x_BT, lb, pos = 4, cex = .8)
}
