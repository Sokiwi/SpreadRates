library(cartography)  # carto.pal
library(RColorBrewer)  # for making plots of diffusion rates in different biomes
library(ddalpha)

# the function as() looks at periods of Y years from "lowest" (default is
# 6000 BP to present)
# within each period it checks the speed of active diffusion events
# and then takes means and medians, outputting a file with that and N
# Can plot make plots for areas and biomes

# areas are: Africa, Australia, Austronesian, Eurasia, North America, Papuan, South America

# plots and analyses
# mean rates across areas (Fig. 1)
as(Y=1, lowest=7000, loess.span=.15, upper=2, present=300, what="areas")
# Mean diffusion across biomes (Fig. 2)
as(what="biomes")
# ratio of EW to SW movements across areas (Fig. 3)
as(what="bearings_areas")
# mean rates across subsistence patterns (Fig. 4)
as(Y=1, lowest=5000, loess.span=.2, upper=0.9, present=300, what="subsistence")
# Ratios of EW to NS movements for AGR vs. HG languages (Fig. 5)
as(what="bearings_subsistence")
# print Pearson and Spearman correlations between speed and 
# altitude, rugosity, and npp
as(what="continuous_variables")

# as stands for average speed
as <- function(Y=200, lowest=7000, loess.span=.5, what="world", upper=10, present=300, ...) {

  ### data preparation ###
  if (exists("Y") == FALSE) {
    Y <- 1
  }
	filename <- paste("diffusion_averages_", Y, "_year_period.txt", sep="")
	cat("lower_bound\taverage_speed\tmedian_speed\tN\n", file=filename)
	de <- read.table(file="diffusion_events_annotated.txt", header=TRUE, 
	  sep="\t", quote="", na.strings="", comment.char="")
	
	# correct a typo
	w_evegreen <- grep("Evegreen", de$biome4_names)
	if (length(w_evegreen) > 0) {
	  de$biome4_names[w_evegreen] <- "Evergreen taiga/montane forest"
	}
  # turn columns with numbers into numerics
  de$perc_b4 <- as.numeric(de$perc_b4)
  de$altitude <- as.numeric(de$altitude)
  de$rugosity <- as.numeric(de$rugosity)
  de$npp <- as.numeric(de$npp)
	# get rid of cases of negative ages,
	# which represent very recent splitting events, estimate to have happened
	# in the future
	w_neg_age <- union(which(de$BP_mo < 0), which(de$BP_da < 0))
	if (length(w_neg_age) > 0 ) { 
		de <- de[-w_neg_age,]
	}
	# get get rid of cases where a daughter is estimated to have split
	# before or at the same time as a mother, which may involve problems 
	# with the classification,
	# misassignment of ISO-code or other sources of error
	# This is not really necessary since it is already done earlier in the pipeline
	w_old_da <- which(de$BP_M <= de$BP_D)
	if ( length(w_old_da) > 0 ) {
		de <- de[-w_old_da,]
	}
	# get get rid of cases where age estimates for
	# mothers or daughters are NA
	w_na <- unique(union(which(de$BP_M=="NA"),which(de$BP_D=="NA")))
	if ( length(w_na) > 0 ) {
		de <- de[-w_na,]
	}

	### function for creating a transparent color when plotting
	create_col <- function(color, percent = 50, name = NULL) {
	  rgb.val <- col2rgb(color)
	  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3], 
	               max = 255,
	               alpha = (100 - percent) * 255 / 100,
	               names = name)
	  return(t.col)
	}
	
	### areas ###
  if ( what=="areas" ) {
    # from outside the function runs as as(Y=1, lowest=7000, loess.span=.15, upper=2, present=300, what="areas")
    areas <- sort(unique(de$area))
    w_uniden <- which(areas=="unidentified")
    if ( length(w_uniden) > 0 ) {
      areas <- areas[-w_uniden]
    }
    for ( a in 1:length(areas) ) {
      filename <- paste("diffusion_averages_", areas[a], "_", Y, "_year_period.txt", sep="")
      cat("lower_bound\tmean_speed\tmean_speed_from\tmean_speed_to\tN\n", file=filename)
      de_sub <- de[which(de$area==areas[a]),]
      lower <- lowest
      while ( lower >= present ) {
        BPM <- as.numeric(de_sub$BP_mo)
        BPD <- as.numeric(de_sub$BP_da)
        window.edge <- lower - Y
        w_in <- which( !((BPM < lower & BPD < window.edge) | (BPM > lower & BPD > window.edge)) )
        if ( length(w_in)==0 ) {
          cat(lower, "\t", NA, "\t", NA, "\t", NA, "\t", length(w_in), "\n", sep="")
          cat(lower, "\t", NA, "\t", NA, "\t", NA, "\t", length(w_in), "\n", sep="", file=filename, append=TRUE)
        } else {
          mean_speed <- mean(as.numeric(de_sub$speed[w_in]))
          mean_speed_from <- mean(as.numeric(de_sub$speed_from[w_in]))
          mean_speed_to <- mean(as.numeric(de_sub$speed_to[w_in]))
          cat(lower, "\t", mean_speed, "\t", mean_speed_from, "\t", mean_speed_to, "\t", length(w_in), "\n", sep="")
          cat(lower, "\t", mean_speed, "\t", mean_speed_from, "\t", mean_speed_to, "\t", length(w_in), "\n", sep="", file=filename, append=TRUE)
        }
        lower <- lower - Y
      }
    }
    nareas <- length(areas)
    colors <- c("blue","red","green3","darkorchid1","orange","maroon1","lightsalmon3")
    linetypes <- c(1,1,1,1,1,1,1)
    plot(100, 100, xlim=c(-lowest,0), ylim=c(0,upper), xlab="BP", ylab="km/year")
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "#EBEBEB")
    abline(v=seq(-lowest,0,200), h=seq(0,upper,0.1), col="white")
    legend((-lowest), upper, areas, cex=1, col=colors, lty=linetypes, lwd=3)
    # plot with lines only
    for ( a in 1:length(areas) ) {
      filename <- paste("diffusion_averages_", areas[a], "_", Y, "_year_period.txt", sep="")
      plot_data <- read.table(file=filename, header=TRUE, sep="\t")
      w_na2 <- which(is.na(plot_data$mean_speed))
      if ( length(w_na2) > 0 ) {
        plot_data <- plot_data[-w_na2,]
      }
      plot_data$lower_bound <- -1*plot_data$lower_bound
      lines(predict(loess(plot_data$mean_speed ~ plot_data$lower_bound, 
            span=loess.span)), x=plot_data$lower_bound, 
            col=colors[a], lty=linetypes[a], lwd=4)
    } # plot with error bars
    for ( a in 1:length(areas) ) {
      filename <- paste("diffusion_averages_", areas[a], "_", Y, "_year_period.txt", sep="")
      plot_data <- read.table(file=filename, header=TRUE, sep="\t")
      w_na2 <- which(is.na(plot_data$mean_speed))
      if ( length(w_na2) > 0 ) {
        plot_data <- plot_data[-w_na2,]
      }
      plot_data$lower_bound <- -1*plot_data$lower_bound
      X <- plot_data$lower_bound
      Yp <- predict(loess(plot_data$mean_speed ~ X, span=loess.span))
      Yp2 <- predict(loess(plot_data$mean_speed_from ~ X, span=loess.span))
      Yp3 <- predict(loess(plot_data$mean_speed_to ~ X, span=loess.span))
      lines(X, Yp, col=colors[a], lty=linetypes[a], lwd=4)
      polygon(c(X, rev(X)), c(Yp2,rev(Yp3)), col=create_col(colors[a], 90), border=NA)
    }

	### biomes ###
	} else if ( what=="biomes") {
	  biome_key <- get_biome_classes("Example")
	  # correcting a typo
	  w_evegreen <- which(biome_key$category=="Evegreen taiga/montane forest")
	  if (length(w_evegreen) > 0) {
	    biome_key$category[w_evegreen] <- "Evergreen taiga/montane forest"
	  }
	  names(biome_key) <- c("key","description")
	  # Sea is not represented in biome4 explicitly, only by NA
	  # so is added here explicitly, given the code 100
	  # It has been verified that all NA cases are in fact Sea in biome98
	  biome_key <- rbind(biome_key, (c(100, "Sea")))

    # define as "Sea" and code 100 biome4 where both biome4 and biome98 are "Sea"
	  # as these values are changed a new data frame is defined
    de2 <- de
	  sea <- intersect(which(de$biome4_names=="Sea"), which(de$biome98_names=="Sea"))
    nas <- setdiff(which(de$biome4_names=="Sea"), which(de$biome98_names=="Sea"))
	  de2$biome4[sea] <- 100
	  de2$biome4_names[sea] <- "Sea"
	  de2$biome4[nas] <- NA
	  de2$biome4_names[nas] <- NA
	  ## find the events where biomes are represented by 90% or more
		## do t-tests for cases where the number of events N >= 10

		# get rid of possible cases where the rate is infinite
	  # (in practice not found)
		is_inf <- which(de2$speed==Inf)
		if ( length(is_inf) > 0 ) {
			de2 <- de2[-is_inf,]
		}

		# get rid of some cases where the rate is NA
		# (in practice not found)
		is_na <- which(is.na(de2$speed))
		if ( length(is_na) > 0 ) {
			de2 <- de2[-is_na,]
		}
		
		rel <- which(as.numeric(de2$perc_b4) >= 90)  # stands for relevant
		x <- de2[rel,]
		
		# focus on biomes that are involved in at least 10 cases
		tb <- table(x$biome4)
		wab <- names(which(tb >= 10))  # well attested biomes
		
		biomes <- biome_key$key[biome_key$key %in% wab]
		designations <- biome_key$description[biome_key$key %in% wab]
		L <- length(biomes)

		cat("biome A\tbiome B\tmean A\tmean B\tp\n", file="t tests biomes.txt")
		for (i in 1:(L-1)) {
			for (j in (i+1):L) {
			  rawA <- as.numeric(x$speed[which(x$biome4==biomes[i])])
				rawB <- as.numeric(x$speed[which(x$biome4==biomes[j])])
				res <- t.test(rawA, rawB, alternative="two.sided")
				p <- res$p.value
				  cat(designations[i], "\t", designations[j], "\t", mean(rawA), "\t", mean(rawB), "\t", p, "\n", file="t tests biomes.txt", append=TRUE)
			}
		}

		# get rates and N for total distances traveled in different biomes 
		# (using the 90% criterion)
		# collect data in vectors to be merged with biome_keys
		N <- c()
		biome_rates <- c()
		biome_rates_from <- c()
		biome_rates_to <- c()
		for ( i in 1:length(biomes) ) {
		  ind <- which(x$biome4==biomes[i])
		  N[i] <- length(ind)
		  biome_rates[i] <- mean(x$speed[ind])
		  biome_rates_from[i] <- mean(x$speed_from[ind])
		  biome_rates_to[i] <- mean(x$speed_to[ind])
		}
		
		biome_rates_ranked <- biome_rates[order(biome_rates, decreasing=TRUE)]
		biome_rates_from_ranked <- biome_rates_from[order(biome_rates, decreasing=TRUE)]
		biome_rates_to_ranked <- biome_rates_to[order(biome_rates, decreasing=TRUE)]
		biomes_ranked <- biomes[order(biome_rates, decreasing=TRUE)]
		biome_descriptions_ranked <- biome_key$description[match(biomes_ranked, biome_key$key)]
		N_ranked <- N[order(biome_rates, decreasing=TRUE)]
		
		biome_results_full <- data.frame(biome_rates_ranked, 
		                                 biome_rates_from_ranked,
		                                 biome_rates_to_ranked,
		                                 biomes_ranked, 
		                                 biome_descriptions_ranked, 
		                                 N_ranked)
		biome_results <- biome_results_full[biome_results_full$N_ranked > 10,]
		
		# plot mean diffusion rates across biomes
		# colors that are just distinct, legend generated through defaults
		pal <- brewer.pal(n = nrow(biome_results), name = "Paired")
		bp <- barplot(biome_results$biome_rates_ranked, ylim=c(0,8), 
		        legend.text=biome_results$biome_descriptions_ranked, 
		        col=pal, ylab="km/year")

		arrows95 <- arrows(x0=bp, x1=bp, y0=biome_results$biome_rates_from_ranked, 
		                   y1=biome_results$biome_rates_to_ranked, angle=90, code=3, length=.15)
		Ns <- paste("N = ", biome_results$N_ranked)
		text(bp, biome_results$biome_rates_to_ranked, Ns, cex=.6, pos=3)

  ### subsistence ###
  } else if ( what=="subsistence") {
    # from outside function run as as(Y=1, lowest=5000, loess.span=.2, upper=.9, present=300, what="subsistence")
    # assign subtypes of hunter-gatherers to HG
		for (i in 1:nrow(de)) {
			if ( de$subsistence[i]=="HG-SED, HG-SAGO" ) { de$subsistence[i] <- "HG" }
			if ( de$subsistence[i]=="HG-SAGO" ) { de$subsistence[i] <- "HG" }
			if ( de$subsistence[i]=="HG-SED" ) { de$subsistence[i] <- "HG" }
			if ( de$subsistence[i]=="HG-SAGO, HG-SED" ) { de$subsistence[i] <- "HG" }
			if ( de$subsistence[i]=="HG-FISH" ) { de$subsistence[i] <- "HG" }
		}
	  # plot mean rates across subsistence patterns
    # prepare plotdata
    subsist <- c("AGR", "HG")
    for ( a in 1:2 ) {
			filename <- paste("diffusion_averages_", subsist[a], "_", Y, "_year_period.txt", sep="")
			cat("lower_bound\taverage_speed\taverage_speed_from\taverage_speed_to\tN\n")
			cat("lower_bound\taverage_speed\taverage_speed_from\taverage_speed_to\tN\n", file=filename)
			de_sub <- de[which(de$subsistence==subsist[a]),]
			lower <- lowest
			while ( lower >= present ) {
				BPM <- as.numeric(de_sub$BP_mo)
				BPD <- as.numeric(de_sub$BP_da)
				window.edge <- lower - Y
				w_in <- which( !((BPM < lower & BPD < window.edge) | (BPM > lower & BPD > window.edge)) )
				if ( length(w_in)==0 ) {
					cat(lower, "\t", NA, "\t", NA, "\t", length(w_in), "\n", sep="")
					cat(lower, "\t", NA, "\t", NA, "\t", length(w_in), "\n", sep="", file=filename, append=TRUE)
				} else {
					speed_average <- mean(as.numeric(de_sub$speed[w_in]))
					speed_average_from <- mean(as.numeric(de_sub$speed_from[w_in]))
					speed_average_to <- mean(as.numeric(de_sub$speed_to[w_in]))
					cat(lower, "\t", speed_average, "\t", speed_average_from, "\t", speed_average_to, "\t", length(w_in), "\n", sep="")
					cat(lower, "\t", speed_average, "\t", speed_average_from, "\t", speed_average_to, "\t", length(w_in), "\n", sep="", file=filename, append=TRUE)
				}
				lower <- lower - Y
			}
    }
    # prepare plot background
    dev.off()
    colors <- c("lightblue", "pink", "blue", "red")
    linetypes <- c(1,1)
    plot(100, 100, xlim=c(-lowest,0), ylim=c(0,upper), xlab="BP", ylab="km/year")
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "#EBEBEB")
    abline(v=seq(-lowest,0,200), h=seq(0,1,0.1), col="white")
    legend(-lowest, upper, subsist, cex=1.3, col=colors[3:4], lty=linetypes, lwd=3)
    # plotting
		for ( a in 1:2 ) {
		  filename <- paste("diffusion_averages_", subsist[a], "_", Y, "_year_period.txt", sep="")
		  plot_data <- read.table(file=filename, header=TRUE, sep="\t")
			w_na2 <- unique(union(which(is.na(plot_data$average_speed)), which(is.na(plot_data$median_speed))))
			if ( length(w_na2) > 0 ) {
				plot_data <- plot_data[-w_na2,]
			}
			plot_data$lower_bound <- -1*plot_data$lower_bound
			
			X <- plot_data$lower_bound
			Y1 <- predict(loess(plot_data$average_speed ~ X, span=loess.span))
      lines(X, Y1, col=colors[a+2], lty=linetypes[a], lwd=4)
      Y2 <- predict(loess(plot_data$average_speed_from ~ X, span=loess.span))
      # in case a line is wanted for the lower bound
      # lines(X, Y2, col=colors[a], lty=linetypes[a], lwd=4)
      Y3 <- predict(loess(plot_data$average_speed_to ~ X, span=loess.span))
      # in case a line is wanted for the upper bound
      # lines(X, Y3, col=colors[a], lty=linetypes[a], lwd=4)
      polygon(c(X, rev(X)), c(Y2,rev(Y3)), col=create_col(colors[a+2], 70), border=NA)
		}

	### bearings and areas ###
	} else if (what=="bearings_areas") {
	  afb <- sort(unique(de$area))  # areas for bearings
	  mab <- matrix(0,ncol=length(afb),nrow=2,dimnames=list(c("EW", "NS"),afb))  # matrix with areas and bearings
	  EW <- grep("EW", de$bearing)
	  NS <- grep("NS", de$bearing)
	  for (i in 1:length(afb)) {
	    w_area <- which(de$area==afb[i])
	    mab[1,i] <- length(intersect(w_area,EW))
	    mab[2,i] <- length(intersect(w_area,NS))
	  }
	  # make barplot similar to the one for bearings_subsistence but 
	  # hardwired colors same as those of diffusion rates
	  # order the matrix by descending proportions
	  ord <- order(mab[1,]/(mab[1,] + mab[2,]), decreasing=TRUE)
    mab_ord <- mab[,ord]
    # prepare data for a column for the World
    World <- apply(mab_ord, 1, sum)
    mab_ord <- as.matrix(data.frame(mab_ord, World))
    bearing_proportions_areas <- 100*(mab_ord[1,]/(mab_ord[1,] + mab_ord[2,]))
    abbr1 <- names(bearing_proportions_areas)
	  abbr2 <- gsub("North\\.", "N. ", abbr1)
	  abbr3 <- gsub("South\\.", "S. ", abbr2)
	  names(bearing_proportions_areas) <- abbr3
	  # calculate margins of error
	  n <- apply(mab_ord, 2, sum)
	  moe95 <- 100*1.96*sqrt(((bearing_proportions_areas/100)*(1-bearing_proportions_areas/100))/n)
	  # make a text vector of sample sizes
	  ss <- paste0("n=",n)
	  colors <- c(c("blue","red","green3","darkorchid1","orange","maroon1","lightsalmon3")[ord],"gray")
	  bp <- barplot(bearing_proportions_areas, ylab="proportion (%) of EW to NS movements", 
	     col=colors, xlim=NULL, ylim=c(0,70))
	  abline(50, 0, col="black", lty=2)
	  arrows95 <- arrows(x0=bp, x1=bp, y0=bearing_proportions_areas-moe95, 
	         y1=bearing_proportions_areas+moe95, angle=90, code=3, length=.15)
    text(bp, 1, ss, cex=1, pos=3)
	  # do pairwise chi-square tests for all pairs of areas
	  ab_chisq <- matrix(100,ncol=length(afb),nrow=length(afb),dimnames=list(afb,afb))
	  for (i in 1:ncol(ab_chisq)) {
	    for (j in 1:nrow(ab_chisq)) {
	      chisq <- chisq.test(mab[,i],mab[,j])
	      p <- chisq$p.value
	      ab_chisq[j,i] <- p
	    }
	  }
	  ab_chisq
	  
	### continuous variables ###
	} else if (what=="continuous_variables") {
	  speed_data <- de$speed
	  log_speed_data <- log(speed_data)
	  w_zero_speed <- which(speed_data == 0)
	  if (length(w_zero_speed) > 0) {
  	  log_speed_data[w_zero_speed] <- log(0.001)
	  }
	  altitude_data <- as.numeric(de$altitude)
	  log_altitude_data <- log(altitude_data)
	  w_zero_altitude <- which(altitude_data == 0)
	  if (length(w_zero_altitude) > 0) {
	    log_altitude_data[w_zero_altitude] <- log(0.001)
	  }
	  rugosity_data <- as.numeric(de$rugosity)
	  log_rugosity_data <- log(rugosity_data)
	  w_zero_rugosity <- which(rugosity_data == 0)
	  if (length(w_zero_rugosity) > 0) {
	    log_rugosity_data[w_zero_rugosity] <- log(0.001)
	  }
	  npp_data <- as.numeric(de$npp)
	  w_sea <- which(de$biome4_names=="Sea" | de$biome98_names=="Sea")
	  if (length(w_sea) > 0) {
	    speed_data[w_sea] <- NA
	    log_speed_data[w_sea] <- NA
	    altitude_data[w_sea] <- NA
	    log_altitude_data[w_sea] <- NA
	    rugosity_data[w_sea] <- NA
	    log_rugosity_data[w_sea] <- NA
	    npp_data[w_sea] <- NA
	  }
	  cor.alt.pearson <- cor.test(altitude_data, log_speed_data)  # this and the two next in case logs are wanted
	  cor.rug.pearson <- cor.test(rugosity_data, log_speed_data)
	  cor.npp.pearson <- cor.test(npp_data, log_speed_data)
	  cor.alt.spearman <- cor.test(altitude_data, speed_data, method="spearman")
	  cor.rug.spearman <- cor.test(rugosity_data, speed_data, method="spearman")
	  cor.npp.spearman <- cor.test(npp_data, speed_data, method="spearman")
	  cat("var1\tvar2\ttest\tcoefficient\tp\n")
	  cat("altitude\tspeed\tpearson\t", round(cor.alt.pearson$estimate, 4), "\t", cor.alt.pearson$p.value, "\n")
	  cat("rugosity\tspeed\tpearson\t", round(cor.rug.pearson$estimate, 4), "\t", cor.rug.pearson$p.value, "\n")
	  cat("npp\tspeed\tpearson\t", round(cor.npp.pearson$estimate, 4), "\t", round(cor.npp.pearson$p.value, 4), "\n")
	  cat("altitude\tspeed\tspearman\t", round(cor.alt.spearman$estimate, 4), "\t", cor.alt.spearman$p.value, "\n")
	  cat("rugosity\tspeed\tspearman\t", round(cor.rug.spearman$estimate, 4), "\t", cor.rug.spearman$p.value, "\n")
	  cat("npp\tspeed\tspearman\t", round(cor.npp.spearman$estimate, 4), "\t", round(cor.npp.spearman$p.value, 4), "\n")
	  
	### bearings and subsistence ###
	} else if (what=="bearings_subsistence") {
	    # assign subtypes of hunter-gatherers to HG
	    for (i in 1:nrow(de)) {
	      if ( de$subsistence[i]=="HG-SED, HG-SAGO" ) { de$subsistence[i] <- "HG" }
	      if ( de$subsistence[i]=="HG-SAGO" ) { de$subsistence[i] <- "HG" }
	      if ( de$subsistence[i]=="HG-SED" ) { de$subsistence[i] <- "HG" }
	      if ( de$subsistence[i]=="HG-SAGO, HG-SED" ) { de$subsistence[i] <- "HG" }
	      if ( de$subsistence[i]=="HG-FISH" ) { de$subsistence[i] <- "HG" }
	    }
	    AGR_EW <- length(which(de$subsistence=="AGR" & de$bearing=="EW"))
	    AGR_NS <- length(which(de$subsistence=="AGR" & de$bearing=="NS"))
	    HG_EW <- length(which(de$subsistence=="HG" & de$bearing=="EW"))
	    HG_NS <- length(which(de$subsistence=="HG" & de$bearing=="NS"))
	    # bearing_ratios <- c(AGR_EW/AGR_NS, HG_EW/HG_NS)
	    bearing_proportions <- 100*c(AGR_EW/(AGR_EW+AGR_NS), HG_EW/(HG_EW+HG_NS))
	    # dev.off()
	    par(mfrow =c(1,2))
	    # bp <- barplot(bearing_ratios, ylab="ratio of EW to NS movements", col=c("lightblue", "pink"), 
	    #         ylim=c(0,1.6), space=c(1,.2), beside=TRUE, width=c(15,15))
	    n <- c(AGR_EW+AGR_NS, HG_EW+HG_NS)
	    moe95 <- 100*1.96*sqrt(((bearing_proportions/100)*(1-bearing_proportions/100))/n)
	    col_AGR <- create_col("blue", percent = 70, name = NULL)
	    col_HG <- create_col("red", percent = 70, name = NULL)
	    bp <- barplot(bearing_proportions, ylab="proportion (%) of EW to NS movements", col=c(col_AGR, col_HG), 
	           ylim=c(0,70), space=c(1,.2), beside=TRUE, width=c(15,15))
	    legend("top", c("AGR" ,"HG"), fill = c(col_AGR,col_HG))

	    # bp <- barplot(bearing_proportions, ylab="proportion (%) of EW to NS movements", col=c("lightblue", "pink"), 
	    #             ylim=c(0,70), space=c(1,.2), beside=TRUE, width=c(15,15))
	    # legend("top", c("AGR" ,"HG"), fill = c("lightblue","pink"))

	    abline(50, 0, col="black", lty=2)
	    arrows95 <- arrows(x0=bp, x1=bp, y0=bearing_proportions-moe95, y1=bearing_proportions+moe95, angle=90, code=3, length=.15)
	    ss <- paste0("n=",n)
	    text(bp, 1, ss, cex=1, pos=3)
	    # resetPar()
	    # do a chi-square test of differences between HG and AGR
	    chisq.test(c(AGR_EW,HG_EW),c(AGR_NS,HG_NS))
	
	} else {
		cat("not specified what to plot\n")
	}
}
