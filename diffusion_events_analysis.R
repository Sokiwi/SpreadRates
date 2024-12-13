library(cartography)  # carto.pal
library(RColorBrewer)  # for making plots of diffusion rates in different biomes
library(ddalpha)

# function looks at periods of Y years from "lowest" (default is
# 6000 BP to present 
# within each period it checks the speed of active diffusion events
# and then takes means and medians, outputting a file with that and N
# Can plot world, areas, biomes

# areas are: Africa, Australia, Austronesian, Eurasia, North America, Papuan, South America

# plots and analyses
# mean rates across areas
as(Y=1, lowest=7000, loess.span=.1, upper=2, present=300, what="areas")  # Fig. 1
# ratio of EW to SW movements across areas
# presently this only works when the corresponding code is run from inside the function
as(what="biomes")  # Fig. 4
as(what="continuous_variables")  # prints correlations between speed and some other variables,
as("bearings_areas")  # Fig. S1
# ratio of EW to SW movements across subsistence
# presently this only works when the corresponding code is run from inside the function
# mean rates across subsistence patterns
as(Y=1, lowest=5000, loess.span=.2, upper=0.8, present=300, what="subsistence")  # Fig. S1
as("bearings_subsistence")  # Fig. S2
# as("latitudes")
# as("area")

# examples of other plots and analyses (not used)
# as(Y=1, lowest=7000, loess.span=.1, what="areas", upper=10, present=0)  # means across world areas
# as(Y=1, lowest=7000, loess.span=.1, what="world")  # mean and median world, not used
# as(300, 5000, 1, "world")  # mean and median world with fitting, not used
# specifically: altitude, rugosity, and net primary productivity
# as(Y=1, lowest=7000, loess.span=.1, upper=5, present=200, what="areas")
# as(Y=1, lowest=7000, loess.span=.2, what="area", upper=1, repetitions=10, area="Africa")
# as(Y=1, lowest=5000, loess.span=.2, what="area", upper=6, repetitions=1, area="Eurasia")
# as(Y=1, lowest=6000, loess.span=.2, what="area", upper=2, repetitions=1, "South America")
# as(Y=1, lowest=6000, loess.span=.2, what="area", upper=2, repetitions=1, "Austronesian")
# as(Y=1, lowest=6000, loess.span=.2, what="area", upper=0.4, repetitions=1, "Australia")
# as(Y=1, lowest=6000, loess.span=.2, what="area", upper=0.4, repetitions=1, "North America")
# as(Y=1, lowest=7000, loess.span=.2, what="area", upper=0.1, repetitions=1, "Papuan")


# as stands for average speed
as <- function(Y=200, lowest=7000, loess.span=.5, what="world", upper=10, present=300, repetitions=10, area, ...) {

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

	### areas ###
  if ( what=="areas" ) {
    areas <- sort(unique(de$area))
    w_uniden <- which(areas=="unidentified")
    if ( length(w_uniden) > 0 ) {
      areas <- areas[-w_uniden]
    }
    nareas <- length(areas)
    colors <- c("blue","red","green3","darkorchid1","orange","maroon1","lightsalmon3")
    linetypes <- c(1,1,1,1,1,1,1)
    plot(100, 100, xlim=c(-lowest,0), ylim=c(0,upper), xlab="BP", ylab="km/year")
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "#EBEBEB")
    abline(v=seq(-lowest,0,200), h=seq(0,upper,0.1), col="white")
    legend((-lowest), upper, areas, cex=1, col=colors, lty=linetypes, lwd=3)
    for ( a in 1:length(areas) ) {
      filename <- paste("diffusion_averages_", areas[a], "_", Y, "_year_period.txt", sep="")
      cat("lower_bound\taverage_speed\tmedian_speed\tN\n", file=filename)
      de_sub <- de[which(de$area==areas[a]),]
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
          speed_median <- median(as.numeric(de_sub$speed[w_in]))
          cat(lower, "\t", speed_average, "\t", speed_median, "\t", length(w_in), "\n", sep="")
          cat(lower, "\t", speed_average, "\t", speed_median, "\t", length(w_in), "\n", sep="", file=filename, append=TRUE)
        }
        lower <- lower - Y
      }
      plot_data <- read.table(file=filename, header=TRUE, sep="\t")
      w_na2 <- unique(union(which(is.na(plot_data$average_speed)), which(is.na(plot_data$median_speed))))
      if ( length(w_na2) > 0 ) {
        plot_data <- plot_data[-w_na2,]
      }
      plot_data$lower_bound <- -1*plot_data$lower_bound
      lines(predict(loess(plot_data$average_speed ~ plot_data$lower_bound, span=loess.span)), x=plot_data$lower_bound, col=colors[a], lty=linetypes[a], lwd=4)
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
	  # as these values are changed a new dataframe is defined
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
			  # w_rawA_zero <- which(rawA==0)
			  # rawApruned <- rawA
			  # if (length(w_rawA_zero) > 0) {
			  #   rawApruned[w_rawA_zero] <- 0.001
			  # }
				# A <- log(rawApruned)
				
				rawB <- as.numeric(x$speed[which(x$biome4==biomes[j])])
				# w_rawB_zero <- which(rawB==0)
				# rawBpruned <- rawB
				# if (length(w_rawB_zero) > 0) {
				#   rawBpruned[w_rawB_zero] <- 0.001
				# }
				# B <- log(rawBpruned)

				res <- t.test(rawA, rawB, alternative="two.sided")
				p <- res$p.value
				  cat(designations[i], "\t", designations[j], "\t", mean(rawA), "\t", mean(rawB), "\t", p, "\n", file="t tests biomes.txt", append=TRUE)
			}
		}

		# get rates and N for total distances traveled in different biomes 
		# (using the 90% criterion)
		# collect data in vectors to be merged with biome_keys
		N <- c()
		km <- c()
		years <- c()
		biome_rates <- c()
		for ( i in 1:length(biomes) ) {
		  ind <- which(x$biome4==biomes[i])
		  N[i] <- length(ind)
		  km[i] <- sum(x$distance[ind])
		  years[i] <- sum(x$BP_mo[ind] - x$BP_da[ind])
		  biome_rates[i] <- km[i]/years[i]
		}
		
		biome_rates_ranked <- biome_rates[order(biome_rates, decreasing=TRUE)]
		biomes_ranked <- biomes[order(biome_rates, decreasing=TRUE)]
		biome_descriptions_ranked <- biome_key$description[match(biomes_ranked, biome_key$key)]
		N_ranked <- N[order(biome_rates, decreasing=TRUE)]
		
		biome_results_full <- data.frame(biome_rates_ranked, biomes_ranked, biome_descriptions_ranked, N_ranked)
    biome_results <- biome_results_full[biome_results_full$N_ranked > 10,]
				
		# plot mean diffusion rates across biomes
		# colors that are just distinct, legend generated through defaults
		pal <- brewer.pal(n = nrow(biome_results), name = "Paired")
		barplot(biome_results$biome_rates_ranked, ylim=c(0,1), 
		        legend.text=biome_results$biome_descriptions_ranked, 
		        col=pal, ylab="km/year")
		
  ### subsistence ###
  } else if ( what=="subsistence") {
		# assign subtypes of hunter-gatherers to HG
		for (i in 1:nrow(de)) {
			if ( de$subsistence[i]=="HG-SED, HG-SAGO" ) { de$subsistence[i] <- "HG" }
			if ( de$subsistence[i]=="HG-SAGO" ) { de$subsistence[i] <- "HG" }
			if ( de$subsistence[i]=="HG-SED" ) { de$subsistence[i] <- "HG" }
			if ( de$subsistence[i]=="HG-SAGO, HG-SED" ) { de$subsistence[i] <- "HG" }
			if ( de$subsistence[i]=="HG-FISH" ) { de$subsistence[i] <- "HG" }
		}
	  # plot mean rates across subsistence patterns
	  dev.off()
		colors <- c("lightblue", "pink")
		linetypes <- c(1,1)
		subsist <- c("AGR", "HG")
		plot(100, 100, xlim=c(-lowest,0), ylim=c(0,upper), xlab="BP", ylab="km/year")
		rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "#EBEBEB")
		abline(v=seq(-lowest,0,200), h=seq(0,upper,0.1), col="white")
		legend(-lowest, upper, subsist, cex=1.3, col=colors, lty=linetypes, lwd=3)
		for ( a in 1:2 ) {
			filename <- paste("diffusion_averages_", subsist[a], "_", Y, "_year_period.txt", sep="")
			cat("lower_bound\taverage_speed\tmedian_speed\tN\n", file=filename)
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
					speed_median <- median(as.numeric(de_sub$speed[w_in]))
					cat(lower, "\t", speed_average, "\t", speed_median, "\t", length(w_in), "\n", sep="")
					cat(lower, "\t", speed_average, "\t", speed_median, "\t", length(w_in), "\n", sep="", file=filename, append=TRUE)
				}
				lower <- lower - Y
			}
			plot_data <- read.table(file=filename, header=TRUE, sep="\t")
			w_na2 <- unique(union(which(is.na(plot_data$average_speed)), which(is.na(plot_data$median_speed))))
			if ( length(w_na2) > 0 ) {
				plot_data <- plot_data[-w_na2,]
			}
			plot_data$lower_bound <- -1*plot_data$lower_bound
			lines(predict(loess(plot_data$average_speed ~ plot_data$lower_bound, span=loess.span)), 
			      x=plot_data$lower_bound, col=colors[a], lty=linetypes[a], lwd=4)
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
    bearing_ratios_areas <- mab[1,]/mab[2,]
    print(bearing_ratios_areas)
    # make barplot similar to the one for bearings_subsistence but 
    # hardwired colors same as those of diffusion rates
    abbr1 <- names(bearing_ratios_areas)
    abbr2 <- gsub("North", "N.", abbr1)
    abbr3 <- gsub("South", "S.", abbr2)
    names(bearing_ratios_areas) <- abbr3
    colors <- c("blue","red","green3","darkorchid1","orange","maroon1","lightsalmon3")
	  barplot(bearing_ratios_areas, ylab="ratio of EW to NS movements", 
	          col=colors, xlim=NULL, ylim=c(0,1.6))
	  abline(1, 0, col="black", lty=2)
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
	    bearing_ratios <- c(AGR_EW/AGR_NS, HG_EW/HG_NS)
	    # dev.off()
	    par(mfrow =c(1,2))
	    barplot(bearing_ratios, ylab="ratio of EW to NS movements", col=c("lightblue", "pink"), 
	            ylim=c(0,1.6), space=c(1,.2), beside=TRUE, width=c(15,15))
	    # text(bearing_ratios + .08, paste("N =", c(sum(AGR_EW, AGR_NS), sum(HG_EW, HG_NS))), cex=1) 
	    legend("top", c("AGR" ,"HG"), fill = c("lightblue","pink"))
	    abline(1, 0, col="black", lty=2)
	    resetPar()
	    # do a chi-square test of differences between HG and AGR
	    chisq.test(c(AGR_EW,HG_EW),c(AGR_NS,HG_NS))
	
	    ### area ###
	} else if ( what=="area") {
	  de_sub <- de[which(de$area==area),]
	  plotfile.name <- paste0(area, "_with_", repetitions, "_calibration_noise_curves.pdf")
	  pdf(plotfile.name)
	  plot(100, 100, xlim=c(-lowest,min(de_sub$BP_da)), ylim=c(0,upper), main=area, xlab="BP", ylab="km/year")
	  # the following plots lines for dates sampled within their calibration error ranges
	  for (r in 1:repetitions) {
	    filename <- paste("diffusion_averages_", area, "_", Y, "_year_period.txt", sep="")
	    cat("lower_bound\taverage_speed\tmedian_speed\tN\n", file=filename)
	    lower <- lowest
	    while ( lower >= Y ) {
	      BPM <- as.numeric(de_sub$BP_mo)
	      BPD <- as.numeric(de_sub$BP_da)
	      BPM_new <- rep(NA, length(BPM))
	      BPD_new <- rep(NA, length(BPD))
	      new <- function(M, D) {
	        M_low <- round(M - .29*M)
	        M_hi <- round(M + .29*M)
	        M_new <- sample(M_low:M_hi, 1)
	        D_low <- round(D - .29*D)
	        D_hi <- round(D + .29*D)
	        if (D_hi >= M_new) {
	          D_hi <- M_new - 1
	        }
	        D_new <- sample(D_low:D_hi, 1)
	        return(c(M_new, D_new))
	      }
	      for (i in 1:length(BPM)) {
	        new_out <- new(BPM[i],BPD[i])
	        BPM_new[i] <- new_out[1]
	        BPD_new[i] <- new_out[2]
	      }
	      window.edge <- lower - Y
	      w_in <- which( !((BPM_new < lower & BPD_new < window.edge) | (BPM_new > lower & BPD_new > window.edge)) )
	      if ( length(w_in)==0 ) {
	        cat(lower, "\t", NA, "\t", NA, "\t", length(w_in), "\n", sep="")
	        cat(lower, "\t", NA, "\t", NA, "\t", length(w_in), "\n", sep="", file=filename, append=TRUE)
	      } else {
	        speed_average <- round(mean(as.numeric(de_sub$distance[w_in])/(BPM_new[w_in] - BPD_new[w_in])),3)
	        # speed_average <- mean(as.numeric(de_sub$speed[w_in]))
	        speed_median <- round(median(as.numeric(de_sub$distance[w_in])/(BPM_new[w_in] - BPD_new[w_in])),3)
	        cat(lower, "\t", speed_average, "\t", speed_median, "\t", length(w_in), "\n", sep="", file=filename, append=TRUE)
	      }
	      lower <- lower - Y
	    }
	    plot_data <- read.table(file=filename, header=TRUE, sep="\t")
	    w_na2 <- unique(union(which(is.na(plot_data$average_speed)), which(is.na(plot_data$median_speed))))
	    if ( length(w_na2) > 0 ) {
	      plot_data <- plot_data[-w_na2,]
	    }
	    plot_data$lower_bound <- -1*plot_data$lower_bound
	    lines(predict(loess(plot_data$average_speed ~ plot_data$lower_bound, span=loess.span)), x=plot_data$lower_bound, lty=1, lwd=1)
	    par(new=TRUE)
	  }
	  # the following plots the original, unsampled data
	  filename <- paste("diffusion_averages_", area, "_", Y, "_year_period.txt", sep="")
	  cat("lower_bound\taverage_speed\tmedian_speed\tN\n", file=filename)
	  de_sub <- de[which(de$area==area),]
	  lower <- lowest
	  while ( lower >= Y ) {
	    BPM <- as.numeric(de_sub$BP_mo)
	    BPD <- as.numeric(de_sub$BP_da)
	    window.edge <- lower - Y
	    w_in <- which( !((BPM < lower & BPD < window.edge) | (BPM > lower & BPD > window.edge)) )
	    if ( length(w_in)==0 ) {
	      cat(lower, "\t", NA, "\t", NA, "\t", length(w_in), "\n", sep="")
	      cat(lower, "\t", NA, "\t", NA, "\t", length(w_in), "\n", sep="", file=filename, append=TRUE)
	    } else {
	      speed_average <- mean(as.numeric(de_sub$speed[w_in]))
	      speed_median <- median(as.numeric(de_sub$speed[w_in]))
	      cat(lower, "\t", speed_average, "\t", speed_median, "\t", length(w_in), "\n", sep="", file=filename, append=TRUE)
	    }
	    lower <- lower - Y
	  }
	  plot_data <- read.table(file=filename, header=TRUE, sep="\t")
	  w_na2 <- unique(union(which(is.na(plot_data$average_speed)), which(is.na(plot_data$median_speed))))
	  if ( length(w_na2) > 0 ) {
	    plot_data <- plot_data[-w_na2,]
	  }
	  plot_data$lower_bound <- -1*plot_data$lower_bound
	  lines(predict(loess(plot_data$average_speed ~ plot_data$lower_bound, span=loess.span)), x=plot_data$lower_bound, lty=1, lwd=3, col="red")
	  grid(nx = NA,
	       ny = NULL,
	       lty = 2, col = "gray", lwd = 2)
	  rug(-1*de_sub$BP_mo)
	  dev.off()
	  
	### world ###
	} else if ( what=="world" ) {
	      lower <- lowest
	      while ( lower >= Y ) {
	        # w_in <- which(as.numeric(de$BP_M) < lower & as.numeric(de$BP_D) > (lower - Y))
	        BPM <- as.numeric(de$BP_mo)
	        BPD <- as.numeric(de$BP_da)
	        window.edge <- lower - Y
	        w_in <- which( !((BPM < lower & BPD < window.edge) | (BPM > lower & BPD > window.edge)) )
	        if ( length(w_in)==0 ) {
	          # cat(lower, "\t", NA, "\t", NA, "\t", length(w_in), "\n", sep="")
	          cat(lower, "\t", NA, "\t", NA, "\t", length(w_in), "\n", sep="", file=filename, append=TRUE)
	        } else {
	          speed_average <- mean(as.numeric(de$speed[w_in]))
	          speed_median <- median(as.numeric(de$speed[w_in]))
	          # cat(lower, "\t", speed_average, "\t", speed_median, "\t", length(w_in), "\n", sep="")
	          cat(lower, "\t", speed_average, "\t", speed_median, "\t", length(w_in), "\n", sep="", file=filename, append=TRUE)
	        }
	        lower <- lower - Y
	      }
	      plot_data <- read.table(file=filename, header=TRUE, sep="\t")
	      w_na2 <- unique(union(which(is.na(plot_data$average_speed)), which(is.na(plot_data$median_speed))))
	      if ( length(w_na2) > 0 ) {
	        plot_data <- plot_data[-w_na2,]
	      }
	      plot_data$lower_bound <- -1*plot_data$lower_bound
	      par(mfrow = c(2,1))
	      plot(plot_data$lower_bound, plot_data$average_speed, main="Mean", xlab="BP", ylab="km/year")
	      lines(predict(loess(plot_data$average_speed ~ plot_data$lower_bound, span=loess.span)), x=plot_data$lower_bound, col='red')
	      abline(h = 1, col="blue", lty=2)
	      plot(plot_data$lower_bound, plot_data$median_speed, main="Median", xlab="BP", ylab="km/year")
	      lines(predict(loess(plot_data$median_speed ~ plot_data$lower_bound, span=loess.span)), x=plot_data$lower_bound, col='red')
	      abline(h = 1, col="blue", lty=2)
	      de <- de[order(de$BP_mo, decreasing=TRUE),]
	      write.table(de, sep=",", file="migration_events_data.csv")
	      
  ### latitudes ###
	} else if (what=="latitudes") {
	  # check if diffusion rates are related to latitude
	  touch <- function(latw1, latw2, latt1, latt2) {
	    if (latt1 < latw1 & latt2 < latw1) {return(FALSE)}
	    else if (latt1 > latw2 & latt2 > latw2) {return(FALSE)}
	    else {return(TRUE)}
	  }
	  max_lat <- ceiling(max(c(de$lat_da, de$lat_mo)))
	  min_lat <- floor(min(c(de$lat_da, de$lat_mo)))
	  lat <- c(min_lat:(max_lat-1))
	  L <- length(lat)
	  speed <- rep(NA, L)
    speed_table <- data.frame(lat, speed)
	  speed_list <- list()
	  counter <- 0
	  for (i in min_lat:(max_lat-1)) {
	    counter <- counter + 1
	    if (counter %% 10 == 0) {
	      cat("doing", counter, "out of", max_lat - min_lat, "\n")
	    }
	    latw1 <- i
	    latw2 <- i+1
	    for (j in 1:nrow(de)) {
	      lat_mo <- de$lat_mo[j]
	      lat_da <- de$lat_da[j]
        sorted <- sort(c(lat_mo, lat_da))
	      latt1 <- sorted[1]  # latt stands for target latitude
	      latt2 <- sorted[2]
	      touches <- touch(latw1, latw2, latt2, latt2)
	      if (touches) {
	        if (length(speed_list) >= counter) {
	          speed_list[[counter]] <- c(speed_list[[counter]], de$speed[j])
	        } else {
	          speed_list[[counter]] <- de$speed[j]
	        }
	      }
	    }
	  }
	  cat("\nnow averaging rates\n\n")
	  for (k in 1:length(speed_list)) {
	    if (length(speed_list[[k]]) > 0) {
	      speed_table$speed[k] <- mean(speed_list[[k]])
	    }
	  }
	  southern_lats <- which(speed_table$lat < 0)
	  northern_lats <- which(speed_table$lat >= 0)
	  southern <- speed_table[southern_lats,]
	  southern$lat <- southern$lat * -1
	  northern <- speed_table[northern_lats,]
	  plot(northern$lat, northern$speed, xlab="latitudes north", ylab="km/year")
	  abline(lm(northern$speed ~ northern$lat))
	  plot(southern$lat, southern$speed, xlab="latitudes south", ylab="km/year")
	  abline(lm(southern$speed ~ southern$lat))
	} else {
		cat("not specified what to plot\n")
	}
}
