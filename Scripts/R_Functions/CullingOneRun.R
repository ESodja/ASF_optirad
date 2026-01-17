
CullingOneRun <- function(pop, idNEW, idZONE, Intensity, alphaC, centroids, Rad, inc, i, POSlive, POSdead, POSlive_locs, POSdead_locs, NEGlive, NEGdead, DetP, cullstyle){
	cells <- nrow(centroids)

	######################
	###### Get Zone ######
	######################

	#if there were infected pigs detected in the last time step
	if(length(idNEW)>0){
		#initiate vector to store pairs of each infected grid cell ID with every other cell ID
		#nrow for each combo, 3 cols
		#col 1-each grid cell ID with detected infection
		#col 2-each paired grid cell
		#col 3-distance between infected grid cell ID with detection and each grid cell
        ## i *think* this is causing a problem with integer overflow when things get too large (nrow is doing that somewhere in this file)
		pairedIDs <- matrix(nrow=cells*length(idNEW),ncol=3)
		#get matrix of each infected grid cell ID paired with every other cell ID
		if(length(idNEW)==1){
			pairedIDs[1:cells,1] <- idNEW
			pairedIDs[1:cells,2] <- 1:cells
		} else if(length(idNEW)>1){
			for(j in 1:length(idNEW)){
				if(j==1){
					pairedIDs[1:cells,1] <- idNEW[j]
					pairedIDs[1:cells,2] <- 1:cells
				}
				if(j>1){
					cells=nrow(centroids)
					pairedIDs[((1)+((j-1)*cells)):((j*cells)),1] <- idNEW[j]
					pairedIDs[((1)+((j-1)*cells)):((j*cells)),2] <- 1:cells
				}
			}

		}
		#get distance between each infected grid cell ID paired with every other cell ID
		#below results in matrix with one column of each grid cell id and the second column the distance to each infected
		## I think this could be handled with a cached matrix of distances between cells that is calculated once
		pairedIDs[,3] <- sqrt((centroids[pairedIDs[,1],1]-centroids[pairedIDs[,2],1])^2 + (centroids[pairedIDs[,1],2]-centroids[pairedIDs[,2],2])^2)
		#get all grid cells where the distance between that grid cell and an infected detected pig is less than the pre-determined radius
		idout <- pairedIDs[pairedIDs[,3]<=Rad,]

	} else { #else, if no infections detected in last time step...no new grid cells added
		idout <- matrix(nrow=0,ncol=3)
	}

	#remove NAs

	#fullZONE contains all grid cells with detected infections, from last time step and all before
	# first column is where it is detected, 2nd column is distance from detected point by each other cell
	fullZONE <- rbind(idZONE,idout)

	#get all unique grid cells in the zone
	allINzone <- unique(fullZONE[,2])
	allINzone <- allINzone[!is.na(allINzone)] ## clean out the rows of NA values, ain't nobody got time for that

	#get total number of grid cells in the zone
	Uall <- length(allINzone)

	#get total area of the zone
	ZONEkm2 <- Uall*inc^2

	#get which rows in the zone
	#total number of sounders=soundINzone
	soundINzone <- which(pop[,3]%in%allINzone)

	#get total number of pigs (live and dead) in zone
	pigsinzone <- sum(pop[soundINzone,1],pop[soundINzone,12],pop[soundINzone,13])

	#total number of exposed/infected/infectious carcass pigs in zone
	EICinzone <- sum(pop[soundINzone,9],pop[soundINzone,10],pop[soundINzone,12])

	#get total number of pigs outside the zone
	pigsoutzone <- sum(pop[-soundINzone,1],pop[-soundINzone,12],pop[-soundINzone,13])

	#get total number of infected pigs outside the zone
	EICoutzone <- sum(pop[-soundINzone,9],pop[-soundINzone,10],pop[-soundINzone,12])

	#get total number of  individuals (inside and outside zone)
	totalpigs <- sum(pop[,1],pop[,12],pop[,13])

	#get total number of infected individuals (inside and outside zone)
	totalEIC <- EICinzone+EICoutzone

	#For output to get landscape-level effective removal rate
	Ct <- pigsinzone/totalpigs

	#####################################
	###### Begin Culling Algorithm ######
	#####################################
	#if there are pigs to cull...
	# if(is.na(pigsinzone)){pigsinzone=0} ## if there are na's here, it means there are issues further up that probably shouldn't be glossed over...
	if(pigsinzone > 0){
		#get number of pigs for each grid cell in zone
		#and get their status, SEIRCZ
		#initiate empty matrix, nrow for each grid cell, 7 for each of SEIRCZ
		#tic()
	#     SEIRCZpigs<-matrix(0,nrow=nrow(fullZONE),ncol=7)
	#     fullZONEpigs<-cbind(fullZONE,SEIRCZpigs)
		popINzone <- pop[soundINzone,,drop=FALSE]
	#     for(u in 1:nrow(popINzone)){
	#         ## this doesn't work because popINzone can have multiple of the same value in the 'cell' column because things moved
	#         u_row<-which(fullZONEpigs[,2]==popINzone[u,3])
	#         fullZONEpigs[u_row,4]<-popINzone[u,1] #total number of pigs
	#         fullZONEpigs[u_row,5]<-popINzone[u,8] #total number susceptible pigs
	#         fullZONEpigs[u_row,6]<-popINzone[u,9] #total number exposted pigs
	#         fullZONEpigs[u_row,7]<-popINzone[u,10] #total number infected pigs
	#         fullZONEpigs[u_row,8]<-popINzone[u,11] #total number recovered pigs
	#         fullZONEpigs[u_row,9]<-popINzone[u,12] #total number infected carcasses
	#         fullZONEpigs[u_row,10]<-popINzone[u,13] #total number uninfected carcasses
	#     }
		## do we actually want the multiple values from fullZONE that indicate pigs from multiple cells converging on a single destination? (mostly a rhetorical question, unless the answer is yes)
		## The way it is written, a u_row that has two values will be written twice, and both rows will end up with the second value.
		## This results from having sounders in two locations end up in one location
		## The problem is that whatever sounder happens to be listed second is duplicated by the for loop
		## a cleaner/probably faster/more accurate way to do this would be to
		# make usable
        fullZONE.dt <- as.data.table(fullZONE)
        # clean out in-zone points that have multiple detection cells; keep the id and distance of the closest detection
        fullZONE.dt <- unique(fullZONE.dt[fullZONE.dt[,min(V3),by=V2], on=.(V2,V3=V1)])
        setnames(fullZONE.dt, c('V1','V2','V3'), c('det.cell','cell','dist'))
        # convert popINzone to a data table and aggregate state variables by summation
        popINzone.agg <- as.data.table(popINzone[,c(3,8:13)])[,lapply(.SD, sum), by=cell, .SDcols=2:7]
        # connect popINzone table with cell distances; clean out rows by cell id which have multiple detection points at the same distance (just keep one, doesn't matter which for this)
        fullZONEpigs <- fullZONE.dt[popINzone.agg, on=.(cell = cell)][,.SD[1],by=cell]
		## all the rows without pigs would be removed anyway below, so don't need all.x=TRUE
		## need to consolidate destination cells to one row (no repeat cell values in fullZONEpigs col. 1)


		#remove rows from fullZONEpigs without pigs
		## these should be redundant now given the changes above
	#     fullZONEpigs<-fullZONEpigs[complete.cases(fullZONEpigs),,drop=FALSE]
	#     fullZONEpigs<-fullZONEpigs[fullZONEpigs[,4]>0,,drop=FALSE]

		if (cullstyle == "startIN"){
			#Cullstyle start in, start with closest pigs from detections
		#     fullZONEpigs<-as.matrix(arrange(as.data.frame(fullZONEpigs),fullZONEpigs[,3]))
# 			fullZONEpigs <- fullZONEpigs[order(fullZONEpigs[,3]),,drop=FALSE] ## avoids copying dataframe
# 		} else if (cullstyle == "startOUT"){
			#Cullstyle start out, start with furthest pigs from detections
# 			fullZONEpigs <- fullZONEpigs[order(-fullZONEpigs[,3]),,drop=FALSE]
            fullZONEpigs[,cell.total := S+E+I+R+C+Z]
            fullZONEpigs[,cprob := 1 - 1/((1 + ..alphaC)^(cell.total/(..inc^2)))]
            norm.val <- TwoDt(0, 3, Rad/2)
            fullZONEpigs[,dist.weight := 0.5*TwoDt(dist, 3, ..Rad/2)/..norm.val]
        #     fullZONEpigs <- cbind(fullZONEpigs, prob.cull = Intensity*TwoDt(fullZONEpigs[,'dist'], 3, Rad/2)/norm.val)
#             fullZONEpigs <- cbind(fullZONEpigs, cprob, dist.weight, rowSums(fullZONEpigs[,4:9]))
#             names(fullZONEpigs)[ncol(fullZONEpigs)] <- 'cell.total'

            fullZONEpigs[,cull.cell := rbinom(nrow(fullZONEpigs), cell.total, cprob * dist.weight)]
            setorder(fullZONEpigs, dist, -cell.total)
            fullZONEpigs[cull.cell > 0,cum.cull := round(cumsum(cell.total)/fullZONEpigs[,sum(cell.total)],2)]
            fullZONEpigs[cull.cell > 0 & cum.cull <= 0.05, actual.cull := 1]
		}
	#     fullZONEpigs<-fullZONEpigs[complete.cases(fullZONEpigs),,drop=FALSE] ## probably don't need this unless bugs are making NA's somewhere

		#%density of all live and dead pigs in the zone
	#     Dr=pigsinzone/ZONEkm2
	#     print(paste('Dr', Dr))
		#determine density-dependent capture probability in this radius
	#     cprob=1-(1/(1+alphaC)^Dr) ## this is stupid
# 		pigdensity <-
# 		browser()
		#get total number culled/removed/sampled in the zone
	#     numb <- rbinom(pigsinzone,1,cprob*Intensity) ## this is stupid
# 		browser()
# 		fullZONEpigs.agg <- cbind(fullZONEpigs.agg, cprob*dist.weight, rbinom(nrow(fullZONEpigs.agg), 1, cprob*Intensity*dist.weight)) ## this is stupid
# 		fullZONEpigs <- merge(fullZONEpigs, fullZONEpigs.agg[,c(1:3,10:13)], by.x=c('V2','V1','V3'), by.y=c('loc.cell','detect.cell','dist'))
		## there is an issue here -- as the densisty of pigs in the zone goes down, the number of pigs you will cull also goes down, which makes sense. BUT because the area is expanding as you initially find more, the density is also decreasing. Eventually you've killed all the pigs that are infected within your zone more or less, but you also have live ones moving in. Eventually you start killing uninfected pigs moving in to the region near your detection points, while not getting any more detections because when you start nearest the detection and work outwards you end up with a limited number of target culled individuals because the area is so big and the population is so small, and the disease becomes so rare in that area that you stop detecting any more. maybe.
		## I would say it would make more sense to base the probability of detection off of distance from previous detection points, because the lowered density (between pop decline and area increase) reduces the number of individuals you are culling faster than the number of individuals are declining. Instead, have each sounder's distance from previous detection and, possibly, the size of the sounder, indicate its probability of detection and if that happens, it gets culled.
		## Even if you just randomly ran around in the zone and detected sounders independent of their exact distance from detection, you could have a better result because you would prbably still have positive detections at least occasionally...
		## applying the probability to the density of everything in the radius is kindof... silly...
		## I think the units are kindof screwey, like why does alphaC or intensity fit in the way that they do?


		#get cumulative sum of targeted pigs
		#cpigs = cumsum(tpigs);
	#     cpigs=sum(numb)
	#     print(paste('num culled', cpigs))

		#determine how far down the list to remove pigs from cells
		# removals=0 #total number of removals, go through loop until first time it is equal to or greater than cpigs
		## This adds an additional row (because incr=incr+1 is at the end of the loop)
		# incr=1 #row number where culling stops
		# incr=0 ## if we use cumsum structure below, this doesn't matter
		# if(cpigs>0&nrow(fullZONEpigs)>0){
		#     while(removals<cpigs&incr<nrow(fullZONEpigs)){
		#         removals<-removals+fullZONEpigs[incr,4]
		#         incr=incr+1
		#     }
# 		fullZONEpigs <- cbind(fullZONEpigs, cull.cell = rbinom(fullZONEpigs[,'prob.cull'], 1, fullZONEpigs[,'prob.cull']))
		## Problem: No randomness in selection and bias toward higher? index numbers drifts detection to the south (or north, however plotted)
		if (sum(fullZONEpigs[,actual.cull], na.rm=TRUE) > 0 & nrow(fullZONEpigs) > 0){
            print('some culling and pigs in zone')
			## use < instead of <= so if they hit cpigs exactly, they will stop there
	#         cull.index <- c(1L, which(cumsum(fullZONEpigs[,4]) < cpigs) + 1L) ## go to the next row b/c that's when they would stop
	#         cull.index <- cull.index[seq(min(length(cull.index), nrow(fullZONEpigs)))] ## handles cases with more culling than pigs
			## randomizes which points are selected, not just the n closest/furthest ones with arbitrary ordering (i.e. sorted by cell id number) within the same distance group, which potentially caused bias in culling
			## also with more individuals in that distance than what would normally be culled, allows for missing some sounders which seems more realistic in some ways...
	#         cull.distance <- max(fullZONEpigs[cull.index,3])
	#         print(paste('cull dist', cull.distance))
	#         browser()
	#         if (cullstyle == "startIN") {
	#             cull.index2 <- c(1L, which(fullZONEpigs[,3] <= cull.distance) + 1L)
	#         } else {
	#             cull.index2 <- c(1L, which(fullZONEpigs[,3] >= cull.distance) + 1L)
	#         }
	#         culled <- sum(fullZONEpigs[cull.index2,4])
			culled <- sum(fullZONEpigs[actual.cull == 1, cell.total])
		#determine which pigs culled
		#     culled=removals[[1]] ## ??? shouldn't be a list...
		#% list of cells that pigs will be eliminated from (column index was 1 in old version)

            removalcells <- as.matrix(fullZONEpigs[actual.cull == 1, cell])
		#get which pigs culled
            removalpigs <- as.matrix(fullZONEpigs[actual.cull == 1,])

		######################################
		###### Update surveillance data ######
		######################################

		#POSlive_i is a matrix with a row for each timestep
		#column one of poslive is the number of exposed/infected pigs detected at that timestep
		#sum removalpigs column 6,7
			POSlive_i <- sum(removalpigs[,7],removalpigs[,6])
			if(!is.null(DetP) & DetP != 1){
	#             POSlive_i_sel=rbinom(POSlive_i,1,DetP)
				POSlive_i_sel <- rbinom(nrow(removalpigs),rowSums(removalpigs[,6:7]), DetP)
				## this is used below to choose rows of removalpigs for where positive locs are, which isn't what POSlive_i_sel is actually counting
				POSlive_i <- sum(removalpigs[POSlive_i_sel,6:7])
			} else if(DetP == 1){
				POSlive_i_sel <- as.numeric(which(rowSums(removalpigs[,6:7,drop=FALSE]) > 0))
			}
		#POSdead
		#POSdead is a matrix with a row for each timestep
		#column one of poslive is the number of infected carcasses detected at that timestep
		#sum removalpigs column 9
			POSdead_i <- sum(removalpigs[,9])
			if(!is.null(DetP) & DetP != 1){
				POSdead_i_sel <- rbinom(nrow(removalpigs), removalpigs[,9], DetP)
				POSdead_i <- sum(removalpigs[POSdead_i_sel,9])
			} else if (DetP==1){
				POSdead_i_sel  <-  as.numeric(which(removalpigs[,9] > 0))
			}
		#list of length thyme, each timestep is vector of grid cell locations where live infected pigs detected at that ts
		#removalpigs col 2 where column 6 or 7 > 0
			if(POSlive_i>0){
	#             lll<-length(removalpigs[removalpigs[,6]>0|removalpigs[,7]>0,2])
		#POSlive_locs_i=vector(mode="integer",length=lll)
		#POSlive_locs[[i]]<-removalpigs[removalpigs[,6]>0|removalpigs[,7]>0,2]
				POSlive_locs_i <- removalpigs[POSlive_i_sel,1]
	#             if(!is.null(DetP)){
	#                 POSlive_locs_i=removalpigs[POSlive_i_sel==1,2]
	# #                 POSlive_locs_i=POSlive_locs_i[POSlive_i_sel==1]
	#             }

			} else {POSlive_locs_i <- 0}

			if(POSdead_i>0){
				POSdead_locs_i <- removalpigs[POSdead_i_sel,1]

	#             if(!is.null(DetP)){
	#                 POSdead_locs_i=POSdead_locs_i[POSdead_i_sel==1]
	#             }
	#
			} else {POSdead_locs_i <- 0}

			#vector of nrow timestop, count of total SR removed
			#sum removalpigs column 5,8
			NEGlive_i <- sum(removalpigs[,5],removalpigs[,8])
			if(!is.null(DetP) & DetP != 1){
				NEGlive_i_missed <- length(POSlive_i_sel[POSlive_i_sel==0])
				NEGlive_i <- NEGlive_i+NEGlive_i_missed
			}
			#vector of nrow timestep, count of total Z removed
			#sum removalpigs column 10
			NEGdead_i <- sum(removalpigs[,10])
			if(!is.null(DetP) & DetP != 1){
				NEGdead_i_missed <- length(POSdead_i_sel[POSdead_i_sel==0])
				NEGdead_i <- NEGdead_i+NEGdead_i_missed
			}

			#remove removed sounders from pop
			removalrows <- which(pop[,3] %in% removalcells)
			removedpop <- pop[-removalrows,,drop=FALSE]

		} else{
            ## this should never happen?
			POSlive_i <- 0
			POSdead_i <- 0
			POSlive_locs_i <- NA
			POSdead_locs_i <- NA
			NEGlive_i <- 0
			NEGdead_i <- 0
			culled <- 0
# 			removedpop <- NA
			Ct <- 0
			removedpop <- pop
		}
        #idZONE:
        #grid cell ids that had a positive detection, grid cell ids that are within the zone, distance
        idZONE <- fullZONE

		#send updated objects to output list
		output.list <- vector(mode="list",length=13)
		output.list[[1]] <- POSlive_i
		output.list[[2]] <- POSdead_i
		output.list[[3]] <- POSlive_locs_i
		output.list[[4]] <- POSdead_locs_i
		output.list[[5]] <- NEGlive_i
		output.list[[6]] <- NEGdead_i
		output.list[[7]] <- idZONE
		output.list[[8]] <- removalcells
		output.list[[9]] <- culled
		output.list[[10]] <- ZONEkm2
		output.list[[11]] <- removedpop
		output.list[[12]] <- Ct
		output.list[[13]] <- allINzone
	} else {
		output.list <- vector(mode="list",length=13)
		output.list[[1]] <- 0
		output.list[[2]] <- 0
		output.list[[3]] <- 0
		output.list[[4]] <- 0
		output.list[[5]] <- 0
		output.list[[6]] <- 0
		output.list[[7]] <- idZONE
		output.list[[8]] <- 0
		output.list[[9]] <- 0
		output.list[[10]] <- 0
		output.list[[11]] <- pop
		output.list[[12]] <- 0
		output.list[[13]] <- allINzone
	}

	return(output.list)

} #function closing bracket

