
### Initialize grid(s): ---------------
    #Method for class 'list'
        #Initialize_Grids(object, parameters=parameters, movement=c(parameters$shape,parameters$rate))
    #Arguments
				#movement
					#default is vector with gamma distribution shape and rate fed from parameters file

      #grid.opts- 
        #"homogenous" or "heterogeneous", default for class numeric is "homogenous" and default for type SpatRaster and SpatRasterCollection is "ras"
          #ras- 
            #use input raster to create grid
          #homogeneous-
            #creates grid with 7 columns, no land class variables
          #heterogeneous-
            #creates a neutral random landscape model with X lc variables
    #Value
      #a nested list of grid parameters
RunSimulationReplicates<-function(land_grid_list, parameters, variables, cpp_functions, reps){

print('runsimreps')
# options(warn=2)
	## Filters variables (parameters with >1 value) out of parameters
	## selecting variables is done in SetVarParms.R
	parameters <- parameters[names(parameters) %in% names(variables) == FALSE]

	#Pull needed parms from parameters for all reps
	list2env(parameters, .GlobalEnv)
	#Need nested loops:
	#1, loop through all landscapes
	#2, loop through all parameter settings
for(v in 1:nrow(variables)){
	print(paste0("variable ",v))
	vars=variables[v,]
# 	names(vars)[3]="dens"
	names(vars)[names(vars) == "density"]="dens" ## more robust, especially as variables defined as parameters with multiple values
	vars=as.list(vars)
	list2env(vars, .GlobalEnv)

	#calc vals based on variables
	## could set these as variables/parameter options
	N0=dens*area
	K=N0*1.5

	#loop through landscapes
	for(l in 1:length(land_grid_list)){
		print(paste0("landscape: ",l))
		centroids=land_grid_list[[l]]$centroids
		grid=land_grid_list[[l]]$grid

		pop=InitializeSounders(centroids,grid,c(N0,ss),pop_init_grid_opts)
		outputs=Initialize_Outputs(parameters)
		pop=InitializeInfection(pop,centroids,grid,parameters)

		for(r in 1:reps){
			print(r)
			#Do simulations
			out.list=SimulateOneRun(outputs,pop,centroids,grid,parameters,cpp_functions,K)
			#Handle outputs

			## test these outputs if out.opts doesn't include them
			## id by variable combination, landscape, rep, and timestep for one row per timestep data
			end.tm = out.list$endtime
			id.r = c(v,l,r)
			tm.mat.r = cbind(matrix(id.r, nrow=end.tm, ncol=3, byrow=TRUE), seq(end.tm))
			colnames(tm.mat.r) = c("var","land","rep","timestep")

			#Handle effective removal rate (timestep output)
# 			Ct.r=out.list$Ct
# 			Ct.r=cbind(matrix(id.r, ncol=3, nrow=nrow(Ct.r), byrow=TRUE), 1:thyme, Ct.r)
# 			Ct.r=cbind(rep(r,times=nrow(Ct.r)),Ct.r)
# 			Ct.r=cbind(rep(l,times=nrow(Ct.r)),Ct.r)
# 			Ct.r=cbind(rep(v,times=nrow(Ct.r)),Ct.r)
			tm.mat.r <- cbind(tm.mat.r, out.list$Ct[seq(end.tm)])
			colnames(tm.mat.r)[ncol(tm.mat.r)]="Ct"

			# Births
# 			BB.r = out.list$BB
			tm.mat.r = cbind(tm.mat.r, out.list$BB[seq(end.tm)])
			colnames(tm.mat.r)[ncol(tm.mat.r)] <- 'BB'

			#Handle sounderlocs (optional...)
			## seems like other things were supposed to happen in sounderlocsSummarize, if we want to use those this will have to change
# 			browser()
			if ("sounderlocs" %in% out.opts){
				solocs.r=sounderlocsSummarize(out.list$sounderlocs,r)[[1]]
	# 			solocs.r=solocs.r$SEIRCZ_total
	# 			solocs.r=cbind(rep(l,times=nrow(solocs.r)),solocs.r)
	# 			solocs.r=cbind(rep(v,times=nrow(solocs.r)),solocs.r)
				tm.mat.r = cbind(tm.mat.r, solocs.r[3:8])
			}

			# detections (optional...)
			## has a row for each timestep AND detection type, with timestep, code (1=live,0=dead), number of individuals detected, and position
			## still need to test sample = 1
			## might match with incidence?
			if (end.tm > detectday & 'alldetections' %in% out.opts){
				detections.r = out.list$alldetections
				detections.r <- detections.r[detections.r[,1] <= end.tm & detections.r[,1] >= detectday,]
				n.det = nrow(detections.r)
				detections.r = cbind(matrix(id.r, ncol=3, nrow=n.det, byrow=TRUE), detections.r)
				colnames(detections.r) <- c('var','land','rep','timestep','code','detected','loc')
			}


			# incidence -- more rows than timesteps, separate output (optional...)
			if ('incidence' %in% out.opts){
				incidence.r = out.list$incidence
				n.inc = nrow(incidence.r)
				incidence.r = cbind(matrix(id.r, ncol=3, nrow=n.inc, byrow=TRUE), incidence.r)
				colnames(incidence.r) <- c('var','land','rep','timestep','state','loc')
			}

			# idzone -- not sure what it looks like yet
			# columms "cell" and "timestep"
# 			if ('idzone' %in% out.opts){
				# do idzone things, probably
# 			}

			# single value per vlr combination outputs
			# Tinc
			# sumTculled
			# Mspread
			# IConDD
			# ICatDD
			# TincToDD
			# TincFromDD
			# DET (total detections)
			#
			summ.vals.r = matrix(c(id.r, unlist(out.list[c(1,2,4:9,length(out.list))])), nrow=1)
			colnames(summ.vals.r) <- c('var','land','rep',names(out.list[c(1,2,4:9,length(out.list))]))

			# Put new results in output objects OR add new results to existing output objects
			if(r==1&l==1&v==1){
				tm.mat = tm.mat.r
# 				solocs=solocs.r
# 				Cto=Ct.r
# 				BB=BB.r
# 				idzone = idzone.r ## not done yet
				summ.vals = summ.vals.r
				if ("incidence" %in% out.opts) incidence = incidence.r
				if (end.tm > detectday & "alldetections" %in% out.opts) detections = detections.r

			} else {
				tm.mat = rbind(tm.mat, tm.mat.r)
# 				solocs = rbind(solocs, solocs.r)
# 				Cto = rbind(Cto, Ct.r)
# 				BB = rbind(BB, BB.r)
# 				idzone = rbind(idzone, idzone.r) ## not done yet
				summ.vals = rbind(summ.vals, summ.vals.r)
				if ("incidence" %in% out.opts) incidence = rbind(incidence, incidence.r)
				if (end.tm > detectday & "alldetections" %in% out.opts) detections = rbind(detections, detections.r)
			}
		}
	}
}

return(list('tm.mat' = tm.mat, 'summ.vals' = summ.vals, 'incidence' = incidence, 'detections' = detections))
# return(list("solocs"=solocs,"Cto"=Cto))
}










