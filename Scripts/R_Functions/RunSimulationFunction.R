
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
RunSimulationReplicates <- function(land_grid_list, parameters, variables, cpp_functions, reps, burn.list){

# burn.list: pop, BB, loc.list
#     # Pull values from variables for mort_val parameter
#     burn.vars <- as.data.table(matrix(as.numeric(unlist(burn.list[18,])), nrow=length(burn.list[18,]), byrow=TRUE))
#     burn.vars[,2] <- as.data.table(matrix(unlist(burn.list[18,]), nrow=length(burn.list[18,]), byrow=TRUE))[,2]
#     colnames(burn.vars) <- c(names(burn.list[18,][[1]]))
#     setnames(burn.vars, ncol(burn.vars), 'mort_val')
# burn.pops <- rbindlist(lapply(seq(burn.list), function(x) as.data.table(burn.list[[x]][[1]])))
# BB <- rbindlist(lapply(seq(burn.list), function(x) as.data.table(burn.list[[x]][[2]])))
# loc.list <- rbindlist(lapply(seq(burn.list), function(x) as.data.table(lapply(seq(length(burn.list[[x]][[3]])), function(y) {print(as.data.table(burn.list[[x]][[3]][[y]]))}))))
#
# compress.locs <- function(loc.lists, v, l){
#     lapply(loc.lists, function(x) x[[1]])
# }
# burn.pops <- rbindlist(lapply(seq(burn.list), function(x) as.data.table(burn.list[[x]][[1]])))

    names(variables)[names(variables) == "density"] <- "dens"

#     variables <- as.data.table(variables)[burn.vars, on=.NATURAL]

    # Filter variables (parameters with >1 value) out of parameters
    ## selecting variables is done in SetVarParms.R
    parameters <- parameters[names(parameters) %in% names(variables) == FALSE]

    #Pull needed parms from parameters for all reps
    list2env(parameters, .GlobalEnv)
    #Need nested loops:
    #1, loop through all landscapes
    #2, loop through all parameter settings

    # looping table for mapply
    lvtable <- expand.grid(vars = seq(nrow(variables)), land = seq(length(land_grid_list)), rep=seq(reps))
    lvtable$iterate <- rep(seq(length(burn.list)))

    rep.list <- mapply(function(v.val, l.val, r.val, i.val){

# browser()
        vars <- variables[v.val,]
        vars <- as.list(vars)
        list2env(vars, .GlobalEnv)

        #calc vals based on variables
        N0=dens*area
        K=N0*1.5

        burn.input <- burn.list[[i.val]]

#         burn.times <- unlist(burn.list[17,])
        #loop through landscapes
        centroids <- land_grid_list[[l.val]]$centroids
        grid <- land_grid_list[[l.val]]$grid

#             pop <- InitializeSounders(centroids, grid, c(N0, ss), pop_init_grid_opts)
        print(i.val)
        pop <- as.matrix(burn.input[[1]][,-c(1,2)])
        pop <- InitializeInfection(pop, centroids, grid, parameters)
        outputs <- Initialize_Outputs(parameters)
        outputs$BB <- burn.input[[2]]
#         outputs$Incidence <- burn.input[[i.val]][[3]]
#         outputs$Tculled <- burn.input[[i.val]][[4]]
#         outputs$ICtrue <- burn.input[[i.val]][[5]]
#         outputs$out <- burn.input[[i.val]][[6]]
#         outputs$detectday <- burn.input[[i.val]][[7]]
#         outputs$Ct <- burn.input[[i.val]][[8]]
        outputs$loc.list <- burn.input[[3]]
#         outputs$POSlive <- burn.input[[i.val]][[11]]
#         outputs$POSdead <- burn.input[[i.val]][[12]]
#         outputs$POSlive_locs <- burn.input[[i.val]][[13]]
#         outputs$POSdead_locs <- burn.input[[i.val]][[14]]
#         outputs$allzone <- burn.input[[i.val]][[15]]
#         outputs$incidence <- cbind(matrix(NA, nrow=0, ncol=3), burn.input[[i.val]][[16]])
        burn.time.vl <- burn.input[[4]]

        # each rep starts in the same post burn-in condition

#     browser()
        #Do simulations
        out.list <- SimulateOneRun(outputs, pop, centroids, grid, parameters, cpp_functions, K, v.val, l.val, r.val, burn.time.vl)
        #Handle outputs
        rep.out <- rep_outputs(out.list, v.val, l.val, r.val, parameters, out.opts)
        return(rep.out)
    },
    v.val=lvtable[,1], l.val=lvtable[,2], r.val = lvtable[,3], i.val=lvtable[,4])

    tm.mat <- rbindlist(rep.list[1,])
    summ.vals <- rbindlist(lapply(rep.list[2,], as.data.table))
    incidence <- rbindlist(rep.list[3,which(lapply(rep.list[3,], ncol) > 1)])
    detections <- rbindlist(lapply(rep.list[4,][!is.na(rep.list[4,])], as.data.table))
    allzone <- rbindlist(lapply(rep.list[5,][!is.na(rep.list[5,])], as.data.table))
    solocs.all <- rbindlist(rep.list[6,])

    return(list('tm.mat' = tm.mat, 'summ.vals' = summ.vals, 'incidence' = incidence, 'detections' = detections, 'allzone' = allzone, 'solocs.all' = solocs.all))
}










