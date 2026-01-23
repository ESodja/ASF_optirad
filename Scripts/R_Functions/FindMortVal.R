## burn-in controller

# Essentially runs a burn-in to find mort_val that gives target density for passing to RunBurnIn.R and RunSimulationFunction.R
# uses the same basic structure as RunSimulationFunction.R


FindMortVal <- function(land_grid_list, parameters, variables, cpp_functions){

    ## Filters variables (parameters with >1 value) out of parameters
    ## selecting variables is done in SetVarParms.R
    ## this might be done in a subsection of the landscape to make it faster?
    parameters <- parameters[names(parameters) %in% names(variables) == FALSE]

    #Pull needed parms from parameters for all reps
    list2env(parameters, .GlobalEnv)
    #Need nested loops:
    #1, loop through all landscapes
    #2, loop through all parameter settings

    setDT(variables)
    vars.in <- copy(variables)

    # clean out variables that only affect detection and/or culling, because they affect nothing on this part of the model
    cull.params.vars <- c('Rad','Intensity','cullstyle','alphaC','sample','DetP','detectday')
    if(length(which(names(variables) %in% cull.params.vars)) > 0){
        cull.vars.i <- cull.params.vars[which(cull.params.vars %in% names(variables))]
        variables <- unique(variables[,-..cull.vars.i])
    }

    vl.list <- expand.grid(v = 1:nrow(variables), l = 1:length(land_grid_list))
    burn.out <- mapply(function(v, l){
#     for(v in 1:nrow(variables)){
        vars <- variables[v,]
        names(vars)[names(vars) == "density"] <- "dens"
        vars <- as.list(vars)
        list2env(vars, .GlobalEnv)

        #calc vals based on variables
        N0=dens*area ## actually need this for initial sounder distribution
        K=N0*1.5

        centroids <- land_grid_list[[l]]$centroids
        grid <- land_grid_list[[l]]$grid

        pop <- InitializeSounders(centroids, grid, c(N0, ss), pop_init_grid_opts)
        outputs <- Initialize_Outputs(parameters)

        # Run a bunch of less-targeted burn-ins to estimate what mortality parameter values give target density
        ## can use a more intelligent algorithm to find this faster
        mortparm <- lapply(seq(mort_val_test[1], mort_val_test[2], by=1), function(x) {
            # temporarily change mortality parameter to test value x
            parameters$mort_val <- x
            # run a mini burn-in to get a range of densities corresponding to the mortality parameter value
            out.vals <- BurnIn_paramfind(outputs, pop, centroids, grid, parameters, cpp_functions, K, v, l)
            return(out.vals)
        })
        # find which output is closest to the target density, and keep that mortality parameter in the parameters for v and l values
        ## probably should record that piece as well.
        # convert output to matrix format
        mortparm <- matrix(unlist(mortparm), ncol=3, byrow=TRUE)
        # linear approximate the target density from linear interpolation between the points (not sexy, but not bad if there are enough points)
        apx.target <- approx(mortparm[,3], mortparm[,1], xout=dens)
        # error check if target density is outside the range of tested values
        if(is.na(apx.target$y)) {
            print(mortparm)
            stop('Target density outside of range estimate; use the above table to adjust the range of mort_val_test range')
        }

#         parameters$mort_val <- apx.target$y
#
#         # Burn-in of pig population (similar to SimulateOneRun.R, but no infection)
#         out.burn <- BurnIn(outputs, pop, centroids, grid, parameters, cpp_functions, K, v, l)
# #         browser()
# #         rep.out <- rep_outputs(out.burn, v, l, 0, parameters, out.opts)
# #         rep.out[[6]] <- cbind(v, l, out.burn$pop)
# #         rep.out[[7]] <- parameters$mort_val
#         out.burn <- append(out.burn, list(append(vars, parameters$mort_val)))
        return(list(v, l, apx.target$y))
    },
    v = vl.list[,1],
    l = vl.list[,2]
    )

    mv_found <- as.data.table(matrix(unlist(burn.out), ncol=3, byrow=TRUE))
    setnames(mv_found, c('var','land','mort_val'))

    ## right??
    variables <- cbind(variables, mv_found)

#     tm.mat <- rbindlist(burn.out[1,])[!is.na(var),]
#     summ.vals <- rbindlist(lapply(burn.out[2,], as.data.table))[!is.na(var),]
#     incidence <- burn.out[3,][[1]]
#     detections <- burn.out[4,][[1]]
#     allzone <- burn.out[5,][[1]]

#     pops <- rbindlist(lapply(burn.out[1,], function(x) as.data.table(x)))
#     BB <- burn.out[2,]
#     Incidence <- burn.out[3,]
#     Tculled <- burn.out[4,]
#     ICtrue <- burn.out[5,]
#     out <- burn.out[6,]
#     detectday <- burn.out[7,]
#     Ct <- burn.out[8,]
#     out.opts <- burn.out[9,]
#     input.opts <- burn.out[10,] # out.opts, loc.list, idzone.mat, POSlive, POSdead, POSlive_locs, POSdead_locs, allzonecells, incidence
#     loc.list <- input.opts[[1]][[2]]
#     idzone.mat <- input.opts[[1]][[3]]
#     out.opts <- burn.out[11,]
#     mortparams <- burn.out[12,]
    variables <- vars.in[variables, on=.NATURAL]

    return(variables)
#     return(list('tm.mat' = tm.mat, 'summ.vals' = summ.vals, 'incidence' = incidence, 'detections' = detections, 'allzone' = allzone, 'pops' = pops, 'mortparms' = unlist(burn.out[7,])))
}










