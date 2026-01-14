## burn-in controller

# Essentially runs a burn-in in the targets framework for passing to RunSimulationFunction.R
# uses the same basic structure as RunSimulationFunction.R


RunBurnIn <- function(land_grid_list, parameters, variables, cpp_functions){

    ## Filters variables (parameters with >1 value) out of parameters
    ## selecting variables is done in SetVarParms.R
    parameters <- parameters[names(parameters) %in% names(variables) == FALSE]

    #Pull needed parms from parameters for all reps
    list2env(parameters, .GlobalEnv)
    #Need nested loops:
    #1, loop through all landscapes
    #2, loop through all parameter settings
    vl.list <- expand.grid(v = 1:nrow(variables), l = 1:length(land_grid_list))
    burn.out <- mapply(function(v, l){
#     for(v in 1:nrow(variables)){
        vars <- variables[v,]
        names(vars)[names(vars) == "density"] <- "dens"
        vars <- as.list(vars)
        list2env(vars, .GlobalEnv)

        #calc vals based on variables
        N0=dens*area
        K=N0*1.5

        centroids <- land_grid_list[[l]]$centroids
        grid <- land_grid_list[[l]]$grid

        pop <- InitializeSounders(centroids, grid, c(N0, ss), pop_init_grid_opts)
        outputs <- Initialize_Outputs(parameters)
        # Burn-in of pig population (similar to SimulateOneRun.R, but no infection)
        out.burn <- BurnIn(outputs, pop, centroids, grid, parameters, cpp_functions, K, v, l, r)
        rep.out <- rep_outputs(out.burn, v, l, 0, parameters, out.opts)
        rep.out[[6]] <- cbind(v, l, out.burn$pop)
        return(rep.out)
    },
    v = vl.list[,1],
    l = vl.list[,2]
    )

    tm.mat <- rbindlist(burn.out[1,])[!is.na(var),]
    summ.vals <- rbindlist(lapply(burn.out[2,], as.data.table))[!is.na(var),]
    incidence <- burn.out[3,][[1]]
    detections <- burn.out[4,][[1]]
    allzone <- burn.out[5,][[1]]
    pops <- rbindlist(lapply(burn.out[6,], function(x) as.data.table(x)))

    return(list('tm.mat' = tm.mat, 'summ.vals' = summ.vals, 'incidence' = incidence, 'detections' = detections, 'allzone' = allzone, 'pops' = pops))
}










