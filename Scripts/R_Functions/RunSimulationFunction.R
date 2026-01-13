
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

    ## Filters variables (parameters with >1 value) out of parameters
    ## selecting variables is done in SetVarParms.R
    parameters <- parameters[names(parameters) %in% names(variables) == FALSE]

    #Pull needed parms from parameters for all reps
    list2env(parameters, .GlobalEnv)
    #Need nested loops:
    #1, loop through all landscapes
    #2, loop through all parameter settings
    for(v in 1:nrow(variables)){
        vars <- variables[v,]
        names(vars)[names(vars) == "density"] <- "dens"
        vars <- as.list(vars)
        list2env(vars, .GlobalEnv)

        #calc vals based on variables
        ## could set these as variables/parameter options
        N0=dens*area
        K=N0*1.5

        #loop through landscapes
        for(l in 1:length(land_grid_list)){
            centroids <- land_grid_list[[l]]$centroids
            grid <- land_grid_list[[l]]$grid

            pop <- InitializeSounders(centroids, grid, c(N0, ss), pop_init_grid_opts)
            outputs <- Initialize_Outputs(parameters)
            # Burn-in of pig population (similar to SimulateOneRun.R, but no infection)
            out.burn <- BurnIn(outputs, pop, centroids, grid, parameters, cpp_functions, K, v, l, r)
            if (v==1 & l==1){
                # the first burn-in
                rep.out <- rep_outputs(out.burn, v, l, 0, parameters, out.opts)
            } else {
                # subsequent burn-ins
                rep.out <- rep_outputs(out.burn, v, l, 0, parameters, out.opts, rep.out)
            }
            pop <- InitializeInfection(out.burn$pop, centroids, grid, parameters)

            for(r in 1:reps){
                # each rep starts in the same post burn-in condition

                #Do simulations
                out.list <- SimulateOneRun(outputs, pop, centroids, grid, parameters, cpp_functions, K, v, l, r)
                #Handle outputs
                print('repoutputs')
                rep.out <- rep_outputs(out.list, v, l, r, parameters, out.opts, rep.out)
            } # end rep loop
        } # end landscape loop
    } # end variable loop

    return(list('tm.mat' = rep.out[[1]], 'summ.vals' = rep.out[[2]], 'incidence' = rep.out[[3]], 'detections' = rep.out[[4]], 'allzone' = rep.out[[5]]))
}










