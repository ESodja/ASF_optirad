## Runs a burn-in for a given parameter set

BurnIn_paramfind <- function(outputs, pop, centroids, grid, parameters, cpp_functions, K, v, l){
    require(dplyr)
    for(i in 1:length(cpp_functions)){
        print(paste0("sourcing ",cpp_functions[[i]]))
        Rcpp::sourceCpp(cpp_functions[[i]])
    }

######## Release parameters to function environment ########
    list2env(parameters, .GlobalEnv)

######## Initialize Output Objects ########
    list2env(outputs, .GlobalEnv)

######## Start simulation ########

    i <- 0
    ## need to base condition on population stability, e.g. when density for the last several timesteps varies by less than some given amount
    ## which means we need to measure density throughout...
    pop.list <- sum(pop[,1])
    pop.var <- 99
    dens.list <- c(sum(pop[,1])/(nrow(centroids) * inc^2))
    dens.var <- 99
    pop.recent <- c(1,1)

#     while(i < burn_weeks){
    while((i < 25 | dens.var > 0.05 ) & i < burn_weeks){
#     while((i < 25 | dens.var > 0.05 | (diff(range(pop.recent))/max(pop.recent)) > 0.1 ) & i < burn_weeks){
#     while((i < 25 | dens.var > 0.05 | pop.var > 1000) & i < burn_weeks){
        i <- i+1
        print(paste0("timestep: ",i))
        print(colSums(pop[,8:13]))

#         if("sounderlocs" %in% out.opts){
#             loc.list[[i]] <- pop[, c(3,8:13)]
#         }

######## Track I/C locations ########
#
#         Isums[i] <- 0
#         Csums[i] <- 0
#
#         I_locs[[i]] <- pop[pop[, 10] > 0, 3]
#         C_locs[[i]] <- pop[pop[, 12] > 0, 3]
#         I_locs[[i]] <- 0
#         C_locs[[i]] <- 0

######## Movement ########

        pop <- FastMovement(pop, centroids, shape, rate, inc, mv_pref)

        ## Most of these next few things are included only for continuity with post-infection simualtions
######## State Changes ########
#births, natural deaths, disease state changes (exposure, infection, recovery, death), carcass decay

        st.list <- StateChanges(pop, centroids, nrow(centroids), parameters, Incidence, BB, i)
#
#         Eep_mat[i,] <- st.list$Eep
#         Sdpb_mat[i,] <- st.list$Sdpb
#         Sdpd_mat[i,] <- st.list$Sdpd
#         Rep_mat[i,] <- st.list$Rep
#         Cep_mat[i,] <- st.list$Cep
#         Rdpd_mat[i,] <- st.list$Rdpd
#         Ccd_mat[i,] <- st.list$Ccd
#         Zcd_mat[i,] <- st.list$Zcd
#         Iep_mat[i,] <- st.list$Iep

        pop <- st.list[[1]]
#         Incidence <- st.list[[2]]
        BB <- st.list[[3]]

#         if("incidence" %in% out.opts){
#             inc.mat <- matrix(nrow=0, ncol=3)
#         }


######## Initiate Culling Zone ########

#         Ct[i, 1] <- 0

####Track true spatial spread

#         out[i,] <- c(0,0,0)

####Summarize infections
#         ICtrue[i] <- 0

#### Update population matrix
#Remove rows in pop with 0 pigs
        pigcols <- c(1, 8:13)
        pop <- pop[which(rowSums(pop[, pigcols, drop=FALSE]) != 0),, drop=FALSE]

        pop.list <- c(pop.list, sum(pop[,1]))
        pop.recent <- pop.list[max(1, (length(pop.list)-10)):length(dens.list)]
        pop.var <- var(pop.recent)
        dens.list <- c(dens.list, sum(pop[,1])/(nrow(centroids) * inc^2))
        dens.recent <- dens.list[max(1,(length(dens.list)-10)):length(dens.list)]
        dens.var <- var(dens.recent)
        print(paste('pop variance', pop.var))
        print(paste('dens variance', dens.var))
    } # end while loop of timesteps
    dens.mean <- mean(dens.recent)

#     print('finish')
#     browser()

#     input.opts <- vector(mode="list", length=1)
#     input.opts[[1]] <- out.opts
#     names(input.opts)[1] <- "out.opts"
#     if("sounderlocs" %in% out.opts){
#         templist <- vector(mode="list",length=1)
#         templist[[1]] <- loc.list
#         input.opts <- append(input.opts,templist)
#         names(input.opts)[length(input.opts)] <- "loc.list"
#     }
#
#     ## many of these could be cleaned out as well
#     if("idzone"%in%out.opts){
#         templist <- vector(mode="list",length=1)
#         templist[[1]] <- matrix(rep(NA, 3), ncol=3)
#         input.opts <- append(input.opts,templist)
#         names(input.opts)[length(input.opts)] <- "idzone.mat"
#     }
#
#     if("alldetections"%in%out.opts){
#         templist <- list(POSlive)    # directly create a list with POSlive
#         input.opts <- append(input.opts, templist)
#         names(input.opts)[length(input.opts)] <- "POSlive"
#
#         templist <- list(POSdead)
#         input.opts <- append(input.opts, templist)
#         names(input.opts)[length(input.opts)] <- "POSdead"
#
#         templist <- list(POSlive_locs)
#         input.opts <- append(input.opts, templist)
#         names(input.opts)[length(input.opts)] <- "POSlive_locs"
#
#         templist <- list(POSdead_locs)
#         input.opts <- append(input.opts, templist)
#         names(input.opts)[length(input.opts)] <- "POSdead_locs"
#
#         templist <- list(allzone)
#         input.opts <- append(input.opts, templist)
#         names(input.opts)[length(input.opts)] <- "allzonecells"
#     }
#
#     if("incidence" %in% out.opts){
#         templist <- vector(mode="list",length=1)
#         templist[[1]] <- inc.mat
#         input.opts <- append(input.opts,templist)
#         names(input.opts)[length(input.opts)] <- "incidence"
#     }
#     }

#     list.all <- GetOutputs(pop, centroids, BB, Incidence, Tculled, ICtrue, out, detectday, Ct, out.opts, input.opts, is.burn=TRUE)

## want the end time as output to trim matrices
#     list.all <- append(list.all, i)
#     names(list.all)[length(list.all)] <- 'endtime'

    return(list(mort_val, dens.var, dens.mean))
    # mort_val is defined in input parameters as the value to be tested
    # dens.var gives an idea of how stable the population is at that level
    # dens.mean gives the output density for the tested mort_val
} #function closing bracket
