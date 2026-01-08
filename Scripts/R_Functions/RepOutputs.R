## repackage outputs after simulation/burn-in run

rep_outputs <- function(out.list, v, l, r, parameters, out.opts, prevrep.in = as.list(rep(NA, 5))){
    list2env(parameters, .GlobalEnv)

    # if it is the burnin, the r=0 below will overwrite NA defaults in this irrelevant
    # if it is after the burn-in, these will have existing values that new values will be added to
    tm.mat <- prevrep.in[[1]]
    summ.vals <- prevrep.in[[2]]
    incidence <- prevrep.in[[3]]
    detections <- prevrep.in[[4]]
    allzone <- prevrep.in[[5]]

    ## test these outputs if out.opts doesn't include them
    ## id by variable combination, landscape, rep, and timestep for one row per timestep data
    end.tm <- out.list$endtime
    id.r <- c(v,l,r)
    tm.mat.r <- cbind(matrix(id.r, nrow=end.tm, ncol=3, byrow=TRUE), seq(end.tm))
    colnames(tm.mat.r) <- c("var","land","rep","timestep")

    #Handle effective removal rate (timestep output)
    tm.mat.r <- cbind(tm.mat.r, out.list$Ct[seq(end.tm)])
    colnames(tm.mat.r)[ncol(tm.mat.r)]="Ct"

    # Births
#             BB.r = out.list$BB
    tm.mat.r <- cbind(tm.mat.r, out.list$BB[seq(end.tm)])
    colnames(tm.mat.r)[ncol(tm.mat.r)] <- 'BB'

    #Handle sounderlocs (optional...)
    ## seems like other things were supposed to happen in sounderlocsSummarize, if we want to use those this will have to change
    if ("sounderlocs" %in% out.opts){
        solocs.r <- sounderlocsSummarize(out.list$sounderlocs,r)[[1]]
        tm.mat.r <- cbind(tm.mat.r[complete.cases(tm.mat.r),], solocs.r[3:8])
    }

    # detections (optional...)
    ## has a row for each timestep AND detection type, with timestep, code (1=live,0=dead), number of individuals detected, and position
    ## still need to test sample = 1
    if (end.tm > detectday & 'alldetections' %in% out.opts){
        detections.r <- out.list$alldetections
        detections.r <- detections.r[detections.r[,1] <= end.tm & detections.r[,1] >= detectday,]
        n.det <- nrow(detections.r)
        detections.r <- cbind(matrix(id.r, ncol=3, nrow=n.det, byrow=TRUE), detections.r)
        colnames(detections.r) <- c('var','land','rep','timestep','code','detected','loc')
        allzone.r <- out.list$allzonecells
        colnames(allzone.r) <- c('var','land','rep','timestep','loc')
    }

    # incidence -- more rows than timesteps, separate output (optional...)
    if ('incidence' %in% out.opts){
        incidence.r <- out.list$incidence
        n.inc = nrow(incidence.r)
        incidence.r <- cbind(matrix(id.r, ncol=3, nrow=n.inc, byrow=TRUE), incidence.r)
        colnames(incidence.r) <- c('var','land','rep','timestep','state','loc')
    }

    # single value per vlr combination outputs
    # Tinc # sumTculled # Mspread # IConDD # ICatDD # TincToDD # TincFromDD # DET (total detections) #
    summ.vals.r <- matrix(c(id.r, unlist(out.list[c(1,2,4:9,length(out.list))])), nrow=1)

    colnames(summ.vals.r) <- c('var','land','rep',names(out.list[c(1,2,4:9,length(out.list))]))

    # Put new results in output objects OR add new results to existing output objects
    if(r==0 & l==1 & v==1){
        tm.mat <- tm.mat.r
        summ.vals <- summ.vals.r
        if ("incidence" %in% out.opts) incidence <- incidence.r
        if (end.tm > detectday & "alldetections" %in% out.opts) {
            detections <- detections.r
            allzone <- allzone.r
            colnames(allzone) <- c('var','land','rep','timestep','loc')
        }

    } else {
        tm.mat <- rbind(tm.mat, tm.mat.r)
        summ.vals <- rbind(summ.vals, summ.vals.r)
        if ("incidence" %in% out.opts) incidence <- rbind(incidence, incidence.r)
        if (end.tm > detectday & "alldetections" %in% out.opts) {
            detections <- rbind(detections, detections.r)
            if (is.null(nrow(allzone))) {
                # if it hasn't been established yet
                allzone <- allzone.r
            } else {
                allzone <- rbind(allzone, allzone.r)
            }
            colnames(allzone) <- c('var','land','rep','timestep','loc')
        }
    }

    return(list(tm.mat, summ.vals, incidence, detections, allzone))
}
