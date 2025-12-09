
VisualOutputs <- function(out.list, variables, land_grid_list){
        # Generates output plots for examples
        library(data.table)

        # Grab output data
#         output.data <- tar_read(out.list)
        output.data <- out.list
#         var.list <- tar_read(variables)
        var.list <- variables
        setDT(var.list)
#         grid.data <- tar_read(land_grid_list)
        grid.data <- land_grid_list

        # put outputs in variables
        lapply(names(output.data), function(x) eval(parse(text=paste0(x, '= as.data.table(output.data["',x,'"])'))))
        lapply(names(output.data), function(x) paste0(x, '= as.data.table(output.data["',x,'"])'))
        lapply(output.data, function(x) print(names(x)))

        # make things into data.tables and name columns for sanity
        tm.mat= as.data.table(output.data["tm.mat"])
        setnames(tm.mat, unlist(lapply(strsplit(names(tm.mat), 'tm.mat.'), function(x) unlist(x)[2])))
        summ.vals= as.data.table(output.data["summ.vals"])
        setnames(summ.vals, unlist(lapply(strsplit(names(summ.vals), 'summ.vals.'), function(x) unlist(x)[2])))
        incidence= as.data.table(output.data["incidence"])
        setnames(incidence, unlist(lapply(strsplit(names(incidence), 'incidence.'), function(x) unlist(x)[2])))
        incidence[,max.time := max(timestep), by=.(var, rep, land)]
        detections= as.data.table(output.data["detections"])
        setnames(detections, unlist(lapply(strsplit(names(detections), 'detections.'), function(x) unlist(x)[2])))
        if(nrow(detections) > 0) {
                detections[,max.time := max(timestep), by=.(var, rep, land)]
                unq.det <- unique(detections[,.(var,land,rep,timestep,loc,max.time,code,detected)])
        }

        # generate an image for environment quality grid
        ## only works for single grid in land_grid_list
#         grid.key = as.data.table(grid.data[[1]][[2]])
#         setnames(grid.key, c('cell','tlX','tlY','trX','trY','ctX','ctY','rasval')[1:ncol(grid.key)]) ## see Make_Grid.R
#         grid.centers = grid.key[,.(cell, ctX, ctY, rasval)]


        ## get some plots of temporal data
        # time-based matrix outputs
        tm.pop.unq <- unique(tm.mat[,.(var, land)])
        png('./test_outputs/tm.plots.png', width=1600, height=1000)
        ## make plot dimensions dynamics based on parameter inputs
        pdef <- par(mfrow=c(2,3), oma=c(8,0.2,0.2,0.2),lwd=2.3, cex.lab=1.8, cex.axis=1.6, cex.main=2, cex.sub=1.4)
        # par(mfrow=c(1,1), oma=c(8,0.2,0.2,0.2),lwd=2.3, cex.lab=1.8, cex.axis=1.6, cex.main=2, cex.sub=1.4)
        # par(mfcol=c(length(unique(tm.mat[,var])), length(unique(tm.mat[,land]))))
        # set plot ranges based on maximum of everything that will be on the multiplot figure
        xrng = range(tm.mat[,timestep])
        yrng = log1p(range(tm.mat[,.(BB,S,E,I,R,C,Z)]))
        # loop over the parameter combinations for each plot panel
        seirczbb.temporal <- function(v, l, dat=tm.mat, vlist = var.list){
                vardat <- paste(vlist[v,], collapse=' ') # quick and dirty parameter inclusion
                # create a blank plot
                plot(0,0,xlim=xrng, ylim=yrng, xlab='time (weeks)', ylab='log(1+individuals)', col=NULL, main=paste('pop dynamics: vars',v,'| land',l), sub=vardat)
                # subset with things that have only the land tile and variable combination
                subdat <- dat[var==v & land==l,]
                # count up reps
                reps = unique(subdat[,rep])
                # draw a line of each type for each rep
                lapply(reps, function(r) lines(log1p(S) ~ timestep, data=subdat[rep==r,], col='olivedrab', lty=1))
                lapply(reps, function(r) lines(log1p(E) ~ timestep, data=subdat[rep==r,], col='orange', lty=1))
                lapply(reps, function(r) lines(log1p(I) ~ timestep, data=subdat[rep==r,], col='red', lty=1))
                lapply(reps, function(r) lines(log1p(R) ~ timestep, data=subdat[rep==r,], col='blue', lty=1))
                lapply(reps, function(r) lines(log1p(C) ~ timestep, data=subdat[rep==r,], col='purple', lty=1))
                lapply(reps, function(r) lines(log1p(Z) ~ timestep, data=subdat[rep==r,], col='black', lty=1))
                lapply(reps, function(r) lines(log1p(BB) ~ timestep, data=subdat[rep==r,], col='pink', lty=1))
                lapply(reps, function(r) abline(v=subdat[rep==r,max(timestep)], col='grey'))
        }
        mapply(seirczbb.temporal, v=tm.pop.unq[,var], l=tm.pop.unq[,land])
        # create a legend under the panels
        par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)
        plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
        legend("bottom", legend=c('S','E','I','R','C','Z','Births'),#'rep1','rep2'),
                lty=c(rep(1,8), 2),
                col=c('olivedrab','orange','red','blue','purple','black','pink'),#'grey','grey'),
                horiz=TRUE, bty='n', cex=2.3, lwd=3)
        dev.off()
        par(pdef)


        ## plots of the spatial data
        ### spatial habitat quality
        rasvalmat <- as.matrix(dcast(grid.centers, ctX ~ ctY, value.var = c('rasval')))
        plot(range(grid.centers[,ctX]), range(grid.centers[,ctY]), type='n')
        ## will need to set cex for points to have even squares based on output plot size
        points(ctY ~ ctX, data=grid.centers, col=terrain.colors(10, rev=TRUE)[grid.centers[,rasval*10]], pch=15)

        ### timing and extent of incidence/detections
        # connect incidence with cell locations
        unq.eic <- unique(incidence[,.(var,land,rep,timestep,loc,max.time)])
        mapply(function(v, l){
                input.dat <- unq.eic[var==v & land==l & max.time >= 10,]
                if(nrow(input.dat) != 0) {
#                 dev.new()
                png(paste0('./test_outputs/incidence_v', v, 'l', l, '.png'), width=1000, height=1000)
        #         if(nrow(input.dat) == 0) {browser(); input.dat <- unq.eic[var==v & land==l,]}
                        defpar <- par(mfrow=c(length(unique(input.dat[,rep])),2))

                        lapply(unique(input.dat[,rep]), function(x){
                                inc.cell <- input.dat[rep==x,][grid.centers, on=.(loc = cell)][!is.na(timestep),]
                                setorder(inc.cell, timestep)
                                inc.cell.early <- inc.cell[, .SD[1], by=.(var, land, rep, loc)]
                                plot(range(grid.centers[,ctX]), range(grid.centers[,ctY]), type='n', main=paste('incidence: vars',v,'| land',l))
                                points(ctY ~ ctX, data=inc.cell.early, col=rainbow(max(timestep))[timestep], pch=16,cex=0.4)
                                points(ctY ~ ctX, data=inc.cell.early[timestep==1,], pch='X', cex=2)
                                hist(inc.cell[,timestep], breaks=max(inc.cell[,timestep]), col=rainbow(max(inc.cell[,timestep])))
                        })
                        par(defpar)
                        dev.off()
                }
                },
        v=tm.pop.unq[,var],
        l=tm.pop.unq[,land]
        )

        # connect detections with cell locations
        unq.eic <- unique(detections[!is.na(loc) & detected != 0,.(var,land,rep,timestep,loc,max.time,code,detected)])
        mapply(function(v, l){
                input.dat <- unq.eic[var==v & land==l & max.time >= 10,]
                if(nrow(input.dat) != 0) {
#                 dev.new()
                png(paste0('./test_outputs/detection_v', v, 'l', l, '.png'), width=1000, height=1000)
        #         if(nrow(input.dat) == 0) {browser(); input.dat <- unq.eic[var==v & land==l,]}
                        defpar <- par(mfrow=c(length(unique(input.dat[,rep])),2))

                        lapply(unique(input.dat[,rep]), function(x){
                                inc.cell <- input.dat[rep==x,][grid.centers, on=.(loc = cell)][!is.na(timestep),]
                                setorder(inc.cell, timestep)
        #                         inc.cell.early <- inc.cell[, .SD[1], by=.(var, land, rep, loc)]
                                plot(range(grid.centers[,ctX]), range(grid.centers[,ctY]), type='n', main=paste('detections: vars',v,'| land',l))
                                points(ctY ~ ctX, data=inc.cell, col=rainbow(max(timestep))[timestep], pch=16,cex=0.4)
                                points(ctY ~ ctX, data=inc.cell[timestep==min(inc.cell[,timestep]),], pch='X', cex=2)
                                hist(inc.cell[,timestep], breaks=max(inc.cell[,timestep]), col=rainbow(max(inc.cell[,timestep])))
                        })
                        par(defpar)
                        dev.off()
                }
                },
        v=tm.pop.unq[,var],
        l=tm.pop.unq[,land]
        )

}
