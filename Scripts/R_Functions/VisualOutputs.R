
VisualOutputs <- function(out.list, variables, land_grid_list, parameters){

    land_grid_list <- land_grid_list[1]
    # Generates output plots for examples
    library(data.table)

    setDT(variables)

    # put outputs in variables
    lapply(names(out.list), function(x) eval(parse(text=paste0(x, '= as.data.table(out.list["',x,'"])'))))
    lapply(names(out.list), function(x) paste0(x, '= as.data.table(out.list["',x,'"])'))
    lapply(out.list, function(x) print(names(x)))

    # make things into data.tables and name columns for sanity
    tm.mat= as.data.table(out.list["tm.mat"])
    setnames(tm.mat, unlist(lapply(strsplit(names(tm.mat), 'tm.mat.'), function(x) unlist(x)[2])))
    tm.mat <- tm.mat[!is.na(timestep),]
    summ.vals= as.data.table(out.list["summ.vals"])
    setnames(summ.vals, unlist(lapply(strsplit(names(summ.vals), 'summ.vals.'), function(x) unlist(x)[2])))
    summ.vals <- summ.vals[!is.na(endtime)]
    incidence= as.data.table(out.list["incidence"])
    setnames(incidence, unlist(lapply(strsplit(names(incidence), 'incidence.'), function(x) unlist(x)[2])))
    incidence[,max.time := max(timestep), by=.(var, rep, land)]
    detections= as.data.table(out.list["detections"])
    setnames(detections, unlist(lapply(strsplit(names(detections), 'detections.'), function(x) unlist(x)[2])))
    if(nrow(detections) > 0) {
        detections[,max.time := max(timestep), by=.(var, rep, land)]
        unq.det <- unique(detections[,.(var,land,rep,timestep,loc,max.time,code,detected)])
    }
    allzones <- as.data.table(out.list["allzone"])
    setnames(allzones, unlist(lapply(strsplit(names(allzones), 'allzone.'), function(x) unlist(x)[2])))
    solocs.all <- as.data.table(out.list['solocs.all'])
    setnames(solocs.all, unlist(lapply(strsplit(names(solocs.all), 'solocs.all.'), function(x) unlist(x)[2])))
    setnames(solocs.all, c('x','y'), c('ctX','ctY'))
    solocs.all[,nlive := S+E+I+R]

    browser()
    tm.mat <- tm.mat[land==1,]
    summ.vals <- summ.vals[land==1,]
    incidence <- incidence[land==1,]
    detections <- detections[land==1,]
    allzones <- allzones[land==1,]
    solocs.all <- solocs.all[l==1,]


    # generate an image for environment quality grid
    ## only works for single grid in land_grid_list
    grid.key = as.data.table(land_grid_list[[1]][[2]])
    setnames(grid.key, c('cell','tlX','tlY','trX','trY','ctX','ctY','rasval')[1:ncol(grid.key)]) ## see Make_Grid.R
    grid.centers = grid.key[,.(cell, ctX, ctY, rasval)]

    ## get some plots of temporal data
    # time-based matrix outputs
    tm.pop.unq <- unique(tm.mat[,.(var, land)])
    rows.plt <- round(sqrt(length(unique(tm.pop.unq[,var]))))
    cols.plt <- ceiling(sqrt(length(unique(tm.pop.unq[,var]))))
    png('./test_outputs/tm.plots.png', width=450*cols.plt, height=500*rows.plt)
    ## make plot dimensions dynamics based on parameter inputs
    pdef <- par(mfrow=c(rows.plt, cols.plt),
            oma=c(8,1.2,0.2,0.2),
            lwd=2.3, cex.lab=1.8, cex.axis=1.6, cex.main=2, cex.sub=1.4)
    # par(mfrow=c(1,1), oma=c(8,0.2,0.2,0.2),lwd=2.3, cex.lab=1.8, cex.axis=1.6, cex.main=2, cex.sub=1.4)
    # par(mfcol=c(length(unique(tm.mat[,var])), length(unique(tm.mat[,land]))))
    # set plot ranges based on maximum of everything that will be on the multiplot figure
    xrng = range(tm.mat[,timestep])
    yrng = log1p(range(tm.mat[,.(BB,S,E,I,R,C,Z)]))
    # loop over the parameter combinations for each plot panel
    seirczbb.temporal <- function(v, l, plt.i, rows.plt, cols.plt, dat=tm.mat, vlist = variables){
        vardat <- paste(vlist[v,], collapse=' ') # quick and dirty parameter inclusion
        # xlabel default
        xlabi <- ''
        # ylabel default
        ylabi <- ''
        # default plot panel margin sizes
        bmar <- lmar <- 2
        if (plt.i / rows.plt >= rows.plt){ # first row of panels
            xlabi <- 'time (weeks)'
            bmar <- 4
        }
        if (plt.i %% cols.plt == 1){ # first column of panels
            ylabi <- 'log(1+individuals)'
            lmar <- 4
        }
        pdef2 <- par(mar = c(bmar, lmar, 1, 1))
        # create a blank plot
#         plot(0,0,xlim=xrng, ylim=yrng, xlab=xlabi, ylab=ylabi, col=NULL, main=paste('pop dynamics: vars',v,'| land',l), sub=vardat)
        plot(0,0,xlim=xrng, ylim=yrng, col=NULL, ann=FALSE)
        mtext(xlabi, 1, line=2.5, cex=1.8)
        mtext(ylabi, 2, line=2.5, cex=1.8)
        mtext(paste('vars',v,'| land',l), 3, line=-1.5, cex=1.7, font=2)
        mtext(vardat, 3, line=-2.7)
        # subset with things that have only the land tile and variable combination
        subdat <- dat[var==v & land==l,]
        # count up reps
        reps = unique(subdat[,rep])
        # draw a line of each type for each rep
        lapply(reps, function(r) lines(log1p(S) ~ timestep, data=subdat[rep==r | rep == 0,], col='olivedrab', lty=1))
        lapply(reps, function(r) lines(log1p(E) ~ timestep, data=subdat[rep==r | rep == 0,], col='orange', lty=1))
        lapply(reps, function(r) lines(log1p(I) ~ timestep, data=subdat[rep==r | rep == 0,], col='red', lty=1))
        lapply(reps, function(r) lines(log1p(R) ~ timestep, data=subdat[rep==r | rep == 0,], col='blue', lty=1))
        lapply(reps, function(r) lines(log1p(C) ~ timestep, data=subdat[rep==r | rep == 0,], col='purple', lty=1))
        lapply(reps, function(r) lines(log1p(Z) ~ timestep, data=subdat[rep==r | rep == 0,], col='black', lty=1))
        lapply(reps, function(r) lines(log1p(BB) ~ timestep, data=subdat[rep==r | rep == 0,], col='pink', lty=1))
        lapply(reps, function(r) abline(v=subdat[rep==r,max(timestep)], col='grey'))
        lapply(reps, function(r) abline(v=subdat[rep==0,max(timestep)], col='black', lty=3))
        par(pdef2)
    }
    mapply(seirczbb.temporal, v=tm.pop.unq[,var], l=tm.pop.unq[,land], plt.i=seq(nrow(tm.pop.unq)), MoreArgs = list(rows.plt=rows.plt, cols.plt=cols.plt))
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
    png('./test_outputs/landscape.png')
    rasvalmat <- as.matrix(dcast(grid.centers, ctX ~ ctY, value.var = c('rasval')))
    plot(range(grid.centers[,ctX]), range(grid.centers[,ctY]), type='n')
    ## will need to set cex for points to have even squares based on output plot size
    points(ctY ~ ctX, data=grid.centers, col=terrain.colors(10, rev=TRUE)[grid.centers[,rasval*10]], pch=15)
    dev.off()

    ### timing and extent of incidence/detections
    # connect incidence with cell locations
    unq.eic <- unique(incidence[,.(var,land,rep,timestep,loc,max.time)])
#     mapply(function(v, l){
#         input.dat <- unq.eic[var==v & land==l & max.time >= 10,]
#         if(nrow(input.dat) != 0) {
# #         dev.new()
# #     browser()
#         usable.reps <- unique(input.dat[,rep])
#         png(paste0('./test_outputs/incidence_v', v, 'l', l, '.png'), width=1000, height=500*length(usable.reps))
#     #     if(nrow(input.dat) == 0) {browser(); input.dat <- unq.eic[var==v & land==l,]}
#             defpar <- par(mfrow=c(length(unique(input.dat[,rep])),2), oma=c(8,1.2,0.2,0.2))
#
#             lapply(usable.reps, function(x){
#                 inc.cell <- input.dat[rep==x,][grid.centers, on=.(loc = cell)][!is.na(timestep),]
#                 setorder(inc.cell, timestep)
#                 inc.cell.early <- inc.cell[, .SD[1], by=.(var, land, rep, loc)]
#                 plot(range(grid.centers[,ctX]), range(grid.centers[,ctY]), type='n', main=paste('incidence: vars',v,'| land',l), cex.lab = 1.4, xlab='x coordinate (Km)', ylab='y coordinate (Km)')
#                 points(ctY ~ ctX, data=inc.cell.early, col=rainbow(max(timestep))[timestep], pch=16,cex=0.4)
#                 points(ctY ~ ctX, data=inc.cell.early[timestep==1,], pch='X', cex=2)
# #                 browser()
#                 hist(inc.cell[,timestep], breaks=max(inc.cell[,timestep]), col=rainbow(max(inc.cell[,timestep])), main='Timing of Incidence', ylab='# New Incidence', xlab='Time (weeks)', cex.lab = 1.4)
#                 # Add line to indicate when first detection began
#                 abline(v=parameters['detectday'], lwd=3, lty=2, col='black')
#             })
#             # add legend
#             # create a legend under the panels
#             par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)
#             plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
#             legend("bottom", legend=c('Intro point','First detect'),#'rep1','rep2'),
#                 lty=c(NA, 2),
#                 col=c('black','black'),#'grey','grey'),
#                 pch=c('X',NA),
#                 horiz=TRUE, bty='n', cex=2.3, lwd=3)
#             par(defpar)
#             dev.off()
#         }
#         },
#     v=tm.pop.unq[,var],
#     l=tm.pop.unq[,land]
#     )

    # connect detections with cell locations
#     unq.eic <- unique(detections[!is.na(loc) & detected != 0 & max.time >= 10,.(var,land,rep,timestep,loc,max.time,code,detected)])
#     mapply(function(v, l){
#         input.dat <- unq.eic[var==v & land==l & max.time >= 10,]
#         if(nrow(input.dat) != 0) {
#             input.dat <- aggregate(input.dat[,detected], list(rep=input.dat[,rep], timestep=input.dat[,timestep], loc=input.dat[,loc]), sum)
#             setDT(input.dat)
#             setnames(input.dat, 'x', 'detected')
#             setorderv(input.dat, cols = c('rep','timestep','loc'))
#         #         dev.new()
#             usable.reps <- unique(input.dat[,rep])
#             png(paste0('./test_outputs/detection_v', v, 'l', l, '.png'), width=1000, height=500*length(usable.reps))
#             defpar <- par(mfrow=c(length(unique(input.dat[,rep])),2))
#
#             lapply(unique(input.dat[,rep]), function(x){
#                 inc.cell <- input.dat[rep==x,][grid.centers, on=.(loc = cell)][!is.na(timestep),]
#                 setorder(inc.cell, timestep)
#     #             inc.cell.early <- inc.cell[, .SD[1], by=.(var, land, rep, loc)]
# #                 browser()
#                 plot(range(grid.centers[,ctX]), range(grid.centers[,ctY]), type='n', main=paste('detections: vars',v,'| land',l), cex.lab = 1.4, xlab='x coordinate (Km)', ylab='y coordinate (Km)')
#                 points(ctY ~ ctX, data=inc.cell, col=rainbow(max(timestep))[timestep], pch=16,cex=1.4)
#                 points(ctY ~ ctX, data=inc.cell[timestep==min(inc.cell[,timestep]),], pch='X', cex=2)
#                 # setup for a stacked barplot
#                 bp <- dcast(unq.eic[var==v & land==l & max.time >= 10 & rep==x,][grid.centers, on=.(loc = cell), nomatch=NULL], timestep ~ code, value.var='detected', fun.aggregate=sum)[!is.na(timestep),]
# #                 setnames(bp, c('0','1'), c('C','IE'))
# #                 bp[,tot := C + IE]
#                 col.raw = rainbow(max(bp[,timestep]))[bp[,timestep]]
#                 if (nrow(as.matrix(bp)) > 1){
#                     barplot(t(as.matrix(bp))[c(2,3),], space=0, xlab='timestep', col='white', ylab='detected (dark=dead)', main = x, names.arg=t(as.matrix(bp))[1,])
#                     # dark portion of bars indicates how many detected were found dead vs. culled
#                     # makes two-tone stacked barplots (https://stackoverflow.com/a/59411350)
#                     # uses darken() from 'colorspace' package; probably is a non-package way to shift color hex codes
#                     for (i in 1:ncol(t(as.matrix(bp)))){
#                         xx = t(as.matrix(bp))[2:3,]
#                         xx[,-i] = NA
#                         colnames(xx)[-i] = NA
#                         barplot(xx,col=c(darken(col.raw[i],0.4),col.raw[i]), add=T, axes=F, space=0)
#                     }
#                 }
#
# #                 plot(detected ~ timestep,data=inc.cell, main='Timing of Detections', ylab='# New Detections', xlab='Time (weeks)', cex.lab = 1.4, type='b', pch=16)
#             })
#             par(defpar)
#             dev.off()
#         }
#     },
#     v=tm.pop.unq[,var],
#     l=tm.pop.unq[,land]
#     )

    radius = parameters['Rad']
#     unq.combos <- unique(unq.eic[max.time > 10,.(var, land, rep)])
    if (radius != 0){
        ## plotting incidence with zone of control over time for a bunch of plots
        # detection locations
        unq.detections <- unique(detections[,.(var,land,rep,timestep,loc,max.time,detected)])

        unq.incidence <- unique(incidence[,.(var,land,rep,timestep,loc,max.time)])
        unq.incidence[,is.inf := 1]
        unq.incidence[,loc.min := min(timestep), by=.(var, land, rep, loc)]
        unq.incidence[,loc.max := max(timestep), by=.(var, land, rep, loc)]
        unq.join <- CJ(timestep=seq(parameters[['thyme']]),
                 var=unique(unq.incidence[,var]),
                 land=unique(unq.incidence[,land]),
                 rep=unique(unq.incidence[,rep]),
                 loc=unique(unq.incidence[,loc]))
        unq.join <- unq.join[!unq.incidence, on=.(var, land, rep, timestep, loc)]
        unq.join[, is.inf := 0]
        unq.join <- unq.join[unique(unq.incidence[,.(var, land, rep, loc, loc.min, loc.max, max.time)]), on=.(var, land, rep, loc)]
        unq.join <- unq.join[timestep >= loc.min,]
#         unq.join <- unq.join[timestep <= loc.max & timestep >= loc.min,]
        unq.incidence <- rbind(unq.incidence, unq.join)

        unq.incidence[,loc.min := NULL]
        unq.incidence[,loc.max := NULL]
        setorder(unq.incidence, var, land, rep, timestep, loc)
        unq.combos <- unique(unq.incidence[,sum(is.inf), by=.(var, land, rep, timestep)][V1>1,.(var,land,rep)])

        browser()
#         unq.combos <- unq.combos[var == 7 & rep == 2,]
        mapply(function(i, j, k){
            sub.incidence <- unique(unq.incidence[var == i & land == j & rep == k & loc != 0,.(timestep, loc, is.inf)])[grid.centers, on=.(loc = cell), nomatch=NULL] # gets rid of code column
            sub.detections <- unique(unq.detections[var == i & land == j & rep == k & detected != 0,.(timestep, loc)])[grid.centers, on=.(loc = cell), nomatch=NULL] # gets rid of code column
            sub.zones <- allzones[var == i & land == j & rep == k & loc != 0,][grid.centers, on=.(loc=cell), nomatch=NULL]
            sub.solocs <- solocs.all[v == i & l == j & r == k,]
            paneldim <- ceiling(sqrt(max(sub.incidence[,timestep])))
            png(paste0('./test_outputs/sptmplot_', i, j, k,'.png'), width=2000, height=2000)
            par(mfrow = c(paneldim, paneldim), oma=c(0,0,0,0), mar=c(0,0,0,0))
            lapply(seq(max(sub.incidence[,timestep])), function(x){
                plot(ctY ~ ctX, data=sub.solocs[timestep == x,], pch='.', col='gray', xlim=c(0,100), ylim=c(0,100))
                points(ctY ~ ctX, data=sub.zones[timestep <= x,], pch='.', col='yellow')#, cex=1.2, xlim=c(0,100), ylim=c(0,100))
#                 points(y ~ x, data=sub.solocs[timestep == x,], pch=1, col='gray')
#                 points(ctY ~ ctX, data=sub.incidence[timestep <= x & is.inf == 0,], pch='.', col='red')
                points(ctY ~ ctX, data=sub.incidence[timestep == x & is.inf == 1,], pch=3, col='red')
                if (x > detectday) points(ctY ~ ctX, data=sub.detections[timestep <= x,], pch=2, col='blue')
            })
            dev.off()

            library(animation)
            plot.step <- function(x){
#         browser()
                panels <- layout(matrix(c(1,1,1,2,3,4), nrow=3,ncol=2), widths=c(3,1), heights=c(1,1,1), respect=FALSE)
                plot(ctY ~ ctX, data=sub.solocs[timestep == x,], pch='.', col='gray', xlim=c(0,100), ylim=c(0,100), main=paste('week', x), cex=log1p(sub.solocs[timestep==x, nlive]))
                points(ctY ~ ctX, data=sub.zones[timestep <= x,], pch=3, col='yellow')
#                 points(ctY ~ ctX, data=sub.incidence[timestep <= x & is.inf == 0,], pch='.', col='red')
                points(ctY ~ ctX, data=sub.incidence[timestep == x & is.inf == 1,], pch=3, col='red')
                if (x > detectday) points(ctY ~ ctX, data=sub.detections[timestep <= x,], pch=2, col='blue')
                plot(allnlive ~ timestep, data=sub.solocs[, allnlive := sum(nlive), by=timestep][order(timestep)], main='Live pigs', type='l')
                points(allnlive ~ timestep, data=sub.solocs[timestep==x,], pch=1, cex=2)
                plot(allnlive ~ timestep, data=unique(sub.incidence[, allnlive := sum(is.inf), by=timestep][,.(timestep, allnlive)])[order(timestep)], main='Infected cells', type='l', col='red')
                points(allnlive ~ timestep, data=sub.incidence[timestep==x,], pch=3, col=2, cex=2)
                if (nrow(sub.detections) > 0){
                    temp.detect <- unique(sub.detections[, allnlive := length(unique(loc)), by=timestep][,.(timestep, allnlive)])[order(timestep)][,cs:=cumsum(allnlive)]
                    plot(cs ~ timestep, data=temp.detect, main='Detections', type='l', col='blue')
                    points(cs ~ timestep, data=temp.detect[timestep==x, ], pch=2, col='blue', cex=2)
                }
            }
#             lapply(seq(min(sub.incidence[,timestep]), max(sub.incidence[,timestep])), plot.step)
            out.gif <- paste0('testgif',i,j,k,'.gif')
            saveGIF({
#                 lapply(seq(parameters$burn_weeks, max(sub.incidence[,timestep])), plot.step)
#                 lapply(seq(max(sub.solocs[,timestep])), plot.step)
                lapply(seq(min(sub.incidence[,timestep]) - 10, max(sub.incidence[,timestep])), plot.step)
            }, movie.name = out.gif, ani.width=850, ani.height=600, interval=0.1, imgdir='./test_outputs')
        }, i=unq.combos[,var], j=unq.combos[,land], k=unq.combos[,rep])

        #lapply(list.files(pattern='*.gif'), function(x) file.rename(x, paste0('./test_outputs/',x)))
    }
}
