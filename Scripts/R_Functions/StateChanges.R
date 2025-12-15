#outputs: pop,Incidence,BB
StateChanges<-function(pop,centroids,cells,parameters,Incidence,BB,i){
# StateChanges<-function(pop,centroids,cells,Pbd,B1,B2,F1,F2_int,F2_B,F2i_int,F2i_B,K,death,Pcr,Pir,Incidence,BB,i){
####################################################################
########### Initialize state change probability matrices ###########
####################################################################
list2env(parameters, .GlobalEnv)


# should have list2env for parameters/outputs here for consistency (instead of 1K+ inputs)
#births
Sdpb=matrix(nrow=nrow(pop),ncol=1)
Sdpb[,1]=0
	
#natural deaths
Sdpd<-matrix(nrow=nrow(pop),ncol=1)
Edpd<-matrix(nrow=nrow(pop),ncol=1)
Idpd<-matrix(nrow=nrow(pop),ncol=1)
Rdpd<-matrix(nrow=nrow(pop),ncol=1)
Sdpd[,1]=0
Edpd[,1]=0
Idpd[,1]=0
Rdpd[,1]=0

#disease state change recording
Eep=matrix(nrow=nrow(pop),ncol=1)
Eep[,1]=0
Iep=matrix(nrow=nrow(pop),ncol=1)
Iep[,1]=0
Rep=matrix(nrow=nrow(pop),ncol=1)
Rep[,1]=0
Cep=matrix(nrow=nrow(pop),ncol=1)
Cep[,1]=0

#Carcass decay recording
Ccd=matrix(nrow=nrow(pop),ncol=1)
Ccd[,1]=0
Zcd=matrix(nrow=nrow(pop),ncol=1)
Zcd[,1]=0

	
########################################
########### Determine Births ########### 
########################################

#subset sounder sets with live, uninfected individuals
idN=pop[pop[,8,drop=FALSE]>0|pop[,9,drop=FALSE]>0|pop[,11,drop=FALSE]>0,,drop=FALSE]

#Number of live, uninfected individuals
liveind<-sum(colSums(pop)[c(8,9,11)])

#get row indices of live individuals
liverows<-which(pop[,8,drop=FALSE]>0|pop[,9,drop=FALSE]>0|pop[,11,drop=FALSE]>0) #rownums with live indiv

#density-dependent birth rate
# Brate=Pbd*liveind*(1-liveind/K) ## K is a parameter that may be good to tie to the environment
## I have questions about where the value of Pbd came from -- look up at some point...
## ALSO: technically the baseline growth rate in the logistic equation includes overall population growth, not just reproduction
## (i.e. you could replace Pbd with birth - death) meaning that the stochastic non-infected death rate below might be kind of redundant?
## you can't just throw a logistic equation in there and say 'now we have density dependence', and actually carrying capacity is stupid.
## density dependent growth can be matched with density dependent death (one process of which is disease transmission) in a mechanistic way
## so that you don't have goofy arbitrary parameters that don't reflect real processes.
## What I would do is select locations with at least 2 pigs to have offspring based on per-capita annual reproductive rates, and if there
## are more pigs they are more likely to have offspring. When they do have offspring, they have several. Then, I would have density-dependent
## mortality included separately because it's not like having a lot of pigs is going to make reproduction attempts slow down (the opposite, probably)
## instead we would have higher mortality. Because birth and mortality are separate in the stochastic functions below, we should really have
## separate birth and death rates. The allee effects represented by the toe-slope of the logistic function are also probably not that useful for this
## context as well, because pigs are good at finding each other based on smell, and these are applied in relatively small spatial scales within
## dedicated groups that are very social. (We might even go for a clustering function, as pigs from decimated populations are supposed to go back to their
## original family groups or some such because they are so social, which could have implications for disease spread if these are exposed individuals
## later in the infection for a given location; i.e. anywhere with only 1 or 2 pigs join the nearest group of >2 pigs in neighboring cells at random)
## young pigs also are supposed to respond differently to the virus than older pigs (read that somewhere), so spatiotemporal age dynamics might also be important
## i.e. if you have a bunch of piglets in one spot, they might all catch the virus but also quickly die
reprod.fem = ceiling(rowSums(pop[,c(8,9,11)])) ## reproductive females (assuming 50/50 M/F split in any location; assume half pigs are female)
# reprod.fem = ceiling(rowSums(pop[,c(8,9,11)])/2) ## reproductive females (assuming 50/50 M/F split in any location; assume half pigs are female)
piglets = pmin(10*reprod.fem,rnbinom(nrow(pop), mu=reprod.fem*Pbd*5, size=0.01)) ## more realistic per-cap reproduction by cell
## caps litters at 10 per female, average is 5, times weekly growth rate from 2022 paper. 0.01 gives a good clumping measure. mu is probably wrong-- check mean of nbinom definitions in help file.

#get total births, using Brate as mean in a poisson
# Tbirths=rpois(1,Brate)
#print(paste0("Tbirths: ", Tbirths))
#record total births this time step
# BB[i]=Tbirths
BB[i] = sum(piglets)

#Pick out enough numbers that sum to Tbirths - this will determine how many cells get births
# id<-0
# n=1
# while (sum(id) < Tbirths) { ## this seems like a slow way to do this -- just run rpois on each location
#     birthset_i <- round(runif(1,min=0, max=min(Tbirths,10)))
#  ## also, this forced pop cap wouldn't be necessary, we could just depend on the model's logistic growth rate parameter (would be more defensible/accepted for pubs)
#  ## because then the local groups/populations would be acting according to their own local density dependence, which keeps birth rates aligned with local conditions
#  ## also, isn't there a general litter size for pigs? seems like they usually have a lot (never just one), so a non-parametric birth distribution might be better? e.g. neg.binomial
#  ## in other words, we have the capacity for fine-scale population dynamics, why not use them
#  ## though it depends on the independence of sounders from each other -- maybe they mix a lot more than I am thinking
#     id[n]<-birthset_i
#     n=n+1
# }
# ## my alternative would look like this: Tbirths = rpois(idN, Pbd*rowSums(idN[,c(8,9,11)])*(1-rowSums(idN[,c(8,9,11)])/K))
# ## with K modified to match something like sounder maximum sustainable size
# ## maybe a probability of having a litter based on gestation period, and then given reproduction for an individual, nbinom number of offspring
# ## that would get at the clumpy reproductive dynamics where just scattering piglets across the land would not...
# ## Maybe we should do a density-dependent DEATH rate so that the classes below with natural death (SER) can be responsive to density but not have everyone doubling up on logistic growth death and "out in the world" death
#
#
# #if there are some births
# if(length(id)>1){
#
# #if there are more live sounders than births needed
#   if(length(id)<=nrow(idN)){
#
#
# #%pick which cells with pigs will get the births
#     id2=sample(1:nrow(idN),length(id))
#
#   } else {
# #%if there are more births than cells only add births in cells where the pigs are (so fewer births will be happening)
#     id2=sample(1:length(liverows),length(liverows))
#   }
#
#
#   #assign births to Sdpb
#   for(j in 1:length(id)){ ## for loops are slow in R
#     Sdpb[id2[j],1]<-id[j]
#   }
#
# }


##############################################################
######## Determine disease state change probabilities ######## 
##############################################################

## Susceptible to Exposed?
#Pse<-FOI(pop,centroids,cells,B1,B2,F1,F2) #force of infection #R version
Pse<-FOI_R(pop,centroids,cells,B1,B2,F1,F2_int,F2_B,F2i_int,F2i_B) #cpp parallel version, 22x faster than R version
#Pse<-FOIParallelFull(pop,centroids,cells,B1,B2,F1,F2_int,F2_B,F2i_int,F2i_B) #cpp parallel version, 22x faster than R version
## gives negative probabilities for cells without infected ? Might be because I have the wrong values for F2_int, etc.
Pse[Pse < 0] <- 0 ## assuming these are actually probabilities

## Exposed to infected
Pei=1-exp(-1/(rpois(cells,4)/7)) #transitions exposure to infected
## why rpois(cells, 4)? and /7? if they are exposed, isn't there a probability of infection?
## Or is this based on R0... somehow... ?? (same question for Pic)
## Ah, trying to incorporate incubation and infectious periods (as per 2022 paper) into a weekly probability
## ...seems weird to do this with variation per spatial cell (also, numbers should be parameters to be clear when the model is reviewed...)
## again, even if you want to have variation in these parameters, Poisson is NOT the function I would choose because
## the variance is waaaaaay too high (it will have a lot of very short incubations and infectious periods if the average is 4 or 5)
## this would also be an easy spot to build out some caches to pre-calculate probabilities based on rpois outputs for some speed gains (maybe)
## OHHHHHHH.... This is the probability that the transition happens on a given week, so depending on what day of the week the infection happens
## it has a different probability of actually moving from exposed to infected. But then... why not just sample seq(7) and assign probabilities
## based on the day of the week? (Again, cacheing)

## Infected to contageous corpse or recovered, depending on Pir (recovery rate)
Pic=1-exp(-1/(rpois(cells,5)/7)) #transitions infected to either dead or recovered

# print(nrow(pop))
# if (Sys.time()-state.time0 > 2) browser()
######For troubleshooting ML/R infection process diffs######
#if(!dir.exists(
#  paste0("/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/ASF_optimal_radius/Output/Pse_Compare_Trblsht/temp/",i))){
#  dir.create(paste0("/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/ASF_optimal_radius/Output/Pse_Compare_Trblsht/temp/",i))
#}

#saveRDS(Pse,paste0("/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/ASF_optimal_radius/Output/Pse_Compare_Trblsht/temp/",i,"/Pse.rds"))
#######

###############################################
######## Conduct the State Transitions ######## 
###############################################
## the problem with this structure is you can have inter-dimensional swine-shifting as they spontaneously appear and disappear from the temporal plane

## density dependent death
## change 25's to modifiable parameter
death=exp(-25+((1+death)*rowSums(pop[,c(8:11)])))/(1+exp(-25+((1+death)*rowSums(pop[,c(8:11)]))))
## when population is too high, top and bottom of the above function become Inf, leading to NaN values which cause errors in rbinom()
## to prevent errors, change NaN values to 1 because this is a probability
death[is.nan(death)] <- 1
#susceptible state changes
Sdpd=rbinom(nrow(pop), pop[,8],death) ## death is given as a rate, not a probability?
# Eep=rbinom(nrow(pop), pop[,8], Pse[pop[,3]]) ## generates NA when Pse has negative values; error when Pse has NA
Eep=rbinom(nrow(pop), pop[,8] - Sdpd, Pse[pop[,3]])
## could suppress warnings (see fix above), would rather have an internally consistent probability function.
## Eep = rbinom(nrow(pop), pop[,8]-Sdpd, Pse[pop[,3]]) if we have death separate from birth (i.e. ditch logistic growth equation)

#exposed state changes
Edpd=rbinom(nrow(pop), pop[,9], death) ## exposed natural death
# Iep=rbinom(nrow(pop), pop[,9], Pei[pop[,3]]) ## exposed -> infected
Iep=rbinom(nrow(pop), pop[,9] - Edpd, Pei[pop[,3]]) ## exposed -> infected
## this is where the sequential calculations could help to avoid the negatives -- set order of events
## such as reproduction, natural death, infection spread, infection death/recovery, infection from corpses
## The choice depends on what you want to emphasize in the model
## (in this case, I would go with a spread rate after reproduction and before death for live-live transitions,
## and then dead-live at the end or beginning of the next timestep to maximize the effect of disease spread
## since that is the purpose of the model)
## Not to be blunt, but this would be rejected by any reviewer worth their salt (but YMMV)
## just set Iep = rbinom(nrow(pop), pop[,9] - Edpd, Pei[pop[,3]])
## assuming some will be left in exposed but not infected? Otherwise do Iep = pop[,9] - Edpd
## A more elegant way would be to use rmultinom to do all the rbinoms simultaneously, e.g.
## rmultinom(nrow(pop), pop[,9], c(0, 0, Pei[pop[,3]], 0, death)) for none staying in E for more than a timestep
## (rmultinom automatically adjusts supplied probabilities to sum to 1);
## rmultinom(nrow(pop), pop[,9], c(0, ((1-Pei[pop[,3]])*(1-death) - Pei[pop[,3]*death]), Pei[pop[,3]], 0, death))
## for some remaining in E beyond one timestep
## (double check; something like that. It's been a while since probability theory; also, have to make sure the function would work on a matrix/df as written)


#infected state changes
Rep=rbinom(nrow(pop),pop[,10],Pir*Pic[pop[,3]]) ## infected -> recovery
# Cep=rbinom(nrow(pop),pop[,10],(1-Pir)*(Pic[pop[,3]])) ## infected -> contageous carcass
Cep = pop[,10] - Rep ## assuming we are expecting all infecteds to die or recover by the end of the week (sounds right)
## add stochasticity in where the division is between the two outcomes for the population, not the numbers of individuals (avoids overlap and the negatives issues below)
## suggest: Cep = pop[,10] - Rep

#recovered state changes
Rdpd=rbinom(nrow(pop),pop[,11],death) ## recovered -> natural death


#infected carcass state changes
Ccd=rbinom(nrow(pop),pop[,12],Pcr) ## infectious carcass removal (natural) -- depends on rate of contact with other pigs at all?

#uninfected carcass state changes
## do we actually care about the location and duration of Z carcasses? do they attract pigs (I know they can be cannibalistic sometimes)
## or is this really to keep track of the accounting to make sure things add up (which hasn't really worked up to this point; see above)
Zcd=rbinom(nrow(pop),pop[,13],Pcr) ## non-infectious carcass removal (natural)

### okay, actually I think we need to do a flat reproductive rate (or the clustered reproduction I talked about... somewhere...) and density dependent death of all individuals
### after (or before) everything else. Maybe not infecteds, but that way you could impose density dependent death on everyone equally instead of all this double-counting garbage.
### you could avoid the stupid carrying capacity thing, and have a function that saturates at whatever level you want (say, 2x the normal size for a family group)
### then you also have put the model focus on infection (i.e. letting it spread as far as possible) to get the outside range of dispersal since there's a lot of unknowns about
### how all these other factors influence the population dynamics, movement, etc.

#collect incidence info
Incidence[i]<-Incidence[i]+sum(Eep)

###################################
######## Update pop matrix ######## 
###################################

#update states in pop matrix
###################################
pop[,8]=pop[,8]-Eep+piglets-Sdpd #S
# pop[,8]=pop[,8]-Eep+Sdpb-Sdpd #S
pop[,9]=pop[,9]-Iep+Eep-Edpd #E
pop[,10]=pop[,10]-Rep-Cep+Iep#-Idpd#I
pop[,11]=pop[,11]+Rep-Rdpd #R
pop[,12]=pop[,12]+Cep-Ccd#+Idpd #C
pop[,13]=pop[,13]+Sdpd+Rdpd+Edpd-Zcd #Z

###ID where there are negative numbers
##then, check to see if part of duplicated cell number
#if(any((pop[,8]<0|pop[,9]<0|pop[,10]<0|pop[,11]<0|pop[,12]<0|pop[,13]<0))){
#  neg.index=which(pop[,8]<0|pop[,9]<0|pop[,10]<0|pop[,11]<0|pop[,12]<0|pop[,13]<0)
#  cellnums=pop[neg.index,3]
#  pop.cellnums=pop[pop[,3]%in%cellnums,,drop=FALSE]
#  if(any(duplicated(pop.cellnums[,3]))){
#    pop.dups=pop[which(pop[,8]<0|pop[,9]<0|pop[,10]<0|pop[,11]<0|pop[,12]<0|pop[,13]<0)&((duplicated(pop[,3],fromLast=FALSE))|duplicated(pop[,3],fromLast=TRUE)),]
#    print(pop.dups[order(pop.dups[,3]),])
#    print(pop.dups[,12])
#    print(pop.dups[,13])
#  }
  
#}
#pop[which(pop[,8]<0|pop[,9]<0|pop[,10]<0|pop[,11]<0|pop[,12]<0|pop[,13]<0),]

#sometimes end up with negative numbers 
#(i.e. all pigs in sounders chosen for natural mort and disease mort)
#just set anything below zero to zero
## ah.... I don't love that for consistency/integrity/realism... what does it mean when they go negative? population crash momentum or some stochastic piece going beyond its boundaries?
## it's generally better to apply those things one at a time to get e.g. number who die naturally (other causes than disease/culling), and from that the number who die from disease
## not necessary if we make the changes above
# pop[which(pop[,8]<0),8]<-0
# pop[which(pop[,9]<0),9]<-0
# pop[which(pop[,10]<0),10]<-0
# pop[which(pop[,11]<0),11]<-0
# pop[which(pop[,12]<0),12]<-0
# pop[which(pop[,13]<0),13]<-0

#move dead individuals (C or Z) into their own rows
#pop[,12] and pop[,13] > 0
deadguys<-pop[pop[,12]>0|pop[,13]>0,,drop=FALSE]

#if there are deadguys....
if(nrow(deadguys)!=0){
#remove abundance and all live guy counts from deadguy set
deadguys[,1]=0
deadguys[,8]=0
deadguys[,9]=0
deadguys[,10]=0
deadguys[,11]=0

#set all deadguys in pop rows to zero
# pop[which(pop[,12]>0),12]<-0
pop[,12][pop[,12] > 0] <- 0 ## as in other places, this handles the 1 row matrix issue
# pop[which(pop[,13]>0),13]<-0
pop[,13][pop[,13] > 0] <- 0

#add deadguys to pop matrix
pop<-rbind(pop,deadguys) ## does this make duplicate rows in the pop matrix?
## yes? I guess the columns add up to the right numbers... but ...why?

}

#Update abundance numbers (live individuals only count in abundance)
pop[,1]=rowSums(pop[,8:11,drop=FALSE]) ## drop=FALSE to handle single row matrix case
# })

return(list(pop,Incidence,BB,"Eep"=sum(Eep),"Sdpb"=sum(Sdpb),"Sdpd"=sum(Sdpd),"Iep"=sum(Iep),"Edp"=sum(Edpd),"Rep"=sum(Rep),"Cep"=sum(Cep),"Rdpd"=sum(Rdpd),"Ccd"=sum(Ccd),"Zcd"=sum(Zcd)))

}
