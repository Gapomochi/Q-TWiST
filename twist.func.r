
#-----------------------------------------------------------------------------------------------------------------------------------------------#
#												qtwist											#
#-----------------------------------------------------------------------------------------------------------------------------------------------#
#																								#
# INPUT																							#
# visit : matrix visit (patnum/delay_tox/grade_tox)																#
# fu : matrix last news (patnum/delay_rfs/status_rfs/delay_death/status_death/group)										#
#	 for "group" : 0=control group and 1=experimental group														#
# gm : integer >0 representing the maximal toxicity grade taken into account 												#
# tp : integer >0 representing the time were estimates are given in the output		 									#
# utility : vector of length 3 representing the utility score for TOX, TWiST and REL in this order								#
# fix_utility : integer between 1 and 3 representing health state for which the utility score is fixed for the threshold 				#
#		    analysis (1=TOX, 2=TWiST, 3=REL)																#
# 		    The value of the utility score will be the same as the value given in the vector 'utility'							# 
# nb_boot : integer >0 representing the number of bootstrap sample used for the estimates										#
#																								#
# OUTPUT																							#
# qtwist : list of length 2 :																				#
#		- survival : list of length 6 containing survival data for TOX, RFS and OS for group 0 and 1							#
#		- threshold : values used for the threshold analysis														#
# Q-TWiST.PDF : PDF file containing results and graphs of the Q-TWiST analysis ; exported in the work directory						#
#-----------------------------------------------------------------------------------------------------------------------------------------------#

qtwist <- function(visit, fu, gm, tmax, tp, utility, fix_utility, nb_boot) {


	#-----------------------------------------------------------------------#
	#		donnee.func (function which fine shape data)			#
	#												#
	# INPUT											#
	# - visit (patnum/delay_tox/grade_tox)						#
	# - fu patnum/delay_rfs/status_rfs/delay_death/status_death/group)	#
	#												#
	# OUTPUT											#
	# - jdd : data with 8 columns (patnum/del_tox/tox/	 			#
	#	    del_rfs/rfs/del_death/death/gpe)					#
	#-----------------------------------------------------------------------#


		donnee.func <- function(visit, fu) {
			#--- column name
				dimnames(visit)[[2]][1] = "patnum"
				dimnames(visit)[[2]][2] = "del_tox"
				dimnames(visit)[[2]][3] = "gr_tox"

				dimnames(fu)[[2]][1] = "patnum"
				dimnames(fu)[[2]][2] = "del_rfs"
				dimnames(fu)[[2]][3] = "rfs"
				dimnames(fu)[[2]][4] = "del_death"
				dimnames(fu)[[2]][5] = "death"
				dimnames(fu)[[2]][6] = "gpe"


			#--- status 0/1 for TOX
				tox = matrix(0,dim(visit)[1],1)
				for ( i in 1:dim(visit)[1] ) {
					if (visit[i,3] >= gm ) {tox[i,1] = 1 }
					else { tox[i,1] = 0 }
				}
				visit = cbind (visit, tox)

			#--- data with delay and status according to the health state (TOX or FU=Follow up)
				visit2 = cbind( visit[,c(1,2,4)], data.frame(type=rep("tox", dim(visit)[1])) )
				colnames(visit2)=c("patnum", "del", "stat", "type")
	
				fu2 = cbind( fu[,1], del=rep(NA, dim(fu)[1]), stat=rep(NA, dim(fu)[1]), data.frame(rep("fu", dim(fu)[1])))
				colnames(fu2)=c("patnum", "del", "stat", "type")
					for ( i in 1:dim(fu)[1] ) {
					if (fu[i,3]==0 & fu[i,5]==0) { fu2[i,2]=max(fu[i,2], fu[i,4]) ; fu2[i,3]=0 }
					if (fu[i,3]==1 & fu[i,5]==1) { fu2[i,2]=min(fu[i,2], fu[i,4]) ; fu2[i,3]=1 }
					if (fu[i,3]==0 & fu[i,5]==1) { fu2[i,2]=fu[i,4] ; fu2[i,3]=1 }
					if (fu[i,3]==1 & fu[i,5]==0) { fu2[i,2]=fu[i,2] ; fu2[i,3]=1 }
				}

				visit3=rbind(visit2, fu2)
				visit3=visit3[order(visit3[,2]),] 
				visit3=visit3[order(visit3[,1]),] 

				# Delete all rows > fu
				visit3=cbind(visit3, pb=matrix(0,dim(visit3)[1],1))
				for (i in 1:(dim(visit3)[1]-1) ) {
					if (visit3[i,1]==visit3[i+1,1] & visit3[i,4]=="fu") {visit3[i,5]=1}
					if (visit3[i,5]==1 & visit3[i,1]==visit3[i+1,1]) {visit3[i+1,5]=1}
				}
				if(length(which(visit3$type=="tox" & visit3$pb==1))>0) {visit3=visit3[-which(visit3$type=="tox" & visit3$pb==1),]}
				visit3=visit3[,-5]
			
			#--- Duration spent in TOX per patient
				visit3 = cbind (visit3, dur_tox=matrix(0, dim(visit3)[1], 1))
					for ( i in 1:(dim(visit3)[1]-1) ) {
					if (visit3[i,1] == visit3[i+1,1] ) { visit3[i,5] = visit3[i+1,2] - visit3[i,2] }
					if (visit3[i,4] == "tox" & visit3[i+1,4] == "fu" &  visit3[i,2]==visit3[i+1,2]) { visit3[i,5] = 1 }
				}

				visit4 = visit3[-which(visit3$stat==0 | visit3$type=="fu"),] # delete all rows in which tox = 0 or type = fu

			#--- if visit4 empty : no toxicity
				if ( dim(visit4)[1]==0 ) {
					jdd = cbind(fu, rep(0,dim(fu)[1]), rep(1,dim(fu)[1]) ) # tox=1 because if tox=0 survival curve for TOX is 1
					jdd = as.data.frame( jdd[,c(1,7:8,2:5,6)] )
					dimnames(jdd)[[2]][2] = "del_tox"
					dimnames(jdd)[[2]][3] = "tox"
				}

			#--- if visit4 not empty : there are toxicities
				if ( dim(visit4)[1]>0 ) {
					res=aggregate(visit4$dur_tox, list(visit4$patnum), sum) # add per patient "del" variable corresponding to the duration spent in TOX
					res=cbind(res,rep(1,dim(res)[1]))
					colnames(res)<-c("patnum", "del_tox", "tox")

					patnum=data.frame(fu[,1])
					colnames(patnum)=c("patnum")
					visit5 = merge(patnum, res, all.x=T)
					visit5[c(which(is.na(visit5$tox))),3] = 1 
					visit5[c(which(is.na(visit5$del_tox))),2] = 0

					#--- merge fu file
						jdd = merge(fu, visit5)
						jdd = jdd[,c(1,7:8,2:5,6)]
				}


			#--- replace delay so that RFS corresponds to survival until relapse or death
				for ( i in 1:dim(jdd)[1] ) {
					if(jdd[i,5]==1 & jdd[i,7]==1) {jdd[i,4] = min(jdd[i,4], jdd[i,6])} # si rfs=os=0 alors del_rfs = min(del_os, del_rfs)
					if(jdd[i,5]==0 & jdd[i,7]==0) {jdd[i,4] = max(jdd[i,4], jdd[i,6])} # si rfs=os=1 alors del_rfs = max(del_os, del_rfs)
					if(jdd[i,5]==0 & jdd[i,7]==1) {jdd[i,4] = jdd[i,6]} # si rfs=0 & os=1 alors del_rfs = del_os
					if(jdd[i,5]==1 & jdd[i,7]==0) {jdd[i,4] = jdd[i,4]} # si rfs=1 & os=0alors del_rfs = del_rfs
 				}


			return(jdd)
		}




	#-----------------------------------------------------------------------#
	# 	surv.func (calculate survival according to Kaplan Meier)		#
	#												#
	# INPUT 											#
	# - time : delay 										#
	# - event : event (indicator 0-1)							#
	#												#
	# OUTPUT											#
	# - surv : Kaplan Meier	survival model						#
	#-----------------------------------------------------------------------# 

		surv.func <- function(time, event) {
			surv = survfit(Surv(time, event)~1)
			return(surv)
		}



	#-----------------------------------------------------------------------------------------------------------#
	#							boot.func (bootstrap sample)							#
	# INPUT																	#
	# - nb_boot : integer representing the number of boostrap sample 								#
	# - jdd : data  (column name : patnum/del_tox/tox/del_rfs/rfs/del_death/death)					#
	# - timepoint : number corresponding to the timepoint										#
	#																		#
	# OUTPUT																	#
	# - rboot : restricted means for TOX, RFS and OS for each nb_boot sample (dim = nb_boot 3)			#
	#-----------------------------------------------------------------------------------------------------------#

		boot.func <- function (nb_boot, jdd, timepoint) {
			rboot = NULL
			for (i in 1:nb_boot) {
				# Data bootstraped from jdd
				sample = sample(jdd$patnum, replace=T)
				mat=NULL
				for (j in sample) {
					k = which(jdd$patnum == j)
					mat = rbind(mat, jdd[k,])
				}

				# calculate rmean for TOX, RFS and OS of this bootstraped data 
					rmean = NULL
					rmean = cbind( rmean, summary( survfit(Surv(mat$del_tox, mat$tox)~1) , rmean = timepoint)$table[5] )
					rmean = cbind( rmean, summary( survfit(Surv(mat$del_rfs, mat$rfs)~1) , rmean = timepoint)$table[5] )
					rmean = cbind( rmean, summary( survfit(Surv(mat$del_death, mat$death)~1) , rmean = timepoint)$table[5] )
					colnames(rmean) = c("TOX", "RFS", "OS")
					rownames(rmean) = i

				# Matrix with nb_boot rmean bootstraped for TOX, RFS and OS
				rboot = rbind(rboot, rmean)	
			}
			return(list(nb_boot=nb_boot, rboot=rboot))
		}



	#-----------------------------------------------------------------------------------------------------------------#
	#						threshold.func (threshold analysis)								#
	# 																			#
	# INPUT																		#
	# - fix_utility : number representing health state for which score is fixed 							#
	# - mean : vector of length 3 (restricted mean of time spent in TOX, TWiST and REL)						#
	# - dur_boot : matrix containing bootstraped restricted means of time spent in TOX, TWiST and PROG 			#
	# - step : number >0 and <1 representing the step of variation for the 2 utility scores					#
	# 																			#
	# OUTPUT																		#
	# - threshold : matrix containing 2 utility scores varying, quality-adjusted survival mean,				#
	# 		    standard error, 95% confidence interval and p.value for all combination for the 2 utility scores	#
	#		    which vary from 0 to 1 and according to the step parameter  							#
	#-----------------------------------------------------------------------------------------------------------------#

		threshold.func <- function (fix_utility, mean, dur_boot, step) {
			threshold = matrix(0,(1/step+1)^2,7)
			if (fix_utility == 1 ) {
				colnames(threshold) = c("utility_twist", "utility_rel", "mean", "se", "CIlow", "CIupp", "p.value")
				nb = 1
				for (i in seq(0,1,by=step)) { 
					for (j in seq(0,1,by=step)) { 
						threshold[nb,c(1,2)] =  round(c(i,j),5)					
						threshold[nb,3] = round( utility[1]*mean[1] + i*mean[2] + j*mean[3],3)
						threshold[nb,4] = round( apply(matrix(utility[1]*dur_boot[,1]+i*dur_boot[,2]+j*dur_boot[,3]), 2, sd) ,3)
						threshold[nb,5] = round( threshold[nb,3]-1.96*(threshold[nb,4]) ,3)			
						threshold[nb,6] = round( threshold[nb,3]+1.96*(threshold[nb,4]) ,3)
 						threshold[nb,7] = round(2*(1-pnorm(abs(threshold[nb,3])/threshold[nb,4])),3)
						nb = nb + 1
					}	
				}
			}
			if (fix_utility == 2 ) {
				colnames(threshold) = c("utility_tox", "utility_rel", "mean", "se", "CIlow", "CIupp", "p.value")
				nb = 1
				for (i in seq(0,1,by=step)) { 
					for (j in seq(0,1,by=step)) { 
						threshold[nb,c(1,2)] =  round(c(i,j),5)							
						threshold[nb,3] = round( i*mean[1] + utility[2]*mean[2] + j*mean[3],3)
						threshold[nb,4] = round( apply(matrix(i*dur_boot[,1]+utility[2]*dur_boot[,2]+j*dur_boot[,3]), 2, sd) ,3)
						threshold[nb,5] = round( threshold[nb,3]-1.96*(threshold[nb,4]) ,3)			
						threshold[nb,6] = round( threshold[nb,3]+1.96*(threshold[nb,4]) ,3)
 						threshold[nb,7] = round(2*(1-pnorm(abs(threshold[nb,3])/threshold[nb,4])),3)
						nb = nb + 1
					}	
				}
			}
			if (fix_utility == 3 ) {
				colnames(threshold) = c("utility_tox", "utility_twist", "mean", "se", "CIlow", "CIupp", "p.value")
				nb = 1
				for (i in seq(0,1,by=step)) { 
					for (j in seq(0,1,by=step)) { 
						threshold[nb,c(1,2)] = round(c(i,j),5)									
						threshold[nb,3] = round( i*mean[1] + j*mean[2] + utility[3]*mean[3],3)
						threshold[nb,4] = round( apply(matrix(i*dur_bootc[,1]+j*dur_bootc[,2]+utility[3]*dur_bootc[,3]), 2, sd) ,3)
						threshold[nb,5] = round( threshold[nb,3]-1.96*(threshold[nb,4]) ,3)			
						threshold[nb,6] = round( threshold[nb,3]+1.96*(threshold[nb,4]) ,3)
 						threshold[nb,7] = round(2*(1-pnorm(abs(threshold[nb,3])/threshold[nb,4])),3)
						nb = nb + 1
					}	
				}
			}
			return(threshold)
		}



	#-----------------------------------------------------------------------#
	#			  		DATA CHECK				 	  	#
	#-----------------------------------------------------------------------#

	cat("+------------------------------------------------------------------+\n")
	cat("|                         Q-TWiST ANALYSIS                         |\n") 
	cat("|                      (E. Bogart & A. Kramar)                     |\n")
	cat("|------------------------------------------------------------------|\n")
	cat("| Centre Oscar Lambret                                             |\n")
	cat("| 3 rue Frédéric Combemale                                         |\n")
	cat("| 59000 LILLE                                                      |\n")
	cat("| France                                                           |\n")
	cat("+------------------------------------------------------------------+\n")
	cat(" \n")


	cat("+--------------------------------------------+\n")
	cat("|                    Checking                |\n")
	cat("+--------------------------------------------+\n")
	cat(" \n")

	#----- check of fu
		if(is.matrix(fu)==FALSE) stop("'fu' must be a matrix !")
		if(dim(fu)[2] != 6) stop("'fu' not appropriate !")

		for (i in 1:6) {
			if(is.numeric(fu[,i])==FALSE) stop("'fu' must be numeric !")
		}

		for (i in 1:6) {
			for (j in 1:dim(fu)[1]) {
				if(is.na(fu[j,i])==TRUE) stop("There are missing values in 'fu' file !")
			}
		}

		for(i in 1:length(fu[,1])){
				if(fu[i,1] <= 0) stop("ID patient must be strictly positive !")
				if(fu[i,2] < 0) stop("Duration need to be positive !")
				if(fu[i,3] != 0 & fu[i,3] !=1) stop("Event must be 0 or 1 !")
				if(fu[i,4] < 0) stop("Duration need to be positive !")
				if(fu[i,5] != 0 & fu[i,5] !=1) stop("Event must be 0 or 1 !")
				if(fu[i,6] !=0 & fu[i,6] !=1) stop("Number of group need to be 0 or 1 !")
				if(fu[i,2] > fu[i,4]) stop("time until relapse > time until death !")
		}
		
		d.list.pat <- fu[,1]
		if(length(d.list.pat) != length(unique(d.list.pat))) stop("Patients must not have the same ID !") 
		if(max(abs(diff(fu[,6]))) > 1) stop("Number of group must be consecutive integer(s) !")

		# we need 2 groups of patients, no more and no less
		group=levels(as.factor(fu[,6]))
		if( length(group)!=2 ) stop("Number of groups must be 2 !")
		if( group[1]!=0 | group[2]!=1 ) stop("Values for group value must be 0 or 1 !")

	#----- check of visit
		if(is.matrix(visit)==FALSE) stop("'visit' must be a matrix !")
		if(dim(visit)[2] != 3) stop("'visit' not appropriate !")

		for (i in 1:3) {
			if(is.numeric(visit[,i])==FALSE) stop("'visit' must be numeric !")
		}

		for (i in 1:3) {
			for (j in 1:dim(visit)[1]) {
				if(is.na(visit[j,i])==TRUE) stop("There are missing values in 'visit' file !")
			}
		}

		s.list.pat <- unique(visit[,1])
		for(i in 1:length(s.list.pat)){
				if(length(which(s.list.pat[i]==d.list.pat)) < 1) stop("ID patient in 'visit' file does not exist in 'fu' file !")
		}
		e <- 0
		for(i in 1:length(visit[,1])){
				if(visit[i,2] < 0) stop("Durations need to be positive !")
				if(visit[i,2] > fu[fu[,1]==visit[i,1],2]){
					#visit[i,1] <- 0 
					e <- e+1
				}
		}

		for(i in s.list.pat){ # it is not possible to have twice the same delay in visit
			x = visit[visit[,1]==i,2]
			if (length(x)>1) {
				x = sort(x)
				for (j in 1:(length(x)-1)) {
					if ( x[j]>=x[j+1] ) { stop("Delays for the same patient in 'visit' file must be different !") }
				}
			}
		}

		visit_check <- visit[visit[,1] > 0,] 
		for(i in s.list.pat){
				x <- visit_check[visit_check[,1]==i,2:3]
				if(is.matrix(x)) visit_check[visit_check[,1]==i,2:3] <- as.matrix(merge(sort(x[,1]),x,by=1))			
		}
		if(e > 0) cat("Visits in 'visit' file after last visit in 'fu' file are ignored\n")

	#----- check of gm
		if(gm!=1 & gm!=2 & gm!=3 & gm!=4 & gm!=5) stop ("Parameter 'gm' not appropriate !")
		#if (max(visit[,3])<gm) stop("gm must be <= to the maximal grade of toxicity in 'visit' file !")

	#----- check of tp
		if(tp <= 0) stop ("Parameter 'tp' not appropriate")

	#----- check of utility
		if(length(utility)!=3) stop ("Parameter 'utility' not appropriate !")
		for (i in 1:3) {if(utility[i]<0 | utility[i]>1) stop ("Parameter 'utility' not appropriate !")}
		
	#----- check of fix_utility
		if(fix_utility!=1 & fix_utility!=2 & fix_utility!=3) stop ("Parameter 'fix_utility' not appropriate !")

	#----- check of nb_boot
		if(missing(nb_boot)) (nb_boot=10)
		if(nb_boot<10) stop ("Parameter 'nb_boot' must be >10 !")
		if(nb_boot != round(nb_boot,0)) stop ("Parameter 'nb_boot' must be an integer !")
		if(nb_boot<500) cat("Parameter 'nb_boot' < 500 so estimations not precise\n")


	cat("Data OK\n")
	cat(" \n")





	#-----------------------------------------------------------------------#
	#			  		EVENT SUMMARY		 	  		#
	#-----------------------------------------------------------------------#

	jdd=donnee.func(visit, fu)

		summary_tox = summary_rfs = summary_os = matrix(0,3,3)

		if ( sum(jdd[,2])==0 ) { # if there in no toxcity grade >=gm in the visit file
		summary_tox[1,1] = length(which(jdd$tox==1 & jdd$gpe==0))
		summary_tox[2,1] = length(which(jdd$tox==1 & jdd$gpe==1))
		summary_tox[1,2] = summary_tox[2,2] = 0
		}
		if ( sum(jdd[,2])>0 ) {
		summary_tox[1,1] = length(which(jdd$del_tox==0 & jdd$gpe==0))
		summary_tox[1,2] = length(which(jdd$del_tox>0 & jdd$gpe==0))
		summary_tox[2,1] = length(which(jdd$del_tox==0 & jdd$gpe==1))
		summary_tox[2,2] = length(which(jdd$del_tox>0 & jdd$gpe==1))
		}
		summary_tox[1,3] = summary_tox[1,1] + summary_tox[1,2]
		summary_tox[2,3] = summary_tox[2,1] + summary_tox[2,2]
		summary_tox[3,1] = summary_tox[1,1] + summary_tox[2,1]
		summary_tox[3,2] = summary_tox[1,2] + summary_tox[2,2]
		summary_tox[3,3] = summary_tox[1,3] + summary_tox[2,3]

		summary_rfs[1,1] = length(which(jdd$rfs==0 & jdd$gpe==0))
		summary_rfs[1,2] = length(which(jdd$rfs==1 & jdd$gpe==0))
		summary_rfs[2,1] = length(which(jdd$rfs==0 & jdd$gpe==1))
		summary_rfs[2,2] = length(which(jdd$rfs==1 & jdd$gpe==1))
		summary_rfs[1,3] = summary_rfs[1,1] + summary_rfs[1,2]
		summary_rfs[2,3] = summary_rfs[2,1] + summary_rfs[2,2]
		summary_rfs[3,1] = summary_rfs[1,1] + summary_rfs[2,1]
		summary_rfs[3,2] = summary_rfs[1,2] + summary_rfs[2,2]
		summary_rfs[3,3] = summary_rfs[1,3] + summary_rfs[2,3]

		summary_os[1,1] = length(which(jdd$death==0 & jdd$gpe==0))
		summary_os[1,2] = length(which(jdd$death==1 & jdd$gpe==0))
		summary_os[2,1] = length(which(jdd$death==0 & jdd$gpe==1))
		summary_os[2,2] = length(which(jdd$death==1 & jdd$gpe==1))
		summary_os[1,3] = summary_os[1,1] + summary_os[1,2]
		summary_os[2,3] = summary_os[2,1] + summary_os[2,2]
		summary_os[3,1] = summary_os[1,1] + summary_os[2,1]
		summary_os[3,2] = summary_os[1,2] + summary_os[2,2]
		summary_os[3,3] = summary_os[1,3] + summary_os[2,3]

		rownames(summary_tox)=rownames(summary_rfs)=rownames(summary_os)=c("Group 0", "Group 1", "Total")
		colnames(summary_tox)=colnames(summary_rfs)=colnames(summary_os)=c("Censored", "Failed", "Total")


	cat("+--------------------------------------------+\n")
	cat("|                Event summary               |\n")
	cat("+--------------------------------------------+\n")
	cat(" \n")
	cat("Event : TOX\n")
	cat(paste("Minimal grade taken into acount : ", gm))
	cat(" \n")
	cat(" \n")
	print(summary_tox)
	cat(" \n")
	if(summary_tox[1,2]<10) cat("Less than 10 failures in TOX for group 0 : estimations can be biaised for TOX\n")
	if(summary_tox[2,2]<10) cat("Less than 10 failures in TOX for group 1 : estimations can be biaised for TOX\n")
	cat("\n")
	cat("Event : RFS\n")
	cat(" \n")
	print(summary_rfs)
	cat(" \n")
	cat("Event : OS\n")
	cat(" \n")
	print(summary_os)
	cat(" \n")


	#-----------------------------------------------------------------------#
	#			  		RESTRICTED MEANS					#
	#-----------------------------------------------------------------------#

	#--- Restricted mean estimate
		# if there is no toxicity in one group tox=1 for all this group (otherwise survival for TOX = 1)
		if ( sum(jdd[which(jdd[,8]==0),2])==0 & sum(jdd[which(jdd[,8]==0),3])==0 ) { jdd[which(jdd[,8]==0),3]=1} # group 0
		if ( sum(jdd[which(jdd[,8]==1),2])==0 & sum(jdd[which(jdd[,8]==1),3])==0 ) { jdd[which(jdd[,8]==1),3]=1} # group 1

		surv_toxa = surv.func(jdd[(jdd$gpe)==0,2], jdd[(jdd$gpe)==0,3])
		surv_toxb = surv.func(jdd[(jdd$gpe)==1,2], jdd[(jdd$gpe)==1,3])
		surv_osa = surv.func(jdd[(jdd$gpe)==0,6], jdd[(jdd$gpe)==0,7])
		surv_osb = surv.func(jdd[(jdd$gpe)==1,6], jdd[(jdd$gpe)==1,7])
		surv_rfsa = surv.func(jdd[(jdd$gpe)==0,4], jdd[(jdd$gpe)==0,5])
		surv_rfsb = surv.func(jdd[(jdd$gpe)==1,4], jdd[(jdd$gpe)==1,5])

		rm = matrix(0,2,3)
		colnames(rm) = c("TOX", "RFS", "OS")
		rownames(rm) = c("Group 0", "Group 1")
		rm[1,1] = round(summary(surv_toxa, rmean = tp)$table[5],3)
		rm[1,2] = round(summary(surv_rfsa, rmean = tp)$table[5],3)
		rm[1,3] = round(summary(surv_osa, rmean = tp)$table[5],3)
		rm[2,1] = round(summary(surv_toxb, rmean = tp)$table[5],3)
		rm[2,2] = round(summary(surv_rfsb, rmean = tp)$table[5],3)
		rm[2,3] = round(summary(surv_osb, rmean = tp)$table[5],3)


	#--- Restricted mean covariance estimate bootstraped
		rboota = boot.func(nb_boot=nb_boot, jdd=jdd[(jdd$gpe)==0,], timepoint=tp)
		rmboota = NULL
		rmboota = rbind(rmboota, round(cov(rboota$rboot),3) )
		rmboota = rbind( rmboota, round(apply(rboota$rboot, 2, mean),3) ) 
		rmboota = rbind( rmboota, round(apply(rboota$rboot, 2, sd),3) ) 
		rownames(rmboota)[c(4:5)] = c("Restricted mean", "SD")

		rbootb = boot.func(nb_boot=nb_boot, jdd=jdd[(jdd$gpe)==1,], timepoint=tp)
		rmbootb = NULL
		rmbootb = rbind( rmbootb, round(cov(rbootb$rboot),3) )
		rmbootb = rbind( rmbootb, round(apply(rbootb$rboot, 2, mean),3) ) 
		rmbootb = rbind( rmbootb, round(apply(rbootb$rboot, 2, sd),3) ) 
		rownames(rmbootb)[c(4:5)] = c("Restricted mean", "SD")

		rbootc = rbootb$rboot - rboota$rboot
		rmbootc = NULL
		rmbootc = rbind( rmbootc, round(cov(rbootc),3) )
		rmbootc = rbind( rmbootc, round(apply(rbootc, 2, mean),3) ) 
		rmbootc = rbind( rmbootc, round(apply(rbootc, 2, sd),3) ) 
		rownames(rmbootc)[c(4:5)] = c("Restricted mean", "SD")

	rmboot = list(gpe0 = rmboota, gpe1 = rmbootb, gpe1vs0 = rmbootc)


	cat("+--------------------------------------------+\n")
	cat("|           Restricted mean estimates        |\n")
	cat("+--------------------------------------------+\n")
	cat(" \n")
	cat(paste("Time point : ", tp))
	cat(" \n")
	cat(" \n")
	print(rm)
	cat(" \n")

	cat("+--------------------------------------------+\n")
	cat("|    Restricted mean covariance estimates    |\n")
	cat("+--------------------------------------------+\n")
	cat(" \n")
	cat(paste("Time point : ", tp))
	cat(" \n")
	cat("Number of bootstrap : ")
	cat(rboota$nb_boot)
	cat(" \n")
	cat(" \n")
	cat("Group 0\n")
	cat(" \n")
	print(rmboot$gpe0)
	cat(" \n")
	cat("Group 1\n")
	cat(" \n")
	print(rmboot$gpe1)
	cat(" \n")
	cat("Group 1 vs group 0\n")
	cat(" \n")
	print(rmboot$gpe1vs0)
	cat(" \n")


	#-----------------------------------------------------------------------#
	#			  		UTILITY ANALYSIS					#
	#-----------------------------------------------------------------------#

	#--------- Group 0
		# Restricted means
			restricteda = NULL
			restricteda = rbind(restricteda, round(rm[1,],3))
			restricteda = rbind(restricteda,  round(rmboot$gpe0[5,],3))
			restricteda = rbind( restricteda,  round(rm[1,] - 1.96*(rmboot$gpe0[5,]),3) ) 
			restricteda = rbind( restricteda,  round(rm[1,] + 1.96*(rmboot$gpe0[5,]),3) )
			rownames(restricteda) = c("Restricted mean", "SE", "CI 95% lower", "CI 95% upper")
		
		# Restricted mean health state duration
			# Calculate durations spent in each health stat for each bootstrap sample 
				dur_boota = rboota$rboot
				dur_boota[,3] = dur_boota[,3] - dur_boota[,2]
				dur_boota[,2] = dur_boota[,2] - dur_boota[,1]
				colnames(dur_boota) = c("TOX duration", "TWiST duration", "REL duration")

			# Fill the matrix
				restricted_dura = NULL
				restricted_dura = rbind(restricted_dura,  round( c(restricteda[1,1], restricteda[1,2]-restricteda[1,1], 
						     restricteda[1,3]-restricteda[1,2]),3) )
				restricted_dura = rbind( restricted_dura, round(apply(dur_boota,2,sd),3) )
				restricted_dura = rbind( restricted_dura, round(restricted_dura[1,] - 1.96*(restricted_dura[2,]),3) )
				restricted_dura = rbind( restricted_dura, round(restricted_dura[1,] + 1.96*(restricted_dura[2,]),3) )
				rownames(restricted_dura) = c("Restricted mean", "SE", "CI 95% lower", "CI 95% upper")
				colnames(restricted_dura) = c("TOX", "TWiST", "REL")
		
		# Restricted mean quality-adjusted survival
			# Calculate quality-adjusted survival for each bootstrap sample 
				adjust_boota = NULL
				for (i in 1:nb_boot) { 
					adjust_boota = rbind(adjust_boota, utility[1]*dur_boota[i,1] + utility[2]*dur_boota[i,2] + utility[3]*dur_boota[i,3]) 
				}
				colnames(adjust_boota) = c("Quality-adjusted survival")
				
			# Fill the matrix
				qualitya = NULL
				qualitya = rbind(qualitya, round(utility[1]*restricted_dura[1,1] +  
						utility[2]*restricted_dura[1,2] + utility[3]*restricted_dura[1,3],3) )
				qualitya = rbind(qualitya, round(apply(adjust_boota,2,sd),3) )
				qualitya  = rbind( qualitya, round(qualitya[1,] - 1.96*(qualitya[2,]),3) )
				qualitya  = rbind( qualitya, round(qualitya[1,] + 1.96*(qualitya[2,]),3) )
				rownames(qualitya) = c("Restricted mean", "SE", "CI 95% lower", "CI 95% upper")
				colnames(qualitya) = c("Quality-adjusted survival")

		# Gather results for group 0 
			utility_analysisa = list(groupe = 0, restricted_mean = restricteda, restricted_mean_state_duration = restricted_dura, 
						   restricted_mean_adjusted = qualitya)



	#--------- Group 1
		# Restricted means
			restrictedb = NULL
			restrictedb = rbind(restrictedb, rm[2,])
			restrictedb = rbind(restrictedb, rmboot$gpe1[5,])
			restrictedb = rbind( restrictedb, round(rm[2,] - 1.96*(rmboot$gpe1[5,]),3) )
			restrictedb = rbind( restrictedb, round(rm[2,] + 1.96*(rmboot$gpe1[5,]),3) )
			rownames(restrictedb) =  c("Restricted mean", "SE", "CI 95% lower", "CI 95% upper")
		
		# Restricted mean health state duration
			# Calculate durations spent in each health stat for each bootstrap sample 
				dur_bootb = rbootb$rboot
				dur_bootb[,3] = dur_bootb[,3] - dur_bootb[,2]
				dur_bootb[,2] = dur_bootb[,2] - dur_bootb[,1]
				colnames(dur_bootb) =  c("TOX duration", "TWiST duration", "REL duration")

			# Fill the matrix
				restricted_durb = NULL
				restricted_durb = rbind(restricted_durb, round(c(restrictedb[1,1], restrictedb[1,2]-restrictedb[1,1], 
						     restrictedb[1,3]-restrictedb[1,2]),3) )
				restricted_durb = rbind( restricted_durb, round(apply(dur_bootb,2,sd),3) )
				restricted_durb = rbind( restricted_durb, round(restricted_durb[1,] - 1.96*(restricted_durb[2,]),3) )
				restricted_durb = rbind( restricted_durb, round(restricted_durb[1,] + 1.96*(restricted_durb[2,]),3) )
				rownames(restricted_durb) = c("Restricted mean", "SE", "CI 95% lower", "CI 95% upper")
				colnames(restricted_durb) = c("TOX", "TWiST", "REL")
		
		# Restricted mean quality-adjusted survival
			# Calculer les quality-adjusted survival pr chaque échantillon bootstrap
				adjust_bootb = NULL
				for (i in 1:nb_boot) { 
					adjust_bootb = rbind(adjust_bootb, utility[1]*dur_bootb[i,1] + utility[2]*dur_bootb[i,2] + utility[3]*dur_bootb[i,3]) 
				}
				colnames(adjust_bootb) = c("Quality-adjusted survival")
			
			# Fill the matrix
				qualityb = NULL
				qualityb = rbind(qualityb, round(utility[1]*restricted_durb[1,1] +  
						utility[2]*restricted_durb[1,2] + utility[3]*restricted_durb[1,3],3) )
				qualityb = rbind(qualityb, round(apply(adjust_bootb,2,sd),3) )
				qualityb  = rbind( qualityb, round(qualityb[1,] - 1.96*(qualityb[2,]),3) )
				qualityb  = rbind( qualityb, round(qualityb[1,] + 1.96*(qualityb[2,]),3) )
				rownames(qualityb) = c("Restricted mean", "SE", "CI 95% lower", "CI 95% upper")
				colnames(qualityb) = c("Quality-adjusted survival")

		# Gather results for group 1
			utility_analysisb = list(groupe = 1, restricted_mean = restrictedb, restricted_mean_state_duration = restricted_durb, 
						   restricted_mean_adjusted = qualityb)


	#--------- Group 1 vs 0
		# Restricted means
			restrictedc = NULL
			restrictedc = rbind(restrictedc, round(rm[2,] - rm[1,],3))
			restrictedc = rbind(restrictedc, rmboot$gpe1vs0[5,])
			restrictedc = rbind( restrictedc,  round(restrictedc[1,] - 1.96*(restrictedc[2,]),3))
			restrictedc = rbind( restrictedc,  round(restrictedc[1,] + 1.96*(restrictedc[2,]),3))
			restrictedc = rbind( restrictedc,  round(2*(1-pnorm(abs(restrictedc[1,])/restrictedc[2,])),3) )
			rownames(restrictedc) = c("Restricted mean", "SE", "CI 95% lower", "CI 95% upper", "P.value")
		
		# Restricted mean health state duration
			# Calculate mean duration spent in each health state for each bootstrap 
				dur_bootc = dur_bootb - dur_boota

			# Fill the matrix
				restricted_durc = NULL
				restricted_durc = rbind(restricted_durc, round(c(restrictedc[1,1], restrictedc[1,2]-restrictedc[1,1], 
						     restrictedc[1,3]-restrictedc[1,2]),3) )
				restricted_durc = rbind( restricted_durc, round(apply(dur_bootc,2,sd),3) )
				restricted_durc = rbind( restricted_durc, round(restricted_durc[1,] - 1.96*(restricted_durc[2,]),3) )
				restricted_durc = rbind( restricted_durc, round(restricted_durc[1,] + 1.96*(restricted_durc[2,]),3) )
				restricted_durc = rbind( restricted_durc, round(2*(1-pnorm(abs(restricted_durc[1,])/restricted_durc[2,])),3) )
				rownames(restricted_durc) = c("Mean", "SE", "CI 95% lower", "CI 95% upper", "P.value")
				colnames(restricted_durc) = c("TOX", "TWiST", "REL")

			# Restricted mean quality-adjusted survival
				adjust_bootc = adjust_bootb - adjust_boota
				
			# Fill the matrix
				qualityc = NULL
				qualityc = rbind(qualityc, round(qualityb[1,1] - qualitya[1,1], 3) )
				qualityc = rbind(qualityc, round(apply(adjust_bootc,2,sd),3) )
				qualityc  = rbind( qualityc, round(qualityc[1,] - 1.96*(qualityc[2,]),3) )
				qualityc  = rbind( qualityc, round(qualityc[1,] + 1.96*(qualityc[2,]),3) )
				qualityc  = rbind( qualityc  , round(2*(1-pnorm(abs(qualityc[1,])/qualityc[2,])),3) )
				rownames(qualityc) = c("Mean", "SE", "CI 95% lower", "CI 95% upper", "P.value")
				colnames(qualityc) = c("Quality-adjusted survival")

		# Gather results for group 1vs0
			utility_analysisc = list(groupe = "1vs0", restricted_mean = restrictedc, restricted_mean_state_duration = restricted_durc, 
						   restricted_mean_adjusted = qualityc)


	#--------- Gather results for group1, 0 and 1vs0
	utility_analysis = list (gpe0 = utility_analysisa, gpe1 = utility_analysisb, gpe1vs0 = utility_analysisc )

	#--------- Utility values used in this analysis
	score=matrix(utility, 3, 1)
	colnames(score) = c("Utility score")
	rownames(score) = c("TOX", "TWiST", "REL")

	cat("+--------------------------------------------+\n")
	cat("|              Utility analysis              |\n")
	cat("+--------------------------------------------+\n")
	cat(" \n")
	cat(paste("Time point : ", tp))
	cat(" \n")
	cat("Number of bootstrap : ")
	cat(rboota$nb_boot)
	cat(" \n")
	cat(" \n")
	cat("------------------------------------\n")
	cat("Utility values used in this analysis\n")
	cat(" \n")
	print(score)
	cat(" \n")
	cat("------------------------------------\n")
	cat("Results for group 0\n")
	cat(" \n")
	cat("Restricted means\n")
	cat(" \n")
	print(utility_analysis$gpe0$restricted_mean)
	cat(" \n")
	cat("Restricted mean health state duration\n")
	cat(" \n")
	print(utility_analysis$gpe0$restricted_mean_state_duration)
	cat(" \n")
	cat("Restricted mean quality adjusted survival\n")
	cat(" \n")
	print(utility_analysis$gpe0$restricted_mean_adjusted)
	cat(" \n")
	cat("------------------------------------\n")
	cat("Results for group 1\n")
	cat(" \n")
	cat("Restricted means\n")
	cat(" \n")
	print(utility_analysis$gpe1$restricted_mean)
	cat(" \n")
	cat("Restricted mean health state duration\n")
	cat(" \n")
	print(utility_analysis$gpe1$restricted_mean_state_duration)
	cat(" \n")
	cat("Restricted mean quality adjusted survival\n")
	cat(" \n")
	print(utility_analysis$gpe1$restricted_mean_adjusted)
	cat(" \n")
	cat("------------------------------------\n")
	cat("Results for group 1 vs group 0\n")
	cat(" \n")
	cat("Restricted means\n")
	cat(" \n")
	print(utility_analysis$gpe1vs0$restricted_mean)
	cat(" \n")
	cat("Restricted mean health state duration\n")
	cat(" \n")
	print(utility_analysis$gpe1vs0$restricted_mean_state_duration)
	cat(" \n")
	cat("Restricted mean quality adjusted survival\n")
	cat(" \n")
	print(utility_analysis$gpe1vs0$restricted_mean_adjusted)
	cat(" \n")


	#-----------------------------------------------------------------------#
	#			  		GRAPH 						#
	#-----------------------------------------------------------------------#

	#---------- Partitioned survival - group 0
		plot_toxa = NULL
		plot_toxa$time = c(0,surv_toxa$time[which(surv_toxa$time<tp)])
		plot_toxa$surv = c(1,surv_toxa$surv[which(surv_toxa$time<tp)])
		if (plot_toxa$surv[length(plot_toxa$surv)]!=0) {plot_toxa$time[length(plot_toxa$time)] = tp}
		if (plot_toxa$surv[length(plot_toxa$surv)]==1) {plot_toxa$time=cbind(plot_toxa$time, tp) 
									plot_toxa$surv=cbind(plot_toxa$surv,plot_toxa$surv)}
		plot_rfsa = NULL
		plot_rfsa$time = c(0,surv_rfsa$time[which(surv_rfsa$time<tp)])
		plot_rfsa$surv = c(1,surv_rfsa$surv[which(surv_rfsa$time<tp)])
		if (plot_rfsa$surv[length(plot_rfsa$surv)]!=0) {plot_rfsa$time[length(plot_rfsa$time)] = tp}
		if (plot_rfsa$surv[length(plot_rfsa$surv)]==1) {plot_rfsa$time=cbind(plot_rfsa$time, tp) 
									plot_rfsa$surv=cbind(plot_rfsa$surv,plot_rfsa$surv)}
		plot_osa = NULL
		plot_osa$time = c(0,surv_osa$time[which(surv_osa$time<tp)])
		plot_osa$surv = c(1,surv_osa$surv[which(surv_osa$time<tp)])
		if (plot_osa$surv[length(plot_osa$surv)]!=0) {plot_osa$time[length(plot_osa$time)] = tp}
		if (plot_osa$surv[length(plot_osa$surv)]==1) {plot_osa$time=cbind(plot_osa$time, tp) 
									plot_osa$surv=cbind(plot_osa$surv,plot_osa$surv)}

	#---------- Partitioned survival curves - group 1
		plot_toxb = NULL
		plot_toxb$time = c(0,surv_toxb$time[which(surv_toxb$time<tp)])
		plot_toxb$surv = c(1,surv_toxb$surv[which(surv_toxb$time<tp)])
		if (plot_toxb$surv[length(plot_toxb$surv)] != 0) {plot_toxb$time[length(plot_toxb$time)] = tp}
		if (plot_toxb$surv[length(plot_toxb$surv)]==1) {plot_toxb$time=cbind(plot_toxb$time, tp) 
									plot_toxb$surv=cbind(plot_toxb$surv,plot_toxb$surv)}
		plot_rfsb = NULL
		plot_rfsb$time = c(0,surv_rfsb$time[which(surv_rfsb$time<tp)])
		plot_rfsb$surv = c(1,surv_rfsb$surv[which(surv_rfsb$time<tp)])
		if (plot_rfsb$surv[length(plot_rfsb$surv)] != 0) {plot_rfsb$time[length(plot_rfsb$time)] = tp}
		if (plot_rfsb$surv[length(plot_rfsb$surv)]==1) {plot_rfsb$time=cbind(plot_rfsb$time, tp) 
									plot_rfsb$surv=cbind(plot_rfsb$surv,plot_rfsb$surv)}
		plot_osb = NULL
		plot_osb$time = c(0,surv_osb$time[which(surv_osb$time<tp)])
		plot_osb$surv = c(1,surv_osb$surv[which(surv_osb$time<tp)])
		if (plot_osb$surv[length(plot_osb$surv)] != 0) {plot_osb$time[length(plot_osb$time)] = tp}
		if (plot_osb$surv[length(plot_osb$surv)]==1) {plot_osb$time=cbind(plot_osb$time, tp) 
									plot_osb$surv=cbind(plot_osb$surv,plot_osb$surv)}


	#---------- Threshold analysis
		threshold=threshold.func(fix_utility=fix_utility, mean=utility_analysis$gpe1vs0$restricted_mean_state_duration[1,c(1:3)], dur_boot=dur_bootc, step=0.01)

		amplitude = max(threshold[,3]) - min(threshold[,3]) # amplitude of variation of mean threshold
		point = signif( max(threshold[,3])- 3*(amplitude/6) , 2 ) # value of the first line (the middle one, there is a total of 5 lines)
		ecartx = threshold[8151,3] - threshold[8141,3] # distance between means for utility_tox fixed and utility_rel reducing of 0.1
		ecarty = threshold[9161,3] - threshold[8151,3] # distance between means for utility_rel fixed and utility_tox reducing of 0.1

		# Draw the first line
			threshold2 = cbind(threshold, matrix(0,dim(threshold)[1],1))
			for (i in 2:(dim(threshold)[1]-1) ) {
			if (threshold2[i,3]<=point & (threshold2[i+1,3]>point |threshold2[i-1,3]>point) & threshold2[i+1,1]==threshold2[i-1,1]) {threshold2[i,8]=1}
			}
			# coord = coordonate of 1 points of the first line 
				num=which(threshold2[,8]==1)
				coord = matrix(0,2,3)
				colnames(coord) = c("utility_tox", "utility_prog", "mean")
					# premier point
					coord[1,1]=threshold2[num[1],1]
					coord[1,2]= threshold2[num[1],2] - (point-threshold2[num[1],3])*-0.1/ecartx
					coord[1,3]=point
					# deuxième point
					if(length(num)!=0) (coord[2,1]=threshold2[num[length(num)],1])
					if(length(num)!=0) (coord[2,2]= threshold2[num[length(num)],2] - (point-threshold2[num[length(num)],3])*-0.1/ecartx)
					coord[2,3]=point
				pente = (coord[1,1]-coord[2,1])/(coord[1,2]-coord[2,2])
				origine1 = coord[1,1]-pente*coord[1,2]

		# Draw the other lines : 2 above et 2 under the first line
			origine2 = origine1 + (((amplitude/6)*0.1)/ecarty)
			origine3 = origine1 + 2*(amplitude/6*0.1/ecarty)
			origine4 = origine1 - (amplitude/6*0.1/ecarty)
			origine5 = origine1 - 2*(amplitude/6*0.1/ecarty)

		# Coordonate of the middle line
			if (pente <=0 & is.na(pente)==F ) { # then middle = intersection line y=x
				milieu1 = c( origine1/(1-pente) , origine1/(1-pente) )
				milieu2 = c( origine2/(1-pente) , origine2/(1-pente) )
				milieu3 = c( origine3/(1-pente) , origine3/(1-pente) )
				milieu4 = c( origine4/(1-pente) , origine4/(1-pente) )
				milieu5 = c( origine5/(1-pente) , origine5/(1-pente) ) }
			if (pente>0 & is.na(pente)==F ) { # then middle = intersection with y=-x+1
				milieu1 = c( (1-origine1)/(pente+1) , (origine1-1)/(pente+1)+1 )
				milieu2 = c( (1-origine2)/(pente+1) , (origine2-1)/(pente+1)+1 )
				milieu3 = c( (1-origine3)/(pente+1) , (origine3-1)/(pente+1)+1 )
				milieu4 = c( (1-origine4)/(pente+1) , (origine4-1)/(pente+1)+1 )
				milieu5 = c( (1-origine5)/(pente+1) , (origine5-1)/(pente+1)+1 ) }

	#---------- Q-TWiST gain function
		# Matrix "gain" represents the gain function for group 1 vs group 0
			gain = matrix(0,ceiling((tp)/(tp/1000))+20,4)
			colnames(gain) = c("timepoint", "min", "0.5", "max")
			nnn=1
			for (i in seq(1,tp,by=round(tp/1000,4))) { 
				# Matrix with the time spent in TOX, TWiST and REL according to i (with i=timepoint)
					tttt=matrix(0,1,3)
					colnames(tttt) = c("TOX", "TWiST", "REL")
					rownames(tttt) = c("Group 1vs0")
					tttt[1,1] = round(summary(surv_toxb, rmean = i)$table[5],3) - round(summary(surv_toxa, rmean = i)$table[5] , 3)
					tttt[1,2] = round(summary(surv_rfsb, rmean = i)$table[5] - summary(surv_toxb, rmean = i)$table[5] , 3) - 
							round(summary(surv_rfsa, rmean = i)$table[5] - summary(surv_toxa, rmean = i)$table[5] , 3)
					tttt[1,3] = round(summary(surv_osb, rmean = i)$table[5] - summary(surv_rfsb, rmean = i)$table[5] , 3) - 
							round(summary(surv_osa, rmean = i)$table[5] - summary(surv_rfsa, rmean = i)$table[5] , 3)

				# Applicate threshold.func for determinate min and max
					aaa=threshold.func(fix_utility=fix_utility, mean=tttt, dur_boot=NULL, step=0.1)

				# Fill the matrix gain 
					gain[nnn,1] = i								# value of timepoint
					gain[nnn,2] = min(aaa[,3])						# min of QoL among all values
					if (fix_utility==1) {gain[nnn,3] = aaa[which(aaa[,1]==utility[2] & aaa[,2]==utility[3]) , 3]}	# values of Qol when 2 other scores = the value in utility
					if (fix_utility==2) {gain[nnn,3] = aaa[which(aaa[,1]==utility[1] & aaa[,2]==utility[3]) , 3]}
					if (fix_utility==3) {gain[nnn,3] = aaa[which(aaa[,1]==utility[1] & aaa[,2]==utility[2]) , 3]}
					gain[nnn,4] = max(aaa[,3])						# max of QoL among all values
				nnn=nnn+1
			}
	gain = gain[-which(gain[,1]==0),]


	#-----------------------------------------------------------------------#
	#			  		EXPORT OUTPUTS					#
	#-----------------------------------------------------------------------#

	   pdf("Q-TWiST analysis.pdf")

		#------- Flyleaf
			garde = data.frame(c(gm, tp, nb_boot))
			colnames(garde)=c(" ")
			rownames(garde)=c("Minimal grade of toxicity taken into acount : ", "Time point : ", "Number of bootstrap for the estimates : ")

			entete = data.frame(c("E. Bogart & A. Kramar", " ", "Centre Oscar Lambret", "3 rue Frédéric Combemale", "59000 LILLE"))
			colnames(entete)=c(" ")
			rownames(entete)=c(" ", "  ", "   " ,"    ", "")

			par(mfrow=c(3,1))
			textplot("RESULTS OF THE Q-TWiST ANALYSIS", cex=2, font=2, valign="bottom")
			textplot(garde, valign="top", cex=1.3)
			textplot(entete, valign="bottom", cex=1, halign="right")

		#------- Event summary
			par(mfrow=c(4,2))
			textplot("EVENTS SUMMARY", cex=1.2, font=2)
			textplot(" ")
			textplot("Toxicity", cex=1, font=2, valign="top", halign="right")
			textplot(summary_tox, cex=1, halign="left")
			textplot("Relapse", cex=1, font=2, valign="top", halign="right")
			textplot(summary_rfs, cex=1, halign="left")
			textplot("Overall survival", cex=1, font=2, valign="top", halign="right")
			textplot(summary_os, cex=1, halign="left")

		#------- Restricted mean estimate
			number_boot = paste("Number of bootstrap : ", nb_boot)

			par(mfrow=c(5,2))
			textplot("RESTRICTED MEAN ESTIMATES", cex=1.2, font=2)
			textplot(rm, cex=1, , halign="left", valign="bottom")
			textplot("RESTRICTED MEAN COVARIANCE ESTIMATES", cex=1.2, font=2)
			textplot(" ")
			textplot("Group 0", cex=1, valign="top", halign="right")
			textplot(rmboot$gpe0, cex=1, valign="top", halign="left")
			textplot("Group 1", cex=1, valign="top", halign="right")
			textplot(rmboot$gpe1, cex=1, valign="top", halign="left")
			textplot("Group 1 vs Group 0", cex=1, valign="top", halign="right")
			textplot(rmboot$gpe1vs0, cex=1, valign="top", halign="left")


		#------- Utility analysis
			par(mfrow=c(5,2))
			textplot("UTILITY ANALYSIS", cex=1.5, font=2)
			textplot(score, cex=1, halign="center")
			textplot("Results for group 0", cex=1.3, font=2)
			textplot(" ")
			textplot("Restricted means", cex=1, valign="top", halign="right", font=3)
			textplot(utility_analysis$gpe0$restricted_mean, cex=1, halign="left")
			textplot("Restricted mean health state duration", cex=1, valign="top", halign="right", font=3)
			textplot(utility_analysis$gpe0$restricted_mean_state_duration, cex=1, halign="left")
			textplot("Restricted mean quality adjusted survival", cex=1, valign="top", halign="right", font=3)
			textplot(utility_analysis$gpe0$restricted_mean_adjusted, cex=1, halign="left")

			par(mfrow=c(4,2))
			textplot("Results for group 1", cex=1.3, font=2)
			textplot(" ")
			textplot("Restricted means", cex=1, valign="top", halign="right", font=3)
			textplot(utility_analysis$gpe1$restricted_mean, cex=1, halign="left")
			textplot("Restricted mean health state duration", cex=1, valign="top", halign="right", font=3)
			textplot(utility_analysis$gpe1$restricted_mean_state_duration, cex=1, halign="left")
			textplot("Restricted mean quality adjusted survival", cex=1, valign="top", halign="right", font=3)
			textplot(utility_analysis$gpe1$restricted_mean_adjusted, cex=1, halign="left")

			par(mfrow=c(4,2))
			textplot("Results for group 1 vs group 0", cex=1.3, font=2)
			textplot(" ")
			textplot("Restricted means", cex=1, valign="top", halign="right", font=3)
			textplot(utility_analysis$gpe1vs0$restricted_mean, cex=1, halign="left")
			textplot("Restricted mean health state duration", cex=1, valign="top", halign="right", font=3)
			textplot(utility_analysis$gpe1vs0$restricted_mean_state_duration, cex=1, halign="left")
			textplot("Restricted mean quality adjusted survival", cex=1, valign="top", halign="right", font=3)
			textplot(utility_analysis$gpe1vs0$restricted_mean_adjusted, cex=1, halign="left")

			par(mfrow=c(1,1))

		#------- Partitioned survival curve - group 0 
			par(xpd=TRUE, mar=c(4.5,4.5,5,3))
			plot(plot_osa$time, plot_osa$surv, type="s", col = "black", xlab = "Time", 
				ylab = "Survival", main = "Group 0", ylim=c(0,1), xlim=c(0,tp))
				for ( i in 1:(length(plot_osa$time)-1) ) {
				rect(0, 0, plot_osa$time[i+1], plot_osa$surv[i], col="#CD4F39", border="#CD4F39")
				}
				for ( i in 1:(length(plot_rfsa$time)-1) ) {
				rect(0, 0, plot_rfsa$time[i+1], plot_rfsa$surv[i], col="#87CEFA", border="#87CEFA")
				}
				for ( i in 1:(length(plot_toxa$time)-1) ) {
				rect(0, 0, plot_toxa$time[i+1], plot_toxa$surv[i], col="#FFD700", border="#FFD700")
				}
				lines(plot_osa$time, plot_osa$surv, col = "black", type="s")
				lines(plot_rfsa$time, plot_rfsa$surv, col = "black", type="s")
				lines(plot_toxa$time, plot_toxa$surv, col = "black", type="s")
				rect(0, 0, 0, 1, border="black")
				rect(0, 0, tp, 0, border="black")
				rect(tp, 0, tp, max(plot_osa$surv[length(plot_osa$surv)],plot_rfsa$surv[length(plot_rfsa$surv)], 
					plot_toxa$surv[length(plot_toxa$surv)]) , border="black")
				legend("topright", inset=c(0,-0.17), legend=c("TOX", "TWiST", "REL"),  pch=15, col=c("#FFD700", "#87CEFA", "#CD4F39"))
				mtext(c("Partitioned survival curves"), 3, line=0.5)

		#------- Partitioned survival curve - group 1 
			par(xpd=TRUE, mar=c(4.5,4.5,5,3))
			plot(plot_osb$time, plot_osb$surv, type="s", col = "black", xlab = "Time", 
				ylab = "Survival", main = "Group 1", ylim=c(0,1), xlim=c(0,tp))
				for ( i in 1:(length(plot_osb$time)-1) ) {
				rect(0, 0, plot_osb$time[i+1], plot_osb$surv[i], col="#CD4F39", border="#CD4F39")
				}
				for ( i in 1:(length(plot_rfsb$time)-1) ) {
				rect(0, 0, plot_rfsb$time[i+1], plot_rfsb$surv[i], col="#87CEFA", border="#87CEFA")
				}
				for ( i in 1:(length(plot_toxb$time)-1) ) {
				rect(0, 0, plot_toxb$time[i+1], plot_toxb$surv[i], col="#FFD700", border="#FFD700")
				}
				lines(plot_osb$time, plot_osb$surv, col = "black", type="s")
				lines(plot_rfsb$time, plot_rfsb$surv, col = "black", type="s")
				lines(plot_toxb$time, plot_toxb$surv, col = "black", type="s")

				rect(0, 0, 0, 1, border="black")
				rect(0, 0, tp, 0, border="black")
				rect(tp, 0, tp, max(plot_osb$surv[length(plot_osb$surv)],plot_rfsb$surv[length(plot_rfsb$surv)], 
					plot_toxb$surv[length(plot_toxb$surv)]) , border="black")
				legend("topright", inset=c(0,-0.17), legend=c("TOX", "TWiST", "REL"),  pch=15, col=c("#FFD700", "#87CEFA", "#CD4F39"))
				mtext(c("Partitioned survival curves"), 3, line=0.5)

		#------- Threshold utility

			#---------- If slope abnormal (that is  0 or -Inf or +Inf or Na) : no possible to have a graph
			  if (pente==0 | pente==-Inf | pente==+Inf | is.na(pente)==T ) {
				par(mfrow=c(3,1))
				textplot("Threshold utility analysis", font=2, cex=3.5, valign="top", halign="center")
				zzz = c("Not enough variability according to the utility score value",
						"in the Q-TWiST difference between group 1 vs group 0", 
						"to refresh a threshold analysis plot")
				textplot(zzz, font=2, cex=1.6, valign="top", halign="center")
				textplot(" ", font=2, cex=1.8, valign="top", halign="center")
			  }
				
			#----------If slope normal (that is neither 0 nor -Inf nor +Inf nor Na) : possible to have a graph
			  if (pente!=0 & pente!=-Inf & pente!=+Inf & is.na(pente)!=T ) {
				par(mfrow=c(1,1))
				par(xpd=TRUE, mar=c(7,4.5,7,3))
				if (fix_utility == 1) {plot(threshold[,2]~threshold[,1], xlab="Utility score for relapse", ylab="Utility score for TWiST", 
					type="n", xaxp=c(0,1,10), yaxp=c(0,1,10))}
				if (fix_utility == 2) {plot(threshold[,2]~threshold[,1], xlab="Utility score for relapse", ylab="Utility score for toxicity", 
					type="n", xaxp=c(0,1,10), yaxp=c(0,1,10))}
				if (fix_utility == 3) {plot(threshold[,2]~threshold[,1], xlab="Utility score for TWiST", ylab="Utility score for toxicity", 
					type="n", xaxp=c(0,1,10), yaxp=c(0,1,10))}
		
				text("Threshold utility analysis", x=0.5, y=1.32, font=2, cex=1.8)

				if (fix_utility == 1) { text(paste("Utility score for toxicity =" , utility[1]), x=0, y=-0.27, font=3, cex=0.9) }
				if (fix_utility == 2) { text(paste("Utility score for TWiST =" , utility[2]), x=0, y=-0.27, font=3, cex=0.9) }
				if (fix_utility == 3) { text(paste("Utility score for relapse =" , utility[3]), x=0, y=-0.27, font=3, cex=0.9) }

				text(paste("Time point =" , tp), x=0, y=-0.32, font=3, cex=0.9)
				legend("topright", inset=c(0.2,-0.21), legend=c("Treatment for group 1 significantly better",
					"Treatment for group 0 significantly better", "Non significant"),  pch=15, col=c("255", "256", "honeydew2"), cex=1)

					clip(0,1,0,1) # nothing outisde x=(0,1) and y=(0,1)
					# Color according to the p value
						for (i in 1:(dim(threshold)[1]-1)) {
							if (threshold[i,7]<0.05 & threshold[i,3]>=0 & is.na(threshold[i,7])!="TRUE") 
							{ rect(threshold[i,2], threshold[i,1], threshold[i,2]+0.01, threshold[i,1], border="255", col="255") }
							if (threshold[i,7]<0.05 & threshold[i,3]<0 & is.na(threshold[i,7])!="TRUE") 
							{ rect(threshold[i,2], threshold[i,1], threshold[i,2]+0.01, threshold[i,1], border="256", col="256") }
							if (threshold[i,7]>=0.05 & is.na(threshold[i,7])!="TRUE") 
							{ rect(threshold[i,2], threshold[i,1], threshold[i,2]+0.01, threshold[i,1], border="honeydew2", col="honeydew2") }
						}
					# Value of the difference between the 2 groups
						text(point, x=milieu1[1], y=milieu1[2], col="#CD4F39", pos=3, lwd=3)
						text(signif(point+amplitude/6 , 2), x=milieu2[1], y=milieu2[2], col="#CD4F39", pos=3, lwd=3)
						text(signif(point+2*amplitude/6 , 2), x=milieu3[1], y=milieu3[2], col="#CD4F39", pos=3, lwd=3)
						text(signif(point-amplitude/6 , 2), x=milieu4[1], y=milieu4[2], col="#CD4F39", pos=3, lwd=3)
						text(signif(point-2*amplitude/6 , 2), x=milieu5[1], y=milieu5[2], col="#CD4F39", pos=3, lwd=3)
					# Line for the difference between the 2 groups 
						abline(a=origine1 , b=pente, col="#CD4F39", lwd=3)
						abline(a=origine2 , b=pente, col="#CD4F39", lwd=3)
						abline(a=origine3 , b=pente, col="#CD4F39", lwd=3)
						abline(a=origine4 , b=pente, col="#CD4F39", lwd=3)
						abline(a=origine5 , b=pente, col="#CD4F39", lwd=3)
					rect(0,0,1,1, border="black")

			  }


		#------- Q-TWiST gain function
			par(mfrow=c(1,1))
			par(xpd=TRUE, mar=c(4.5,4.5,6,3))
			plot(gain[,1], gain[,3], type="s", xlab = "Time from study entry", ylab = "Time gained", 
				ylim=c(min(gain[,2]),max(gain[,4])), xlim=c(0,tp))

			title(main=paste("Q-TWiST gain function","\n",sep=""),cex.main=1.5) 
			title(main=paste("\n","for group 1 compared with group 0",sep=""),cex.main=1)

			if (fix_utility == 1) { mtext(paste("Utility score for toxicity =" , utility[1]), 3, line=1, cex=0.9, font=3) }
			if (fix_utility == 2) { mtext(paste("Utility score for TWiST =" , utility[2]), 3, line=1, cex=0.9, font=3) }
			if (fix_utility == 3) { mtext(paste("Utility score for relapse =" , utility[3]), 3, line=1, cex=0.9, font=3) }
			
			if (fix_utility == 1) { mtext(paste("The red middle line corresponds to utility score for TWiST = " , utility[2],  
							" and for relapse = ", utility[3]), 3, line=0.3, cex=0.8, font=3) }
			if (fix_utility == 2) { mtext(paste("The red middle line corresponds to utility score for toxicity = " , utility[1],  
							" and for relapse = ", utility[3]), 3, line=0, cex=0.8, font=3) }
			if (fix_utility == 3) { mtext(paste("The red middle line corresponds to utility score for toxicity = " , utility[1],  
							" and for TWiST = ", utility[2]), 3, line=0, cex=0.8, font=3) }


			for ( i in 1:(dim(gain)[1]-1) ) {
			rect(gain[i,1], gain[i,2], gain[i+1,1], gain[i,4], col="darksalmon", border="darksalmon")
			}
			lines(c(0,tp), c(0,0), lty="dashed")
			lines(gain[,1], gain[,3], type="s", col = "lightsalmon4", lwd=2)
			lines(gain[,1], gain[,2], type="s", col = "lightsalmon4")
			lines(gain[,1], gain[,4], type="s", col = "lightsalmon4")



	   fermer=dev.off()

	cat("+------------------------------------+\n")
	cat("| Send feed back via e-mail:         |\n")
	cat("|     - e-bogart@o-lambret.fr        |\n")
	cat("|     - a-kramar@o-lambret.fr        |\n")
	cat("+------------------------------------+\n")



	return(list( survival=list(tox_gpe0=summary(surv_toxa), tox_gpe1=summary(surv_toxb), rfs_gpe0=summary(surv_rfsa), rfs_gpe1=summary(surv_rfsb),
					   os_gpe0=summary(surv_osa), os_gpe1=summary(surv_osb)),
			 threshold=threshold ))
	

}


