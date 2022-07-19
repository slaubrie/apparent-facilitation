# this code is written to run a model testing for the effect of shared pollinators on species coexistence in flowering plant communities where the storage effect operates. 
# date last modified 15 January 2020
## code for ch 3, storage effect and pollinator support
# treatments will be: (1) two species, with variable correlation + covariance, (2) one species
# need to run for all germinations

setwd('~/Dropbox/Dissertation/CH3_AmNat_StorageEffectandFacilitation/Simulations/Data/')

#remember to set seed before analyzing
set.seed(23)

# load parameter sets 
paramsets<-read.csv('paramsets.csv', header=T) #this is currently parameterized for high competition (all=1)

# load germination rates 
# there are 5 mean germination rates: 0.02 [x], 0.16 [x], 0.5 [x], 0.84 [x], 0.98 [] and single species
germ1<-readRDS("germ1_mu_0.98")  
germ2<-readRDS("germ2_mu_0.98")  

# number of germination correlation levels 
l<-nrow(germ1)

# number of parameters 
k<-nrow(paramsets)

#number of timesteps per timeseries
ts<-5000

#libraries
require(MASS)
require(ggplot2)
require(reshape2)
require(Rmisc)
require(boot)
require(tidyr)
require(dplyr)
require(plyr)
require(progress)



##########################################
##########################################
################# TIMESERIES  ############
##########################################

#want to store information for six state variables from the k timeseries

#initializing most state variable values with a number requires adding a column to fill in, which is reflected in the "ts+1".

#three species 
B<-array(data=0, c(l, k, (ts+1)))
E<-array(data=0, c(l, k, (ts+1)))
P1<-array(data=0, c(l, k, (ts+1)))
P2<-array(data=0, c(l, k, (ts+1)))
S1<-array(data=0, c(l, k, (ts+1)))
S2<-array(data=0, c(l, k, (ts+1)))
O1<-array(data=0, c(l, k, (ts)))
O2<-array(data=0, c(l, k, (ts)))
F1<-array(data=0, c(l, k, (ts)))
F2<-array(data=0, c(l, k, (ts)))
CompS1<-array(data=0, c(l,k,(ts)))
CompS2<-array(data=0, c(l,k,(ts)))
x<-array(data=0, c(l, k, (ts)))
#initialize 

# Variables: Bees, Eggs, Plants, Seeds
B[,,1]<-10
E[,,1]<-10
P1[,,1]<-10
P2[,,1]<-10
S1[,,1]<-10
S2[,,1]<-10 


###############################################################
######### model for all three species - p1, p2, B #############
###############################################################
pb1<-txtProgressBar(0,l, style=3)
for (l in 1:l){
	for (j in 1:k){
		for(i in 1:ts){


#number of plants 
# P_{i}(t) = g_{i}S_{i}(t) 
P1[l,j,i+1]<-germ1[l,j,i]*S1[l,j,i]

P2[l,j,i+1]<-germ2[l,j,i]*S2[l,j,i]

#number of bees, must be greater than two to continue to reproduce
if(B[l,j,i]<1.0e+10){
	if(B[l,j,i]>2){
	B[l,j,i+1]<-paramsets$p[j]*E[l,j,i] 
	} else {
	B[l,j,i+1]=0
	}
	} else{
	B[l,j,i]=1.0e+10
	}
	

# number of ovules produced in competition 
# O_{i}(t)=\frac{\lambda_{i}P_{i}(t)}{1+\sum_{j}{\alpha_{ij}P_{j}}}
O1[l,j,i]<-(paramsets$l1[j]*P1[l,j,i])/(1+(paramsets$a_ii[j]*P1[l,j,i])+(paramsets$a_ij[j]*P2[l,j,i]))

O2[l,j,i]<-(paramsets$l2[j]*P2[l,j,i])/(1+(paramsets$a_jj[j]*P2[l,j,i])+(paramsets$a_ji[j]*P1[l,j,i]))

#ratio of bees to ovules 
x[l,j,i]<-B[l,j,i]/(O1[l,j,i]+O2[l,j,i])

#eggs produced
# E(t+1) = \lambda_{B}B(t)(1-e^{-\frac{M}{x[l,j,i]}}) 
E[l,j,i+1]<-paramsets$lambdab[j]*B[l,j,i]*(1-exp(-paramsets$M[j]/x[l,j,i]))

#seeds produced 
# S_{i}(t+1) = O_{i}(t)(1-e^{-\frac{x[l,j,i]}{L}})

# keep interdependence 
S1[l,j,i+1]<-(O1[l,j,i]*(1-exp(-x[l,j,i]/paramsets$L[j])))+(S1[l,j,i]*(paramsets$s1[j]*(1-germ1[l,j,i])))

S2[l,j,i+1]<-(O2[l,j,i]*(1-exp(-x[l,j,i]/paramsets$L[j])))+(S2[l,j,i]*(paramsets$s2[j]*(1-germ2[l,j,i]))) 

# get rid of interdependence
#S1[l,j,i+1]<-(O1[l,j,i]*(0.5))+(S1[l,j,i]*(paramsets$s1[j]*(1-germ1[l,j,i])))

#S2[l,j,i+1]<-(O2[l,j,i]*(0.5))+(S2[l,j,i]*(paramsets$s2[j]*(1-germ2[l,j,i]))) 


#seeds produced by pollination only
F1[l,j,i]<-(O1[l,j,i]*(1-exp(-x[l,j,i]/paramsets$L[j])))
F2[l,j,i]<-(O2[l,j,i]*(1-exp(-x[l,j,i]/paramsets$L[j])))

## competition experienced 
CompS1[l,j,i]<-((paramsets$a_ii[j]*P1[l,j,i])+(paramsets$a_ij[j]*P2[l,j,i]))
CompS2[l,j,i]<-((paramsets$a_jj[j]*P2[l,j,i])+(paramsets$a_ji[j]*P1[l,j,i]))


		}
	}
setTxtProgressBar(pb1,l)}


#############################################
#############################################
########## working the data  ################
#############################################
#############################################
# the size of melted data for state variables B, E, P1, P2, S1, S2 should be nrow = l * k * ts = 3 * 100 * 5001 (!!!) = 150300; the rest should be 150300
####  plant species 1 
# melt data to give variable 1, which is parameter set, and variable 2, which is time step. 

P1melt1<-melt(P1[1,,]) # negative correlation
P1melt2<-melt(P1[2,,])
P1melt3<-melt(P1[3,,]) # zero 
P1melt4<-melt(P1[4,,])
P1melt5<-melt(P1[5,,]) # positive correlation 

# rename the columns and rows. "ps" stands for "parameter set" and "ts" stands for "time step". I will be using this set up to collapse variation introduced by different parameter sets (taken as replicates) into mean and sd of all parameter sets for each time step.

colnames(P1melt1)<-c("ps", "ts", "nplants")
P1melt1$cor<-rep("lo", nrow(P1melt1))
P1melt1$sp<-rep("sp1", nrow(P1melt1))

colnames(P1melt2)<-c("ps", "ts", "nplants")
P1melt2$cor<-rep("medlo", nrow(P1melt2))
P1melt2$sp<-rep("sp1", nrow(P1melt2))

colnames(P1melt3)<-c("ps", "ts", "nplants")
P1melt3$cor<-rep("zero", nrow(P1melt3))
P1melt3$sp<-rep("sp1", nrow(P1melt3))

colnames(P1melt4)<-c("ps", "ts", "nplants")
P1melt4$cor<-rep("medhi", nrow(P1melt4))
P1melt4$sp<-rep("sp1", nrow(P1melt4))

colnames(P1melt5)<-c("ps", "ts", "nplants")
P1melt5$cor<-rep("hi", nrow(P1melt5))
P1melt5$sp<-rep("sp1", nrow(P1melt5))

### plant species 2 
# melt data to give variable 1, which is parameter set, and variable 2, which is time step. 

P2melt1<-melt(P2[1,,]) # negative correlation
P2melt2<-melt(P2[2,,])
P2melt3<-melt(P2[3,,]) # zero 
P2melt4<-melt(P2[4,,])
P2melt5<-melt(P2[5,,]) # positive correlation 

# rename the columns and rows. "ps" stands for "parameter set" and "ts" stands for "time step". I will be using this set up to collapse variation introduced by different parameter sets (taken as replicates) into mean and sd of all parameter sets for each time step. 

colnames(P2melt1)<-c("ps", "ts", "nplants")
P2melt1$cor<-rep("lo", nrow(P2melt1))
P2melt1$sp<-rep("sp2", nrow(P2melt1))

colnames(P2melt2)<-c("ps", "ts", "nplants")
P2melt2$cor<-rep("medlo", nrow(P2melt2))
P2melt2$sp<-rep("sp2", nrow(P2melt2))

colnames(P2melt3)<-c("ps", "ts", "nplants")
P2melt3$cor<-rep("zero", nrow(P2melt3))
P2melt3$sp<-rep("sp2", nrow(P2melt3))

colnames(P2melt4)<-c("ps", "ts", "nplants")
P2melt4$cor<-rep("medhi", nrow(P2melt4))
P2melt4$sp<-rep("sp2", nrow(P2melt4))

colnames(P2melt5)<-c("ps", "ts", "nplants")
P2melt5$cor<-rep("hi", nrow(P2melt5))
P2melt5$sp<-rep("sp2", nrow(P2melt5))

# keeps plant species separate to look at if both survive 
#try2<-rbind(P2melt1, P2melt2,P2melt3, P1melt1,P1melt2,P1melt3)
#try.sumry<-as.data.frame(try2 %>% group_by(ps, cov, sp) %>% summarize(meanplants=mean(nplants, na.rm=T)))

# add together plants in diff covariance scenarios, plant species 1 + plant species 2 
P1melt1$np<-P1melt1$nplants+P2melt1$nplants
P1melt2$np<-P1melt2$nplants+P2melt2$nplants
P1melt3$np<-P1melt3$nplants+P2melt3$nplants
P1melt4$np<-P1melt4$nplants+P2melt4$nplants
P1melt5$np<-P1melt5$nplants+P2melt5$nplants

# bind together and plot
PTS<-rbind(P1melt1, P1melt2, P1melt3, P1melt4, P1melt5)

#### bees 
# melt data to give variable 1, which is parameter set, and variable 2, which is time step.
 
B1<-melt(B[1,,]) #low (-1) correlation
B2<-melt(B[2,,]) 
B3<-melt(B[3,,]) #medium (0) correlation 
B4<-melt(B[4,,])
B5<-melt(B[5,,]) # high (1) correlation 

#rename the columns and rows of  "ps" stands for "parameter set" and "ts" stands for "time step". I will be using this set up to collapse variation introduced by different parameter sets (taken as replicates) into mean and sd of all parameter sets for each time step.

colnames(B1)<-c("ps", "ts", "nbees")
colnames(B2)<-c("ps", "ts", "nbees")
colnames(B3)<-c("ps", "ts", "nbees")
colnames(B4)<-c("ps", "ts", "nbees")
colnames(B5)<-c("ps", "ts", "nbees")

# bees 
B1$cor<-rep("lo", nrow(B1))
B2$cor<-rep("medlo", nrow(B2))
B3$cor<-rep("zero", nrow(B3))
B4$cor<-rep("medhi", nrow(B4))
B5$cor<-rep("hi", nrow(B5))

#bind for comparison
BTS<-rbind(B1, B2, B3, B4, B5)

#############################################
#############################################
###### summarizing the data #################
#############################################
#############################################

# goal: to average the count data from the last half of the timesteps of each parameter set. this can be used as a set of response variables that are associated with parameter sets. each response variable can be regressed against all parameters in the model to ask how variation in particular parameters account for variation in the state variables. 

### bee summaries: mean 
aa<-BTS[which(BTS$ts>(3999)),]
aa<-aa[which(aa$ts!=ts+1),]
aa$nbees[which(aa$nbees<1)]<-0 #replace averages less than 1 with 0
ba<-as.data.frame(aa %>% group_by(ps, cor) %>% summarize(mean_bees=mean(nbees, na.rm=T)))
ca<-dcast(ba, ps~cor, value.var="mean_bees") # for the same parameter set at different correlations 
names(ca)<-c("ps", "hicor", "locor", "medhicor", "medlocor", "zerocor")
ca$meas<-rep("mean", nrow(ca))

### bee summaries: variance
b.a<-as.data.frame(aa %>% group_by(ps, cor) %>% summarize(var_bees=var(nbees, na.rm=T)))
c.a<-dcast(b.a, ps~cor, value.var="var_bees")
names(c.a)<-c("ps", "hicor", "locor", "medhicor", "medlocor", "zerocor")
c.a$meas<-rep("var", nrow(c.a))

# smoosh together mean and vars 
paramsets.bmean<-cbind(paramsets, ca[,c(2:7)])
paramsets.bvar<-cbind(paramsets, c.a[,c(2:7)])

# alltogether
paramsets.b<-rbind(paramsets.bmean,paramsets.bvar)

# melt
da<-melt(paramsets.b, id.vars=c('nps', 'p','M','s1','s2','l1','l2', 'a_ii', 'a_jj', 'a_ji', 'a_ij', "lambdab", "L", "meas"), measure.vars=c("locor", "medlocor", "zerocor", "medhicor", "hicor"))
ea<-plyr::rename(da, replace=c("variable"="cor","value"="bees"))

#### plant summaries: means
ab<-PTS[which(PTS$ts>(3999)),]
ab<-ab[which(ab$ts!=ts+1),]
ab$np[which(is.nan(ab$np))]<-0 #replace 0/0 = NaN with 0
ab$np[which(ab$np<1)]<-0 #replace numbers less than 1 with 0. There are many examples of this happening
bb<-as.data.frame(ab %>% group_by(ps, cor) %>% summarize(mean_plants=mean(np, na.rm=T)))
cb<-dcast(bb, ps~cor, value.var="mean_plants")
names(cb)<-c("ps", "hicor", "locor", "medhicor", "medlocor", "zerocor")
cb$meas<-rep("mean", nrow(cb))

### plant summaries: variance
b.b<-as.data.frame(ab %>% group_by(ps, cor) %>% summarize(var_plants=var(nplants, na.rm=T)))
c.b<-dcast(b.b, ps~cor, value.var="var_plants")
names(c.b)<-c("ps", "hicor", "locor", "medhicor", "medlocor", "zerocor")
c.b$meas<-rep("var", nrow(c.b))

# smoosh together mean and vars
paramsets.pmean<-cbind(paramsets, cb[,c(2:7)])
paramsets.pvar<-cbind(paramsets, c.b[,c(2:7)])

paramsets.c<-rbind(paramsets.pmean,paramsets.pvar)

# melt 
db<-melt(paramsets.c, id.vars=c('nps', 'p','M','s1','s2','l1','l2', 'a_ii', 'a_jj', 'a_ji', 'a_ij', "lambdab", "L", "meas"), measure.vars=c("locor", "medlocor", "zerocor", "medhicor", "hicor"))

## put the plants and bees together 
dat1<-as.data.frame(cbind(ea, db$value))
dat<-plyr::rename(dat1, replace=c("db$value"="plants"))


##############################################################################################################################
##############################################################################################################################
########################################## end of two plant species code #####################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
#################################### beginning of one plant species code #####################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################


##########################################
################# TIMESERIES  ############
##########################################

#want to store information for six state variables from the k timeseries

#initializing most state variable values with a number requires adding a column to fill in, which is reflected in the "ts+1".

#three species 
B.1<-array(data=0, c(l, k, (ts+1)))
E.1<-array(data=0, c(l, k, (ts+1)))
P1.1<-array(data=0, c(l, k, (ts+1)))
S1.1<-array(data=0, c(l, k, (ts+1)))
O1.1<-array(data=0, c(l, k, (ts)))
F1.1<-array(data=0, c(l, k, (ts)))
CompS1.1<-array(data=0, c(l,k,(ts)))
x.1<-array(data=0, c(l, k, (ts)))
#initialize 

# Variables: B.1ees, Eggs, Plants, Seeds
B.1[,,1]<-10
E.1[,,1]<-10
P1.1[,,1]<-10
S1.1[,,1]<-10


###############################################################
######### model for two species - P1.1, B.1 ###################
###############################################################
pb2<-txtProgressBar(0,l, style=3)
for (l in 1:l){
	for (j in 1:k){
		for(i in 1:ts){


#number of plants 
# P_{i}(t) = g_{i}S_{i}(t) 
P1.1[l,j,i+1]<-germ1[l,j,i]*S1.1[l,j,i]

#number of bees, must be greater than two to continue to reproduce

if(B.1[l,j,i]<1.0e+10){
	if(B.1[l,j,i]>2){
	B.1[l,j,i+1]<-paramsets$p[j]*E.1[l,j,i] 
	} else {
	B.1[l,j,i+1]=0
	}
	} else{
	B.1[l,j,i]=1.0e+10
	}
	
# number of ovules produced in competition 
# O_{i}(t)=\frac{\lambda_{i}P_{i}(t)}{1+\sum_{j}{\alpha_{ij}P_{j}}}
O1.1[l,j,i]<-(paramsets$l1[j]*P1.1[l,j,i])/(1+(paramsets$a_ii[j]*P1.1[l,j,i]))

# ratio of bees to ovules
x.1[l,j,i]<-B.1[l,j,i]/O1.1[l,j,i]

#eggs produced
# E(t+1) = \lambda_{B.1}B.1(t)(1-e^{-\frac{M}{x[l,j,i]}}) 
E.1[l,j,i+1]<-paramsets$lambdab[j]*B.1[l,j,i]*(1-exp(-paramsets$M[j]/x.1[l,j,i]))

#seeds produced 
# S_{i}(t+1) = O_{i}(t)(1-e^{-\frac{x[l,j,i]}{L}})

# keep interdependence 
S1.1[l,j,i+1]<-(O1.1[l,j,i]*(1-exp(-x.1[l,j,i]/paramsets$L[j])))+(S1.1[l,j,i]*(paramsets$s1[j]*(1-germ1[l,j,i])))

# get rid of interdependence 
#S1.1[l,j,i+1]<-(O1[l,j,i]*(0.5))+(S1.1[l,j,i]*(paramsets$s1[j]*(1-germ1[l,j,i])))

#seeds produced by pollination only
F1.1[l,j,i]<-(O1.1[l,j,i]*(1-exp(-x.1[l,j,i]/paramsets$L[j])))

## competition experienced 
CompS1.1[l,j,i]<-((paramsets$a_ii[j]*P1.1[l,j,i]))

		}
	}
setTxtProgressBar(pb2,l)}


#############################################
#############################################
########## working the data  ################
#############################################
#############################################

####  plant species 1 
# melt data to give variable 1, which is parameter set, and variable 2, which is time step. 

P1.1melt1<-melt(P1.1[1,,]) # negative covariance
P1.1melt2<-melt(P1.1[2,,])
P1.1melt3<-melt(P1.1[3,,]) 
P1.1melt4<-melt(P1.1[4,,])
P1.1melt5<-melt(P1.1[5,,]) # positive covariance 

# rename the columns and rows. "ps" stands for "parameter set" and "ts" stands for "time step". I will be using this set up to collapse variation introduced by different parameter sets (taken as replicates) into mean and sd of all parameter sets for each time step.

colnames(P1.1melt1)<-c("ps", "ts", "nplants")
P1.1melt1$cor<-rep("zerocor", nrow(P1.1melt1))
P1.1melt1$sp<-rep("sp", nrow(P1.1melt1))

colnames(P1.1melt2)<-c("ps", "ts", "nplants")
P1.1melt2$cor<-rep("medlocor", nrow(P1.1melt2))
P1.1melt2$sp<-rep("sp", nrow(P1.1melt2))

colnames(P1.1melt3)<-c("ps", "ts", "nplants")
P1.1melt3$cor<-rep("locor", nrow(P1.1melt3))
P1.1melt3$sp<-rep("sp", nrow(P1.1melt3))

colnames(P1.1melt4)<-c("ps", "ts", "nplants")
P1.1melt4$cor<-rep("medhicor", nrow(P1.1melt4))
P1.1melt4$sp<-rep("sp", nrow(P1.1melt4))

colnames(P1.1melt5)<-c("ps", "ts", "nplants")
P1.1melt5$cor<-rep("hicor", nrow(P1.1melt5))
P1.1melt5$sp<-rep("sp", nrow(P1.1melt5))

# add together plants in diff covariance scenarios 
P1.1melt1$np<-P1.1melt1$nplants
P1.1melt2$np<-P1.1melt2$nplants
P1.1melt3$np<-P1.1melt3$nplants
P1.1melt4$np<-P1.1melt4$nplants
P1.1melt5$np<-P1.1melt5$nplants

# bind together and plot
PTS1<-rbind(P1.1melt1, P1.1melt2, P1.1melt3, P1.1melt4, P1.1melt5)

#### bees 
# melt data to give variable 1, which is parameter set, and variable 2, which is time step.
 
B.11<-melt(B.1[1,,]) #low (-1) covariance
B.12<-melt(B.1[2,,]) 
B.13<-melt(B.1[3,,]) 
B.14<-melt(B.1[4,,]) 
B.15<-melt(B.1[5,,]) 

#rename the columns and rows of  "ps" stands for "parameter set" and "ts" stands for "time step". I will be using this set up to collapse variation introduced by different parameter sets (taken as replicates) into mean and sd of all parameter sets for each time step.

colnames(B.11)<-c("ps", "ts", "nbees")
colnames(B.12)<-c("ps", "ts", "nbees")
colnames(B.13)<-c("ps", "ts", "nbees")
colnames(B.14)<-c("ps", "ts", "nbees")
colnames(B.15)<-c("ps", "ts", "nbees")

# bees 
B.11$cor<-rep("locor", nrow(B.11))
B.12$cor<-rep("medcor", nrow(B.12))
B.13$cor<-rep("zerocor", nrow(B.13))
B.14$cor<-rep("medhicor", nrow(B.14))
B.15$cor<-rep("hicor", nrow(B.15))

#bind for comparison
B.1TS<-rbind(B.11, B.12, B.13, B.14, B.15)

#############################################
#############################################
###### summarizing the data #################
#############################################
#############################################

# goal: to average the count data from the last 1000 of the timesteps of each parameter set. this can be used as a set of response variables that are associated with parameter sets. each response variable can be regressed against all parameters in the model to ask how variation in particular parameters account for variation in the state variables. 

### bee summaries : mean
aa.1<-B.1TS[which(B.1TS$ts>(3999)),] 
aa.1<-aa.1[which(aa.1$ts!=ts+1),]
aa.1$nbees[which(aa.1$nbees<1)]<-0
ba.1<-as.data.frame(aa.1 %>% group_by(ps, cor) %>% summarize(mean_bees=mean(nbees, na.rm=T)))
ca.1<-dcast(ba.1, ps~cor, value.var="mean_bees") # for the same parameter set at different covariances 
names(ca.1)<-c("ps", "hicor", "locor", "medhicor", "medlocor", "zerocor")
ca.1$meas<-rep("mean", nrow(ca.1))

### bee summaries: variance
b.a.1<-as.data.frame(aa.1 %>% group_by(ps, cor) %>% summarize(var_bees=var(nbees, na.rm=T)))
c.a.1<-dcast(b.a.1, ps~cor, value.var="var_bees")
names(c.a.1)<-c("ps", "hicor", "locor", "medhicor", "medlocor", "zerocor")
c.a.1$meas<-rep("var", nrow(c.a.1))

# smoosh together mean and vars
paramsets.bmean1<-cbind(paramsets, ca.1[,c(2:7)])
paramsets.bvar1<-cbind(paramsets, c.a.1[,c(2:7)])

# bee df 
paramsets.b1<-rbind(paramsets.bmean1, paramsets.bvar1)

# melt
da.1<-melt(paramsets.b1, id.vars=c('nps', 'p','M','s1','s2','l1','l2', 'a_ii', 'a_jj', 'a_ji', 'a_ij', "lambdab", "L", "meas"), measure.vars=c("locor", "medlocor", "zerocor", "medhicor", "hicor"))
ea.1<-plyr::rename(da.1, replace=c("variable"="cor","value"="bees"))

### plant summaries : mean
ab.1<-PTS1[which(PTS1$ts>(3999)),]
ab.1<-ab.1[which(ab.1$ts!=ts+1),]
ab.1$nplants[which(is.nan(ab.1$nplants))]<-0  # replace 0/0 = NaN to 0
ab.1$nplants[which(ab.1$nplants<1)]<-0 # replace averages of less than one plant with zero
bb.1<-as.data.frame(ab.1 %>% group_by(ps, cor) %>% summarize(var_plants=var(nplants, na.rm=T), mean_plants=mean(nplants, na.rm=T)))
cb.1<-dcast(bb.1, ps~cor, value.var="mean_plants")
names(cb.1)<-c("ps", "hicor", "locor", "medhicor", "medlocor", "zerocor")
cb.1$meas<-rep("mean", nrow(cb.1))

### plant summaries : variance 
b.b.1<-as.data.frame(ab.1 %>% group_by(ps, cor) %>% summarize(var_plants=var(nplants, na.rm=T)))
c.b.1<-dcast(b.b.1, ps~cor, value.var="var_plants")
names(c.b.1)<-c("ps", "hicor", "locor", "medhicor", "medlocor", "zerocor")
c.b.1$meas<-rep("var", nrow(c.b.1))

# smoosh together mean and var 
paramsets.pmean1<-cbind(paramsets, cb.1[,c(2:7)])
paramsets.pvar1<-cbind(paramsets, c.b.1[,c(2:7)])

# plant df 
paramsets.c1<-cbind(paramsets, cb.1[,c(2:7)])

# melt 
db.1<-melt(paramsets.c1, id.vars=c('nps', 'p','M','s1','s2','l1','l2', 'a_ii', 'a_jj', 'a_ji', 'a_ij', "lambdab", "L", "meas"), measure.vars=c("locor", "medlocor", "zerocor", "medhicor", "hicor"))

## put the plants and bees together 
dat1.1<-as.data.frame(cbind(ea.1, db.1$value))
dat.p1<-plyr::rename(dat1.1, replace=c("db.1$value"="plants"))

## write csv 
    write.csv(dat, "twoSp_0.98Germ.survsims.csv")
    write.csv(dat.p1, "oneSp_0.98Germ.survsims.csv")


