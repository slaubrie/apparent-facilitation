## code for ch 3, storage effect and pollinator support
# parameter generation
# after the jump: germination values 
setwd('~/Dropbox/Dissertation/CH3_AmNat_StorageEffectandFacilitation/Simulations/Data/')

#libraries
require(MASS)
require(ggplot2)
require(reshape2)
require(Rmisc)
require(boot)
require(tidyr)
require(dplyr)
require(plyr)

# the goal is to make a whole bunch of simulations where i 
# have multiple parameter sets with the same covariance strucutre

############################################
############################################
############# Parameters  ##################
############################################
############################################

set.seed(12)

## want these to be sampled anew for each set of timeseries 
# so I know how the behavior of timeseries changes 
# with various parameter values .

## number of parameter sets 
k<-1000

#matrix to store time-independent parameters
paramsets<-matrix(data=NA,nrow=k, ncol=13)
paramsets<-as.data.frame(paramsets)
colnames(paramsets)<-c('nps', 'p','M','s1','s2','l1','l2', 'a_ii', 'a_jj', 'a_ji', 'a_ij', "lambdab", "L")
#number of parameter sets; egg -> bee transition; floral resource limitation to bee reproduction ;seed survival of species 1; seed survival of species 2; growth rate of species 1; growth rate of species 2; intraspecific interactions of species 1; intraspecific interactions of species 2; interspecific effect of species 1 on species 2; interspecific effect of species 2 on species 1.

# row to keep track of what timeseries it is 
nps<-seq(1:length(k), length.out=k)
paramsets$nps<-nps

# rate of bee egg developing into adult bee 
p<-sample(seq(0.0001,1, length.out=k))
paramsets$p<-p

#maximum rate of growth for bees, in eggs per bee
# mean number of eggs per bee: 0-80 https://xerces.org/pollinator-conservation/native-bees/
lambdab<-sample(seq(0,5, length.out=k))
paramsets$lambdab<-lambdab

#saturation constant for floral resource limitation to bee reproduction 
M<-sample(seq(1e-4,1, length.out=k)) # use the scale of variation possible in the term L (1e-7, 1e7) 
paramsets$M<-M

#survival rate of seeds in species 1 and 2 
#this is meaningful in the seed-related rates. 
#s1<-sample(seq(0.0001,1, length.out=k))
#s2<-sample(seq(0.0001,1, length.out=k))

## zero survival
s1<-sample(rep(0, length.out=k))
s2<-sample(rep(0, length.out=k))

paramsets$s1<-s1
paramsets$s2<-s1

##growth rates (max seed production, lambda) for both species 
l1<-sample(seq(1,1000, length.out=k))
l2<-sample(seq(1,1000, length.out=k))
paramsets$l1<-l1
paramsets$l2<-l2

#competition of plant species j on plant species j - per capita
#reference http://www.stat.umn.edu/geyer/old/5101/rlook.html
 a_ii<-rep(1,k)
 a_jj<-rep(1,k)

# for the storage effect to be the only thing that promotes coexistence
a_ji<-rep(1,k)
a_ij<-rep(1,k)

#another option seq(0,1,length.out=k)
paramsets$a_ii<-a_ii
paramsets$a_ij<-a_ij

paramsets$a_jj<-a_jj
paramsets$a_ji<-a_ji

#per-bee number of visits made that result in pollen deposition modified by a saturation constant

## first have to define saturation constant 
chi<-sample(seq(0.0001,1, length.out=k))

####### 
#then define per bee # visits that result in pollen deposition (idk how many lol????). 
gamma<-sample(seq(0.0001,1, length.out=k))

########## 
#then define L
L<-chi/gamma
paramsets$L<-L

### then save params 
paramsets.csv<-write.csv(paramsets, file='paramsets.csv')

# ##########################################
# ##########################################
# ################ GERMINATION  ############
# ##########################################
# ##########################################
# #number of timesteps per timeseries
ts<-5000
k<-1000
# #maximum ovule production of plant species 1 and 2 
# # testing using code example from Koons et al 2008
# # comments look the same as their comments because its le same
# #documentation: http://stat.ethz.ch/R-manual/R-devel/library/MASS/html/mvrnorm.html
# ###################################################
cor_cov_mean<-data.frame()

# # correlation through time 
rho=seq(-1,1, by=0.5) 

# # this is the length of values i will loop through correlation info
l=length(rho)

# # where i will store germination values
# #species 1
germ1<-array(0,c(l,k,ts)) 

# #species 2
germ2<-array(0,c(l,k,ts)) 

muval=2 ## i'm going to do (-2 [x], -1 [X],  0 [x],  1 [x],  2 [])

for (l in 1:l){  
	for (j in 1:k){  

# mvrnorm with correlation rho
Z=mvrnorm(ts, mu=c(muval,muval), Sigma=matrix(c(1,rho[l], rho[l], 1),2,2)) 
# more between -2 and 2 

#normal function to Z to get data uniform on interval [0,1] but correlated how i want B) 
# this acutally doesn't need to be fucked with anymore because it's on the interval 0,1
U<-pnorm(Z)

germ1[l,j,]<-U[,1]
germ2[l,j,]<-U[,2]
		}
		# where i will store covariance values - these will be different for every mean X rho 
		# rho = 11. for every muval there should only be 11 values, which correspond to the 11 correlation 				values. 
		covy<-cov(U)
		cov.extract<-covy[1,2]
		mu<-pnorm(muval)
		out<-cbind(cov.extract, mu, rho[l])
		cor_cov_mean<-rbind(cor_cov_mean, out)
		# store out somewhere!

}
colnames(cor_cov_mean)<-c("covariance", "mu", "correlation")
cor_cov_mean


saveRDS(germ1, file=paste("germ1_mu_",round(pnorm(muval), digits=2),sep=""))
saveRDS(germ2, file=paste("germ2_mu_",round(pnorm(muval), digits=2),sep=""))
write.csv(cor_cov_mean, file=paste("cor_cov_mean_mu_",round(pnorm(muval), digits=2),".csv",sep=""))

