## code for ch 3, storage effect and pollinator support
setwd('~/Dropbox/Dissertation/CH3_AmNat_StorageEffectandFacilitation/Simulations/Data/')

#last edit 24 October
rm(list=ls())

#libraries
require(MASS)
require(ggplot2)
require(reshape2)
require(Rmisc)
require(boot)
require(tidyr)
require(dplyr)
require(plyr)
require(ppcor)
require(MuMIn)
require(ggpubr)
require(emmeans)
require(logistf)


##########################################
####### load response sims ###############
##########################################
# two species 
g2_hi<-read.csv('twoSp_0.98Germ.sims.csv', header=T)
g2_hi$gmu<-rep(0.98, nrow(g2_hi))
g2_hi$nsp<-rep(2, nrow(g2_hi))

g2_medhi<-read.csv('twoSp_0.84Germ.sims.csv', header=T)
g2_medhi$gmu<-rep(0.84, nrow(g2_medhi))
g2_medhi$nsp<-rep(2, nrow(g2_medhi))

g2_med<-read.csv('twoSp_0.5Germ.sims.csv', header=T)
g2_med$gmu<-rep(0.5, nrow(g2_med))
g2_med$nsp<-rep(2, nrow(g2_med))

g2_medlo<-read.csv('twoSp_0.16Germ.sims.csv', header=T)
g2_medlo$gmu<-rep(0.16, nrow(g2_medlo))
g2_medlo$nsp<-rep(2, nrow(g2_medlo))

g2_lo<-read.csv('twoSp_0.02Germ.sims.csv', header=T)
g2_lo$gmu<-rep(0.02, nrow(g2_lo))
g2_lo$nsp<-rep(2, nrow(g2_lo))

# one species 
g1_hi<-read.csv('oneSp_0.98Germ.sims.csv', header=T)
g1_hi$gmu<-rep(0.98, nrow(g1_hi))
g1_hi$nsp<-rep(1, nrow(g1_hi))

g1_medhi<-read.csv('oneSp_0.84Germ.sims.csv', header=T)
g1_medhi$gmu<-rep(0.84, nrow(g1_medhi))
g1_medhi$nsp<-rep(1, nrow(g1_medhi))

g1_med<-read.csv('oneSp_0.5Germ.sims.csv', header=T)
g1_med$gmu<-rep(0.5, nrow(g1_med))
g1_med$nsp<-rep(1, nrow(g1_med))

g1_medlo<-read.csv('oneSp_0.16Germ.sims.csv', header=T)
g1_medlo$gmu<-rep(0.16, nrow(g1_medlo))
g1_medlo$nsp<-rep(1, nrow(g1_medlo))

g1_lo<-read.csv('oneSp_0.02Germ.sims.csv', header=T)
g1_lo$gmu<-rep(0.02, nrow(g1_lo))
g1_lo$nsp<-rep(1, nrow(g1_lo))

# single species bind 
cordata1<-rbind(g1_hi, g1_medhi, g1_med, g1_medlo, g1_lo)
cordata1$gmu<-as.factor(cordata1$gmu)
cordata1$nsp<-as.factor(cordata1$nsp)
cordata1<-plyr::rename(cordata1, replace=c("gmu"="Germination"))
cordata1[is.na(cordata1)]<-0

## means 
mu.cor1<-cordata1[which(cordata1$meas=="mean"),] # this is the est. mean (samples=timesteps) in each timeseries
mu.cor1<-plyr::rename(mu.cor1, replace=c("nsp"="N.Species"))

# var 
sig.cor1<-cordata1[which(cordata1$meas=="var"),]
sig.cor1<-plyr::rename(sig.cor1, replace=c("nsp"="N.Species"))

cordat1<-cbind(mu.cor1, sig.cor1$bees, sig.cor1$plants)
cordat1<-plyr::rename(cordat1, replace=c("sig.cor1$bees"="bvar", "sig.cor1$plants"="pvar"))

cordat1<-cordat1[which(cordat1$cor=='locor'),]
cordat1$Treatment <-rep('onesp', nrow(cordat1))

#####################
###### 2 species bind
#####################
cordata<-rbind(g2_hi, g2_medhi, g2_med, g2_medlo, g2_lo)
cordata$gmu<-as.factor(cordata$gmu)
cordata$nsp<-as.factor(cordata$nsp)
cordata<-plyr::rename(cordata, replace=c("gmu"="Germination"))
cordata[is.na(cordata)]<-0

## means 
mu.cor<-cordata[which(cordata$meas=="mean"),] # this is the est. mean (samples=timesteps) in each timeseries
mu.cor<-plyr::rename(mu.cor, replace=c("nsp"="N.Species"))

# var 
sig.cor<-cordata[which(cordata $meas=="var"),]
sig.cor<-plyr::rename(sig.cor, replace=c("nsp"="N.Species"))

# 
cordat<-cbind(mu.cor, sig.cor$bees, sig.cor$plants)
cordat<-plyr::rename(cordat, replace=c("sig.cor$bees"="bvar", "sig.cor$plants"="pvar"))
cordat$Treatment<-cordat$cor



###############################################
############ full dataset #####################
###############################################
fulldat<-rbind(cordat,cordat1)

## calculate CVs 
fulldat$b.cv<-((sqrt(fulldat$bvar))/fulldat$bees)*100
fulldat$p.cv<-((sqrt(fulldat$pvar))/fulldat$plants)*100

# sometimes bees do not survive or plants do not survive, and the others continue. this is due to seed survival, I think - there are just cases where germiantion and seed survival as high, and the seed bank hasn't been exhausted yet. maybe run for longer timesteps...
## here are the cases (only occur in single species communities)
nobee<-fulldat[which(fulldat$bees==0 & fulldat$plants>0),] # do exist
noplant<-sort(fulldat$nps[which(fulldat$plants==0 & fulldat$bees>0)]) #don't exist

# we could leave these out of the analysis, but discuss them in the paper. here is the code for that: 
#mu.data<-mu.data[which(mu.data$bees>0 & mu.data$plants>0|mu.data$bees==0 & mu.data$plants==0),]
# including them or not does not change the main results of the paper.

## logs
fulldat$lbee<-log(fulldat$bees) 
fulldat$lbee[!is.finite(fulldat$lbee)]<-0
fulldat$lplant<-log(fulldat$plants)
fulldat$lplant[!is.finite(fulldat$lplant)]<-0



########### add the covariances 
# each of these has a different mean value
# for each mean value, there are treatments, which are already in the datasets. however, each of those germination X treatment values carries with it a covariance of germination as well, which is where the following comes in. 

logerm_cov<-read.csv("cor_cov_mean_mu_0.02.csv", header=T)
lgc<-c(rep(logerm_cov$covariance[1], 100),rep(logerm_cov$covariance[2], 100),rep(logerm_cov$covariance[3], 100),rep(logerm_cov$covariance[4], 100), rep(logerm_cov$covariance[5], 100), rep(NA, 100))

medlogerm_cov<-read.csv("cor_cov_mean_mu_0.16.csv", header=T)
mlgc<-c(rep(medlogerm_cov$covariance[1], 100),rep(medlogerm_cov$covariance[2], 100),rep(medlogerm_cov$covariance[3], 100),rep(medlogerm_cov$covariance[4], 100), rep(medlogerm_cov$covariance[5], 100), rep(NA, 100))


medgerm_cov<-read.csv("cor_cov_mean_mu_0.5.csv", header=T)
mgc<-c(rep(medgerm_cov$covariance[1], 100),rep(medgerm_cov$covariance[2], 100),rep(medgerm_cov$covariance[3], 100),rep(medgerm_cov$covariance[4], 100), rep(medgerm_cov$covariance[5], 100), rep(NA, 100))

medhigerm_cov<-read.csv("cor_cov_mean_mu_0.84.csv", header=T)
mhgc<-c(rep(medhigerm_cov$covariance[1], 100),rep(medhigerm_cov$covariance[2], 100),rep(medhigerm_cov$covariance[3], 100),rep(medhigerm_cov$covariance[4], 100), rep(medhigerm_cov$covariance[5], 100), rep(NA, 100))

higerm_cov<-read.csv("cor_cov_mean_mu_0.98.csv", header=T)
hgc<-c(rep(higerm_cov$covariance[1], 100),rep(higerm_cov$covariance[2], 100),rep(higerm_cov$covariance[3], 100),rep(higerm_cov$covariance[4], 100), rep(higerm_cov$covariance[5], 100), rep(NA, 100))


### sort fulldat by Germination 
fulldat1<-fulldat[with(fulldat, order(Germination)),]
fulldat1$covs<-c(lgc, mlgc, mgc, mhgc, hgc)

### analysis 
# choose continuous variables to use based on rho 
# choose germination and nsp automatically - bc they are categorical 

### data - numeric 

## look at relationshipe of means to the rest of the variables 
numdat1<-(fulldat1[,c(3:5,7, 13:14,26:27)]) #should i add the variance in these two?
numdat.mu1<-cbind(numdat1)
numdat.mu<-numdat.mu1[which(complete.cases(numdat.mu1)==T),]
names(numdat.mu)

## look at relationshipe of CV to the rest of the variables 
numdat2<-(fulldat1[,c(3:5,7, 13:14, 24:25)]) #should i add the variance in these two?
numdat.mu2<-cbind(numdat2)
numdat.cv<-numdat.mu2[which(complete.cases(numdat.mu2)==T),]
names(numdat.cv)
numdat.cv<-numdat.cv[which(complete.cases(numdat.cv)==T),]
numdat.cv$p.cv[which(is.infinite(numdat.cv$p.cv))]<-0


######################### look at rank correlations ######################### 
cor.dat<-pcor(numdat.mu, method=c('pearson'))
partialCorrs.mu<-cor.dat$estimate[,c(7:8)]
partialCorrs.mu

cor.dat2<-pcor(numdat.cv, method=c('pearson'))
partialCorrs.cv<-cor.dat2$estimate[,c(7:8)]
partialCorrs.cv

corrsdata<-cbind(partialCorrs.mu, partialCorrs.cv)


###################### plot plant and bee relationship#######################
cor(fulldat1$lbee, fulldat1$lplant, method='pearson')

ggplot(fulldat1, aes(lplant, lbee))+geom_jitter(width=0.1, height=0,col='gray77', alpha=0.3, size=4)+theme_classic()+geom_smooth(method='lm', lwd=2.5, se=F, col='black')+labs(x="log(Plants)", y="log(Bees)")+theme(axis.title=element_text(size=14), axis.text=element_text(size=14))+annotate("text", label="Pearson's r = 0.96", x=2, y=7, size=4)

ggplot(fulldat1, aes(Germination, lplant))+stat_summary()+theme_classic()+labs(y="log(Plants)", x="Mean Germination")

#######################
## crashing analysis ## 
#######################
gcolz<-c("#D99693", "#268184","#93C7DA", "#555B6E", "#BEE3DB")

# zero columns
fulldat1$bzero<-ifelse(fulldat1$lbee>0,1,0)
fulldat1$pzero<-ifelse(fulldat1$lplant>0,1,0)

# look at if things survived without each other 
sort(fulldat1$nps[which(fulldat1 $bees==0 & fulldat1 $plants>0)])
sort(fulldat1 $nps[which(fulldat1 $plants==0 & fulldat1 $bees>0)])


### bees 
# look at how many bee communities survive with different treatments 
positions<-c("locor", "medlocor", "zerocor", "medhicor", "hicor", "onesp") #fill in treatments 

## plot cors and covs 
ggplot(fulldat1, aes(Treatment, covs))+geom_hline(yintercept=0, lwd=1, lty=2, col='gray88')+stat_summary(fun.y='mean', geom='point', size=4)+scale_x_discrete(limits=positions)+theme_classic()+labs(x="Treatment", y="Average germination covariance")+theme(axis.title=element_text(size=13), axis.text=element_text(size=19))

## bee population survival: does it matter according to how correlated the germination rates are? 
cors<-as.data.frame(fulldat1%>%group_by(Treatment, Germination)%>%summarize(bzc=sum(bzero)))

bcorplot<-ggplot(cors, aes(Treatment, bzc))+geom_bar(stat="identity", position='stack')+scale_x_discrete(limits=positions)+theme_classic()+labs(y="Total Bee Population Successes",x="")+scale_fill_manual()+theme(axis.title=element_text(size=13), axis.text=element_text(size=20), legend.title=element_text(size=20), legend.text=element_text(size=20))
bcorplot


bcorplot2<-ggplot(cors, aes(Treatment, bzc, fill=Germination))+geom_bar(stat="identity", position='dodge2')+scale_x_discrete(limits=positions)+theme_bw()+labs(y="Bee Population Successes",x="")+scale_fill_manual(values=gcolz)+theme(axis.title=element_text(size=13), axis.text=element_text(size=20), legend.title=element_text(size=20), legend.text=element_text(size=20))

bcorplot2
 
## bee population survival: does it matter according to covariance in germination rates?
covs<-as.data.frame(fulldat1%>%group_by(covs, Treatment)%>%summarize(bzc=sum(bzero)))
ggplot(covs, aes(covs, bzc))+geom_point()+theme_bw()+ylim(0,50)+stat_smooth(method=lm, se=F)

## plant survival: does it matter according to how correlated the germination rates are? 
pcors<-as.data.frame(fulldat1%>%group_by(Treatment, Germination)%>%summarize(pzc=sum(pzero)))

pcorplot<-ggplot(pcors, aes(Treatment, pzc))+geom_bar(stat="identity")+scale_x_discrete(limits=positions)+theme_classic()+labs(y="Total Plant Community Successes",x="")+theme(axis.title=element_text(size=13), axis.text=element_text(size=20), legend.title=element_text(size=20), legend.text=element_text(size=20))

pcorplot

pcorplot2<-ggplot(pcors, aes(Treatment, pzc, fill=Germination))+geom_bar(stat="identity", position='dodge2')+scale_x_discrete(limits=positions)+theme_bw()+labs(y="Plant Community Successes",x="")+scale_fill_manual(values=gcolz)+theme(axis.title=element_text(size=13), axis.text=element_text(size=20), legend.title=element_text(size=20), legend.text=element_text(size=20))

pcorplot2

pcov<-as.data.frame(fulldat1%>%group_by(covs)%>%summarize(pzc=sum(pzero)))
ggplot(pcov, aes(covs, pzc))+geom_point()+theme_bw()+stat_smooth(method=lm, se=F)

plot1<-ggarrange(bcorplot, pcorplot, bcorplot2, pcorplot2, common.legend=T, legend="right", labels=c("A","B","C","D"))
annotate_figure(plot1, bottom=text_grob("Treatment", vjust=-0.5, hjust=0.75, size=30))
plot1

# plot
# pdf(file='barSuccess.pdf', width=10, height=8)
# plot(success)
# dev.off()


###### calculate CV and crashing ######
# each treatment: 500 times per treatment
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3444174/

#####################################
########### Risk ratios #############
#####################################

########  bees ######## 

## Total Successes
h<-fulldat1%>%group_by(Treatment)%>%summarize(bzc=sum(bzero))
h<-as.data.frame(h)
h$percent<-(h$bzc/500)*100
h

## risk ratio 
rr.b.tot<-h[,3]/h[2,3]
rr.b.tot

## Successes by germination 
# data 
h2<-fulldat1%>%group_by(Treatment, Germination)%>%summarize(bzc=sum(bzero))
h2<-as.data.frame(h2)
h2$percent<-(h2$bzc/100)*100

# 0.02 germ
hlo<-h2[which(h2$Germination==0.02),]
hlo$rr<-hlo[2,4]/hlo[,4]

# 0.16 germ
hmlo<-h2[which(h2$Germination==0.16),]
hmlo$rr<-hmlo[2,4]/hmlo[,4]

# 0.5 germ
hm<-h2[which(h2$Germination==0.5),]
hm$rr<-hm[2,4]/hm[,4]

# 0.84 germ
hmhi<-h2[which(h2$Germination==0.84),]
hmhi$rr<-hmhi[2,4]/hmhi[,4]

# 0.98 germ
hhi<-h2[which(h2$Germination==0.98),]
hhi$rr<-hhi[2,4]/hhi[,4]

## bind 
b.rr.dat<-rbind(hlo, hmlo,hm,hmhi,hhi)

########  plants ######## 
i<-fulldat1%>%group_by(Treatment)%>%summarize(bzc=sum(bzero))
i<-as.data.frame(i)
i$percent<-(i$bzc/500)*100
i

## risk ratio 
rr.b.tot<-i[,3]/i[2,3]
rr.b.tot

## Successes by germination 
# data 
i2<-fulldat1%>%group_by(Treatment, Germination)%>%summarize(pzc=sum(pzero))
i2<-as.data.frame(i2)
i2$percent<-(i2$pzc/100)*100

# 0.02 germ
ilo<-i2[which(i2$Germination==0.02),]
ilo$rr<-ilo[2,4]/ilo[,4]

# 0.16 germ
imlo<-i2[which(i2$Germination==0.16),]
imlo$rr<-imlo[2,4]/imlo[,4]

# 0.5 germ
im<-i2[which(i2$Germination==0.5),]
im$rr<-im[2,4]/im[,4]

# 0.84 germ
imii<-i2[which(i2$Germination==0.84),]
imii$rr<-imii[2,4]/imii[,4]

# 0.98 germ
iii<-i2[which(i2$Germination==0.98),]
iii$rr<-iii[2,4]/iii[,4]

## bind 
p.rr.dat<-rbind(ilo, imlo,im,imii,iii)

#####################################
############## CVs ##################
#####################################

## get mean and sd of coefficients of variation in bee communities 
j<-fulldat1%>%group_by(Germination, Treatment, bzero)%>%summarize(bcv=mean(b.cv), bcvs=sqrt(var(b.cv)))
j<-as.data.frame(j)
j<-j[which(complete.cases(j)==T),]


## plants
# choose just the survivors
k<-fulldat1%>%group_by(Germination, Treatment, pzero)%>%summarize(pcv=mean(p.cv), pcvs=sqrt(var(p.cv)))
k<-as.data.frame(k)
k<-k[which(complete.cases(k)==T),]

################
### plot CVs ###
################

bcv<-ggplot(data=j, aes(Treatment, bcv, color=Germination, group=Germination))+geom_point(size=3)+theme_bw()+scale_x_discrete(limits=positions)+labs(x="", y="Bee Population CV")+scale_color_manual(values=gcolz)+lims(y=c(0,150))+theme(axis.title=element_text(size=13), axis.text=element_text(size=20), legend.title=element_text(size=20), legend.text=element_text(size=20))
pcv
bcv

pcv<-ggplot(data=k, aes(Treatment,pcv, color=Germination, group=Germination))+geom_point(size=3)+theme_bw()+scale_x_discrete(limits=positions)+labs(x="", y="Plant Community CV")+scale_color_manual(values=gcolz)+lims(y=c(0,150))+theme(axis.title=element_text(size=13), axis.text=element_text(size=20), legend.title=element_text(size=20), legend.text=element_text(size=20))
pcv

plot2<-ggarrange(pcv, bcv, common.legend=T, legend='right', labels=c("A","B"), ncol=2)
annotate_figure(plot2, bottom=text_grob("Treatment", vjust=-1, hjust=1, size=14))


## cv together
cvdat<-rbind(j,k)

# data tables 
#write.csv(corrsdata, file="corrdat.csv")
# write.csv(k, file="plantcohensDs.csv")
# write.csv(j, file="beecohensDs.csv")
# write.csv(p.rr.dat, file="plantRR.csv")
# write.csv(b.rr.dat, file="beeRR.csv")


