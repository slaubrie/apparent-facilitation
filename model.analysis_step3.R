# draft manuscript title: Pollinator support drives inter-annual facilitation in a model of annual plant communities 
# code for analyzing and plotting model output
# last edit July 2022
# run this after running steps 1 and 2 


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
require(rstatix)
require(grDevices)


####### load response sims ###############

### two species 

## diff blocks are different mean germination rates 
g2_hi<-read.csv('twoSp_0.98Germ.sims.csv', header=T) #read file 
g2_hi$gmu<-rep(0.98, nrow(g2_hi))  # add treatment 
g2_hi$nsp<-rep(2, nrow(g2_hi)) # add number of plant species

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

####### data shaping ###############


###### 1 species bind

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


###### 2 species bind


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

cordat<-cbind(mu.cor, sig.cor$bees, sig.cor$plants)
cordat<-plyr::rename(cordat, replace=c("sig.cor$bees"="bvar", "sig.cor$plants"="pvar"))
cordat$Treatment<-cordat$cor


############ full dataset 


fulldat<-rbind(cordat,cordat1)

## calculate CVs 
fulldat$b.cv<-((sqrt(fulldat$bvar))/fulldat$bees)*100
fulldat$p.cv<-((sqrt(fulldat$pvar))/fulldat$plants)*100

# sometimes bees do not survive or plants do not survive, and the others continue. 
# can look at that here 
## here are the cases (only occur in single species communities)
# nobee<-fulldat[which(fulldat$bees==0 & fulldat$plants>0),] 
# noplant<-sort(fulldat$nps[which(fulldat$plants==0 & fulldat$bees>0)]) 

# we could leave these out of the analysis, but discuss them in the paper. here is the code for that: 
# fulldat<-fulldat[which(fulldat$bees>0 & fulldat$plants>0|fulldat$bees==0 & fulldat$plants==0),]
# including them or not does not change the main results of the paper.

## logs of abundance
fulldat$lbee<-log(fulldat$bees) 
fulldat$lbee[!is.finite(fulldat$lbee)]<-0
fulldat$lplant<-log(fulldat$plants)
fulldat$lplant[!is.finite(fulldat$lplant)]<-0

######## add covariances 
 
# each of these has a different mean value
# for each mean value, there are treatments 
# each of those germination rate X germination correlation treatment carries with it the 
# covariance of germination as well
# each value needs to be repeated X the number of parametersets 

logerm_cov<-read.csv("cor_cov_mean_mu_0.02.csv", header=T)
lgc<-c(rep(logerm_cov$covariance[1], 1000),rep(logerm_cov$covariance[2], 1000),rep(logerm_cov$covariance[3], 1000),rep(logerm_cov$covariance[4], 1000), rep(logerm_cov$covariance[5], 1000), rep(NA, 1000))

medlogerm_cov<-read.csv("cor_cov_mean_mu_0.16.csv", header=T)
mlgc<-c(rep(medlogerm_cov$covariance[1], 1000),rep(medlogerm_cov$covariance[2], 1000),rep(medlogerm_cov$covariance[3], 1000),rep(medlogerm_cov$covariance[4], 1000), rep(medlogerm_cov$covariance[5], 1000), rep(NA, 1000))


medgerm_cov<-read.csv("cor_cov_mean_mu_0.5.csv", header=T)
mgc<-c(rep(medgerm_cov$covariance[1], 1000),rep(medgerm_cov$covariance[2], 1000),rep(medgerm_cov$covariance[3], 1000),rep(medgerm_cov$covariance[4], 1000), rep(medgerm_cov$covariance[5], 1000), rep(NA, 1000))

medhigerm_cov<-read.csv("cor_cov_mean_mu_0.84.csv", header=T)
mhgc<-c(rep(medhigerm_cov$covariance[1], 1000),rep(medhigerm_cov$covariance[2], 1000),rep(medhigerm_cov$covariance[3], 1000),rep(medhigerm_cov$covariance[4], 1000), rep(medhigerm_cov$covariance[5], 1000), rep(NA, 1000))

higerm_cov<-read.csv("cor_cov_mean_mu_0.98.csv", header=T)
hgc<-c(rep(higerm_cov$covariance[1], 1000),rep(higerm_cov$covariance[2], 1000),rep(higerm_cov$covariance[3], 1000),rep(higerm_cov$covariance[4], 1000), rep(higerm_cov$covariance[5], 1000), rep(NA, 1000))

### sort fulldat by Germination 
fulldat1<-fulldat[with(fulldat, order(Germination)),]
fulldat1$covs<-c(lgc, mlgc, mgc, mhgc, hgc)

####### analysis ###############

### data - numeric 

## look at relationship of means to the rest of the variables 
numdat1<-(fulldat1[,c(3:5,7, 13:14,26:27)]) #should i add the variance in these two?
numdat.mu1<-cbind(numdat1)
numdat.mu<-numdat.mu1[which(complete.cases(numdat.mu1)==T),]
names(numdat.mu)

## look at relationships of CV to the rest of the variables 
numdat2<-(fulldat1[,c(3:5,7, 13:14, 24:25)]) #should i add the variance in these two?
numdat.mu2<-cbind(numdat2)
numdat.cv<-numdat.mu2[which(complete.cases(numdat.mu2)==T),]
names(numdat.cv)
numdat.cv<-numdat.cv[which(complete.cases(numdat.cv)==T),]
numdat.cv$p.cv[which(is.infinite(numdat.cv$p.cv))]<-0


######################### look at rank correlations
cor.dat<-pcor(numdat.mu, method=c('pearson'))
partialCorrs.mu<-cor.dat$estimate[,c(7:8)]
partialCorrs.mu

cor.dat2<-pcor(numdat.cv, method=c('pearson'))
partialCorrs.cv<-cor.dat2$estimate[,c(7:8)]
partialCorrs.cv

corrsdata<-cbind(partialCorrs.mu, partialCorrs.cv)


###################### plot plant and bee relationship
cor(fulldat1$lbee, fulldat1$lplant, method='pearson')

ggplot(fulldat1, aes(lplant, lbee))+geom_jitter(width=0.1, height=0,col='gray77', alpha=0.3, size=4)+theme_classic()+geom_smooth(method='lm', lwd=2.5, se=F, col='black')+labs(x="log(Plants)", y="log(Bees)")+theme(axis.title=element_text(size=14), axis.text=element_text(size=14))+annotate("text", label="Pearson's r = 0.96", x=2, y=7, size=4)

ggplot(fulldat1, aes(Germination, lplant))+stat_summary()+theme_classic()+labs(y="log(Plants)", x="Mean Germination")


#### crashing analysis #### 


gcolz<-c("#D99693", "#268184","#67B0CB", "#555B6E", "#D59B2F") #colors 

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
corcovplot<-ggplot(fulldat1, aes(Treatment, covs))+geom_hline(yintercept=0, lwd=1, lty=2, col='gray88')+stat_summary(fun.y='mean', geom='point', size=4)+scale_x_discrete(limits=positions[1:5], labels=c("-1", "-0.5", "0", "0.5", "1"))+theme_classic()+labs(x="Germination correlation treatment", y="Average germination covariance")+theme(axis.title=element_text(size=20), axis.text=element_text(size=19))

## bee population survival: does it matter according to how correlated the germination rates are? 
cors<-as.data.frame(fulldat1%>%group_by(Treatment, Germination)%>%summarize(bzc=sum(bzero)))

bcorplot<-ggplot(cors, aes(Treatment, bzc))+
  geom_bar(stat="identity", position='stack')+
  scale_x_discrete(limits=positions, 
                   label=c("-1", "-0.5", "0", "0.5", "1", '1sp'))+
  theme_classic()+labs(y="Bee population",x="")+
  scale_fill_manual()+
  theme(axis.title=element_text(size=15), 
        axis.text=element_text(size=15), 
        legend.title=element_text(size=15), 
        legend.text=element_text(size=15))


bcorplot2<-ggplot(cors, aes(Treatment, bzc, fill=Germination))+
  geom_bar(stat="identity", position='dodge2')+
  scale_x_discrete(limits=positions, 
                   label=c("-1", "-0.5", "0", "0.5", "1", '1sp'))+
  theme_bw()+
  labs(y="",x="")+
  scale_fill_manual(values=gcolz)+
  theme(axis.title=element_text(size=15), 
        axis.text=element_text(size=15), 
        legend.title=element_text(size=15), 
        legend.text=element_text(size=15))


## bee population survival: does it matter according to covariance in germination rates?
covs<-as.data.frame(fulldat1%>%group_by(covs, Treatment)%>%summarize(bzc=sum(bzero)))
ggplot(covs, aes(covs, bzc))+geom_point()+theme_bw()+ylim(0,500)+stat_smooth(method=lm, se=F)

## plant survival: does it matter according to how correlated the germination rates are? 
pcors<-as.data.frame(fulldat1%>%group_by(Treatment, Germination)%>%summarize(pzc=sum(pzero)))

pcorplot<-ggplot(pcors, aes(Treatment, pzc))+
  geom_bar(stat="identity")+
  scale_x_discrete(limits=positions, label=c("-1", "-0.5", "0", "0.5", "1", '1sp'))+
  theme_classic()+labs(y="Plant community",x="")+
  theme(axis.title=element_text(size=15), 
        axis.text=element_text(size=15), 
        legend.title=element_text(size=15), 
        legend.text=element_text(size=15))

pcorplot2<-ggplot(pcors, aes(Treatment, pzc, fill=Germination))+
  geom_bar(stat="identity", position='dodge2')+
  scale_x_discrete(limits=positions, label=c("-1", "-0.5", "0", "0.5", "1", '1sp'))+
  theme_bw()+
  labs(y="",x="")+
  scale_fill_manual(values=gcolz)+
  theme(axis.title=element_text(size=15), 
        axis.text=element_text(size=15), 
        legend.title=element_text(size=15), 
        legend.text=element_text(size=15))


pcov<-as.data.frame(fulldat1%>%group_by(covs)%>%summarize(pzc=sum(pzero)))
ggplot(pcov, aes(covs, pzc))+geom_point()+theme_bw()+stat_smooth(method=lm, se=F)+ylim(0,500)

plot1<-ggarrange(pcorplot, pcorplot2, bcorplot, bcorplot2, common.legend=T, legend="right", labels=c("","","",""))
annotate_figure(plot1, left=text_grob("Persistence", rot=90, size=20), bottom=text_grob("Germination correlation treatment", size=20, vjust=-1))

########### Risk ratios #############

# each treatment: 5000 times per treatment (check this by looking at table(fulldat1$Treatment))
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3444174/

# the way i'm doing this is clunky, but whatever 

########  bees

## Total Successes
h<-fulldat1%>%group_by(Treatment)%>%summarize(bzc=sum(bzero))
h<-as.data.frame(h)
h$percent<-(h$bzc/5000)*100
h

## risk ratio 
rr.b.tot<-h[,3]/h[2,3]
rr.b.tot

## Successes by germination 
# data 
h2<-fulldat1%>%group_by(Treatment, Germination)%>%summarize(bzc=sum(bzero))
h2<-as.data.frame(h2)
h2$percent<-(h2$bzc/1000)*100

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

########  plants
i<-fulldat1%>%group_by(Treatment)%>%summarize(bzc=sum(bzero))
i<-as.data.frame(i)
i$percent<-(i$bzc/5000)*100
i

## risk ratio 
rr.b.tot<-i[,3]/i[2,3]
rr.b.tot

## Successes by germination 
# data 
i2<-fulldat1%>%group_by(Treatment, Germination)%>%summarize(pzc=sum(pzero))
i2<-as.data.frame(i2)
i2$percent<-(i2$pzc/1000)*100

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

########### CVs #############

## plants
# choose the survivors
k<-fulldat1%>%group_by(Germination, Treatment, pzero)%>%summarize(pcv=mean(p.cv), pcvs=sqrt(var(p.cv)))
k<-as.data.frame(k)
k<-k[which(complete.cases(k)==T),]

k1<-k%>%filter(Treatment %in% 'onesp')
k2<-k[which(k$Treatment!='onesp'),]

# bees 
# choose the survivors
j<-fulldat1%>%group_by(Germination, Treatment, bzero)%>%summarize(bcv=mean(b.cv), pcvs=sqrt(var(b.cv)))
j<-as.data.frame(j)
j<-j[which(complete.cases(j)==T),]

j1<-j%>%filter(Treatment %in% 'onesp')
j2<-j[which(k$Treatment!='onesp'),]

########### Cohen's d #############

dsubdat<-fulldat1[c(19,23:25)]
dsubdat<-dsubdat[which(dsubdat$b.cv>0),]

grp_j<-as.data.frame(dsubdat%>%pivot_longer(c(b.cv, p.cv), values_to='value')%>% filter(name=="b.cv")) 
grp_j$trt<-paste(grp_j$Germination, grp_j$Treatment, sep="")
tryj<-as.data.frame(grp_j %>% group_by(Germination) %>% cohens_d(value~Treatment))%>%filter(group1=="locor"|group2=="locor")

tryj2<-tryj[c(3:5,8)]
tryj2$group2[which(tryj2$group2=="locor")]<-"hicor"

# this is for plotting 
threshhold1<-cbind.data.frame(rep(0, 5), levels(k$Germination))
names(threshhold1)<-c("thresh", 'Germination')

bee_cohens<-ggplot(tryj2, aes(group2, effsize, group=Germination))+
  geom_hline(data=threshhold1,aes(yintercept=thresh), color='darkgray', lty=2)+
  annotate("rect", xmin=0, xmax=6, ymin=-0.1999,ymax=0.1999, alpha=0.15)+
  ylim(-0.5,0.5)+
  geom_point(aes(color=Germination), size=2.5)+
  facet_grid(cols=vars(Germination))+
  scale_x_discrete(limits=positions[2:6], labels=c("-0.5","0","0.5","1", "1sp"))+
  theme_bw()+
  scale_color_manual(values=gcolz)+
  theme(axis.title.x=element_blank(), 
        axis.title.y=element_text(size=13),
        axis.text=element_text(size=12), 
        axis.text.y=element_text(angle=25),
        legend.title=element_text(size=15), 
        legend.text=element_text(size=12),
        strip.text.x = element_blank(),
  )+
  labs(y="Bee population")+
  scale_y_continuous(
    labels = scales::label_number(accuracy = 0.1))
bee_cohens

grp_k<-as.data.frame(dsubdat%>%pivot_longer(c(b.cv, p.cv), values_to='value')%>% filter(name=="p.cv")) 
grp_k$trt<-paste(grp_k$Germination, grp_k$Treatment, sep="")
tryk<-as.data.frame(grp_k%>%group_by(Germination)%>%cohens_d(value~Treatment))%>%filter(group1=="locor"|group2=="locor")

tryk2<-tryk[c(3:5,8)]
tryk2$group2[which(tryk2$group2=="locor")]<-"hicor"

plant_cohens<-ggplot(tryk2, aes(group2, effsize, group=Germination))+
  geom_hline(data=threshhold1,aes(yintercept=thresh), color='darkgray', lty=2)+
  annotate("rect", xmin=0, xmax=6, ymin=-0.1999,ymax=0.1999, alpha=0.15)+
  ylim(-0.5,1.5)+
  geom_point(aes(color=Germination), size=2.5)+
  facet_grid(cols=vars(Germination))+
  scale_x_discrete(limits=positions[2:6], 
                   labels=c("-0.5","0","0.5","1", "1sp"))+
  theme_bw()+
  scale_color_manual(values=gcolz)+
  theme(axis.title.x=element_blank(), 
        axis.title.y=element_text(size=13),
        axis.text=element_text(size=12), 
        axis.text.y=element_text(angle=25),
        legend.title=element_text(size=15), 
        legend.text=element_text(size=12),
        strip.text.x = element_blank(),
  )+
  labs(y="Plant community")+
  scale_y_continuous(
    labels = scales::label_number(accuracy = 0.1))
  
plant_cohens

plot4<-ggarrange(plant_cohens, 
                 bee_cohens, 
                 common.legend=T, 
                 legend='right', nrow=2)
annotate_figure(plot4, 
                left=text_grob("Cohen's d", rot=90, size=20), 
                bottom=text_grob("Germination correlation treatment", 
                                 size=20))


### plot CVs


bcv_2sp<-ggplot(data=j2, aes(Treatment,bcv, 
                             color=Germination, 
                             group=Germination))+
  geom_point(size=3)+theme_bw()+
  scale_x_discrete(limits=positions[1:5], labels=c("-1","-0.5","0","0.5","1"))+
  labs(x="", y="Bee population")+
  scale_color_manual(values=gcolz)+
  lims(y=c(0,125))+
  theme(axis.title=element_text(size=15), 
        axis.text=element_text(size=15), 
        axis.text.x=element_text(size=20), 
        legend.title=element_text(size=15), 
        legend.text=element_text(size=15))+
  geom_line()
bcv_2sp

bcv_1sp<-ggplot(data=j1, aes(Treatment,bcv, 
                             color=Germination, 
                             group=Germination))+
  geom_point(size=3)+theme_bw()+
  scale_x_discrete(limits=positions[6], labels=c("1sp"))+
  labs(x="", y="")+
  scale_color_manual(values=gcolz[2:5])+
  lims(y=c(0,125))+
  theme(axis.title=element_text(size=18), 
        axis.ticks.y=element_blank(), 
        axis.text.y=element_blank(),
        axis.text=element_text(size=20), 
        legend.title=element_text(size=15), 
        legend.text=element_text(size=15))
bcv_1sp

pcv_2sp<-ggplot(data=k2, aes(Treatment,pcv, 
                             color=Germination, 
                             group=Germination))+
  geom_point(size=3)+theme_bw()+
  scale_x_discrete(limits=positions[1:5], 
                   labels=c("-1","-0.5","0","0.5","1"))+
  labs(x="", y="Plant community")+
  scale_color_manual(values=gcolz)+
  lims(y=c(0,125))+
  theme(axis.title=element_text(size=15), 
        axis.text=element_text(size=15), 
        axis.text.x=element_text(size=20), 
        legend.title=element_text(size=15), 
        legend.text=element_text(size=15))+
  geom_line()
pcv_2sp

pcv_1sp<-ggplot(data=k1, aes(Treatment,pcv, 
                             color=Germination, 
                             group=Germination))+
  geom_point(size=3)+theme_bw()+
  scale_x_discrete(limits=positions[6], labels=c("1sp"))+
  labs(x="", y="")+
  scale_color_manual(values=gcolz[2:5])+
  lims(y=c(0,125))+
  theme(axis.title=element_text(size=18), 
        axis.ticks.y=element_blank(), 
        axis.text.y=element_blank(),
        axis.text=element_text(size=20), 
        legend.title=element_text(size=15), 
        legend.text=element_text(size=15))
pcv_1sp

plot2<-ggarrange(pcv_2sp, pcv_1sp, bcv_2sp, bcv_1sp, common.legend=T, legend='top', ncol=2, nrow=2, widths=c(5,1))
annotate_figure(plot2, left=text_grob("Mean Coefficient of Variation", rot=90, size=20), bottom=text_grob("Germination correlation treatment", hjust=0.55, vjust=0, size=20))

## cv together
cvdat<-rbind(j,k)

## risk ratio figure 
p.rr.dat$rr[!is.finite(p.rr.dat$rr)]<-NA
b.rr.dat$rr[!is.finite(b.rr.dat$rr)]<-NA

# for adding a line at "2" 
threshhold<-cbind.data.frame(rep(2, 5), levels(p.rr.dat$Germination))
names(threshhold)<-c("thresh", 'Germination')

# plot 
plant.rr.plot<-ggplot(data=p.rr.dat, aes(x=Treatment, 
                                         y=rr, 
                                         group=Germination))+
  geom_hline(data=threshhold,aes(yintercept=thresh), 
             color='darkgray', lty=2)+
  geom_point(aes(color=Germination), size=2)+
  geom_line(aes(color=Germination))+
  scale_color_manual(values=gcolz)+
  facet_grid(rows=vars(Germination), scales='free')+
  scale_x_discrete(limits=positions, 
                   label=c("-1","-0.5","0","0.5","1", "1sp"))+
  theme_bw()+
  theme(axis.title=element_text(size=15), 
        axis.text=element_text(size=12), 
        axis.text.y=element_text(angle=25),
        legend.title=element_text(size=15), 
        legend.text=element_text(size=12),
        strip.text.y = element_blank(),
        )+
  labs(x="Plant community", y="Risk Ratio")
plant.rr.plot


bee.rr.plot<-ggplot(data=b.rr.dat, aes(x=Treatment, y=rr, group=Germination))+
  geom_hline(data=threshhold,aes(yintercept=thresh), color='darkgray', lty=2)+
  geom_point(aes(color=Germination), size=2)+
  geom_line(aes(color=Germination), lty=1)+
  scale_color_manual(values=gcolz)+
  facet_grid(rows=vars(Germination), scales='free')+
  scale_x_discrete(limits=positions, label=c("-1","-0.5","0","0.5","1", "1sp"))+
  theme_bw()+
  theme(axis.title=element_text(size=15), 
        axis.text=element_text(size=12), 
        axis.text.y=element_text(angle=25),
        legend.title=element_text(size=15), 
        legend.text=element_text(size=12),
        strip.text.y = element_blank(),
  )+
  labs(x="Bee population", y="")
bee.rr.plot

plot3<-ggarrange(plant.rr.plot, bee.rr.plot, common.legend=T, legend='top', ncol=2)
annotate_figure(plot3, bottom=text_grob("Germination correlation treatment", size=20))


# data tables 
# write.csv(corrsdata, file="corrdat.csv")


# ####### plots to keep - as of 20 July 2022########
# # show the correlation covariance relationship
# corcovplot

# ######## PERSISTENCE
# ## show the persistence of different treatments
# annotate_figure(plot1, left=text_grob("Persistence", rot=90, size=20), bottom=text_grob("Germination correlation treatment", size=20, vjust=-1))
# 
# ## Show the risk ratios
# annotate_figure(plot3, bottom=text_grob("Germination correlation treatment", size=20, hjust=0.45))
# 
# ######## STABILITY
# ## Show the coefficients of variation
# annotate_figure(plot2, left=text_grob("Mean Coefficient of Variation", rot=90, size=20), bottom=text_grob("Germination correlation treatment", hjust=0.55, vjust=0, size=20))
# 
# ## Show the cohen's d for CV
# annotate_figure(plot4, left=text_grob("Cohen's d", rot=90, size=20), bottom=text_grob("Germination correlation treatment", size=20))
# 
