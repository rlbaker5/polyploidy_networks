# goal: generate a single script that encompasses all the analyses for an upcoming manuscript.
# note: these analyses were performed years ago and this script is a compilation of multiple separate files. The script is provided as an example and not intended to be stand-alone reproducible code. The original scripts were subject editing including removing output, removing analyses that did not pan out, editing documentation, and obscuring file systems. There is no guarantee that all analyses from the manuscript are included in the script and there is no guarantee that all analyses included in the script are in the manuscript.


#AverageBMT.csv
#This dataset has had outliers removed.
#There are 3 lines of data for each plant - one for the Base, Middle, and Apex ("tip") - BMT.  All data should be the same for each plant except for venation data.
#The "areole area" value is an average of all areole areas for that part of the leav (B/M/T). 
#analysing areole area will likely require using a different dataset that includes each areole's area and slightly more complicated models (see APPS 2020 paper)

rm(list=ls())
gc()
setwd("~")

dat<-read.csv(file="AveragedBMT.csv")

dat1<-dat[,-c(1, 12,14,15,20,21,22,23,24,25,26,27,29, 37,38,46,47,48,49)]

names(dat1)
plot(dat1$Lf4_wetMg, dat1$Lf4_dryMg) #this plot looks really good. There is a nice (close to) 1:1 ratio of wet to dry.

Lf4<-lm(Lf4_dryMg~Lf4_wetMg, data=dat1)

# Call:
#   lm(formula = Lf4_dryMg ~ Lf4_wetMg, data = dat1)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -45.002 -12.322  -1.305  11.996  45.184 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 25.810377   1.906392   13.54   <2e-16 ***
#   Lf4_wetMg    0.059662   0.001056   56.50   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 19.63 on 303 degrees of freedom
# (17 observations deleted due to missingness)
# Multiple R-squared:  0.9133,	Adjusted R-squared:  0.913 
# F-statistic:  3192 on 1 and 303 DF,  p-value: < 2.2e-16

#R-squared = 0.913. That's pretty freaking good.

hist(Lf4$residuals) # looks super normal and centered around zero means this is a good fit.

library(ggplot2)

#https://www.dataquest.io/blog/statistical-learning-for-predictive-modeling-r/
#seems to be a very thorough explanation, unfortunatey all in very obtuse ggplot

pdf(file="Supplemental1_Lf4wet_vs_dry_mass.pdf")
ggplot(data=dat1, aes(x=Lf4_wetMg, y=Lf4_dryMg))+
  geom_point() +
  stat_smooth(method=lm, col="dodgerblue3")+
  theme(panel.background = element_rect(fill="white"), 
        axis.line.x=element_line(),
        axis.line.y=element_line())+
  geom_text(x=1000, y=325, label="Rsquared=0.9133", size=5)+
  ggtitle("Linear model fitted to data")
dev.off()

Lf4_wetMg<-dat1$Lf3_wetMg

Lf3_dryMg_pred<-predict(Lf4, data.frame(Lf4_wetMg))

plot(dat1$Lf3_wetMg, Lf3_dryMg_pred)

min(dat1$Lf4_wetMg, na.rm=TRUE) 88.9
max(dat1$Lf4_wetMg, na.rm=TRUE) 5933.1
min(dat1$Lf3_wetMg, na.rm=TRUE) 352.8
max(dat1$Lf3_wetMg, na.rm=TRUE) 5710.8
min(dat1$Lf4_dryMg, na.rm=TRUE) 12.1
max(dat1$Lf4_dryMg, na.rm=TRUE) 365.6
min(Lf3_dryMg_pred, na.rm=TRUE) 46.85898
max(Lf3_dryMg_pred, na.rm=TRUE) 366.6256

dat1$Lf3_dryMg_pred<-dat1$Lf3_dryMg_pred

#table2.R
#28 May 2020
#Test for the effect of species or parent/hybrid on venation and stomatal traits

#dat: 3 rows per plant: different data for B/M/T, all other data the same. Use for veination except areole area
#dat1: keeps only one of the 3 lines from dat ("Tip"). Used for non-veination analyses
#dat2: 56,000+ lines with many many observations of areole area from each B/M/T location on each plant. Used for areole area analyses.

rm(list=ls())
gc()

library(lme4)
library(emmeans)
install.packages("lsr") #for effect sizes

#Calculate effect size using etaSquared: SStreatment/SStotal

#########################################################
#### load datasets and get rid of NAs as necessary:  ####
#########################################################

dat<-read.csv(file="avgBMT2020.csv") #venation except areole area

nrow(dat) #319
sum(is.na(dat$species)) #2
dat<-dat[-c(which(is.na(dat$species))),]
sum(is.na(dat$species)) 
nrow(dat) #317
319-317 #2

dat1<-dat[which(dat$loc=="T"),] #keep just one copy of each stomatal measurement
nrow(dat1) #107 Good! #LMA and stomata

dat2<-read.csv(file="all_veination_recalib.csv") #areole area
nrow(dat2) #56401
sum(is.na(dat2$species)) #191
dat2<-dat2[-c(which(is.na(dat2$species))),]
nrow(dat2) #56210
56401-56210 #191

############################################
######### Stomatal Traits ##################
############################################

hist(dat1$adaxial_density) #very normal
hist(dat1$abaxial_density) #also quite normal
hist(dat1$ab_ad_density) # a bit left-skewed a fair amount, may want to transform
hist(dat1$abaxial_density-dat$adaxial_density) #pretty normal

####Contrasts #########
str(dat1$species)
levels(dat1$species)
#[1] "carinata" "juncea"   "napus"    "nigra"    "oleracea" "rapa"
#levels 1, 2, 3 (carinata, juncea, napus) are allotetraploid hybrids
#levels 4, 5, 6 (nigra, oleracea, rapa) are diploid parents

#set up contrasts
c1<-c(0,0,0,1,1,1)
mat<-cbind(c1)
contrasts(dat1$species)<-mat

#Run models:
model1<-aov(adaxial_density~species, data=dat1)
#             Df    Sum Sq   Mean Sq F value  Pr(>F)   
#species       5 2.550e-08 5.100e-09   4.355 0.00126 **
#Residuals   100 1.171e-07 1.171e-09                   
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#1 observation deleted due to missingness

#summary(model1, split=list(species=list("Hybrids vs. Parents"=1)))
#                                Df    Sum Sq   Mean Sq F value  Pr(>F)   
#species                          5 2.550e-08 5.100e-09   4.355 0.00126 **
#  species: Hybrids vs. Parents   1 2.920e-09 2.925e-09   2.498 0.11718   
#Residuals                      100 1.171e-07 1.171e-09                   
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#1 observation deleted due to missingness

#EtaSquared.species<-2.550e-8/(2.550e-08+2.92e-09+1.171e-07)
#[1] 0.1751975
#etaSquared(model1) #0.1788194 #from lsr package. Darn close
#EtaSq.cont<-2.920e-09/(2.550e-08+2.92e-09+1.171e-07)
#[1] 0.02006597

model2<-aov(abaxial_density ~ species, data=dat1)
summary(model2)
#             Df    Sum Sq   Mean Sq F value Pr(>F)  
#species       5 2.586e-08 5.173e-09    2.82   0.02 *
#Residuals   100 1.835e-07 1.835e-09                 
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#1 observation deleted due to missingness

summary(model2, split=list(species=list("Hybrids vs Parents"=1)))
#                               Df    Sum Sq   Mean Sq F value Pr(>F)  
#species                         5 2.586e-08 5.173e-09   2.820 0.0200 *
#  species: Hybrids vs Parents   1 5.710e-09 5.709e-09   3.112 0.0808 .
#Residuals                     100 1.835e-07 1.835e-09                 
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#1 observation deleted due to missingness

model3<-aov(ab_ad_density ~ species, data=dat1)
summary(model3)
#             Df Sum Sq Mean Sq F value Pr(>F)  
#species       5  2.615  0.5231   2.735 0.0233 *
#Residuals   100 19.123  0.1912                 
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#1 observation deleted due to missingness

summary(model3, split=list(species=list("Hybrids vs Parents"=1)))
#                               Df Sum Sq Mean Sq F value Pr(>F)  
#species                         5  2.615  0.5231   2.735 0.0233 *
#  species: Hybrids vs Parents   1  0.133  0.1331   0.696 0.4061  
#Residuals                     100 19.123  0.1912                 
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#1 observation deleted due to missingness

model4<-aov(abaxial_density-adaxial_density~ species, data=dat1)
summary(model4)
#             Df    Sum Sq   Mean Sq F value Pr(>F)
#species       5 3.770e-09 7.537e-10   0.795  0.556
#Residuals   100 9.481e-08 9.481e-10               
#1 observation deleted due to missingness

summary(model4, split=list(species=list("Hybrids vs Parents"=1)))
#                               Df    Sum Sq   Mean Sq F value Pr(>F)
#species                         5 3.770e-09 7.537e-10   0.795  0.556
#  species: Hybrids vs Parents   1 4.600e-10 4.612e-10   0.486  0.487
#Residuals                     100 9.481e-08 9.481e-10               
#1 observation deleted due to missingness

model5<-aov(LMA ~ species, data=dat1)
summary(model5)
#            Df Sum Sq Mean Sq F value Pr(>F)
#species      5  3.125  0.6249   1.887  0.104
#Residuals   94 31.124  0.3311               
#7 observations deleted due to missingness

summary(model5, split=list(species=list("Hybrids vs Parents"=1)))
#                              Df Sum Sq Mean Sq F value Pr(>F)
#species                        5  3.125  0.6249   1.887  0.104
#  species: Hybrids vs Parents  1  0.684  0.6839   2.066  0.154
#Residuals                     94 31.124  0.3311               
#7 observations deleted due to missingness

datB<-dat[which(dat$loc=="B"),]
sum(is.na(datB$species)) #0

datM<-dat[which(dat$loc=="M"),]
sum(is.na(datM$species)) #0

datT<-dat[which(dat$loc=="T"),]
sum(is.na(datT$species)) #0

####################BASE
hist(datB$branch_points) #right skewed a bit.
hist(datB$end_points) #better
hist(datB$areole_num) #also right skewed a bit
hist(datB$skel_length_new_mm) #not bad.
hist(datB$vein_dens_mm2) #also right skewed
#None of these are wonderfully normally distributed.  Consider sqrt transformations.

################### MID
hist(datM$branch_points) #good
hist(datM$end_points) #good
hist(datM$areole_num) #good
hist(datM$skel_length_new_mm) #nice
hist(datM$vein_dens_mm2) #right skewed a bit

################### Tip
hist(datT$branch_points) #good
hist(datT$end_points) #good
hist(datT$areole_num) #good
hist(datT$skel_length_new_mm) #nice
hist(datT$vein_dens_mm2) #right skewed a bit

#################################################
################### BASE MODELS #################
#################################################
c1<-c(0,0,0,1,1,1)
mat<-cbind(c1)
contrasts(datB$species)<-mat

modB_vd<-aov(vein_dens_mm2~species, data=datB) 
summary(modB_vd)
#            Df Sum Sq Mean Sq F value Pr(>F)
#species      5  149.7   29.94   1.508  0.194
#Residuals   98 1945.9   19.86               
#2 observations deleted due to missingness

summary(modB_vd, split=list(species=list("Hybrids vs Parents"=1)))
#                              Df Sum Sq Mean Sq F value Pr(>F)  
#species                        5  149.7   29.94   1.508 0.1943  
#  species: Hybrids vs Parents  1   59.1   59.10   2.976 0.0877 .
#Residuals                     98 1945.9   19.86                 
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#2 observations deleted due to missingness

ModB_sl<-aov(skel_length_new_mm ~ species, data=datB)
summary(ModB_sl)
#             Df Sum Sq Mean Sq F value Pr(>F)
#species       5  39460    7892   1.434  0.219
#Residuals   100 550233    5502               

summary(ModB_sl, split=list(species=list("Hybrids vs Parents"=1)))
#                               Df Sum Sq Mean Sq F value Pr(>F)
#species                         5  39460    7892   1.434  0.219
#  species: Hybrids vs Parents   1    333     333   0.061  0.806
#Residuals                     100 550233    5502               

ModB_an<-aov(areole_num ~ species, data=datB)
summary(ModB_an)
#            Df  Sum Sq Mean Sq F value Pr(>F)  
#species      5  177823   35565   2.962 0.0156 *
#Residuals   99 1188668   12007                 
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#1 observation deleted due to missingness

summary(ModB_an, split=list(species=list("Hybrids vs Parents"=1)))
#                               Df  Sum Sq Mean Sq F value Pr(>F)  
# species                        5  177823   35565   2.962 0.0156 *
#   species: Hybrids vs Parents  1    3142    3142   0.262 0.6101  
# Residuals                     99 1188668   12007                 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 1 observation deleted due to missingness

ModB_bp<-aov(branch_points ~ species, data=datB)
summary(ModB_bp)
#              Df  Sum Sq Mean Sq F value  Pr(>F)   
# species       5 1595435  319087   4.054 0.00216 **
# Residuals   100 7870891   78709                   
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(ModB_bp, split=list(species=list("Hybrids vs Parents"=1)))
#                                Df  Sum Sq Mean Sq F value  Pr(>F)   
# species                         5 1595435  319087   4.054 0.00216 **
#   species: Hybrids vs Parents   1   52645   52645   0.669 0.41540   
# Residuals                     100 7870891   78709                   
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

ModB_ep<-aov(end_points ~ species, data=datB)
summary(ModB_ep)
#              Df Sum Sq Mean Sq F value Pr(>F)  
# species       5 108856   21771   3.178 0.0105 *
# Residuals   100 684954    6850                 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(ModB_ep, split=list(species=list("Hybrids vs Parents"=1)))
#                                Df Sum Sq Mean Sq F value Pr(>F)  
# species                         5 108856   21771   3.178 0.0105 *
#   species: Hybrids vs Parents   1   2948    2948   0.430 0.5133  
# Residuals                     100 684954    6850                 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#################################################
################### Middle MODELS #################
#################################################
c1<-c(0,0,0,1,1,1)
mat<-cbind(c1)
contrasts(datM$species)<-mat

ModM_vd<-aov(vein_dens_mm2~species, data=datM) 
#             Df Sum Sq Mean Sq F value  Pr(>F)   
# species      5  219.7   43.95   3.923 0.00281 **
# Residuals   95 1064.3   11.20                   
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 3 observations deleted due to missingness

summary(ModM_vd, split=list(species=list("Hybrids vs Parents"=1)))
#                               Df Sum Sq Mean Sq F value  Pr(>F)   
# species                        5  219.7   43.95   3.923 0.00281 **
#   species: Hybrids vs Parents  1   20.5   20.46   1.826 0.17979   
# Residuals                     95 1064.3   11.20                   
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 3 observations deleted due to missingness

ModM_sl<-aov(skel_length_new_mm ~ species, data=datM)
# summary(ModM_sl)
#             Df Sum Sq Mean Sq F value Pr(>F)
# species      5  29237    5847   1.582  0.172
# Residuals   98 362124    3695               

summary(ModM_sl, split=list(species=list("Hybrids vs Parents"=1)))
#                               Df Sum Sq Mean Sq F value Pr(>F)
# species                        5  29237    5847   1.582  0.172
#   species: Hybrids vs Parents  1   1468    1468   0.397  0.530
# Residuals                     98 362124    3695               

ModM_an<-aov(areole_num ~ species, data=datM)
#             Df Sum Sq Mean Sq F value  Pr(>F)   
# species      5 109932   21986   3.239 0.00957 **
# Residuals   96 651700    6789                   
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 2 observations deleted due to missingness

summary(ModM_an, split=list(species=list("Hybrids vs Parents"=1)))
#                               Df Sum Sq Mean Sq F value  Pr(>F)   
# species                        5 109932   21986   3.239 0.00957 **
#   species: Hybrids vs Parents  1    179     179   0.026 0.87134   
# Residuals                     96 651700    6789                   
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 2 observations deleted due to missingness

ModM_bp<-aov(branch_points ~ species, data=datM)
summary(ModM_bp)
#             Df  Sum Sq Mean Sq F value Pr(>F)   
# species      5  897677  179535   3.753 0.0038 **
# Residuals   96 4592514   47839                  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 2 observations deleted due to missingness

summary(ModM_bp, split=list(species=list("Hybrids vs Parents"=1)))
#                               Df  Sum Sq Mean Sq F value Pr(>F)   
# species                        5  897677  179535   3.753 0.0038 **
#   species: Hybrids vs Parents  1    2451    2451   0.051 0.8214   
# Residuals                     96 4592514   47839                  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 2 observations deleted due to missingness

ModM_ep<-aov(end_points ~ species, data=datM)
# summary(ModM_ep)
#             Df Sum Sq Mean Sq F value   Pr(>F)    
# species      5 132591   26518   4.556 0.000895 ***
# Residuals   97 564558    5820                     
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 1 observation deleted due to missingness

summary(ModM_ep, split=list(species=list("Hybrids vs Parents"=1)))
#                               Df Sum Sq Mean Sq F value   Pr(>F)    
# species                        5 132591   26518   4.556 0.000895 ***
#   species: Hybrids vs Parents  1      3       3   0.001 0.981734    
# Residuals                     97 564558    5820                     
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 1 observation deleted due to missingness

#######################################################
################### Tip (Apex) MODELS #################
#######################################################
c1<-c(0,0,0,1,1,1)
mat<-cbind(c1)
contrasts(datT$species)<-mat

ModT_vd<-aov(vein_dens_mm2~species, data=datT) 
summary(ModT_vd)
#             Df Sum Sq Mean Sq F value Pr(>F)  
# species      5   61.6  12.314   1.964 0.0906 .
# Residuals   99  620.6   6.269                 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 2 observations deleted due to missingness

summary(ModT_vd, split=list(species=list("Hybrids vs Parents"=1)))
#                               Df Sum Sq Mean Sq F value Pr(>F)  
# species                        5   61.6  12.314   1.964 0.0906 .
#   species: Hybrids vs Parents  1    1.7   1.651   0.263 0.6090  
# Residuals                     99  620.6   6.269                 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 2 observations deleted due to missingness

ModT_sl<-aov(skel_length_new_mm ~ species, data=datT)
summary(ModT_sl)
#              Df Sum Sq Mean Sq F value Pr(>F)  
# species       5  28101    5620    2.36 0.0454 *
# Residuals   100 238161    2382                 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 1 observation deleted due to missingness

summary(ModT_sl, split=list(species=list("Hybrids vs Parents"=1)))
#                                Df Sum Sq Mean Sq F value Pr(>F)  
# species                         5  28101    5620   2.360 0.0454 *
#   species: Hybrids vs Parents   1   3391    3391   1.424 0.2356  
# Residuals                     100 238161    2382                 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 1 observation deleted due to missingness


ModT_an<-aov(areole_num ~ species, data=datT)
summary(ModT_an)
#              Df Sum Sq Mean Sq F value Pr(>F)   
# species       5  75343   15069   3.891 0.0029 **
# Residuals   100 387284    3873                  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 1 observation deleted due to missingness

summary(ModT_an, split=list(species=list("Hybrids vs Parents"=1)))
#                                Df Sum Sq Mean Sq F value Pr(>F)   
# species                         5  75343   15069   3.891 0.0029 **
#   species: Hybrids vs Parents   1   7920    7920   2.045 0.1558   
# Residuals                     100 387284    3873                  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 1 observation deleted due to missingness


ModT_bp<-aov(branch_points ~ species, data=datT)
summary(ModT_bp)

summary(ModT_bp, split=list(species=list("Hybrids vs Parents"=1)))
#                                Df  Sum Sq Mean Sq F value   Pr(>F)    
# species                         5  796978  159396   5.100 0.000328 ***
#   species: Hybrids vs Parents   1   76589   76589   2.451 0.120596    
# Residuals                     101 3156359   31251                     
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

ModT_ep<-aov(end_points ~ species, data=datT)
summary(ModT_ep)

summary(ModT_ep, split=list(species=list("Hybrids vs Parents"=1)))
#                                Df Sum Sq Mean Sq F value  Pr(>F)   
# species                         5  88337   17667   3.823 0.00328 **
#   species: Hybrids vs Parents   1    448     448   0.097 0.75631   
# Residuals                     100 462186    4622                   
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 1 observation deleted due to missingness

############################################
########## Areole Area (uses dat2)   #######
############################################

#for reference in table 1 used sqrt transformation:
hist(dat2$areole_area_new_mmsq)
hist(sqrt(dat2$areole_area_new_mmsq))
hist(log(dat2$areole_area_new_mmsq))

#subset data for individual analyses:
dat2B<-dat2[which(dat2$loc=="B"),]
sum(is.na(dat2B$species)) #0
nrow(dat2B) #21307

dat2M<-dat2[which(dat2$loc=="M"),]
sum(is.na(dat2M$species)) #0
nrow(dat2M) #19448

dat2T<-dat2[which(dat2$loc=="T"),]
sum(is.na(dat2T$species)) #0
nrow(dat2T) #15455

#histograms of subset data
hist(dat2B$areole_area_new_mmsq) #holy right skew
hist(sqrt(dat2B$areole_area_new_mmsq)) #better but...
hist(log(dat2B$areole_area_new_mmsq)) #much better!

hist(dat2M$areole_area_new_mmsq) #crazy right skew
hist(sqrt(dat2M$areole_area_new_mmsq)) #much better
hist(log(dat2M$areole_area_new_mmsq)) #better than raw, but sqrt might be best.

hist(dat2T$areole_area_new_mmsq) #pretty right-skewed as well
hist(sqrt(dat2T$areole_area_new_mmsq)) #better. Improves right skew
hist(log(dat2T$areole_area_new_mmsq)) #hard to say.. left skewed now

#set up contrasts
c1<-c(0,0,0,1,1,1)
mat<-cbind(c1)
contrasts(dat2B$species)<-mat

modB1<-aov(sqrt(areole_area_new_mmsq) ~ species, data=dat2B)
# summary(modB1, split=list(species=list("Hybrids vs Parents"=1)))
#                                  Df Sum Sq Mean Sq F value Pr(>F)    
# species                           5   32.2   6.442   213.7 <2e-16 ***
#   species: Hybrids vs Parents     1    5.3   5.346   177.3 <2e-16 ***
# Residuals                     21301  642.2   0.030                   
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

c1<-c(0,0,0,1,1,1)
mat<-cbind(c1)
contrasts(dat2M$species)<-mat

dat2M$sqrt_aa<-sqrt(dat2M$areole_area_new_mmsq)
modM1<-aov(sqrt_aa ~ species, data=dat2M)
summary(modM1, split=list(species=list("Hybrids vs Parents"=1)))
#                                  Df Sum Sq Mean Sq F value   Pr(>F)    
# species                           5   40.7   8.131  244.43  < 2e-16 ***
#   species: Hybrids vs Parents     1    1.0   0.987   29.68 5.17e-08 ***
# Residuals                     19442  646.8   0.033                     
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

etaSquared(modM1)
#             eta.sq eta.sq.part
# species 0.05914316  0.05914316

c1<-c(0,0,0,1,1,1)
mat<-cbind(c1)
contrasts(dat2T$species)<-mat

modT1<-aov(sqrt(areole_area_new_mmsq) ~ species, data=dat2T)
# summary(modT1, split=list(species=list("Hybrids vs Parents"=1)))
#                                  Df Sum Sq Mean Sq F value  Pr(>F)    
# species                           5   32.7   6.548  159.87 < 2e-16 ***
#   species: Hybrids vs Parents     1    0.4   0.420   10.27 0.00136 ** 
# Residuals                     15449  632.8   0.041                    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

###########

rm(list=ls())
gc()
setwd("~")
library(psych)
library(GGally)

#data set with one line for each of B/M/T. Areole area is an average.
dat<-read.csv(file="AveragedBMT.csv") #has areole area (averaged for each location/leaf combo) but not LMA

dat2020<-read.csv("avgBMT2020.csv") #somehow areole area was removed from this one; probably because it's (typically) an unifomormative average. Has LMA

#remove some useless columns
dat<-dat[,c(2:11, 16:18, 28,29,31,35,39:45, 50:65)]
#add LMA
dat$LMA<-dat2020$LMA

#data from just leaf bases
datB<-dat[which(dat$loc=="B"),]
nrow(datB) #106
write.csv(file="DatB.csv", datB, row.names=FALSE)

#data from just leaf middles
datM<-dat[which(dat$loc=="M"),]
rownames(datM)<-c(1:nrow(datM))
write.csv(file="DatM.csv", datM, row.names=FALSE)
nrow(datM) #106

#data from just leaf apices (tips)
datT<-dat[which(dat$loc=="T"),]
rownames(datT)<-c(1:nrow(datT))
nrow(datT) #107
write.csv(file="DatT.csv", datT, row.names=FALSE)

#note that dataframes have different numbers of rows because of some missing data. These dataframes were recombined manually to account for missing data to generate CorsDat.csv

#load dataset for correlation matrices
cordat<-read.csv(file="CorsDat.csv")

nrow(cordat) #110
sum(is.na(cordat$species)) #2
cordat<-cordat[-c(which(is.na(cordat$species))),]
sum(is.na(cordat$species)) #0
nrow(cordat) #108
                                         
#reduce dataframe to only those variables that had either species or ploidy level effects in one-way ANOVAS (current MS, see above, or Baker et al 2017):

cordatred<-cordat[,c(1:19,22:28, 30:33,36,37,41,44,45,52:54)]

full<-corr.test(cordatred[,8:38])
print(full, short=FALSE) # lots of output; run this at your own risk.

corCI.out<-corCi(cordatred[,8:37], p=0.05)

corPlotUpperLowerCi(corCI.out)

cor.plot(cordatred[,8:37], scale=FALSE, cex=0.8, stars=TRUE, diag=FALSE)

cordatred$Category<-as.factor(cordatred$Category)
levels(cordatred$Category)<-strtrim(levels(cordatred$Category), 1)

cordatred$species<-as.factor(cordatred$species)
levels(cordatred$species)<-strtrim(levels(cordatred$species), 3)

cordatred$Category<-as.factor(cordatred$Category)
levels(cordatred$Category)<-strtrim(levels(cordatred$Category), 1)

pdf("parent_hybrid_test.pdf", h=15, w=15)
ggpairs(cordatred, columns=c(8:37), ggplot2::aes(colour=Category, alpha=0.5), cardinality_threshold=110, upper=list(continuous = wrap(cor.plot, r=cordatred[,8:37])), lower=list(continuous=wrap("points", alpha=0.3, size=0.1)))
dev.off()

pdf("species_test.pdf", h=15, w=15)
ggpairs(cordatred, columns=c(8:37), ggplot2::aes(colour=species, alpha=0.5), cardinality_threshold=110, upper=list(continuous = wrap("cor", size = 1)), lower=list(continuous=wrap("points", alpha=0.3, size=0.1)))
dev.off()

printVar = function(x,y){
      vals = corr.test(x,y,
      method="spearman")[c("estimate","p.value")]

      vals[[1]]<-round(vals[[1]],2)   
      vals[[2]]<-ifelse(test = vals[[2]]<0.001,"<0.001",ifelse(test=vals[[2]]<0.01,"<0.01",round(vals[[2]],2)))

          names(vals) = c("rho","p")
      paste(names(vals),unlist(vals),collapse="\n")
}

my_fn <- function(data, mapping, ...){
  # takes in x and y for each panel
  xData <- eval_data_col(data, mapping$x)
  yData <- eval_data_col(data, mapping$y)
  colorData <- eval_data_col(data, mapping$colour)

mainCor = printVar(xData,yData)

p <- ggplot(data = data, mapping = mapping) +
annotate(x=0.5,y=0.8,label=mainCor,geom="text",size=3) +
geom_text(data=byGroup,inherit.aes=FALSE,
aes(x=x,y=y,col=col,label=label),size=3)+ 
theme_void() + ylim(c(0,1))
  p
}

ggpairs(df[,-1],columns = 1:ncol(df[,-1]),
mapping=ggplot2::aes(colour = df$Group),
axisLabels = "show", 
upper = list(continuous = my_fn))+
theme(panel.grid.major = element_blank(), 
panel.grid.minor = element_blank(),
axis.line = element_line(colour = "black")) 

attach(cordat)
dev.new()
pdf(file="baseRattempt.pdf", h=15, w=15)
pairs(cordat[,c(8:26)], diag.panel=panel.hist, lower.panel=rlb_plot, upper.panel=panel.cor_text, gap=0.3)
dev.off()

ggpairs(df[,-1],columns = 1:ncol(df[,-1]),
mapping=ggplot2::aes(colour = df$Group),legends = T,axisLabels = "show", 
upper = list(continuous = wrap("cor", method = "spearman", size = 2.5, hjust=0.7)))+ 
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
axis.line = element_line(colour = "black"))

#https://www.thetopsites.net/projects/ggplot2/ggpairs.shtml
#GGally plots:
pdf("ggallytest.pdf", h=15, w=15)
ggpairs(cordat, columns=c(8:26), ggplot2::aes(colour=Category, alpha=0.5), cardinality_threshold=110)
dev.off()

pdf("parent_hybrid_test.pdf", h=15, w=15)
ggpairs(cordat, columns=c(8:26), ggplot2::aes(colour=Category, alpha=0.5), cardinality_threshold=110, upper=list(continuous = wrap("cor", size = 1)))
dev.off()

pdf("parent_hybrid_test1.pdf", h=15, w=15)
ggpairs(cordat, columns=c(8:26) 
, ggplot2::aes(colour=Category, alpha=0.5), cardinality_threshold=110, upper=list(continuous = wrap("cor", title=NULL, size = 2)))
dev.off()

pdf("species_test.pdf", h=15, w=15)
ggpairs(cordat, columns=c(8:26), ggplot2::aes(colour=species, alpha=0.5), cardinality_threshold=110, upper=list(continuous = wrap("cor", size = 1)))
dev.off()

pdf("test.pdf", h=15, w=15)
ggpairs(cordat, columns=c(8:12), ggplot2::aes(colour=species, alpha=0.5), cardinality_threshold=110, upper=list(continuous = wrap("cor", size = 1)), p.adjust="bonferroni")
dev.off()

pdf("test_corrected.pdf", h=15, w=15) #p.adjust not actually doing anything.
ggpairs(cordat, columns=c(8:12), ggplot2::aes(colour=species, alpha=0.5, p.adjust="bonferroni"), cardinality_threshold=110, upper=list(continuous = wrap("cor", size = 1)))
dev.off()

ggpairs(cordat, columns=c(8:12), ggplot2::aes(colour=Category, alpha=0.5), cardinality_threshold=110, upper=list(continuous = wrap("cor", size = 3)))

#custom GGally correlation matrix with corrections for multiple testing.
## makes figure S1. 

#adapted from:
#https://stackoverflow.com/questions/61686171/how-to-add-the-spearman-correlazion-p-value-along-with-correlation-coefficient-t

install.packages("ggExtra")
library(ggplot2)
library(GGally)
library(ggExtra)
library(data.table)

## read in data
cordat<-read.csv(file="CorsDat.csv")


nrow(cordat) #110
sum(is.na(cordat$species)) #2
cordat<-cordat[-c(which(is.na(cordat$species))),]
sum(is.na(cordat$species)) #0
nrow(cordat) #108

## reduce data set to just the pertinent info:
cordatred<-cordat[,c(1:19,22:28, 30:33,36,37,41,44,45,52:54)]

cordatred$Category<-as.factor(cordatred$Category)
levels(cordatred$Category)<-strtrim(levels(cordatred$Category), 1)

cordatred$species<-as.factor(cordatred$species)
levels(cordatred$species)<-strtrim(levels(cordatred$species), 3)

colnames(cordatred)<-c("Plant","concat","loc","Block", "ID","species","Category","adaxial","abaxial","ab:ad","B ends","M ends","A ends","B branch","M branch","A branch","B areole#","M areole#","A areole#","A skel","B areole area","M areole area","A areole area","B density","M density","A density","Photo","Cond","WUE","Fo","Fv","FvFm","Palisade","Spongey","Palisade:spongy","Leaf area","perimeter","dissection")

output<-data.frame()

printVar = function(x,y){
      vals = cor.test(x,y,
      method="pearson")[c("estimate","p.value")]

      vals[[1]]<-round(vals[[1]],2)   
      vals[[2]]<-p.adjust(vals[[2]], method="holm")
      
      vals[[2]]<-ifelse(test = vals[[2]]<0.001,"<0.001",ifelse(test=vals[[2]]<0.01,"<0.01",round(vals[[2]],2)))

          names(vals) = c("r","p")
      paste(names(vals),unlist(vals),collapse="\n")
}

my_fn <- function(data, mapping, ...){
  # takes in x and y for each panel
  xData <- eval_data_col(data, mapping$x)
  yData <- eval_data_col(data, mapping$y)
  colorData <- eval_data_col(data, mapping$colour)

# if you have colors, split according to color group and calculate cor
  byGroup =by(data.frame(xData,yData),colorData,function(i)printVar(i[,1],i[,2]))
  byGroup = data.frame(col=names(byGroup),label=as.character(byGroup))
  byGroup$x = 0.5
  byGroup$y = c(0.5, 0.1)
  #byGroup$y = seq(0.8-0.3,0.2,length.out=nrow(byGroup))

#main correlation
mainCor = printVar(xData,yData)

p <- ggplot(data = data, mapping = mapping) +
annotate(x=0.5,y=0.9,label=mainCor,geom="text",size=1) +
geom_text(data=byGroup,inherit.aes=FALSE, lineheight=1,
aes(x=x,y=y,col=col,label=label),size=1)+ 
theme_void() + ylim(c(0,1))
  p
}

pdf("S1.ParentHybrid_holm2023.pdf", h=15, w=15)
#### Parent is BLUE; Hybrid is RED
ggpairs(cordatred[,8:38],
mapping=ggplot2::aes(colour = cordatred$Category, alpha=0.5), cardinality_threshold=110,
lower = list(continuous=wrap("points", alpha=0.6, size=0.1)),
upper = list(continuous = my_fn))+
rotateTextX(angle = 90, hjust = 1, vjust = 0.5)+

theme(strip.background = element_rect(fill = "white"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line=element_line(colour="black"), strip.text.x = element_text(size = 4), strip.text.y=element_text(size = 4)
)

dev.off()

####################################
##### Extract data from ggplot #####
####################################
g<-ggpairs(cordatred[,8:37],
mapping=ggplot2::aes(colour = cordatred$Category, alpha=0.5), cardinality_threshold=110,
lower = list(continuous=wrap("points", alpha=0.6, size=0.1)),
upper = list(continuous = my_fn))+
rotateTextX(angle = 90, hjust = 1, vjust = 0.5)+

theme(strip.background = element_rect(fill = "white"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line=element_line(colour="black"), strip.text.x = element_text(size = 4), strip.text.y=element_text(size = 4)
)

library(stringr)
maincors<-data.frame(NULL)
for(m in 1:29){
	for(j in (m+1):30){
		l<-layer_data(g[m,j])
		maincors<-rbind(maincors,l)
		}
	}
#results in a dataframe with 435 correlations coefficients & p values

split<-data.frame(str_split_fixed(maincors$label, " ", 3))
maincors<-cbind(maincors, split)
head(maincors$X3)
sum(maincors$X3==1) #0; if >0, will generate NAs in the next step
maincors$X3<-sub('.', '', maincors$X3) #removes first character from a string, essentially removes the "<" symbol (or anything else). So <0.001 becomes 0.001. >0.05 becomes 0.05 and 0.05 becomes .05. This is somewhat problematic as now it is not clear whether it is <0.05 of 0.05. Luckily there are no instances of "<0.05". So in this case, not a problem:

sum(split$X3=="0.05")
#[1] 5
sum(split$X3=="<0.05")
#[1] 0

maincors$X3<-as.numeric(maincors$X3)
sum(maincors$X3<0.05) #215

ph_cors<-data.frame(NULL)
for(m in 1:29){
	for(j in (m+1):30){
		l<-layer_data(g[m,j], i=2L)
		l$yaxis<-g$yAxisLabel[m] #axis labels of correlations
		l$xaxis<-g$xAxisLabel[j] #axis labels of correlations
		ph_cors<-rbind(ph_cors,l)
	}
}	
#results in a dataframe with 870 correlation coefficients & p values
# color #F8766D is the hybrid (pale orange)
# color #00BFC4 is the parent (pale blue)

head(ph_cors)$label # gives r & p values
split2<-data.frame(str_split_fixed(ph_cors$label, " ",3))
sum(split2$X3=="0.05") #15
sum(split2$X3=="<0.05") #0 noice.
sum(split2$X3=="1") #2; these will generate NAs and will need to be dealt with individually.

ph_cors<-cbind(ph_cors, split2)
ph_cors$X3<-sub('.', '', ph_cors$X3) #gets rid of first character in string (often <)
ph_cors$X3<-as.numeric(ph_cors$X3)

sum(is.na(ph_cors$X3)) #2 These need to become "1".
ph_cors[813,19] #NA
ph_cors[813,19]<-1
ph_cors[813,19] #1

which(is.na(ph_cors$X3)) #106

ph_cors[106,19] #NA
ph_cors[106,19] <-1
ph_cors[106,19] #1

sum(is.na(ph_cors$X3)) #0

parents<-ph_cors[which(ph_cors$colour=='#00BFC4'),]
hybrids<-ph_cors[which(ph_cors$colour=='#F8766D'),]
nrow(parents)#435
nrow(hybrids)#435

colnames(parents)[19]<-"parent.p"
colnames(hybrids)[19]<-"hybrid.p"

p.values<-data.frame(parents$yaxis, parents$xaxis, parents$parent.p, hybrids$hybrid.p)
colnames(p.values)<-c("yaxis", "xaxis", "parent.p", "hybrid.p")
p.values$all.p<-maincors$X3

p.values$parent.p<-as.numeric(p.values$parent.p)
p.values$hybrid.p<-as.numeric(p.values$hybrid.p)

sum(p.values$parent<0.05) #181
sum(p.values$hybrid<0.05) #204

p.values$h<-p.values$hybrid<0.05
p.values$p<-p.values$parent<0.05
p.values$a<-p.values$all.p<0.05
p.values$sum<-p.values$h+p.values$p+p.values$a

head(p.values)

sum(p.values$sum==3) #166 - the number of correlations where P&H sig - no change in pheno integration

hsig<-p.values[which(p.values$hybrid.p<0.05 & p.values$parent.p>=0.05),]
nrow(hsig) #38; the number where parent is nonsig and hybrid is sig (gain of pheno integration)
psig<-p.values[which(p.values$hybrid.p>=0.05 & p.values$parent.p<0.05),]
nrow(psig) #15; the number where parent is sig and hybrid is non-sig (loss of pheno integration)

435-38-15 #382 - the number of correlations where there was no state change, irregardless of whehther they were sig or non-sig.


###############################################################
###### determine which (if any) traits are over-represented in correlations where phenotypic integratin changes "status" between diploid parents and allotetraploid hybrids ###########
###############################################################

############### find the frequency that each trait participates in a state-change, convert to Z-scale, and then find sig. difference.

### Issues: 1) clearly non-normal distribution (so did a log transformation, but still not awesome) 2) performed a one-tailed test but not super happy with it.

sigplayers<-c(unlist(hsig$yaxis), unlist(hsig$xaxis), unlist(psig$yaxis), unlist(psig$xaxis)) #list of all traits involved in state changes and the number of times they are listed corresponds to the number of state changes they are involved in.

length(sigplayers) #106
tab<-data.frame(table(sigplayers))
tab
#         sigplayers Freq
# 1    A areole area    2
# 2        A areole#    2
# 3         A branch    2
# 4        A density    5
# 5           A ends    5
# 6           A skel    2
# 7            ab:ad    7
# 8          abaxial    3
# 9          adaxial    2
# 10   B areole area    3
# 11       B areole#    1
# 12        B branch    1
# 13       B density    3
# 14          B ends    1
# 15            Cond    3
# 16              Fo    1
# 17              Fv    1
# 18            FvFm    4
# 19       Leaf area   15
# 20   M areole area    3
# 21       M areole#    4
# 22        M branch    3
# 23       M density    2
# 24          M ends    3
# 25        Palisade    8
# 26 Palisade:spongy    4
# 27       perimeter    2
# 28           Photo    3
# 29         Spongey   10
# 30             WUE    1

# interesting that all traits were involved in at least 1 state change.

hist(tab$Freq) # highly left-skewed
dev.new()
hist(log(tab$Freq))#better
dev.new()
hist(sqrt(tab$Freq))#nope.

tab$log<-log(tab$Freq) #log transform improves normality...still not normal.
sd(tab$log) #[1] 0.7095751 
mean(tab$log) #[1] 1.001822

tab$Z<-( (tab$log-1.001822)/0.7095751 ) #compute Z-score
#anything with Z>1.645 is significant
	#Leaf area and spongy mesophyll (area?) are significantly over represented in state-shifts. 

############ Chi-squared tests ###############

#38+15=53 total state changes. If equally likely expect 26.5 state changes.

#### is one type of state change (loss/gain) favored over another? YES.
statechange<-c(38, 15)
res<-chisq.test(statechange, p=c(1/2, 1/2))
res

# 	Chi-squared test for given probabilities
# 
# data:  statechange
# X-squared = 9.9811, df = 1, p-value = 0.001582

#signfiicant deviation from expected frequency.

###### is stable integration status favored over changes? YES.
216+166 #382 number of integration states stable across polyploidy
15+38 #53 number of integration states that change
int<-c(382,53)
res<-chisq.test(int, p=c(1/2, 1/2))
res

# 	Chi-squared test for given probabilities
# 
# data:  int
# X-squared = 248.83, df = 1, p-value < 2.2e-16

#reject Ho.

### is a given trait over represented in terms of significant correlations that switch phenotypic integration status?
#435 possible correlations
#53 state changes
53/435 #0.1218391 So would expect on average a 12 percent chance of a status change.

#there are 30 traits, each participating in 29 potential bivariate correlations.
29*.12 #3.48. Expect on average that each trait participates in 3.48 correlations that change state.

29-3.48 #25.52 and 25.52 that do not.

##### any trait participating in more than 7 bivariate correlations with state changes can be considered as "over represented"

#generic chi_square test:
overrep<-c(observed_uncorrelated, observed_correlated)
res<-chisq.test(overrep, p=c(382/435, 53/435))
res

#for "adaxial" (stomata) with 2 sig cors.
overrep<-c(27, 2)
res<-chisq.test(overrep, p=c(382/435, 53/435))
res
#X-squared = 0.75773, df = 1, p-value = 0.384

# for "M areole #" (4 sig cors)
overrep<-c(25, 4)
res<-chisq.test(overrep, p=c(382/435, 53/435))
res
#X-squared = 0.070187, df = 1, p-value = 0.7911

#for "A ends" (5 sig cors)
overrep<-c(24, 5)
res<-chisq.test(overrep, p=c(382/435, 53/435))
res
# X-squared = 0.69327, df = 1, p-value = 0.4051

overrep<-c(22, 7)
res<-chisq.test(overrep, p=c(382/435, 53/435))
res

#based on chi-squared tests, the following traits are "over represented": ab:ad axial stomata; Leaf aera; Palisade parenchyma; spongy mesophyll.

#################################################
#################################################

#generates a mutual rank co-expression trait network.
#generates networks based on focal traits; uses both diploid and tetraploid data

library(igraph)
library(tidyverse)
library(netranker)
library(corrr)
library(stringr)

#Network tools: for quantifying differences in networks:
#https://github.com/galanisl/DisparityFilter
#https://github.com/alextkalinka/linkcomm

cordat<-read.csv(file="CorsDat.csv")

nrow(cordat) #110
sum(is.na(cordat$species)) #2
cordat<-cordat[-c(which(is.na(cordat$species))),]
sum(is.na(cordat$species)) #0
nrow(cordat) #108

## reduce data set to just the pertinent info:
cordatred<-cordat[,c(1:19,22:28, 30:33,36,37,41,44,45,52:54)]

cordatred$Category<-as.factor(cordatred$Category)
levels(cordatred$Category)<-strtrim(levels(cordatred$Category), 1)

cordatred$species<-as.factor(cordatred$species)
levels(cordatred$species)<-strtrim(levels(cordatred$species), 3)

colnames(cordatred)<-c("Plant","concat","loc","Block", "ID","species","Category","adaxial","abaxial","ab:ad","B ends","M ends","A ends","B branch","M branch","A branch","B areole#","M areole#","A areole#","A skel","B areole area","M areole area","A areole area","B density","M density","A density","Photo","Cond","WUE","Fo","Fv","FvFm","Palisade","Spongey","Palisade:spongy","Leaf area","perimeter","dissection")

#function to calculate mutual ranks:
mutrank.wrap <- function(data){
    m <- cor(data, use= "pairwise.complete.obs", method = "pearson")
    r <- apply(m,1,rank)
    diag(r) <- 0
    net <- r*t(r)/2
    colnames(net) <- colnames(data)
    rownames(net) <- colnames(data)
    return(net)
}

#calculate mutual ranks:
#all.MR<-mutrank.wrap(cordatred[,8:38])

#rescale mutual ranks:
#MRrescaled<-apply(all.MR, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))

#re-arrange the dataframe, clean up names, etc:
rcors<-cordatred[,8:38]
rcors<-cbind(rcors[,c(3,29,26,27)], rcors[,-c(3,29,26,27)])
colnames(rcors)<-c("ab_ad", "Leaf_area", "Palisade", "Spongey", "adaxial", "abaxial", "B_ends", "M_ends", "A_ends", "B_branch", "M_branch", "A_branch", "B_areole_num", "M_areole_num", "A_areole_num", "A_skel", "B_areole_area", "M_areole_area", "A_areole_area", "B_density", "M_density", "A_density", "Photo", "Cond", "WUE", "Fo", "Fv", "FvFm", "Palisade_spongy", "perimeter", "dissection")      
rcors.MR<-mutrank.wrap(rcors)

#get edgelist from a mutual rank correlation matrix:
g  <- graph.adjacency(rcors.MR,weighted=TRUE)
edges <- get.data.frame(g)
nrow(edges) #930

#Make a nodelist
nodes<-data.frame(colnames(rcors.MR))
nrow(nodes)

#add in variables to nodes (over vs. underrpresented traits))
represented<-c(replicate(4, "tomato"), replicate(27, "grey")) #tomato is over-represented; grey is not over represented
length(represented) #31
nodes<-cbind(nodes, represented)
colnames(nodes)<-c("Trait", "Influence")

nrow(nodes); length(unique(nodes$Trait)) #31, 31. good
nrow(edges); nrow(unique(edges[,c("from", "to")])) #930, 930 - good. If there were more links than unique two-from combinations they wuld need to be collapsed

links<-edges
net<-graph_from_data_frame(d=links, vertices=nodes, directed=F) #F removes arrowheads

colrs <- c("tomato", "grey") #needed for fig. legend.
V(net)$color<-V(net)$Influence

deg<-degree(net, mode="all")
V(net)$size<-deg/5
#V(net)$label<-NA
E(net)$width<-E(net)$weight/100
#E(net)$arrow.size<-.2
E(net)$edge.color<-"gray80"
E(net)$width<-1+E(net)$weight/100

set.seed(1234)
pdf("full_net.pdf")
plot(net)
legend(x=-1.5, y=-1.1, c("Over represented", "Not over represented"), pch=21, col="#777777", pt.bg=colrs, pt.cex=2, cex=0.8, bty="n", ncol=1)
dev.off()

#delete edges with less than mean weights:
set.seed(1234)
pdf("More_than_mean_vertices.pdf")
cut.off <- mean(links$weight) 
net.sp <- delete_edges(net, E(net)[weight<cut.off])

plot(net.sp, layout=layout_with_lgl)
legend(x=-1.5, y=-1.1, c("Over represented", "Not over represented"), pch=21, col="#777777", pt.bg=colrs, pt.cex=2, cex=0.8, bty="n", ncol=1)
dev.off()

#list all edges:
attr(E(net), "vnames") #930 instance
length(attr(E(net), "vnames")) #930 yup.

#get rid of duplicated edges:
to.delete2<-which(duplicated(attr(E(net), "vnames")) %in% "TRUE")
net<-delete.edges(net, to.delete2)
length(attr(E(net), "vnames")) #465. good.

#list all elements to be deleted
named_list<-str_starts(attr(E(net), "vnames"), "ab_ad\\||Leaf_area\\||Palisade\\||Spongy\\|", negate=FALSE)

#get vector locations of elements to be deleted
to.delete<-which(named_list %in% "FALSE")

#delete all non-focal vertices:
focus.net<-delete.edges(net, to.delete)

attr(E(focus.net), "vnames") #87 instances

#graph focal vertices (and all things connected to them):
set.seed(1324)
pdf("focus_all.pdf")
plot(focus.net, layout=layout_with_lgl)
legend(x=-1.5, y=-1.1, c("Over represented", "Not over represented"), pch=21, col="#777777", pt.bg=colrs, pt.cex=2, cex=0.8, bty="n", ncol=1)
dev.off()

#min and max.. not super useful.
max(E(focus.net)$weight) #450
min(E(focus.net)$weight) #0.5

#Find value of top XXXX edges:
E(focus.net)$weight[order(-E(focus.net)$weight)]
#378 is the 10th number

focus.net.10<-delete.edges(focus.net, E(focus.net)[weight<378]) #delete all vertices with wieght less than 378 (just use top 10 vertices)

set.seed(1234)
pdf("focus.net.10.pdf")
plot(delete.vertices(focus.net.10, degree(focus.net.10)==0), layout=layout_with_lgl) # plot and delete orphane nodes.
legend(x=-1.5, y=-1.1, c("Over represented", "Not over represented"), pch=21, col="#777777", pt.bg=colrs, pt.cex=2, cex=0.8, bty="n", ncol=1)
dev.off()

focus.net.20<-delete.edges(focus.net, E(focus.net)[weight<264])

set.seed(1234)
pdf("focus.net.20.pdf")
plot(delete.vertices(focus.net.20, degree(focus.net.20)==0), layout=layout_with_lgl)
legend(x=-1.5, y=-1.1, c("Over represented", "Not over represented"), pch=21, col="#777777", pt.bg=colrs, pt.cex=2, cex=0.8, bty="n", ncol=1)
dev.off()


focus.net.50<-delete.edges(focus.net, E(focus.net)[weight<26])

set.seed(1234)
pdf("focus.net.50.pdf")
plot(delete.vertices(focus.net.50, degree(focus.net.50)==0), layout=layout_with_lgl)
legend(x=-1.5, y=-1.1, c("Over represented", "Not over represented"), pch=21, col="#777777", pt.bg=colrs, pt.cex=2, cex=0.8, bty="n", ncol=1)
dev.off()

set.seed(1234)
pdf("focus.net.50.grid.pdf")
plot(delete.vertices(focus.net.50, degree(focus.net.50)==0), layout=layout_on_grid)
legend(x=-1.5, y=-1.1, c("Over represented", "Not over represented"), pch=21, col="#777777", pt.bg=colrs, pt.cex=2, cex=0.8, bty="n", ncol=1)
dev.off()

focus.connectivity<-vertex_connectivity(focus.net)

#Make network plots separately for diploids and tetraploids:
rm(list=ls())
library(tidyverse) #tidy functions
library(igraph) #graph theory
library(corrr) #correlations
library(netranker) #???
library(stringr) #regular expressions
library(graph4lg) #percolation pruning

#load data
cordat<-read.csv(file="CorsDat.csv")

nrow(cordat) #110
sum(is.na(cordat$species)) #2
cordat<-cordat[-c(which(is.na(cordat$species))),]
sum(is.na(cordat$species)) #0
nrow(cordat) #108


#reduce data frame to only traits that had significant species or ploidy effecgts (plus metadata)
cordatred<-cordat[,c(1:19,22:28, 30:33,36,37,41,44,45,52:54)]

cordatred$Category<-as.factor(cordatred$Category)
levels(cordatred$Category)<-strtrim(levels(cordatred$Category), 1)

cordatred$species<-as.factor(cordatred$species)
levels(cordatred$species)<-strtrim(levels(cordatred$species), 3)

colnames(cordatred)<-c("Plant","concat","loc","Block", "ID","species","Category","adaxial","abaxial","ab_ad","B_ends","M_ends","A_ends","B_branch","M_branch","A_branch","B_areole_num","M_areole_num","A_areole_num","A_skel","B_areole_area","M_areole_area","A_areole_area","B_density","M_density","A_density","Photo","Cond","WUE","Fo","Fv","FvFm","Palisade","Spongy","Palisade_spongy","Leaf_area","perimeter","dissection")

#split into two dataframes: one for diploids and one for tetraploids:
diploid<-cordatred[which(cordatred$Category=="p"),]
head(diploid)

tetraploid<-cordatred[which(cordatred$Category=="h"),]

#strip metadata; re-arrange columns so that the 4 overrepresented traits are at the front:
diploid<-cbind(diploid[,c(10,36,33,34)], diploid[c(8,9,11:32,35,37,38)])

tetraploid<-cbind(tetraploid[,c(10,36,33,34)], tetraploid[c(8,9,11:32,35,37,38)])

#function to calculate mutual rank correlation coefficients (based on Bennot et al. 2015, R package called 'netbenchmark'. Although note that the package did not install/compile so I had to go to github and pull out the code to run it independently.)

#changes to the function: added the parameter "use = "pairwise.complete.obs" to the cor function
mutrank.wrap <- function(data){
    m <- cor(data, use= "pairwise.complete.obs", method = "pearson")
    r <- apply(m,1,rank)
    diag(r) <- 0
    net <- r*t(r)/2
    colnames(net) <- colnames(data)
    rownames(net) <- colnames(data)
    return(net)
}

#generate mutual rank matrices:
dip.MR<-mutrank.wrap(diploid)
tet.MR<-mutrank.wrap(tetraploid)

#get edgelist from MR matrix:
dip.g  <- graph.adjacency(dip.MR,weighted=TRUE)
dip.links <- get.data.frame(dip.g)
nrow(dip.links) #930

tet.g <- graph.adjacency(tet.MR, weighted=TRUE)
tet.links<-get.data.frame(tet.g)
nrow(tet.links) #930

#Make a nodelist
dip.nodes<-data.frame(colnames(dip.MR))
nrow(dip.nodes) #31
tet.nodes<-data.frame(colnames(tet.MR))
nrow(tet.nodes) #31

#make networks:
dip.net<-graph_from_data_frame(d=dip.links, vertices=dip.nodes, directed=F) #directed=F removes arrow heads on edges

tet.net<-graph_from_data_frame(d=tet.links, vertices=tet.nodes, directed=F) #directed=F removes arrow heads on edges

# get rid of duplicated edges
dip.to.delete<-which(duplicated(attr(E(dip.net), "vnames")) %in% "TRUE")
dip.net<-delete.edges(dip.net, dip.to.delete)
length(attr(E(dip.net), "vnames")) #465. good.

tet.to.delete<-which(duplicated(attr(E(tet.net), "vnames")) %in% "TRUE")
tet.net<-delete.edges(tet.net, tet.to.delete)
length(attr(E(tet.net), "vnames")) #465. good.

#plot full networks:
set.seed(123)
dev.new()
pdf("full_dip_network.pdf", h=15, w=15)
plot(dip.net, layout=layout_in_circle, vertex.color="grey")
dev.off()

set.seed(123)
dev.new()
pdf("full_tet_network.pdf", h=15, w=15)
plot(tet.net, layout=layout_in_circle, vertex.color="grey")
dev.off()


#First prune all links > mean weight
mean(E(dip.net)$weight) #148.3366
dip.cut<-mean(E(dip.net)$weight)
dip.sp<-delete_edges(dip.net, E(dip.net)[weight<dip.cut])

tet.cut<-mean(E(dip.net)$weight)
tet.sp<-delete_edges(tet.net, E(tet.net)[weight<tet.cut])

#make adjacency matrices with only values > mean
dip.mean.mat<-ifelse(dip.MR<mean(E(dip.net)$weight), 0, dip.MR) #setting to NA doesn't seem to work
tet.mean.mat<-ifelse(tet.MR<mean(E(tet.net)$weight), 0, tet.MR)

#percol matrices with just MR>mean value
dip.mean.percol<-g_percol(dip.mean.mat)
	#Number of conserved links : 73
	#Maximum weight of the conserved links : 220
plot(dip.mean.percol, layout=layout_with_lgl)
plot(dip.mean.percol, layout=layout_in_circle)

tet.mean.percol<-g_percol(tet.mean.mat)
	#Number of conserved links : 123
	#Maximum weight of the conserved links : 300
plot(tet.mean.percol, layout=layout_with_lgl)
plot(tet.mean.percol, layout=layout_in_circle)


#Edge width based on edge weight:
E(dip.mean.percol)$width<-E(dip.mean.percol)$weight/100
E(dip.mean.percol)$edge.color<-"gray80"
E(dip.mean.percol)$width<-1+E(dip.mean.percol)$weight/75
plot(dip.mean.percol)

#community assessment:
gg.com<-cluster_fast_greedy(dip.mean.percol)
V(dip.mean.percol)$color<-gg.com$membership+1
E(dip.mean.percol)$width<-1+E(dip.mean.percol)$weight/50
plot(dip.mean.percol, layout=layout_in_circle)

gg.com<-cluster_fast_greedy(tet.mean.percol)
V(tet.mean.percol)$color<-gg.com$membership+1
E(tet.mean.percol)$width<-1+E(tet.mean.percol)$weight/50
#plot(dip.mean.percol, layout=layout_in_circle)

set.seed(138)
pdf("dip.mean.percol.pdf", 20,20)
plot(dip.mean.percol, layout=layout_with_lgl, vertex.label.font=2, vertex.label.cex=2)
dev.off()

set.seed(123)
pdf("tet.mean.percol.pdf", 20,20)
plot(tet.mean.percol, layout=layout_with_lgl, vertex.label.font=2, vertex.label.cex=2)
dev.off()


#mean percol closeness:
##### DIP: 29 traits
dip.mean.percol.close<-data.frame(closeness(dip.mean.percol, mode="all"))
dip.mean.percol.close$trait<-rownames(dip.mean.percol.close)
colnames(dip.mean.percol.close)<-c("dip.mean.percol.close", "trait")
#tet.close$rank<-rank(tet.close$closeness) # column of rank order closeness
#tet.close<-tet.close[order(-tet.close$rank),] # reorder based on rank order closeness
dip.mean.percol.close$dip.mean.percol.rank<-rank(dip.mean.percol.close$dip.mean.percol.close)
dip.mean.percol.close<-cbind(dip.mean.percol.close[,c(2,1,3)]) #29 traits

####TET: 31 traits
tet.mean.percol.close<-data.frame(closeness(tet.mean.percol, mode="all"))
tet.mean.percol.close$trait<-rownames(tet.mean.percol.close)
colnames(tet.mean.percol.close)<-c("tet.mean.percol.close", "trait")
tet.mean.percol.close$tet.mean.percol.rank<-rank(tet.mean.percol.close$tet.mean.percol.close)
tet.mean.percol.close<-cbind(tet.mean.percol.close[,c(2,1,3)]) #31 traits

####### Find critical cutoffs for edges:
dist.dip<-(E(dip.net)$weight)
hist(dist.dip) #very right-skewed
hist(sqrt(dist.dip)) #better.
p.dip<-quantile(sqrt(dist.dip), probs=c(0.99, 0.95, 0.90, 0.75))
p.dip*p.dip
#  99%   95%   90%   75% 
#450.0 406.0 362.5 238.0 <-- use these values for critical cutoffs.

dist.tet<-(E(tet.net)$weight)
hist(dist.tet) #looks similar to dist.dip, which is good I guess.
hist(sqrt(dist.tet)) #better. kinda.
p.tet<-quantile(sqrt(dist.tet), probs=c(0.99, 0.95, 0.90, 0.75))
p.tet*p.tet
# 99% 95% 90% 75% 
# 450 405 351 250  <-- use these values for critical cutoffs

###### 95% confidence network ****** DO NOT USE 95 % NETWORKS *********
dip.net.95<-delete.edges(dip.net, E(dip.net)[weight<406]) 
tet.net.95<-delete.edges(tet.net, E(tet.net)[weight<405])

##### retain largest networks
decomp.dip.95<-decompose(dip.net.95) #returns 9 networks
plot(decomp.dip.95[[1]]) # 15 node network interestingly retains all 4 "key traits"
closeness(decomp.dip.95[[1]])

decomp.tet.95<-decompose(tet.net.95) #returns 10 networks
plot(decomp.tet.95[[1]]) #4 nodes includes ab_ad
plot(decomp.tet.95[[2]]) #4 nodes includes Leaf area, spongy, palisade
plot(decomp.tet.95[[3]]) #2 nodes, abaxial and adaxial
plot(decomp.tet.95[[4]]) #3 nodes
plot(decomp.tet.95[[5]]) #3 nodes
plot(decomp.tet.95[[6]]) #3 nodes
plot(decomp.tet.95[[7]]) #1 node
plot(decomp.tet.95[[8]]) #3 nodes
plot(decomp.tet.95[[9]]) #7 nodes
plot(decomp.tet.95[[10]]) #1 node
 
##### 95% confidence networks
dip.net.90<-delete.edges(dip.net, E(dip.net)[weight<362.5]) 
tet.net.90<-delete.edges(tet.net, E(tet.net)[weight<351])

dip.comp.90<-decompose(dip.net.90) #3 networks
plot(dip.comp.90[[1]]) #16 nodes contains all "key traits"
#network 2: 2 traits
#network 3: 13 traits

tet.comp.90<-decompose(tet.net.90) #5 netowrks
plot(tet.comp.90[[1]]) #14 nodes, contains all "key traits"
#network2: 2 traits; network3: 12 traits; network4: 1 trait; network5: 1 trait

#Closenes for 90% networks:
dip.90.close<-data.frame(closeness(dip.comp.90[[1]], mode="all"))
dip.90.close$trait<-rownames(dip.90.close)
colnames(dip.90.close)<-c("dip.90.close", "trait")
dip.90.close$dip.90.close.rank<-rank(dip.90.close$dip.90.close)
dip.90.close<-cbind(dip.90.close[,c(2,1,3)]) #16 traits

tet.90.close<-data.frame(closeness(tet.comp.90[[1]], mode="all"))
tet.90.close$trait<-rownames(tet.90.close)
colnames(tet.90.close)<-c("tet.90.close", "trait")
tet.90.close$tet.90.close.rank<-rank(tet.90.close$tet.90.close)
tet.90.close<-cbind(tet.90.close[,c(2,1,3)]) #16 traits


### community assessment
#Edge width based on edge weight:
dip.90.width<-E(dip.comp.90[[1]])$weight
E(dip.comp.90[[1]])$width<-dip.90.width
E(dip.comp.90[[1]])$edge.color<-"gray80"
E(dip.comp.90[[1]])$width<-1+E(dip.comp.90[[1]])$weight/100
plot(dip.comp.90[[1]])

#community assessment:
gg.com<-cluster_fast_greedy(dip.comp.90[[1]])
V(dip.comp.90[[1]])$color<-gg.com$membership+1

set.seed(126)
plot(dip.comp.90[[1]], layout=layout_with_lgl,vertex.label.font=2, vertex.label.cex=2)

set.seed(127)
pdf("dip.90.1.comm.pdf", 20,20)
plot(dip.comp.90[[1]], layout=layout_with_lgl,vertex.label.font=2, vertex.label.cex=2)
dev.off()

tet.90.width<-E(tet.comp.90[[1]])$weight
E(tet.comp.90[[1]])$width<-tet.90.width
E(tet.comp.90[[1]])$edge.color<-"gray80"
E(tet.comp.90[[1]])$width<-1+E(tet.comp.90[[1]])$weight/100

gg.com<-cluster_fast_greedy(tet.comp.90[[1]])

#re-order community membership to better coincide diploids:
gg.com$membership
#[1] 2 1 2 2 2 2 2 3 3 3 3 3 3 1 1
mem<-c(3,1,3,3,3,3,3,2,2,2,2,2,2,1,1)
gg.com$membership<-mem
V(tet.comp.90[[1]])$color<-gg.com$membership+1

set.seed(123)
pdf("tet.90.1.com.pdf", 20,20)
plot(tet.comp.90[[1]], layout=layout_with_lgl,vertex.label.font=2, vertex.label.cex=2)
dev.off()

#s_percol_function.R:

s_percol <- function(x, val_step = 20){

# Check whether x is a symmetric matrix
  if(!inherits(x, c("matrix", "dist"))){
    stop("'x' must be a matrix or a dist object")
  }

  if(inherits(x, "matrix")){
    if(!isSymmetric(x)){
      stop("The matrix 'x' must be symmetric")
    }
  }
  if(inherits(x, "dist")){
    x <- as.matrix(x)
  }

  # Creation of the complete initial graph
  g1 <- igraph::graph.adjacency(x,
  mode = "undirected",
  weighted = TRUE,
  diag =FALSE)

  # Edge list of the complete graph
  g1_df <- data.frame(igraph::as_edgelist(g1))

  # We add the weight of each link in the edge list data frame
  g1_df$w <- igraph::E(g1)$weight

  # We order the df with the SMALLEST weights first
  g1_df <- g1_df[order(g1_df$w, decreasing=FALSE), ]

  # We calculate the limits of the classes (id_sup)
  val_part <- round(nrow(g1_df)/val_step, digits = 0)

  # We create id_sup. Its first element is1
  id_sup <- c(1)
  for (i in 1:(val_step-1)){
    id_sup[i+1] <- val_part*i
  }
  id_sup[val_step] <- nrow(g1_df)-1

  # The first value of id_sup is 1, the last is the number of edges - 1

  # We calculate the weight of the link at the class limits (threshold)
  threshold <- g1_df[id_sup, 'w']

  # We will find in which class, weights become so small that the graph
  # breaks into two components
  # If we remove the links with a big weight, the graph stays connected
  # logically.
  #Once it breaks into two components, it means that the threshold is
  # in the last class before the value
  comp <- 1
  i <- 1
  while(comp < 2){
    mat_g2 <- x
    mat_g2[mat_g2 > threshold[i]] <- 0
    g2 <- igraph::graph.adjacency(mat_g2,
    mode = "undirected",
    weighted = TRUE,
    diag =FALSE)
    comp <- igraph::components(g2)$no[1]
    i <- i + 1
  }
  t1 <- i - 1

  # t1 is the first iteration for which we get two components
  # So threshold[t1] is lower than the real threshold
  # The real threshold is between the id_sup[t1-1]-th and the id_sup[t1]-th
  # rows of g1_df
  # We reiterate the operation in the identified class
  val_step2 <- val_step/2
  val_part2 <- round(val_part/val_step2, digits = 0)
  id_sup2 <- c(id_sup[t1 - 1])

  for (i in 1:(val_step2-1)){
    id_sup2[i+1] <- id_sup[t1-1]+val_part2*i
  }
  id_sup2[val_step2] <- id_sup[t1]

  # We calculate the weight of the link at the class limits (threshold)
  threshold2 <- g1_df[id_sup2, 'w']

  # We will find in which class, weights become so small that the graph
  # breaks into two components
  # If we remove the links with a big weight, the graph stays connected
  # logically. Once it breaks into two components, it means that the threshold is
  # in the last class before the value
  comp <- 1
  i <- 1
  while(comp < 2){
    mat_g2 <- x
    mat_g2[mat_g2 > threshold2[i]] <- 0
    g2 <- igraph::graph.adjacency(mat_g2,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE)
    comp <- igraph::components(g2)$no[1]
    i <- i + 1
  }
  t2 <- i - 1

# The threshold is the weight of one of the link
# between the rows id_sup2[t2 - 1] and id_sup2[t2]
# of g1_df
  id_sup3 <- c(id_sup2[t2-1]:id_sup2[t2])
  threshold3 <- g1_df[id_sup3, 'w']
  comp <- 1
  i <- 1
  while(comp < 2){
    mat_g2 <- x
    mat_g2[mat_g2 > threshold3[i]] <- 0
    g2 <- igraph::graph.adjacency(mat_g2,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE)
    comp <- igraph::components(g2)$no[1]
    i <- i + 1
  }
  t3 <- i - 1
  thr2 <- g1_df[id_sup3[t3], 'w']
  mat_g2 <- x
  mat_g2[mat_g2 > thr2] <- 0
  g2 <- igraph::graph.adjacency(mat_g2,
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE)
  comp2 <- igraph::components(g2)$no[1]
  thr <- g1_df[id_sup3[t3]-1, 'w']
  mat_g_end <- x
  mat_g_end[mat_g_end > thr] <- 0
  g_end <- igraph::graph.adjacency(mat_g_end,
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE)
  comp1 <- igraph::components(g_end)$no[1]

  if(sum(comp1, comp2) > 3){
    message("There are probably equal link weights.")
  }

  message(paste("Number of conserved links :",
  length(igraph::E(g_end))), sep = "")
  message(paste("Maximum weight of the conserved links :",
  max(igraph::E(g_end)$weight)), sep = "")
  return(g_end)
}

#Make neighborhood networks based on "key traits":

rm(list=ls())
library(tidyverse) #tidy functions
library(igraph) #graph theory
library(corrr) #correlations
library(netranker) #???
library(stringr) #regular expressions
library(graph4lg) #percolation pruning

#load data

cordat<-read.csv(file="CorsDat.csv")

nrow(cordat) #110
sum(is.na(cordat$species)) #2
cordat<-cordat[-c(which(is.na(cordat$species))),]
sum(is.na(cordat$species)) #0
nrow(cordat) #108

#reduce data frame to only traits that had significant species or ploidy effecgts (plus metadata)
cordatred<-cordat[,c(1:19,22:28, 30:33,36,37,41,44,45,52:54)]
cordatred$Category<-as.factor(cordatred$Category)
levels(cordatred$Category)<-strtrim(levels(cordatred$Category), 1)
cordatred$species<-as.factor(cordatred$species)
levels(cordatred$species)<-strtrim(levels(cordatred$species), 3)
colnames(cordatred)<-c("Plant","concat","loc","Block", "ID","species","Category","adaxial","abaxial","ab_ad","B_ends","M_ends","A_ends","B_branch","M_branch","A_branch","B_areole_num","M_areole_num","A_areole_num","A_skel","B_areole_area","M_areole_area","A_areole_area","B_density","M_density","A_density","Photo","Cond","WUE","Fo","Fv","FvFm","Palisade","Spongy","Palisade_spongy","Leaf_area","perimeter","dissection")

#split into two dataframes: one for diploids and one for tetraploids:
diploid<-cordatred[which(cordatred$Category=="p"),]
head(diploid)

tetraploid<-cordatred[which(cordatred$Category=="h"),]

#strip metadata; re-arrange columns so that the 4 overrepresented traits are at the front:
diploid<-cbind(diploid[,c(10,36,33,34)], diploid[c(8,9,11:32,35,37,38)])

tetraploid<-cbind(tetraploid[,c(10,36,33,34)], tetraploid[c(8,9,11:32,35,37,38)])

#function to calculate mutual rank correlation coefficients (based on Bennot et al. 2015 Although note that MutRank did not work so I had to go to github and pull out the code to run it independently.)

mutrank.wrap <- function(data){
    m <- cor(data, use= "pairwise.complete.obs", method = "pearson")
    r <- apply(m,1,rank)
    diag(r) <- 0
    net <- r*t(r)/2
    colnames(net) <- colnames(data)
    rownames(net) <- colnames(data)
    return(net)
}

#generate mutual rank matrices:
dip.MR<-mutrank.wrap(diploid)
tet.MR<-mutrank.wrap(tetraploid)

#get edgelist from MR matrix:
dip.g  <- graph.adjacency(dip.MR,weighted=TRUE)
dip.links <- get.data.frame(dip.g)
nrow(dip.links) #930

tet.g <- graph.adjacency(tet.MR, weighted=TRUE)
tet.links<-get.data.frame(tet.g)
nrow(tet.links) #930

#Make a nodelist
dip.nodes<-data.frame(colnames(dip.MR))
nrow(dip.nodes) #31
tet.nodes<-data.frame(colnames(tet.MR))
nrow(tet.nodes) #31

#make networks:
dip.net<-graph_from_data_frame(d=dip.links, vertices=dip.nodes, directed=F) #directed=F removes arrow heads on edges

tet.net<-graph_from_data_frame(d=tet.links, vertices=tet.nodes, directed=F) #directed=F removes arrow heads on edges

# get rid of duplicated edges
dip.to.delete<-which(duplicated(attr(E(dip.net), "vnames")) %in% "TRUE")
dip.net<-delete.edges(dip.net, dip.to.delete)
length(attr(E(dip.net), "vnames")) #465. good.

tet.to.delete<-which(duplicated(attr(E(tet.net), "vnames")) %in% "TRUE")
tet.net<-delete.edges(tet.net, tet.to.delete)
length(attr(E(tet.net), "vnames")) #465. good.

dist.dip<-(E(dip.net)$weight)
hist(dist.dip) #very right-skewed
hist(sqrt(dist.dip)) #better.
p.dip<-quantile(sqrt(dist.dip), probs=c(0.99, 0.95, 0.90, 0.8, 0.60))
p.dip*p.dip
#  99%   95%   90%   80%   60% 
#450.0 406.0 362.5 266.0 165.0 

plot(dip.net)
sub.dip99<-graph.neighborhood(subgraph.edges(dip.net, eids=which(E(dip.net)$weight>=450)), nodes=c("dissection", "Leaf_area", "Palisade", "Spongy"), order=2)
#error: invaldid vertex names: likely some vertices not connected at 99%!

#generates 4 graphs: each one starts with one of the "nodes" and works out to 2nd order nodes considering only edges in the top 95% of all edges:
sub.dip95<-graph.neighborhood(subgraph.edges(dip.net, eids=which(E(dip.net)$weight>=406)), nodes=c("dissection", "Leaf_area", "Palisade", "Spongy"), order=2)
#No single graph has all 4 "key traits" listed in "nodes"

pdf("sub.dip95.pdf", 20, 20)
par(mfrow=c(2,2))
set.seed(1)
V(sub.dip95[[1]])$color<-"gray"
V(sub.dip95[[1]])["dissection"]$color<-"red"
plot(sub.dip95[[1]], layout=layout_with_lgl, sub="Leaf dissection index", vertex.label.font=2, vertex.label.cex=2)  
V(sub.dip95[[2]])$color<-"gray"
V(sub.dip95[[2]])["Leaf_area"]$color<-"red"
plot(sub.dip95[[2]], layout=layout_with_lgl, sub="Leaf area", vertex.label.font=2, vertex.label.cex=2)  
V(sub.dip95[[3]])$color<-"gray"
V(sub.dip95[[3]])["Palisade"]$color<-"red"
plot(sub.dip95[[3]], layout=layout_with_lgl, sub="Palisde parenchyma", vertex.label.font=2, vertex.label.cex=2)  
V(sub.dip95[[4]])$color<-"gray"
V(sub.dip95[[4]])["Spongy"]$color<-"red"
plot(sub.dip95[[4]], layout=layout_with_lgl, sub="Spongy mesophyll", vertex.label.font=2, vertex.label.cex=2)  
title("Seeded Diploid Networks\n 95% threshold", line = -3, outer = TRUE)
dev.off()

sub.dip90<-graph.neighborhood(subgraph.edges(dip.net, eids=which(E(dip.net)$weight>=362.5)),nodes=c("dissection", "Leaf_area", "Palisade", "Spongy"), order=2)
sub.dip90

dev.new()
pdf("sub.dip90.pdf", 20, 20)
par(mfrow=c(2,2))
set.seed(1)
V(sub.dip90[[1]])$color<-"gray"
V(sub.dip90[[1]])["dissection"]$color<-"red"
plot(sub.dip90[[1]], layout=layout_with_lgl, sub="Leaf dissection index", vertex.label.font=2, vertex.label.cex=2)  
V(sub.dip90[[2]])$color<-"gray"
V(sub.dip90[[2]])["Leaf_area"]$color<-"red"
plot(sub.dip90[[2]], layout=layout_with_lgl, sub="Leaf area", vertex.label.font=2, vertex.label.cex=2)  
V(sub.dip90[[3]])$color<-"gray"
V(sub.dip90[[3]])["Palisade"]$color<-"red"
plot(sub.dip90[[3]], layout=layout_with_lgl, sub="Palisde parenchyma", vertex.label.font=2, vertex.label.cex=2)   #has all 4
V(sub.dip90[[4]])$color<-"gray"
V(sub.dip90[[4]])["Spongy"]$color<-"red"
plot(sub.dip90[[4]], layout=layout_with_lgl, sub="Spongy mesophyll", vertex.label.font=2, vertex.label.cex=2)  
title("Seeded Diploid Networks\n 90% threshold", line = -3, outer = TRUE)
dev.off()

sub.dip80<-graph.neighborhood(subgraph.edges(dip.net, eids=which(E(dip.net)$weight>=266)),nodes=c("dissection", "Leaf_area", "Palisade", "Spongy"), order=2)
dev.new()
pdf("sub.dip80.pdf", 20, 20)
par(mfrow=c(2,2))
set.seed(1)
V(sub.dip80[[1]])$color<-"gray"
V(sub.dip80[[1]])["dissection"]$color<-"red"
plot(sub.dip80[[1]], layout=layout_with_lgl, sub="Leaf dissection index", vertex.label.font=2, vertex.label.cex=2)   #has all 4
V(sub.dip80[[2]])$color<-"gray"
V(sub.dip80[[2]])["Leaf_area"]$color<-"red"
plot(sub.dip80[[2]], layout=layout_with_lgl, sub="Leaf area", vertex.label.font=2, vertex.label.cex=2)   #has all 4
V(sub.dip80[[3]])$color<-"gray"
V(sub.dip80[[3]])["Palisade"]$color<-"red"
plot(sub.dip80[[3]], layout=layout_with_lgl, sub="Palisde parenchyma", vertex.label.font=2, vertex.label.cex=2)    #has all 4
V(sub.dip80[[4]])$color<-"gray"
V(sub.dip80[[4]])["Spongy"]$color<-"red"
plot(sub.dip80[[4]], layout=layout_with_lgl, sub="Spongy mesophyll", vertex.label.font=2, vertex.label.cex=2)   #has all 4
title("Seeded Diploid Networks\n 80% threshold", line = -3, outer = TRUE)
dev.off()

dist.tet<-(E(tet.net)$weight)
p.tet<-quantile(sqrt(dist.tet), probs=c(0.99, 0.95, 0.90, 0.85, 0.8, .6))
p.tet*p.tet
#99% 95% 90% 85% 80% 60% 
#450 405 351 312 276 162 

sub.tet99<-graph.neighborhood(subgraph.edges(tet.net, eids=which(E(tet.net)$weight>=450)),nodes=c("dissection", "Leaf_area", "Palisade", "Spongy"), order=2)
#invalid vertiex: likley some vertices not connected at 99%!

sub.tet95<-graph.neighborhood(subgraph.edges(tet.net, eids=which(E(tet.net)$weight>=405)),nodes=c("dissection", "Leaf_area", "Palisade", "Spongy"), order=2)
#returns 4 disconnected graphs none have all 4 traits. 3 of 4 are just lines!
dev.new()
pdf("sub.tet95.pdf",20,20)
par(mfrow=c(2,2))
set.seed(1)
V(sub.tet95[[1]])$color<-"gray"
V(sub.tet95[[1]])["dissection"]$color<-"red"
plot(sub.tet95[[1]], layout=layout_with_lgl, sub="Leaf dissection index", vertex.label.font=2, vertex.label.cex=2)   
V(sub.tet95[[2]])$color<-"gray"
V(sub.tet95[[2]])["Leaf_area"]$color<-"red"
plot(sub.tet95[[2]], layout=layout_with_lgl, sub="Leaf area", vertex.label.font=2, vertex.label.cex=2)    
V(sub.tet95[[3]])$color<-"gray"
V(sub.tet95[[3]])["Palisade"]$color<-"red"
plot(sub.tet95[[3]], layout=layout_with_lgl, sub="Palisde parenchyma", vertex.label.font=2, vertex.label.cex=2)     
V(sub.tet95[[4]])$color<-"gray"
V(sub.tet95[[4]])["Spongy"]$color<-"red"
plot(sub.tet95[[4]], layout=layout_with_lgl, sub="Spongy mesophyll", vertex.label.font=2, vertex.label.cex=2)  
title("Seeded Tetraploid Networks\n 95% threshold", line = -3, outer = TRUE)
dev.off()

sub.tet90<-graph.neighborhood(subgraph.edges(tet.net, eids=which(E(tet.net)$weight>=351)),nodes=c("dissection", "Leaf_area", "Palisade", "Spongy"), order=2)
#All 4 contain all 4 "key traits"
dev.new()
pdf("sub.tet90.pdf", 20,20)
par(mfrow=c(2,2))
set.seed(1)
V(sub.tet90[[1]])$color<-"gray"
V(sub.tet90[[1]])["dissection"]$color<-"red"
plot(sub.tet90[[1]], layout=layout_with_lgl, sub="Leaf dissection index", vertex.label.font=2, vertex.label.cex=2)    
V(sub.tet90[[2]])$color<-"gray"
V(sub.tet90[[2]])["Leaf_area"]$color<-"red"
plot(sub.tet90[[2]], layout=layout_with_lgl, sub="Leaf area", vertex.label.font=2, vertex.label.cex=2)  
V(sub.tet90[[3]])$color<-"gray"
V(sub.tet90[[3]])["Palisade"]$color<-"red"
plot(sub.tet90[[3]], layout=layout_with_lgl, sub="Palisde parenchyma", vertex.label.font=2, vertex.label.cex=2)  
V(sub.tet90[[4]])$color<-"gray"
V(sub.tet90[[4]])["Spongy"]$color<-"red"
plot(sub.tet90[[4]], layout=layout_with_lgl, sub="Spongy mesophyll", vertex.label.font=2, vertex.label.cex=2)  
title("Seeded Tetraploid Networks\n 90% threshold", line = -3, outer = TRUE)
dev.off()

sub.tet80<-graph.neighborhood(subgraph.edges(tet.net, eids=which(E(tet.net)$weight>=276)),nodes=c("dissection", "Leaf_area", "Palisade", "Spongy"), order=2)
#All 4 contain all 4 "key traits"
dev.new()
pdf("sub.tet80.pdf", h=20, w=20)
par(mfrow=c(2,2))
set.seed(1)
V(sub.tet80[[1]])$color<-"gray"
V(sub.tet80[[1]])["dissection"]$color<-"red"
plot(sub.tet80[[1]], layout=layout_with_lgl, sub="Leaf dissection index", vertex.label.font=2, vertex.label.cex=2)  
V(sub.tet80[[2]])$color<-"gray"
V(sub.tet80[[2]])["Leaf_area"]$color<-"red"
plot(sub.tet80[[2]], layout=layout_with_lgl, sub="Leaf area", vertex.label.font=2, vertex.label.cex=2)       
V(sub.tet80[[3]])$color<-"gray"
V(sub.tet80[[3]])["Palisade"]$color<-"red"
plot(sub.tet80[[3]], layout=layout_with_lgl, sub="Palisde parenchyma", vertex.label.font=2, vertex.label.cex=2)       
V(sub.tet80[[4]])$color<-"gray"
V(sub.tet80[[4]])["Spongy"]$color<-"red"
plot(sub.tet80[[4]], layout=layout_with_lgl, sub="Spongy mesophyll", vertex.label.font=2, vertex.label.cex=2)   
title("Seeded Tetraploid Networks\n 80% threshold", line = -3, outer = TRUE)
dev.off()

diameter(sub.tet90[[1]]) #1637.5 #length of shortest path between two nodes 
get_diameter(sub.tet90[[1]])

#edge densities
Trait<-c("dissection", "Leaf_area", "Palisade", "Spongy")
dip95<-vector("list", 4)
dip95[[1]]<-transitivity(sub.dip95[[1]], type="global", isolates="zero") # no. of edges/all possible no. edges
dip95[[2]]<-transitivity(sub.dip95[[2]], type="global")
dip95[[3]]<-transitivity(sub.dip95[[3]], type="global")
dip95[[4]]<-transitivity(sub.dip95[[4]], type="global")

dip90<-vector("list", 4)
dip90[[1]]<-transitivity(sub.dip95[[1]], type="global") # no. of edges/all possible no. edges
dip90[[2]]<-transitivity(sub.dip90[[2]], type="global")
dip90[[3]]<-transitivity(sub.dip90[[3]], type="global")
dip90[[4]]<-transitivity(sub.dip90[[4]], type="global")

dip80<-vector("list", 4)
dip80[[1]]<-transitivity(sub.dip80[[1]]) # no. of edges/all possible no. edges
dip80[[2]]<-transitivity(sub.dip80[[2]])
dip80[[3]]<-transitivity(sub.dip80[[3]])
dip80[[4]]<-transitivity(sub.dip80[[4]])

tet95<-vector("list", 4)
tet95[[1]]<-transitivity(sub.tet95[[1]]) # no. of edges/all possible no. edges
tet95[[2]]<-transitivity(sub.tet95[[2]])
tet95[[3]]<-transitivity(sub.tet95[[3]])
tet95[[4]]<-transitivity(sub.tet95[[4]])

tet90<-vector("list", 4)
tet90[[1]]<-transitivity(sub.tet90[[1]]) # no. of edges/all possible no. edges
tet90[[2]]<-transitivity(sub.tet90[[2]])
tet90[[3]]<-transitivity(sub.tet90[[3]])
tet90[[4]]<-transitivity(sub.tet90[[4]])

tet80<-vector("list", 4)
tet80[[1]]<-transitivity(sub.tet80[[1]]) # no. of edges/all possible no. edges
tet80[[2]]<-transitivity(sub.tet80[[2]])
tet80[[3]]<-transitivity(sub.tet80[[3]])
tet80[[4]]<-transitivity(sub.tet80[[4]])

e_trans<-data.frame(cbind(Trait, dip95, tet95, dip90, tet90, dip80, tet80))

e_trans<-lapply(e_trans, unlist)
e_trans<-data.frame(e_trans)
head(e_trans)
#      Trait     dip95 tet95     dip90     tet90     dip80     tet80
#1 dissection 0.3333333 0.375 0.3333333 0.6000000 0.4906542 0.6467890
#2  Leaf_area 0.3333333 0.000 0.4883721 0.6315789 0.5217391 0.7241379
#3   Palisade 0.0000000 0.000 0.5121951 0.6500000 0.4906542 0.6498195
#4     Spongy 0.0000000 0.000 0.3870968 0.6101695 0.4716981 0.6410256
 
t.test(e_trans$dip95, e_trans$tet95, paired = TRUE, alternative = "two.sided")

# 	Paired t-test
# 
# data:  e_trans$dip95 and e_trans$tet95
# t = 0.83468, df = 3, p-value = 0.4651
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  -0.2050998  0.3509331
# sample estimates:
# mean of the differences 
#              0.07291667 

t.test(e_trans$dip90, e_trans$tet90, paired = TRUE, alternative = "two.sided")

# 	Paired t-test
# 
# data:  e_trans$dip90 and e_trans$tet90
# t = -6.1301, df = 3, p-value = 0.008729
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  -0.29272251 -0.09265305
# sample estimates:
# mean of the differences 
#              -0.1926878 

t.test(e_trans$dip80, e_trans$tet80, paired = TRUE, alternative = "two.sided")

# 	Paired t-test
# 
# data:  e_trans$dip80 and e_trans$tet80
# t = -16.209, df = 3, p-value = 0.0005109
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  -0.2054794 -0.1380338
# sample estimates:
# mean of the differences 
#              -0.1717566 

Trait<-c("dissection", "Leaf_area", "Palisade", "Spongy")
dip95<-vector("list", 4)
dip95[[1]]<-edge_density(sub.dip95[[1]]) # no. of edges/all possible no. edges
dip95[[2]]<-edge_density(sub.dip95[[2]])
dip95[[3]]<-edge_density(sub.dip95[[3]])
dip95[[4]]<-edge_density(sub.dip95[[4]])

dip90<-NULL
dip90[[1]]<-edge_density(sub.dip95[[1]]) # no. of edges/all possible no. edges
dip90[[2]]<-edge_density(sub.dip90[[2]])
dip90[[3]]<-edge_density(sub.dip90[[3]])
dip90[[4]]<-edge_density(sub.dip90[[4]])

dip80<-vector("list", 4)
dip80[[1]]<-edge_density(sub.dip80[[1]]) # no. of edges/all possible no. edges
dip80[[2]]<-edge_density(sub.dip80[[2]])
dip80[[3]]<-edge_density(sub.dip80[[3]])
dip80[[4]]<-edge_density(sub.dip80[[4]])

tet95<-vector("list", 4)
tet95[[1]]<-edge_density(sub.tet95[[1]]) # no. of edges/all possible no. edges
tet95[[2]]<-edge_density(sub.tet95[[2]])
tet95[[3]]<-edge_density(sub.tet95[[3]])
tet95[[4]]<-edge_density(sub.tet95[[4]])

tet90<-vector("list", 4)
tet90[[1]]<-edge_density(sub.tet90[[1]]) # no. of edges/all possible no. edges
tet90[[2]]<-edge_density(sub.tet90[[2]])
tet90[[3]]<-edge_density(sub.tet90[[3]])
tet90[[4]]<-edge_density(sub.tet90[[4]])

tet80<-vector("list", 4)
tet80[[1]]<-edge_density(sub.tet80[[1]]) # no. of edges/all possible no. edges
tet80[[2]]<-edge_density(sub.tet80[[2]])
tet80[[3]]<-edge_density(sub.tet80[[3]])
tet80[[4]]<-edge_density(sub.tet80[[4]])

e_density<-data.frame(cbind(Trait, dip95, tet95, dip90, tet90, dip80, tet80))
#        Trait     dip95     tet95     dip90     tet90     dip80     tet80
# 1 dissection 0.3333333       0.5 0.3333333 0.3928571 0.3666667 0.3308824
# 2  Leaf_area 0.3333333       0.5 0.3555556 0.4722222 0.4545455       0.5
# 3   Palisade       0.5 0.6666667 0.4166667 0.4222222 0.3666667 0.3676471
# 4     Spongy 0.2857143       0.5 0.3111111 0.4722222 0.4358974 0.4095238

#custom GGally correlation matrix with corrections for multiple testing.
#adapted from:
#https://stackoverflow.com/questions/61686171/how-to-add-the-spearman-correlation-p-value-along-with-correlation-coefficient-t

install.packages("ggExtra")
library(ggplot2)
library(GGally)
library(ggExtra)

## read in data
cordat<-read.csv(file="CorsDat.csv")

nrow(cordat) #110
sum(is.na(cordat$species)) #2
cordat<-cordat[-c(which(is.na(cordat$species))),]
sum(is.na(cordat$species)) #0
nrow(cordat) #108

## reduce data set to just the pertinent info:
cordatred<-cordat[,c(1:19,22:28, 30:33,36,37,41,44,45,52:54)]

cordatred$Category<-as.factor(cordatred$Category)
levels(cordatred$Category)<-strtrim(levels(cordatred$Category), 1)

cordatred$species<-as.factor(cordatred$species)
levels(cordatred$species)<-strtrim(levels(cordatred$species), 3)

colnames(cordatred)<-c("Plant","concat","loc","Block", "ID","species","Category","adaxial","abaxial","ab:ad","B ends","M ends","A ends","B branch","M branch","A branch","B areole#","M areole#","A areole#","A skel","B areole area","M areole area","A areole area","B density","M density","A density","Photo","Cond","WUE","Fo","Fv","FvFm","Palisade","Spongey","Palisade:spongy","Leaf area","perimeter","dissection")

printVar = function(x,y){
      vals = cor.test(x,y,
      method="pearson")[c("estimate","p.value")]

      vals[[1]]<-round(vals[[1]],2)   
      vals[[2]]<-p.adjust(vals[[2]], method="holm")
      
      vals[[2]]<-ifelse(test = vals[[2]]<0.001,"<0.001",ifelse(test=vals[[2]]<0.01,"<0.01",round(vals[[2]],2)))

          names(vals) = c("r","p")
      paste(names(vals),unlist(vals),collapse="\n")
}

#this printVar will return just adjusted p.values, not r values:
printVar = function(x,y){
      vals = cor.test(x,y,
      method="pearson")[c("p.value")]

      vals[[1]]<-p.adjust(vals[[1]], method="holm")
      
      vals[[1]]<-ifelse(test = vals[[1]]<0.001,"<0.001",ifelse(test=vals[[1]]<0.01,"<0.01",round(vals[[1]],2)))

          #names(vals) = c("r","p")
          names(vals) = c("p") #try printing just p-values
      paste(names(vals),unlist(vals),collapse="\n")
}

my_fn <- function(data, mapping, ...){
  # takes in x and y for each panel
  xData <- eval_data_col(data, mapping$x)
  yData <- eval_data_col(data, mapping$y)
  colorData <- eval_data_col(data, mapping$colour)

# if you have colors, split according to color group and calculate cor

  byGroup =by(data.frame(xData,yData),colorData,function(i)printVar(i[,1],i[,2]))
  byGroup = data.frame(col=names(byGroup),label=as.character(byGroup))
  byGroup$x = 0.5
  #byGroup$y = c(0.5, 0.1)
  #byGroup$y = seq(0.8-0.3,0.2,length.out=nrow(byGroup))
	byGroup$y = seq(1.0-0.3,0.2,length.out=nrow(byGroup))

#main correlation
mainCor = printVar(xData,yData)

p <- ggplot(data = data, mapping = mapping) +
annotate(x=0.5,y=0.9,label=mainCor,geom="text",size=1) +
geom_text(data=byGroup,inherit.aes=FALSE, lineheight=1,
aes(x=x,y=y,col=col,label=label),size=1)+ 
theme_void() + ylim(c(0,1))
  p
}

#so this gets pretty ugly. Seems to basically work but because there are 6 r values and 6 p values per cell, they print on top of each other. 

pdf("S2.Species_holm.pdf", h=15, w=15)

ggpairs(cordatred[,8:37],
mapping=ggplot2::aes(colour = cordatred$species, alpha=0.5), cardinality_threshold=110,
lower = list(continuous=wrap("points", alpha=0.6, size=0.1)),
upper = list(continuous = my_fn))+
rotateTextX(angle = 90, hjust = 1, vjust = 0.5)+

theme(strip.background = element_rect(fill = "white"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line=element_line(colour="black"), strip.text.x = element_text(size = 4), strip.text.y=element_text(size = 4))

dev.off()

trip.text.x = element_text(size = 5)

#custom GGally correlation matrix with corrections for multiple testing.

#adapted from:
#https://stackoverflow.com/questions/61686171/how-to-add-the-spearman-correlation-p-value-along-with-correlation-coefficient-t

install.packages("ggExtra")

library(ggplot2)
library(GGally)
library(ggExtra)

## read in data
cordat<-read.csv(file="CorsDat.csv")

nrow(cordat) #110
sum(is.na(cordat$species)) #2
cordat<-cordat[-c(which(is.na(cordat$species))),]
sum(is.na(cordat$species)) #0
nrow(cordat) #108

## reduce data set to just the pertinent info:
cordatred<-cordat[,c(1:19,22:28, 30:33,36,37,41,44,45,52:54)]

cordatred$Category<-as.factor(cordatred$Category)
levels(cordatred$Category)<-strtrim(levels(cordatred$Category), 1)
cordatred$species<-as.factor(cordatred$species)
levels(cordatred$species)<-strtrim(levels(cordatred$species), 3)
colnames(cordatred)<-c("Plant","concat","loc","Block", "ID","species","Category","adaxial","abaxial","ab:ad","B ends","M ends","A ends","B branch","M branch","A branch","B areole#","M areole#","A areole#","A skel","B areole area","M areole area","A areole area","B density","M density","A density","Photo","Cond","WUE","Fo","Fv","FvFm","Palisade","Spongey","Palisade:spongy","Leaf area","perimeter","dissection")

printVar = function(x,y){
      vals = cor.test(x,y,
      method="pearson")[c("estimate","p.value")]

      vals[[1]]<-round(vals[[1]],2)   
      vals[[2]]<-p.adjust(vals[[2]], method="holm")
      
      vals[[2]]<-ifelse(test = vals[[2]]<0.001,"<0.001",ifelse(test=vals[[2]]<0.01,"<0.01",round(vals[[2]],2)))

          names(vals) = c("r","p")
      paste(names(vals),unlist(vals),collapse="\n")
}

my_fn <- function(data, mapping, ...){
  # takes in x and y for each panel
  xData <- eval_data_col(data, mapping$x)
  yData <- eval_data_col(data, mapping$y)
  colorData <- eval_data_col(data, mapping$colour)

# if you have colors, split according to color group and calculate cor

  byGroup =by(data.frame(xData,yData),colorData,function(i)printVar(i[,1],i[,2]))
  byGroup = data.frame(col=names(byGroup),label=as.character(byGroup))
  byGroup$x = 0.5
  byGroup$y = c(0.5, 0.1)
  #byGroup$y = seq(0.8-0.3,0.2,length.out=nrow(byGroup))

#main correlation
mainCor = printVar(xData,yData)

p <- ggplot(data = data, mapping = mapping) +
annotate(x=0.5,y=0.9,label=mainCor,geom="text",size=1) +
geom_text(data=byGroup,inherit.aes=FALSE, lineheight=1,
aes(x=x,y=y,col=col,label=label),size=1)+ 
theme_void() + ylim(c(0,1))
  p
}

pdf("ParentHybrid_holm.pdf", h=15, w=15)

ggpairs(cordatred[,8:37],
mapping=ggplot2::aes(colour = cordatred$Category, alpha=0.5), cardinality_threshold=110,
lower = list(continuous=wrap("points", alpha=0.6, size=0.1)),
upper = list(continuous = my_fn))+
rotateTextX(angle = 90, hjust = 1, vjust = 0.5)+

theme(strip.background = element_rect(fill = "white"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line=element_line(colour="black"), strip.text.x = element_text(size = 4), strip.text.y=element_text(size = 4))

dev.off()
