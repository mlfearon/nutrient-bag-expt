
## Manuscript working title: "Nutrient enrichment, habitat structure, and disease"

## First author: Rachel M. Penczykowski
## Co-authors: Jessica L. Hite, Marta S. Shocket, Spencer R. Hall, and Meghan A. Duffy

## Experiment designed by MAD, SRH, and RMP. Mesocosum sampling performed by JLH, MSS, and RMP.  

## Type of experiment: bags (mesocosums)
## Treatments: Nutrients N & P (high/low), water column mixing (yes/no), Metschnikowia spores (yes/no)
## Location of experiment: University Lake at Indiana University in Bloomington
## Date of experiment: Autumn 2011

## Code by Rachel Penczykowski, with advice from Ben Bolker on prevalence analysis on 22 Aug 2012. 
## Latest code update: 19 Jan 2021

## R Project name: BagExpt2011_Reanalysis2018 

##################################################################################################
## Import data and load packages and a function
##################################################################################################
## Load packages
library(car)
library(ggplot2)
library(gridExtra)
library(lme4)
library(nlme)
library(plyr) ## do not load dplyr, or the SummarySE function won't run due to issue with "rename"
library(scales) # for 'squish'
library(splines)
library(stats)

## Load function
source("SummarySE.R")

## Import the data
mydata <- read.csv("BagExpt2011_SummaryLong_21Jan2021.csv")

##################################################################################################
## Define some variables and data subsets
##################################################################################################

## Exclude these replicates (bags) from all analyses. Bags 3, 9, 14, 20, and 21 were also excluded from Jessica's analyses.
mydata<-subset(mydata,Bag!=3) # LowNutrients/Unmixed/Spores: evidence of nutrient contamination in first week (= huge Daphnia population).
mydata<-subset(mydata,Bag!=9) # HighNutrients/Unmixed/NoSpores: entangled and destroyed by anchor line that broke loose from a raft on 7 Oct.
mydata<-subset(mydata,Bag!=14) # LowNutrients/Unmixed/NoSpores: >10x more Chaoborus than mean of other bags; Daphnia population crashed early.
mydata<-subset(mydata,Bag!=20) # HighNutrients/Mixed/Spores: Daphnia population crashed early.
mydata<-subset(mydata,Bag!=21) # LowNutrients/Unmixed/Spores: Daphnia population crashed early (but bag 21 was included in Rachel's PhD thesis).
## Note: In Rachel's PhD thesis, bag 19 was excluded due to >10x more Chaob than mean of other bags. But Jessica did not exclude bag 19 in her analyses.

# ### TEMPORARY EXCLUSION TO LOOK AT INFLUENCE OF THESE TWO BAGS
# mydata<-subset(mydata,Bag!=23) # LowNutrients/Mixed/Spores: Daphnia population crashed at end
# mydata<-subset(mydata,Bag!=24) # HighNutrients/Mixed/NoSpores: Daphnia population crashed at end

names(mydata) # for variable descriptions and units, see "BagExpt2011_Metadata.xlsx"
## Note: there are a few additional types of data not included in "BagExpt2011_SummaryLong_18May2018.csv" -- see FieldSchedule 

## Define variables
ZoopNetArea<-pi*0.065^2 # Hall lab zoop net radius = 0.065 m (diameter = 0.13 m)
mydata <- within(mydata,
             {
                 Raft <- factor(Raft)
                 Bag <- factor(Bag) 
                 Trt <- factor(paste(mydata$Nutrients,mydata$Mixing,mydata$Spores))
                 TotalCount <- UA+UJ+UM+MA+MJ+MM  
                 AdultCount <- UA+MA
                 UninfCount <- UA+UJ+UM
                 InfCount <- MA+MJ+MM
                 PrevAdult <- MA/(UA+MA)	
                 PrevTotal	<- InfCount/TotalCount
                 PropJuv	<- (UJ+MJ)/TotalCount
                 TotalDensity <- TotalCount/ZoopNetArea # TotalDensity units = number/(m^2)
                 UninfDensity <- UninfCount/ZoopNetArea # UninfDensity units = number/(m^2)
                 InfDensity <- InfCount/ZoopNetArea # InfDensity units = number/(m^2)
             })

mydata[is.na(mydata)] <- NA # replace NaN with NA (in cases where prevalence is NaN because denominator = 0)
str(mydata)

## For Daphnia density data, replace zeroes with areal density equivalent of 1 Daphnia/tow
## (Don't want any zeroes because we will log-transform these variables.)
summary(mydata$TotalDensity)
summary(mydata$UninfDensity)
summary(mydata$InfDensity)

NonzeroTotDens<-subset(mydata,TotalDensity>0)
MinDaphnia<-min(NonzeroTotDens$TotalDensity) # this is areal density equivalent of 1 Daphnia/tow

for (i in 1:nrow(mydata)){
  if(is.na(mydata$TotalDensity[i])) {
    mydata[i,"TotalDensity"]<-NA
  } else if (mydata[i,"TotalDensity"]==0){
    mydata[i,"TotalDensity"]<-MinDaphnia
  } else {
    mydata[i,"TotalDensity"]<-mydata[i,"TotalDensity"]
  } 
  if(is.na(mydata$UninfDensity[i])) {
    mydata[i,"UninfDensity"]<-NA
  } else if (mydata[i,"UninfDensity"]==0){
    mydata[i,"UninfDensity"]<-MinDaphnia
  } else {
    mydata[i,"UninfDensity"]<-mydata[i,"UninfDensity"]
  } 
  if(is.na(mydata$InfDensity[i])) {
    mydata[i,"InfDensity"]<-NA
  } else if (mydata[i,"InfDensity"]==0){
    mydata[i,"InfDensity"]<-MinDaphnia
  } else {
    mydata[i,"InfDensity"]<-mydata[i,"InfDensity"]
  } 
}

## Check that replacements worked
summary(mydata$TotalDensity) 
summary(mydata$UninfDensity)
summary(mydata$InfDensity) 

## Create additional variables
mydata$LogTotDens<-log(mydata$TotalDensity)
mydata$LogUninfDens<-log(mydata$UninfDensity)
mydata$LogInfDens<-log(mydata$InfDensity)

## Check that variables make sense
# str(mydata)
# hist(mydata$TotalDensity)
# hist(mydata$LogTotDens)
# hist(mydata$UninfDensity)
# hist(mydata$LogUninfDens)
# hist(mydata$InfDensity)
# hist(mydata$LogInfDens)

## Create additional variables
summary(mydata$TotChl)
summary(mydata$EdChl)
mydata$LogTotChl<-log(mydata$TotChl)
mydata$LogEdChl<-log(mydata$EdChl)

## I calculate InedChl as the difference between TotChl and EdChl.
## However, this is problematic when EdChl > TotChl.
## (This can happen by chance when we extract chl from separate filtered and unfiltered samples.)
## My solution: when EdChl > TotChl, set InedChl to zero.
for (i in 1:nrow(mydata)){
  if(is.na(mydata$TotChl[i])) {
    mydata[i,"InedChl"]<- NA
  } else if (is.na(mydata$EdChl[i])) {
    mydata[i,"InedChl"]<- NA
  } else if (mydata[i,"EdChl"]>mydata[i,"TotChl"]){
    mydata[i,"InedChl"]<- 0
  } else {
    mydata[i,"InedChl"]<-mydata[i,"TotChl"]-mydata[i,"EdChl"]
  }
}

NonzeroInedChl<-subset(mydata,InedChl>0)
MinInedChl<-min(NonzeroInedChl$InedChl) 

for (i in 1:nrow(mydata)){
  if(is.na(mydata$InedChl[i])) {
    mydata[i,"InedChl"]<- NA
  } else if(mydata[i,"InedChl"]==0){
    mydata[i,"InedChl"]<-MinInedChl/2 
  } else {
    mydata[i,"InedChl"]<-mydata[i,"InedChl"]
  }
}

summary(mydata$InedChl)
mydata$LogInedChl<-log(mydata$InedChl)

mydata$PropInedChl<-mydata$InedChl/mydata$TotChl
summary(mydata$PropInedChl)

## Check that variables make sense
# hist(mydata$LogTotChl)
# hist(mydata$LogEdChl)
# hist(mydata$LogInedChl)
# hist(mydata$PropInedChl)

## Create additional variables
summary(mydata$TP)
summary(mydata$TN)
mydata$LogTP<-log(mydata$TP)
mydata$LogTN<-log(mydata$TN)

## Check that variables make sense
# hist(mydata$LogTP)
# hist(mydata$LogTN)

## Create additional variables
summary(mydata$C)
summary(mydata$N)
summary(mydata$P)
mydata$LogC<-log(mydata$C)
mydata$LogN<-log(mydata$N)
mydata$LogP<-log(mydata$P)

## Check that variables make sense
# hist(mydata$LogC)
# hist(mydata$LogN)
# hist(mydata$LogP)
# hist(mydata$CN)
# hist(mydata$CP)
# hist(mydata$NP)

## Subset data for -spores and +spores bags
NoSpores<-subset(mydata,Spores=="No")
Spores<-subset(mydata,Spores=="Spores")

## export cleaned data set for SEM analysis
write.csv(mydata, "BagExpt_clean.csv", quote = F)

##################################################################################################
## Set standard formatting for figures 
##################################################################################################

SporesNames<-c("No"="-Spores","Spores"="+Spores")
NutColors<-c("red","blue") # nutrient treatment color code:  high = red, low = blue

fav_theme<-theme_bw()+
  theme(panel.border = element_rect(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=10), 
        axis.title=element_text(size=10),
        legend.text=element_text(size=10),
        legend.title=element_text(size=10),
        plot.title=element_text(size=12),
        axis.title.y=element_text(vjust=1.2),
        axis.title.x=element_text(vjust=0.4))

##################################################################################################
## SRP and C:N:P
##################################################################################################
NutData<-mydata[,c("Raft","Nutrients","Mixing","Spores","Trt","Bag","Day","NDay",
                  "SRP","CP","NP","NP","ad320")]
NutData<-na.omit(NutData)
write.csv(NutData,"NutData.csv")
DayNames<-c("31"="Day 31","46"="Day 46")

ggplot(NutData)+
  geom_boxplot(aes(x=Spores,y=SRP,color=Nutrients))+
  geom_point(aes(x=Spores,y=SRP,color=Nutrients),position=position_jitterdodge(),size=0.8,alpha=0.8)+
  facet_wrap(~Day,labeller = labeller(Day = DayNames))+
  scale_color_manual(values=NutColors)+
  ylab("SRP (ug P/L)")+
  fav_theme

SRP.1<-lme(SRP~Nutrients*Mixing*Spores*factor(Day)+Raft,random=~1|Bag,data=NutData)
Anova(SRP.1,type=3)

ggplot(NutData)+
  geom_boxplot(aes(x=Spores,y=CP,color=Nutrients))+
  geom_point(aes(x=Spores,y=CP,color=Nutrients),position=position_jitterdodge(),size=0.8,alpha=0.5)+
  facet_wrap(~Day,labeller = labeller(Day = DayNames))+
  scale_color_manual(values=NutColors)+
  ylab("C:P")+
  fav_theme

CP.1<-lme(CP~Nutrients*Mixing*Spores*factor(Day)+Raft,random=~1|Bag,data=NutData)
Anova(CP.1,type=3)

CP.99<-lme(CP~Nutrients+Spores+factor(Day)+Mixing+Raft,random=~1|Bag,data=NutData)
Anova(CP.99,type=3)
summary(CP.99)

ggplot(NutData)+
  geom_boxplot(aes(x=Spores,y=NP,color=Nutrients))+
  geom_point(aes(x=Spores,y=NP,color=Nutrients),position=position_jitterdodge(),size=0.8,alpha=0.5)+
  facet_wrap(~Day,labeller = labeller(Day = DayNames))+
  scale_color_manual(values=NutColors)+
  ylab("N:P")+
  fav_theme

ggplot(NutData)+
  geom_boxplot(aes(x=Spores,y=ad320,color=Nutrients))+
  geom_point(aes(x=Spores,y=ad320,color=Nutrients),position=position_jitterdodge(),size=0.8,alpha=0.5)+
  facet_wrap(~Day,labeller = labeller(Day = DayNames))+
  scale_color_manual(values=NutColors)+
  ylab("ad 320 (proxy for colored DOC")+
  fav_theme

##################################################################################################
## Egg ratios
##################################################################################################

names(mydata)
ggplot(mydata)+
  geom_boxplot(aes(x=factor(Day),y=UninfEggRatio,color=Nutrients))+
  facet_wrap(~Spores,labeller = labeller(Spores = SporesNames))+
  scale_color_manual(values=NutColors)+
  ylab("Uninfected Egg Ratio")+
  fav_theme

ggplot(mydata)+
  geom_boxplot(aes(x=factor(Day),y=InfEggRatio,color=Nutrients))+
  scale_color_manual(values=NutColors)+
  ylab("Infected Egg Ratio")+
  fav_theme

StartDay29<-subset(mydata,Day>=29) # all bags
ggplot(StartDay29)+
  geom_boxplot(aes(x=Spores,y=UninfEggRatio,color=Nutrients),notch=TRUE)+
  scale_color_manual(values=NutColors)+
  ylab("Uninf egg ratio")+
  labs(title="Day 29-46")+
  fav_theme

StartDay32<-subset(mydata,Day>=32) # all bags
ggplot(StartDay32)+
  geom_boxplot(aes(x=Spores,y=UninfEggRatio,color=Nutrients),notch=TRUE)+
  scale_color_manual(values=NutColors)+
  ylab("Uninf egg ratio")+
  labs(title="Day 32-46")+
  fav_theme

StartDay35<-subset(mydata,Day>=35) # all bags
ggplot(StartDay35)+
  geom_boxplot(aes(x=Spores,y=UninfEggRatio,color=Nutrients),notch=TRUE)+
  geom_point(aes(x=Spores,y=UninfEggRatio,color=Nutrients),position=position_jitterdodge(),size=0.8,alpha=0.8)+
  scale_color_manual(values=NutColors)+
  ylab("Uninf egg ratio")+
  labs(title="Day 35-46")+
  fav_theme

StartDay40<-subset(mydata,Day>=40) # all bags
StartDay40UER<-StartDay40[,c("Spores","Nutrients","Mixing","Raft","Bag","UninfEggRatio")]
StartDay40UER<-na.omit(StartDay40UER)

ggplot(StartDay40UER)+
  geom_boxplot(aes(x=Spores,y=UninfEggRatio,color=Nutrients),notch=TRUE)+
  geom_point(aes(x=Spores,y=UninfEggRatio,color=Nutrients),position=position_jitterdodge(),size=0.8,alpha=0.8)+
  scale_color_manual(values=NutColors)+
  ylab("Uninf egg ratio")+
  labs(title="Day 40-46")+
  fav_theme

ggplot(StartDay40UER)+
  geom_boxplot(aes(x=Spores,y=UninfEggRatio),notch=TRUE)+
  geom_point(aes(x=Spores,y=UninfEggRatio),size=0.8,alpha=0.8)+
  ylab("Uninf egg ratio")+
  labs(title="Day 40-46")+
  fav_theme

m.40.uer.1<-lme(UninfEggRatio~Nutrients*Mixing*Spores+Raft,random=~1|Bag,method="REML",data=StartDay40UER) 
Anova(m.40.uer.1,type=3)

m.40.uer.2<-update(m.40.uer.1,~.-Nutrients:Mixing:Spores)
Anova(m.40.uer.2,type=3)

m.40.uer.3<-update(m.40.uer.2,~.-Mixing:Spores)
Anova(m.40.uer.3,type=3)

m.40.uer.4<-update(m.40.uer.3,~.-Nutrients:Mixing)
Anova(m.40.uer.4,type=3)

m.40.uer.5<-update(m.40.uer.4,~.-Nutrients:Spores)
Anova(m.40.uer.5,type=3)
summary(m.40.uer.5)

m.40.uer.6<-update(m.40.uer.5,~.-Mixing)
Anova(m.40.uer.6,type=3)

m.40.uer.7<-update(m.40.uer.6,~.-Nutrients)
Anova(m.40.uer.7,type=3)


StartDay40<-subset(mydata,Day>=40) # all bags
ggplot(StartDay40)+
  geom_boxplot(aes(x=Spores,y=LogEdChl,color=Nutrients),notch=TRUE)+
  scale_color_manual(values=NutColors)+
  ylab("Log (Edible chl)")+
  labs(title="Day 40-46")+
  fav_theme



## Light
##################################################################################################
LightData<-mydata[,c("Raft","Nutrients","Mixing","Spores","Trt","Bag","Day","NDay",
                       "Light")]
LightData<-na.omit(LightData)
Day32Data<-subset(mydata,Day==32)
Day32DataLight<-merge(LightData,Day32Data,by="Bag")

ggplot(LightData)+
  geom_boxplot(aes(x=Spores,y=Light,color=Nutrients))+
  geom_point(aes(x=Spores,y=Light,color=Nutrients),position=position_jitterdodge(),size=0.8,alpha=0.5)+
  scale_color_manual(values=NutColors)+
  fav_theme

ggplot(Day32DataLight)+
  geom_point(aes(x=LogTotChl,y=Light.x,color=Nutrients.x))+
  stat_smooth(aes(x=LogTotChl,y=Light.x,color=Nutrients.x),method="lm")+
  facet_wrap(~Spores.x,labeller = labeller(Spores.x = SporesNames))+
  scale_color_manual(values=NutColors,name="Nutrients")+
  fav_theme+
  xlab("Log(Total chlorophyll)")+
  ylab("Light extinction, k")

ggplot(Day32DataLight)+
  geom_point(aes(x=LogTotChl,y=Light.x,color=Nutrients.x))+
  stat_smooth(aes(x=LogTotChl,y=Light.x,color=Nutrients.x),method="lm")+
  facet_wrap(~Mixing.x)+
  scale_color_manual(values=NutColors,name="Nutrients")+
  fav_theme+
  xlab("Log(Total chlorophyll)")+
  ylab("Light extinction, k")

m.light.1<-lm(Light.x~Nutrients.x*Mixing.x*Spores.x+Raft.x,data=Day32DataLight)
Anova(m.light.1,type=3)

m.light.2<-update(m.light.1,~.-Nutrients.x:Mixing.x:Spores.x)
Anova(m.light.2,type=3)

m.light.3<-update(m.light.2,~.-Nutrients.x:Spores.x)
Anova(m.light.3,type=3)

m.light.4<-update(m.light.3,~.-Nutrients.x:Mixing.x)
Anova(m.light.4,type=3)

m.light.5<-update(m.light.4,~.-Mixing.x:Spores.x)
Anova(m.light.5,type=3)

ggplot(Day32DataLight)+
  geom_point(aes(x=LogTotChl,y=Light.x,color=Nutrients.x))+
  stat_smooth(aes(x=LogTotChl,y=Light.x,color=Nutrients.x),method="lm")+
  facet_wrap(~Mixing.x)+
  scale_color_manual(values=NutColors,name="Nutrients")+
  fav_theme+
  xlab("Log(Total chlorophyll)")+
  ylab("Light extinction, k")

m.light.chl.1<-lm(Light.x~LogTotChl*Nutrients.x*Mixing.x*Spores.x+Raft.x,data=Day32DataLight)
Anova(m.light.chl.1,type=3)

m.light.chl.2<-update(m.light.chl.1,~.-LogTotChl:Nutrients.x:Mixing.x:Spores.x)
Anova(m.light.chl.2,type=3)

m.light.chl.3<-update(m.light.chl.2,~.-LogTotChl:Mixing.x:Spores.x)
Anova(m.light.chl.3,type=3)

m.light.chl.4<-update(m.light.chl.3,~.-LogTotChl:Nutrients.x:Mixing.x)
Anova(m.light.chl.4,type=3)

m.light.chl.5<-update(m.light.chl.4,~.-LogTotChl:Nutrients.x:Spores.x)
Anova(m.light.chl.5,type=3)

m.light.chl.6<-update(m.light.chl.5,~.-Nutrients.x:Mixing.x:Spores.x)
Anova(m.light.chl.6,type=3)

m.light.chl.7<-update(m.light.chl.6,~.-Mixing.x:Spores.x)
Anova(m.light.chl.7,type=3)

m.light.chl.8<-update(m.light.chl.7,~.-LogTotChl:Nutrients.x)
Anova(m.light.chl.8,type=3)

m.light.chl.9<-update(m.light.chl.8,~.-Nutrients.x:Mixing.x)
Anova(m.light.chl.9,type=3)

m.light.chl.10<-update(m.light.chl.9,~.-Nutrients.x:Spores.x)
Anova(m.light.chl.10,type=3)

m.light.chl.11<-update(m.light.chl.10,~.-LogTotChl:Mixing.x)
Anova(m.light.chl.11,type=3)

m.light.chl.12<-update(m.light.chl.11,~.-LogTotChl:Spores.x)
Anova(m.light.chl.12,type=3)
summary(m.light.chl.12)

##################################################################################################
## Total host density and chl
##################################################################################################
DensChlData<-mydata[,c("Raft","Nutrients","Mixing","Spores","Trt","Bag","Day","NDay",
                    "TotalDensity","UninfDensity","InfDensity","LogTotDens","LogUninfDens","LogInfDens","LogTotChl","LogEdChl","LogInedChl","PropInedChl","UninfEggRatio")]
DensChlData<-na.omit(DensChlData)

ggplot(DensChlData)+
  geom_point(aes(x=LogTotDens,y=LogTotChl,color=Nutrients))+
  scale_color_manual(values=NutColors)+
  facet_wrap(~Spores)

ggplot(DensChlData)+
  geom_point(aes(x=LogTotDens,y=LogEdChl,color=Nutrients))+
  scale_color_manual(values=NutColors)+
  facet_wrap(~Spores)

ggplot(DensChlData)+
  geom_point(aes(x=LogTotDens,y=LogInedChl,color=Nutrients))+
  scale_color_manual(values=NutColors)+
  facet_wrap(~Spores)

ggplot(DensChlData)+
  geom_point(aes(x=LogTotDens,y=PropInedChl,color=Nutrients))+
  scale_color_manual(values=NutColors)+
  facet_wrap(~Spores)

m001<-lme(LogTotChl~LogTotDens*Nutrients*Mixing*Spores,random=~1|Raft/Bag,data=DensChlData)
summary(m001)
Anova(m001,type=3)

m002<-update(m001,~.-LogTotDens:Nutrients:Mixing:Spores)
Anova(m002,type=3)

m003<-update(m002,~.-Mixing:Spores:LogTotDens)
Anova(m003,type=3)

m004<-update(m003,~.-Nutrients:Mixing:LogTotDens)
Anova(m004,type=3)

m005<-update(m004,~.-Nutrients:Spores:LogTotDens)
Anova(m005,type=3)

m006<-update(m005,~.-Nutrients:Spores)
Anova(m006,type=3)

m007<-update(m006,~.-Mixing:LogTotDens)
Anova(m007,type=3)

m008<-update(m007,~.-Nutrients:Mixing:Spores)
Anova(m008,type=3)

m009<-update(m008,~.-Mixing:Spores)
Anova(m009,type=3)

m010<-update(m009,~.-Nutrients:Mixing)
Anova(m010,type=3)

m011<-update(m010,~.-Spores:LogTotDens)
Anova(m011,type=3)
summary(m011)

##################################################################################################
## Change in log host density, log chl, and uninf egg ratios from start to end of experiment
##################################################################################################
DensDataStartEnd<-mydata[,c("Raft","Nutrients","Mixing","Spores","Trt","Bag","Day",
                            "TotalDensity","LogTotDens","LogTotChl","LogEdChl","LogInedChl","UninfEggRatio")]
#DensDataStartEnd<-na.omit(DensDataStartEnd)


DensData.1.12<-subset(DensDataStartEnd,Day==1|Day==12)
DensData.1.12<-reshape(DensData.1.12,direction="wide",idvar = "Bag",timevar="Day")
DensData.1.12$ChangeLogDens<-DensData.1.12$LogTotDens.12-DensData.1.12$LogTotDens.1
DensData.1.12$ChangeLogTotChl<-DensData.1.12$LogTotChl.12-DensData.1.12$LogTotChl.1
DensData.1.12$ChangeLogEdChl<-DensData.1.12$LogEdChl.12-DensData.1.12$LogEdChl.1
DensData.1.12$ChangeLogInedChl<-DensData.1.12$LogInedChl.12-DensData.1.12$LogInedChl.1
DensData.1.12$ChangeUninfEggRatio<-DensData.1.12$UninfEggRatio.12-DensData.1.12$UninfEggRatio.1

DensData.1.46<-subset(DensDataStartEnd,Day==1|Day==46)
DensData.1.46<-reshape(DensData.1.46,direction="wide",idvar = "Bag",timevar="Day")
DensData.1.46$ChangeLogDens<-DensData.1.46$LogTotDens.46-DensData.1.46$LogTotDens.1
DensData.1.46$ChangeLogTotChl<-DensData.1.46$LogTotChl.46-DensData.1.46$LogTotChl.1
DensData.1.46$ChangeLogEdChl<-DensData.1.46$LogEdChl.46-DensData.1.46$LogEdChl.1
DensData.1.46$ChangeLogInedChl<-DensData.1.46$LogInedChl.46-DensData.1.46$LogInedChl.1
DensData.1.46$ChangeUninfEggRatio<-DensData.1.46$UninfEggRatio.46-DensData.1.46$UninfEggRatio.1

DensData.26.46<-subset(DensDataStartEnd,Day==26|Day==46)
DensData.26.46<-reshape(DensData.26.46,direction="wide",idvar = "Bag",timevar="Day")
DensData.26.46$ChangeLogDens<-DensData.26.46$LogTotDens.46-DensData.26.46$LogTotDens.26
DensData.26.46$ChangeLogTotChl<-DensData.26.46$LogTotChl.46-DensData.26.46$LogTotChl.26
DensData.26.46$ChangeLogEdChl<-DensData.26.46$LogEdChl.46-DensData.26.46$LogEdChl.26
DensData.26.46$ChangeLogInedChl<-DensData.26.46$LogInedChl.46-DensData.26.46$LogInedChl.26
DensData.26.46$ChangeUninfEggRatio<-DensData.26.46$UninfEggRatio.46-DensData.26.46$UninfEggRatio.26

DensData.29.46<-subset(DensDataStartEnd,Day==29|Day==46)
DensData.29.46<-reshape(DensData.29.46,direction="wide",idvar = "Bag",timevar="Day")
DensData.29.46$ChangeLogDens<-DensData.29.46$LogTotDens.46-DensData.29.46$LogTotDens.29
DensData.29.46$ChangeLogTotChl<-DensData.29.46$LogTotChl.46-DensData.29.46$LogTotChl.29
DensData.29.46$ChangeLogEdChl<-DensData.29.46$LogEdChl.46-DensData.29.46$LogEdChl.29
DensData.29.46$ChangeLogInedChl<-DensData.29.46$LogInedChl.46-DensData.29.46$LogInedChl.29
DensData.29.46$ChangeUninfEggRatio<-DensData.29.46$UninfEggRatio.46-DensData.29.46$UninfEggRatio.29

DensData.32.46<-subset(DensDataStartEnd,Day==32|Day==46)
DensData.32.46<-reshape(DensData.32.46,direction="wide",idvar = "Bag",timevar="Day")
DensData.32.46$ChangeLogDens<-DensData.32.46$LogTotDens.46-DensData.32.46$LogTotDens.32
DensData.32.46$ChangeLogTotChl<-DensData.32.46$LogTotChl.46-DensData.32.46$LogTotChl.32
DensData.32.46$ChangeLogEdChl<-DensData.32.46$LogEdChl.46-DensData.32.46$LogEdChl.32
DensData.32.46$ChangeLogInedChl<-DensData.32.46$LogInedChl.46-DensData.32.46$LogInedChl.32
DensData.32.46$ChangeUninfEggRatio<-DensData.32.46$UninfEggRatio.46-DensData.32.46$UninfEggRatio.32

ggplot(DensData.1.46)+
  geom_point(aes(x=ChangeLogTotChl,y=ChangeLogEdChl,color=Nutrients.1),size=0.5,alpha=0.5)+ 
  scale_color_manual(values=NutColors,name="Nutrients")+
  stat_smooth(aes(x=ChangeLogTotChl,y=ChangeLogEdChl,color=Nutrients.1),method="lm")+
  facet_wrap(~Spores.1)+
  xlab("Change in Log (Total chl)")+
  ylab("Change in Log (Edible chl)")+
  labs(title="Change from Day 1 to 46")+
  fav_theme

ggplot(DensData.1.46)+
  geom_point(aes(x=ChangeLogEdChl,y=ChangeLogInedChl,color=Nutrients.1),size=0.5,alpha=0.5)+ 
  scale_color_manual(values=NutColors,name="Nutrients")+
  stat_smooth(aes(x=ChangeLogEdChl,y=ChangeLogInedChl,color=Nutrients.1),method="lm")+
  facet_wrap(~Spores.1)+
  xlab("Change in Log (Edible chl)")+
  ylab("Change in Log (Inedible chl)")+
  labs(title="Change from Day 1 to 46")+
  fav_theme

ggplot(DensData.1.46)+
  geom_point(aes(x=ChangeLogEdChl,y=ChangeUninfEggRatio,color=Nutrients.1),size=0.5,alpha=0.5)+ 
  scale_color_manual(values=NutColors,name="Nutrients")+
  stat_smooth(aes(x=ChangeLogEdChl,y=ChangeUninfEggRatio,color=Nutrients.1),method="lm")+
  facet_wrap(~Spores.1)+
  xlab("Change in Log (Edible chl)")+
  ylab("Change in Uninf Egg Ratio")+
  labs(title="Change from Day 1 to 46")+
  fav_theme

ggplot(DensData.1.46)+
  geom_point(aes(x=ChangeLogEdChl,y=ChangeUninfEggRatio,color=Nutrients.1),size=0.5,alpha=0.5)+ 
  scale_color_manual(values=NutColors,name="Nutrients")+
  stat_smooth(aes(x=ChangeLogEdChl,y=ChangeUninfEggRatio,color=Nutrients.1),method="lm")+
  xlab("Change in Log (Edible chl)")+
  ylab("Change in Uninf Egg Ratio")+
  labs(title="Change from Day 1 to 46")+
  fav_theme

ggplot(DensData.1.46)+
  geom_point(aes(x=ChangeLogEdChl,y=ChangeUninfEggRatio,color=Spores.1),size=0.5,alpha=0.5)+ 
  scale_color_manual(values=c("pink","black"),name="Spores")+
  stat_smooth(aes(x=ChangeLogEdChl,y=ChangeUninfEggRatio,color=Spores.1),method="lm")+
  xlab("Change in Log (Edible chl)")+
  ylab("Change in Uninf Egg Ratio")+
  labs(title="Change from Day 1 to 46")+
  fav_theme

ggplot(DensData.1.12)+
  geom_boxplot(aes(x=Mixing.1,y=ChangeUninfEggRatio,color=Nutrients.1))+ 
  scale_color_manual(values=NutColors,name="Nutrients")+
  ylab("Change in Uninf Egg Ratio")+
  labs(title="Change from Day 1 to 12")+
  fav_theme

ggplot(DensData.1.46)+
  geom_boxplot(aes(x=Spores.1,y=ChangeUninfEggRatio,color=Nutrients.1))+ 
  scale_color_manual(values=NutColors,name="Nutrients")+
  xlab("Spores")+
  ylab("Change in Uninf Egg Ratio")+
  labs(title="Change from Day 1 to 46")+
  fav_theme

ggplot(DensData.26.46)+
  geom_boxplot(aes(x=Spores.26,y=ChangeUninfEggRatio,color=Nutrients.26))+ 
  scale_color_manual(values=NutColors,name="Nutrients")+
  xlab("Spores")+
  ylab("Change in Uninf Egg Ratio")+
  labs(title="Change from Day 26 to 46")+
  fav_theme

ggplot(DensData.29.46)+
  geom_boxplot(aes(x=Spores.29,y=ChangeUninfEggRatio,color=Nutrients.29))+ 
  scale_color_manual(values=NutColors,name="Nutrients")+
  xlab("Spores")+
  ylab("Change in Uninf Egg Ratio")+
  labs(title="Change from Day 29 to 46")+
  fav_theme

ggplot(DensData.32.46)+
  geom_boxplot(aes(x=Spores.32,y=ChangeUninfEggRatio,color=Nutrients.32))+ 
  scale_color_manual(values=NutColors,name="Nutrients")+
  xlab("Spores")+
  ylab("Change in Uninf Egg Ratio")+
  labs(title="Change from Day 32 to 46")+
  fav_theme

ggplot(DensData.1.46)+
  geom_point(aes(x=ChangeLogTotChl,y=ChangeUninfEggRatio,color=Nutrients.1),size=0.5,alpha=0.5)+ 
  scale_color_manual(values=NutColors,name="Nutrients")+
  stat_smooth(aes(x=ChangeLogTotChl,y=ChangeUninfEggRatio,color=Nutrients.1),method="lm")+
  facet_wrap(~Spores.1)+
  xlab("Change in Log (Total chl)")+
  ylab("Change in Uninf Egg Ratio")+
  labs(title="Change from Day 1 to 46")+
  fav_theme

ggplot(DensData.1.46)+
  geom_point(aes(x=ChangeLogDens,y=ChangeUninfEggRatio,color=Nutrients.1),size=0.5,alpha=0.5)+ 
  scale_color_manual(values=NutColors,name="Nutrients")+
  stat_smooth(aes(x=ChangeLogDens,y=ChangeUninfEggRatio,color=Nutrients.1),method="lm")+
  facet_wrap(~Spores.1)+
  xlab("Change in Log (Total host density, #/m^2)")+
  ylab("Change in Uninf Egg Ratio")+
  labs(title="Change from Day 1 to 46")+
  fav_theme

ggplot(DensData.1.46)+
  geom_point(aes(x=ChangeLogDens,y=ChangeLogTotChl,color=Nutrients.1),size=0.5,alpha=0.5)+ 
  scale_color_manual(values=NutColors,name="Nutrients")+
  stat_smooth(aes(x=ChangeLogDens,y=ChangeLogTotChl,color=Nutrients.1),method="lm")+
  facet_wrap(~Spores.1)+
  xlab("Change in Log (Total host density, #/m^2)")+
  ylab("Change in Log (Total chlorophyll)")+
  labs(title="Change from Day 1 to 46")+
  fav_theme

ggplot(DensData.1.46)+
  geom_point(aes(x=ChangeLogDens,y=ChangeLogTotChl,color=Nutrients.1),size=0.5,alpha=0.5)+ 
  scale_color_manual(values=NutColors,name="Nutrients")+
  stat_smooth(aes(x=ChangeLogDens,y=ChangeLogTotChl,color=Nutrients.1),method="lm")+
  xlab("Change in Log (Total host density, #/m^2)")+
  ylab("Change in Log (Total chlorophyll)")+
  labs(title="Change from Day 1 to 46")+
  fav_theme

ggplot(DensData.1.46)+
  geom_point(aes(x=ChangeLogDens,y=ChangeLogTotChl),size=0.5,alpha=0.5)+ 
  stat_smooth(aes(x=ChangeLogDens,y=ChangeLogTotChl),method="lm")+
  xlab("Change in Log (Total host density, #/m^2)")+
  ylab("Change in Log (Total chlorophyll)")+
  labs(title="Change from Day 1 to 46")+
  fav_theme

ggplot(DensData.1.46)+
  stat_smooth(aes(x=ChangeLogDens,y=ChangeLogEdChl),color="black",method="lm")+
  geom_point(aes(x=ChangeLogDens,y=ChangeLogEdChl,color=Nutrients.1,shape=Mixing.1,fill=Spores.1),size=1,alpha=1)+ 
  geom_text(aes(x=ChangeLogDens,y=ChangeLogEdChl,label=Bag))+
  scale_color_manual(values=NutColors,name="Nutrients")+
  scale_shape_manual(values=c(21,24),name="Mixing")+
  scale_fill_manual(values=c("white","black"),name="Spores")+
  xlab("Change in Log (Total host density, #/m^2)")+
  ylab("Change in Log (Edible chlorophyll)")+
  labs(title="Change from Day 1 to 46")+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  geom_vline(xintercept=0, linetype="dashed", color = "black")+
  fav_theme

m1<-lm(ChangeLogEdChl~ChangeLogDens*Nutrients.1*Mixing.1*Spores.1+Raft.1,data=DensData.1.46)
summary(m1)
Anova(m1,type=3)

m2<-update(m1,~.-ChangeLogDens:Nutrients.1:Mixing.1:Spores.1) 
Anova(m2,type=3)

m3<-update(m2,~.-Nutrients.1:Mixing.1:Spores.1) 
Anova(m3,type=3)

m4<-update(m3,~.-ChangeLogDens:Nutrients.1:Spores.1) 
Anova(m4,type=3)

m5<-update(m4,~.-ChangeLogDens:Mixing.1:Spores.1) 
Anova(m5,type=3)

m6<-update(m5,~.-ChangeLogDens:Nutrients.1:Mixing.1) 
Anova(m6,type=3)

m7<-update(m6,~.-ChangeLogDens:Spores.1) 
Anova(m7,type=3)

m8<-update(m7,~.-Nutrients.1:Spores.1) 
Anova(m8,type=3)

m9<-update(m8,~.-ChangeLogDens:Mixing.1) 
Anova(m9,type=3)

m10<-update(m9,~.-ChangeLogDens:Nutrients.1) 
Anova(m10,type=3)

m11<-update(m10,~.-Mixing.1:Spores.1) 
Anova(m11,type=3)

m12<-update(m11,~.-Nutrients.1:Mixing.1) 
Anova(m12,type=2)

# m13<-update(m12,~.-Mixing.1) 
# Anova(m13,type=2)
# 
# m14<-update(m13,~.-Spores.1) 
# Anova(m14,type=2)

ggplot(DensData.1.46)+
  geom_boxplot(aes(x=Nutrients.1,y=ChangeLogEdChl,color=Nutrients.1),notch=TRUE)+ 
  scale_color_manual(values=NutColors,name="Nutrients")+
  ylab("Change in Log (Edible chlorophyll)")+
  labs(title="Change from Day 1 to 46")+
  fav_theme

m1<-lm(ChangeLogEdChl~Nutrients.1,data=DensData.1.46)
summary(m1)
Anova(m1,type=3)

cor.test(DensData.1.46$ChangeLogDens,DensData.1.46$ChangeLogEdChl)
ggplot(DensData.26.46)+
  geom_point(aes(x=ChangeLogDens,y=ChangeLogTotChl),size=0.5,alpha=0.5)+ 
  stat_smooth(aes(x=ChangeLogDens,y=ChangeLogTotChl),method="lm")+
  xlab("Change in Log (Total host density, #/m^2)")+
  ylab("Change in Log (Total chlorophyll)")+
  labs(title="Change from Day 26 to 46")+
  fav_theme

ggplot(DensData.29.46)+
  geom_point(aes(x=ChangeLogDens,y=ChangeLogTotChl),size=0.5,alpha=0.5)+ 
  stat_smooth(aes(x=ChangeLogDens,y=ChangeLogTotChl),method="lm")+
  xlab("Change in Log (Total host density, #/m^2)")+
  ylab("Change in Log (Total chlorophyll)")+
  labs(title="Change from Day 29 to 46")+
  fav_theme

ggplot(DensData.29.46)+
  geom_point(aes(x=ChangeLogDens,y=ChangeLogEdChl),size=0.5,alpha=0.5)+ 
  stat_smooth(aes(x=ChangeLogDens,y=ChangeLogEdChl),method="lm")+
  xlab("Change in Log (Total host density, #/m^2)")+
  ylab("Change in Log (Edible chlorophyll)")+
  labs(title="Change from Day 29 to 46")+
  fav_theme

ggplot(DensData.29.46)+
  geom_point(aes(x=ChangeLogEdChl,y=ChangeUninfEggRatio,color=Nutrients.29),size=0.5,alpha=0.5)+ 
  scale_color_manual(values=NutColors,name="Nutrients")+
  stat_smooth(aes(x=ChangeLogEdChl,y=ChangeUninfEggRatio,color=Nutrients.29),method="lm")+
  facet_wrap(~Spores.29)+
  xlab("Change in Log (Edible chl)")+
  ylab("Change in Uninf Egg Ratio")+
  labs(title="Change from Day 29 to 46")+
  fav_theme

ggplot(DensData.26.46)+
  geom_point(aes(x=ChangeLogEdChl,y=ChangeUninfEggRatio,color=Nutrients.26),size=0.5,alpha=0.5)+ 
  scale_color_manual(values=NutColors,name="Nutrients")+
  stat_smooth(aes(x=ChangeLogEdChl,y=ChangeUninfEggRatio,color=Nutrients.26),method="lm")+
  facet_wrap(~Spores.26)+
  xlab("Change in Log (Edible chl)")+
  ylab("Change in Uninf Egg Ratio")+
  labs(title="Change from Day 26 to 46")+
  fav_theme

m.change.1.46.1<-lm(ChangeLogTotChl~ChangeLogDens*Nutrients.1*Mixing.1*Spores.1+Raft.1,data=DensData.1.46)
Anova(m.change.1.46.1)  

m.change.1.46.2<-update(m.change.1.46.1,~.-ChangeLogDens:Nutrients.1:Mixing.1:Spores.1) 
Anova(m.change.1.46.2,type=3)

m.change.1.46.3<-update(m.change.1.46.2,~.-Nutrients.1:Mixing.1:Spores.1) 
Anova(m.change.1.46.3,type=3)

m.change.1.46.4<-update(m.change.1.46.3,~.-ChangeLogDens:Nutrients.1:Spores.1) 
Anova(m.change.1.46.4,type=3)

m.change.1.46.5<-update(m.change.1.46.4,~.-ChangeLogDens:Mixing.1:Spores.1) 
Anova(m.change.1.46.5,type=3)

m.change.1.46.6<-update(m.change.1.46.5,~.-ChangeLogDens:Nutrients.1:Mixing.1) 
Anova(m.change.1.46.6,type=3)

m.change.1.46.7<-update(m.change.1.46.6,~.-ChangeLogDens:Mixing.1) 
Anova(m.change.1.46.7,type=3)

m.change.1.46.8<-update(m.change.1.46.7,~.-Mixing.1:Spores.1) 
Anova(m.change.1.46.8,type=3)

m.change.1.46.9<-update(m.change.1.46.8,~.-Nutrients.1:Spores.1) 
Anova(m.change.1.46.9,type=3)

m.change.1.46.10<-update(m.change.1.46.9,~.-ChangeLogDens:Nutrients.1) 
Anova(m.change.1.46.10,type=3)

m.change.1.46.11<-update(m.change.1.46.10,~.-ChangeLogDens:Spores.1) 
Anova(m.change.1.46.11,type=3)

m.change.1.46.12<-update(m.change.1.46.11,~.-Nutrients.1:Mixing.1) 
Anova(m.change.1.46.12,type=3)
summary(m.change.1.46.12)

ggplot(DensData.1.46)+
  geom_boxplot(aes(x=Spores.1,y=ChangeUninfEggRatio,color=Nutrients.1))+
  geom_point(aes(x=Spores.1,y=ChangeUninfEggRatio,color=Nutrients.1),size=0.5,position = position_jitterdodge(),alpha=0.5)+ 
  scale_color_manual(values=NutColors,name="Nutrients")+
  xlab("Spores")+
  ylab("Change in Uninf Egg Ratio")+
  labs(title="Change from Day 1 to 46")+
  fav_theme

ggplot(DensData.1.46)+
  geom_boxplot(aes(x=Spores.1,y=ChangeLogTotChl,color=Nutrients.1))+
  geom_point(aes(x=Spores.1,y=ChangeLogTotChl,color=Nutrients.1),size=0.5,position = position_jitterdodge(),alpha=0.5)+ 
  scale_color_manual(values=NutColors,name="Nutrients")+
  xlab("Spores")+
  ylab("Change in Log (Total chl)")+
  labs(title="Change from Day 1 to 46")+
  fav_theme

m.1.46.1<-lm(ChangeLogTotChl~Nutrients.1*Mixing.1*Spores.1+Raft.1,data=DensData.1.46)
Anova(m.1.46.1,type=3)
summary(m.1.46.1)

m.1.46.2<-update(m.1.46.1,~.-Nutrients.1:Mixing.1:Spores.1) 
Anova(m.1.46.2,type=3)

m.1.46.3<-update(m.1.46.2,~.-Nutrients.1:Spores.1) 
Anova(m.1.46.3,type=3)

m.1.46.4<-update(m.1.46.3,~.-Mixing.1:Spores.1) 
Anova(m.1.46.4,type=3)

m.1.46.5<-update(m.1.46.4,~.-Nutrients.1:Mixing.1) 
Anova(m.1.46.5,type=3)
summary(m.1.46.5)

DensData.29.46<-subset(DensDataStartEnd,Day==29|Day==46)
DensData.29.46<-reshape(DensData.29.46,direction="wide",idvar = "Bag",timevar="Day")
DensData.29.46$ChangeLogDens<-DensData.29.46$LogTotDens.46-DensData.29.46$LogTotDens.29
DensData.29.46$ChangeLogTotChl<-DensData.29.46$LogTotChl.46-DensData.29.46$LogTotChl.29

ggplot(DensData.29.46)+
  geom_point(aes(x=ChangeLogDens,y=ChangeLogTotChl,color=Nutrients.29),size=0.5,alpha=0.5)+ 
  scale_color_manual(values=NutColors,name="Nutrients")+
  stat_smooth(aes(x=ChangeLogDens,y=ChangeLogTotChl,color=Nutrients.29),method="lm")+
  facet_wrap(~Spores.29)+
  xlab("Change in Log (Total host density, #/m^2)")+
  ylab("Change in Log (Total chlorophyll)")+
  labs(title="Change from Day 29 to 46")+
  fav_theme

ggplot(DensData.29.46)+
  geom_boxplot(aes(x=Spores.29,y=ChangeLogTotChl,color=Nutrients.29))+
  geom_point(aes(x=Spores.29,y=ChangeLogTotChl,color=Nutrients.29),size=0.5,position = position_jitterdodge(),alpha=0.5)+ 
  scale_color_manual(values=NutColors,name="Nutrients")+
  xlab("Spores")+
  ylab("Change in Log (Total chl)")+
  labs(title="Change from Day 29 to 46")+
  fav_theme

##################################################################################################
## Plots of host density over different time ranges
##################################################################################################
## Subset data 
DensData<-mydata[,c("Raft","Nutrients","Mixing","Spores","Trt","Bag","Day","NDay",
                    "TotalDensity","UninfDensity","InfDensity","LogTotDens","LogUninfDens","LogInfDens")]

## Omit missing data
DensData<-na.omit(DensData)

## Subset data for -spores and +spores bags
DensDataNoSpores<-subset(DensData,Spores=="No")
DensDataSpores<-subset(DensData,Spores=="Spores")

DayList<-c(12,15,19,22,26,29,32,35,40,43,46)
for(i in DayList){
  source("SummarySE.R")
  print(i)
  DensData<-subset(DensData,Day>=i)
  DensDataNoSpores<-subset(DensDataNoSpores,Day>=i)
  DensDataSpores<-subset(DensDataSpores,Day>=i)
  
  p1<-(ggplot(DensDataNoSpores,aes(x=Day,y=LogTotDens,group=Bag,color=Nutrients))+
         scale_color_manual(values=NutColors)+
         geom_line(aes(linetype=Mixing))+
         ylim(c(3,14))+
         ylab("Log (Total host density, #/m^2)")+
         xlab("Day of experiment")+
         labs(title=paste("- Spores, Day",i))+
         guides(color="none",linetype="none")+
         fav_theme)
  
  p2<-(ggplot(DensDataSpores,aes(x=Day,y=LogTotDens,group=Bag,color=Nutrients))+
         scale_color_manual(values=NutColors)+
         geom_line(aes(linetype=Mixing))+
         ylim(c(3,14))+
         ylab("Log (Total host density, #/m^2)")+
         xlab("Day of experiment")+
         labs(title=paste("+ Spores, Day",i))+
         guides(color="none",linetype="none")+
         fav_theme)
  
  UninfLogTotDensSummary<-summarySE(DensDataNoSpores, measurevar="LogTotDens", groupvars=c("Nutrients","Mixing","Day"))

  p3<-(ggplot(UninfLogTotDensSummary)+
         geom_point(aes(x=Day, y=LogTotDens, color=Nutrients)) +
         geom_line(aes(x=Day, y=LogTotDens, color=Nutrients, linetype=Mixing))+
         geom_errorbar(aes(x=Day, ymin=LogTotDens-se, ymax=LogTotDens+se,color=Nutrients), width=0) +
         scale_color_manual(values=NutColors)+
         fav_theme+
         ylab("Log (Total host density, #/m^2)")+
         xlab("Day of experiment")+
         labs(title=paste("- Spores, Day",i))+
         ylim(c(3,14)))
  
  InfLogTotDensSummary<-summarySE(DensDataSpores, measurevar="LogTotDens", groupvars=c("Nutrients","Mixing","Day"))
 
  p4<-(ggplot(InfLogTotDensSummary)+
         geom_point(aes(x=Day, y=LogTotDens, color=Nutrients)) +
         geom_line(aes(x=Day, y=LogTotDens, color=Nutrients, linetype=Mixing))+
         geom_errorbar(aes(x=Day, ymin=LogTotDens-se, ymax=LogTotDens+se,color=Nutrients), width=0) +
         scale_color_manual(values=NutColors)+
         fav_theme+
         ylab("Log (Total host density, #/m^2)")+
         xlab("Day of experiment")+
         labs(title=paste("+ Spores, Day",i))+
         ylim(c(3,14)))
  
  grid.arrange(p1,p3,p2,p4, nrow=2,ncol=2)
  
  UninfLogTotDensSummary2<-summarySE(DensDataNoSpores, measurevar="LogTotDens", groupvars=c("Nutrients","Mixing"))
  
  p3.1<-(ggplot(UninfLogTotDensSummary2)+
         geom_point(aes(x=Nutrients, y=LogTotDens, color=Nutrients,shape=Mixing)) +
       #  geom_line(aes(x=Nutrients, y=LogTotDens, color=Nutrients, linetype=Mixing))+
         geom_errorbar(aes(x=Nutrients, ymin=LogTotDens-se, ymax=LogTotDens+se,color=Nutrients), width=0) +
         scale_color_manual(values=NutColors)+
         fav_theme+
         ylab("Log (Total host density, #/m^2)")+
         xlab("Nutrients")+
         labs(title=paste("- Spores, Day",i))+
         ylim(c(3,14)))
  
  InfLogTotDensSummary2<-summarySE(DensDataSpores, measurevar="LogTotDens", groupvars=c("Nutrients","Mixing"))
  
  p4.1<-(ggplot(InfLogTotDensSummary2)+
         geom_point(aes(x=Nutrients, y=LogTotDens, color=Nutrients,shape=Mixing)) +
       #  geom_line(aes(x=Nutrients, y=LogTotDens, color=Nutrients, linetype=Mixing))+
         geom_errorbar(aes(x=Nutrients, ymin=LogTotDens-se, ymax=LogTotDens+se,color=Nutrients), width=0) +
         scale_color_manual(values=NutColors)+
         fav_theme+
         ylab("Log (Total host density, #/m^2)")+
         xlab("Nutrients")+
         labs(title=paste("+ Spores, Day",i))+
         ylim(c(3,14)))
  
  grid.arrange(p3.1,p4.1, nrow=2,ncol=1)

  p5<-(ggplot(DensDataSpores,aes(x=Day,y=LogInfDens,group=Bag,color=Nutrients))+
         scale_color_manual(values=NutColors)+
         geom_line(aes(linetype=Mixing))+
         ylim(c(4,12))+
         ylab("Log (Infected host density, #/m^2)")+
         xlab("Day of experiment")+ 
         labs(title=paste("+ Spores, Day",i))+
         guides(color="none",linetype="none")+  
         fav_theme)
  
  LogInfDensSummary<-summarySE(DensDataSpores, measurevar="LogInfDens", groupvars=c("Nutrients","Mixing","Day"))
  
  p6<-(ggplot(LogInfDensSummary)+
         geom_point(aes(x=Day, y=LogInfDens, color=Nutrients,shape=Mixing)) +  
         geom_line(aes(x=Day, y=LogInfDens, color=Nutrients, linetype=Mixing))+
         geom_errorbar(aes(x=Day, ymin=LogInfDens-se, ymax=LogInfDens+se,color=Nutrients), width=0) + 
         scale_color_manual(values=NutColors)+
         fav_theme+
         ylab("Log (Infected host density, #/m^2)")+
         xlab("Day of experiment")+ 
         labs(title=paste("+ Spores, Day",i))+
         ylim(c(4,12)))
  
  grid.arrange(p5,p6, nrow=1,ncol=2)
  
  LogInfDensSummary2<-summarySE(DensDataSpores, measurevar="LogInfDens", groupvars=c("Nutrients","Mixing"))
  
  p6.1<-(ggplot(LogInfDensSummary2)+
         geom_point(aes(x=Nutrients, y=LogInfDens, color=Nutrients,shape=Mixing)) +  
  #       geom_line(aes(x=Nutrients, y=LogInfDens, color=Nutrients, linetype=Mixing))+
         geom_errorbar(aes(x=Nutrients, ymin=LogInfDens-se, ymax=LogInfDens+se,color=Nutrients), width=0) + 
         scale_color_manual(values=NutColors)+
         fav_theme+
         ylab("Log (Infected host density, #/m^2)")+
         xlab("Day of experiment")+ 
         labs(title=paste("+ Spores, Day",i))+
         ylim(c(4,12)))
  print(p6.1)
  
}

##################################################################################################
## Change in log host density from start to end of experiment
##################################################################################################
DensDataStartEnd<-mydata[,c("Raft","Nutrients","Mixing","Spores","Trt","Bag","Day",
                    "TotalDensity","LogTotDens")]
DensDataStartEnd<-na.omit(DensDataStartEnd)


DensData.1.12<-subset(DensDataStartEnd,Day==1|Day==12)
DensData.1.12<-reshape(DensData.1.12,direction="wide",idvar = "Bag",timevar="Day")
DensData.1.12$ChangeLogDens<-DensData.1.12$LogTotDens.12-DensData.1.12$LogTotDens.1
ggplot(DensData.1.12)+
  geom_boxplot(aes(x=Spores.1,y=ChangeLogDens,color=Nutrients.1),notch=TRUE)+
  geom_point(aes(x=Spores.1,y=ChangeLogDens,color=Nutrients.1),size=0.5,position = position_jitterdodge(),alpha=0.5)+ 
  scale_color_manual(values=NutColors,name="Nutrients")+
  xlab("Spores")+
  ylab("Change in Log (Total host density, #/m^2)")+
  labs(title="Change from Day 1 to 12")+
  fav_theme


DensData.1.15<-subset(DensDataStartEnd,Day==1|Day==15)
DensData.1.15<-reshape(DensData.1.15,direction="wide",idvar = "Bag",timevar="Day")
DensData.1.15$ChangeLogDens<-DensData.1.15$LogTotDens.15-DensData.1.15$LogTotDens.1
ggplot(DensData.1.15)+
  geom_boxplot(aes(x=Spores.1,y=ChangeLogDens,color=Nutrients.1),notch=TRUE)+
  geom_point(aes(x=Spores.1,y=ChangeLogDens,color=Nutrients.1),size=0.5,position = position_jitterdodge(),alpha=0.5)+ 
  scale_color_manual(values=NutColors,name="Nutrients")+
  xlab("Spores")+
  ylab("Change in Log (Total host density, #/m^2)")+
  labs(title="Change from Day 1 to 15")+
  fav_theme

ggplot(DensData.1.15)+
  geom_boxplot(aes(x=Nutrients.1,y=ChangeLogDens,color=Nutrients.1),notch=TRUE)+
  geom_point(aes(x=Nutrients.1,y=ChangeLogDens,color=Nutrients.1),size=0.5,position = position_jitterdodge(),alpha=0.5)+ 
  scale_color_manual(values=NutColors,name="Nutrients")+
  ylab("Change in Log (Total host density, #/m^2)")+
  labs(title="Change from Day 1 to 15")+
  fav_theme

DensData.1.19<-subset(DensDataStartEnd,Day==1|Day==19)
DensData.1.19<-reshape(DensData.1.19,direction="wide",idvar = "Bag",timevar="Day")
DensData.1.19$ChangeLogDens<-DensData.1.19$LogTotDens.19-DensData.1.19$LogTotDens.1
ggplot(DensData.1.19)+
  geom_boxplot(aes(x=Spores.1,y=ChangeLogDens,color=Nutrients.1),notch=TRUE)+
  geom_point(aes(x=Spores.1,y=ChangeLogDens,color=Nutrients.1),size=0.5,position = position_jitterdodge(),alpha=0.5)+ 
  scale_color_manual(values=NutColors,name="Nutrients")+
  xlab("Spores")+
  ylab("Change in Log (Total host density, #/m^2)")+
  labs(title="Change from Day 1 to 19")+
  fav_theme

ggplot(DensData.1.19)+
  geom_boxplot(aes(x=Mixing.1,y=ChangeLogDens,color=Nutrients.1),notch=TRUE)+
  geom_point(aes(x=Mixing.1,y=ChangeLogDens,color=Nutrients.1),size=0.5,position = position_jitterdodge(),alpha=0.5)+ 
  scale_color_manual(values=NutColors,name="Nutrients")+
  ylab("Change in Log (Total host density, #/m^2)")+
  labs(title="Change from Day 1 to 19")+
  fav_theme

DensData.1.26<-subset(DensDataStartEnd,Day==1|Day==26)
DensData.1.26<-reshape(DensData.1.26,direction="wide",idvar = "Bag",timevar="Day")
DensData.1.26$ChangeLogDens<-DensData.1.26$LogTotDens.26-DensData.1.26$LogTotDens.1
ggplot(DensData.1.26)+
  geom_boxplot(aes(x=Spores.1,y=ChangeLogDens,color=Nutrients.1))+
  geom_point(aes(x=Spores.1,y=ChangeLogDens,color=Nutrients.1),size=0.5,position = position_jitterdodge(),alpha=0.5)+ 
  scale_color_manual(values=NutColors,name="Nutrients")+
  xlab("Spores")+
  ylab("Change in Log (Total host density, #/m^2)")+
  labs(title="Change from Day 1 to 26")+
  fav_theme

DensData.1.46<-subset(DensDataStartEnd,Day==1|Day==46)
DensData.1.46<-reshape(DensData.1.46,direction="wide",idvar = "Bag",timevar="Day")
DensData.1.46$ChangeLogDens<-DensData.1.46$LogTotDens.46-DensData.1.46$LogTotDens.1
ggplot(DensData.1.46)+
  geom_boxplot(aes(x=Spores.1,y=ChangeLogDens,color=Nutrients.1))+
  geom_point(aes(x=Spores.1,y=ChangeLogDens,color=Nutrients.1),size=0.5,position = position_jitterdodge(),alpha=0.5)+ 
  scale_color_manual(values=NutColors,name="Nutrients")+
  xlab("Spores")+
  ylab("Change in Log (Total host density, #/m^2)")+
  labs(title="Change from Day 1 to 46")+
  fav_theme

m.1.46.1<-lm(ChangeLogDens~Nutrients.1*Mixing.1*Spores.1+Raft.1,data=DensData.1.46)
Anova(m.1.46.1,type=3)
summary(m.1.46.1)

m.1.46.2<-update(m.1.46.1,~.-Nutrients.1:Mixing.1:Spores.1) 
Anova(m.1.46.2,type=3)

m.1.46.3<-update(m.1.46.2,~.-Nutrients.1:Spores.1) 
Anova(m.1.46.3,type=3)

m.1.46.4<-update(m.1.46.3,~.-Mixing.1:Spores.1) 
Anova(m.1.46.4,type=3)

m.1.46.5<-update(m.1.46.4,~.-Nutrients.1:Mixing.1) 
Anova(m.1.46.5,type=3)

DensData.29.46<-subset(DensDataStartEnd,Day==29|Day==46)
DensData.29.46<-reshape(DensData.29.46,direction="wide",idvar = "Bag",timevar="Day")
DensData.29.46$ChangeLogDens<-DensData.29.46$LogTotDens.46-DensData.29.46$LogTotDens.29
ggplot(DensData.29.46)+
  geom_boxplot(aes(x=Spores.29,y=ChangeLogDens,color=Nutrients.29))+
  geom_point(aes(x=Spores.29,y=ChangeLogDens,color=Nutrients.29),size=0.5,position = position_jitterdodge(),alpha=0.5)+ 
  scale_color_manual(values=NutColors,name="Nutrients")+
  xlab("Spores")+
  ylab("Change in Log (Total host density, #/m^2)")+
  labs(title="Change from Day 29 to 46")+
  fav_theme

m.29.46.1<-lm(ChangeLogDens~Nutrients.29*Mixing.29*Spores.29+Raft.29,data=DensData.29.46)
Anova(m.29.46.1,type=3)
summary(m.29.46.1)

m.29.46.2<-update(m.29.46.1,~.-Nutrients.29:Mixing.29:Spores.29) 
Anova(m.29.46.2,type=3)

m.29.46.3<-update(m.29.46.2,~.-Mixing.29:Spores.29) 
Anova(m.29.46.3,type=3)

m.29.46.4<-update(m.29.46.3,~.-Nutrients.29:Spores.29) 
Anova(m.29.46.4,type=3)

m.29.46.5<-update(m.29.46.4,~.-Nutrients.29:Mixing.29) 
Anova(m.29.46.5,type=3)
summary(m.29.46.5)

DensData.1.43<-subset(DensDataStartEnd,Day==1|Day==43)
DensData.1.43<-reshape(DensData.1.43,direction="wide",idvar = "Bag",timevar="Day")
DensData.1.43$ChangeLogDens<-DensData.1.43$LogTotDens.43-DensData.1.43$LogTotDens.1
ggplot(DensData.1.43)+
  geom_boxplot(aes(x=Spores.1,y=ChangeLogDens,color=Nutrients.1))+
  geom_point(aes(x=Spores.1,y=ChangeLogDens,color=Nutrients.1),size=0.5,position = position_jitterdodge(),alpha=0.5)+ 
  scale_color_manual(values=NutColors,name="Nutrients")+
  xlab("Spores")+
  ylab("Change in Log (Total host density, #/m^2)")+
  labs(title="Change from Day 1 to 43")+
  fav_theme

m.1.43.1<-lm(ChangeLogDens~Nutrients.1*Mixing.1*Spores.1+Raft.1,data=DensData.1.43)
Anova(m.1.43.1,type=3)
summary(m.1.43.1)

m.1.43.2<-update(m.1.43.1,~.-Nutrients.1:Mixing.1:Spores.1) 
Anova(m.1.43.2,type=3)

m.1.43.3<-update(m.1.43.2,~.-Nutrients.1:Mixing.1) 
Anova(m.1.43.3,type=3)

m.1.43.4<-update(m.1.43.3,~.-Mixing.1:Spores.1) 
Anova(m.1.43.4,type=3)

m.1.43.5<-update(m.1.43.4,~.-Nutrients.1:Spores.1) 
Anova(m.1.43.5,type=3)
summary(m.1.43.5)

DensData.29.43<-subset(DensDataStartEnd,Day==29|Day==43)
DensData.29.43<-reshape(DensData.29.43,direction="wide",idvar = "Bag",timevar="Day")
DensData.29.43$ChangeLogDens<-DensData.29.43$LogTotDens.43-DensData.29.43$LogTotDens.29
ggplot(DensData.29.43)+
  geom_boxplot(aes(x=Spores.29,y=ChangeLogDens,color=Nutrients.29),notch=FALSE)+
  geom_point(aes(x=Spores.29,y=ChangeLogDens,color=Nutrients.29),size=0.5,position = position_jitterdodge(),alpha=0.5)+ 
  scale_color_manual(values=NutColors,name="Nutrients")+
  xlab("Spores")+
  ylab("Change in Log (Total host density, #/m^2)")+
  labs(title="Change from Day 29 to 43")+
  fav_theme

m.29.43.29<-lm(ChangeLogDens~Nutrients.29*Mixing.29*Spores.29+Raft.29,data=DensData.29.43)
Anova(m.29.43.29,type=3)
summary(m.29.43.29)

m.29.43.2<-update(m.29.43.29,~.-Nutrients.29:Mixing.29:Spores.29) 
Anova(m.29.43.2,type=3)

m.29.43.3<-update(m.29.43.2,~.-Mixing.29:Spores.29) 
Anova(m.29.43.3,type=3)

m.29.43.4<-update(m.29.43.3,~.-Nutrients.29:Mixing.29) 
Anova(m.29.43.4,type=3)

m.29.43.5<-update(m.29.43.4,~.-Nutrients.29:Spores.29) 
Anova(m.29.43.5,type=3)
summary(m.29.43.5)

DensData.26.46<-subset(DensDataStartEnd,Day==26|Day==46)
DensData.26.46<-reshape(DensData.26.46,direction="wide",idvar = "Bag",timevar="Day")
DensData.26.46$ChangeLogDens<-DensData.26.46$LogTotDens.46-DensData.26.46$LogTotDens.26
ggplot(DensData.26.46)+
  geom_boxplot(aes(x=Spores.26,y=ChangeLogDens,color=Nutrients.26))+
  geom_point(aes(x=Spores.26,y=ChangeLogDens,color=Nutrients.26),size=0.5,position = position_jitterdodge(),alpha=0.5)+ 
  scale_color_manual(values=NutColors,name="Nutrients")+
  xlab("Spores")+
  ylab("Change in Log (Total host density, #/m^2)")+
  labs(title="Change from Day 26 to 46")+
  fav_theme

m.26.46.1<-lm(ChangeLogDens~Nutrients.26*Mixing.26*Spores.26+Raft.26,data=DensData.26.46)
Anova(m.26.46.1,type=3)
summary(m.26.46.1)

m.26.46.2<-update(m.26.46.1,~.-Nutrients.26:Mixing.26:Spores.26) 
Anova(m.26.46.2,type=3)

m.26.46.3<-update(m.26.46.2,~.-Nutrients.26:Spores.26) 
Anova(m.26.46.3,type=3)

m.26.46.4<-update(m.26.46.3,~.-Nutrients.26:Mixing.26) 
Anova(m.26.46.4,type=3)

m.26.46.5<-update(m.26.46.4,~.-Mixing.26:Spores.26) 
Anova(m.26.46.5,type=3)
summary(m.26.46.5)

DensData.26.43<-subset(DensDataStartEnd,Day==26|Day==43)
DensData.26.43<-reshape(DensData.26.43,direction="wide",idvar = "Bag",timevar="Day")
DensData.26.43$ChangeLogDens<-DensData.26.43$LogTotDens.43-DensData.26.43$LogTotDens.26
ggplot(DensData.26.43)+
  geom_boxplot(aes(x=Spores.26,y=ChangeLogDens,color=Nutrients.26))+
  geom_point(aes(x=Spores.26,y=ChangeLogDens,color=Nutrients.26),size=0.5,position = position_jitterdodge(),alpha=0.5)+ 
  scale_color_manual(values=NutColors,name="Nutrients")+
  xlab("Spores")+
  ylab("Change in Log (Total host density, #/m^2)")+
  labs(title="Change from Day 26 to 43")+
  fav_theme

m.26.43.1<-lm(ChangeLogDens~Nutrients.26*Mixing.26*Spores.26+Raft.26,data=DensData.26.43)
Anova(m.26.43.1,type=3)
summary(m.26.43.1)

m.26.43.2<-update(m.26.43.1,~.-Nutrients.26:Mixing.26:Spores.26) 
Anova(m.26.43.2,type=3)

m.26.43.3<-update(m.26.43.2,~.-Nutrients.26:Spores.26) 
Anova(m.26.43.3,type=3)

m.26.43.4<-update(m.26.43.3,~.-Mixing.26:Spores.26) 
Anova(m.26.43.4,type=3)

m.26.43.5<-update(m.26.43.4,~.-Nutrients.26:Mixing.26) 
Anova(m.26.43.5,type=3)
summary(m.26.43.5)

##################################################################################################
## Analyses of total and infected host density, starting from day 1, 29, 32, 35, 40
##################################################################################################
## Re-subset the data here, otherwise we're left with whatever was the final subset in the for-loop above
## Subset data 
DensData<-mydata[,c("Raft","Nutrients","Mixing","Spores","Trt","Bag","Day","NDay",
                    "TotalDensity","UninfDensity","InfDensity","LogTotDens","LogUninfDens","LogInfDens","LogTotChl","LogEdChl","LogTP","LogTN")]
## Omit missing data
DensData<-na.omit(DensData)

## Subset data for -spores and +spores bags
DensDataNoSpores<-subset(DensData,Spores=="No")
DensDataSpores<-subset(DensData,Spores=="Spores")

##################################################################################################
## Full time series (day 1-25)
EndDay25<-subset(DensData,Day<26) # all bags
EndDay25Spores<-subset(DensDataSpores,Day<26) # +spores only

ggplot(EndDay25)+
  geom_boxplot(aes(x=Spores,y=LogTotDens,color=Nutrients),notch=TRUE)+
  geom_point(aes(x=Spores,y=LogTotDens,color=Nutrients),size=0.5,position = position_jitterdodge(),alpha=0.5)+ 
  scale_color_manual(values=NutColors)+
  facet_wrap(~Mixing)+
  ylab("Log (Total host density, #/m^2)")+
  labs(title="Day 1-25")+
  fav_theme

ggplot(EndDay25)+
  geom_boxplot(aes(x=Spores,y=LogTotDens,color=Nutrients),notch=TRUE)+
  geom_point(aes(x=Spores,y=LogTotDens,color=Nutrients),size=0.5,position = position_jitterdodge(),alpha=0.5)+ 
  scale_color_manual(values=NutColors)+
  ylab("Log (Total host density, #/m^2)")+
  labs(title="Day 1-25")+
  fav_theme

##################################################################################################
## Full time series (day 1-46)

ggplot(DensData)+
  geom_boxplot(aes(x=Spores,y=LogTotDens,color=Nutrients),notch=TRUE)+
  geom_point(aes(x=Spores,y=LogTotDens,color=Nutrients),size=0.5,position = position_jitterdodge(),alpha=0.5)+ 
  scale_color_manual(values=NutColors)+
  facet_wrap(~Mixing)+
  ylab("Log (Total host density, #/m^2)")+
  labs(title="Day 1-46")+
  fav_theme

ggplot(DensData)+
  geom_boxplot(aes(x=Spores,y=LogTotDens,color=Nutrients),notch=TRUE)+
  geom_point(aes(x=Spores,y=LogTotDens,color=Nutrients),size=0.5,position = position_jitterdodge(),alpha=0.5)+ 
  scale_color_manual(values=NutColors)+
  ylab("Log (Total host density, #/m^2)")+
  labs(title="Day 1-46")+
  fav_theme

m1.1<-lme(LogTotDens~Nutrients*Mixing*Spores+Raft,random=~1|Bag,method="REML",data=DensData) 
summary(m1.1)
Anova(m1.1,type=3)

m1.2<-update(m1.1,~.-Nutrients:Mixing:Spores)
Anova(m1.2,type=3)

m1.3<-update(m1.2,~.-Mixing:Spores)
Anova(m1.3,type=3)

m1.4<-update(m1.3,~.-Nutrients:Mixing)
Anova(m1.4,type=3)

m1.5<-update(m1.4,~.-Nutrients:Spores)
Anova(m1.5,type=3)
summary(m1.5)

##################################################################################################
## Subset data to begin at Day 26
StartDay26<-subset(DensData,Day>=26) # all bags
StartDay26Spores<-subset(DensDataSpores,Day>=26) # +spores only

ggplot(StartDay26Spores)+
  geom_boxplot(aes(x=Mixing,y=LogInfDens,color=Nutrients),notch=TRUE)+
  geom_point(aes(x=Mixing,y=LogInfDens,color=Nutrients),size=0.5,position = position_jitterdodge(),alpha=0.5)+ 
  scale_color_manual(values=NutColors)+
  ylim(c(0,12))+
  ylab("Log (Infected host density, #/m^2)")+
  labs(title="Day 26-46")+
  fav_theme

StartDayEquals32Spores<-subset(DensDataSpores,Day==32) # +spores only
ggplot(StartDayEquals32Spores)+
  geom_boxplot(aes(x=Mixing,y=LogInfDens,color=Nutrients),notch=FALSE)+
  geom_point(aes(x=Mixing,y=LogInfDens,color=Nutrients),size=0.5,position = position_jitterdodge(),alpha=0.5)+ 
  scale_color_manual(values=NutColors)+
  ylab("Log (Infected host density, #/m^2)")+
  labs(title="Day 32")+
  fav_theme

StartDayEquals35Spores<-subset(DensDataSpores,Day==35) # +spores only
ggplot(StartDayEquals35Spores)+
  geom_boxplot(aes(x=Mixing,y=LogInfDens,color=Nutrients),notch=FALSE)+
  geom_point(aes(x=Mixing,y=LogInfDens,color=Nutrients),size=0.5,position = position_jitterdodge(),alpha=0.5)+ 
  scale_color_manual(values=NutColors)+
  ylab("Log (Infected host density, #/m^2)")+
  labs(title="Day 35")+
  fav_theme

m35.1<-lm(LogInfDens~Nutrients*Mixing+Raft,data=StartDayEquals35Spores) 
summary(m35.1)
Anova(m35.1,type=3)

m35.2<-lm(LogInfDens~Nutrients+Mixing+Raft,data=StartDayEquals35Spores) 
summary(m35.2)
Anova(m35.2,type=3)


##################################################################################################
## Subset data to begin at Day 29
StartDay29<-subset(DensData,Day>=29) # all bags
StartDay29Spores<-subset(DensDataSpores,Day>=29) # +spores only

SporesNames<-c("No"="-Spores","Spores"="+Spores")
ggplot(StartDay29)+
  geom_boxplot(aes(x=factor(Day),y=LogTotDens,color=Nutrients))+
  geom_point(aes(x=factor(Day),y=LogTotDens,color=Nutrients),size=0.5,position = position_jitterdodge(),alpha=0.5)+ 
  facet_wrap(~Spores,labeller = labeller(Spores = SporesNames))+
  scale_color_manual(values=NutColors)+
  ylab("Log (Total host density, #/m^2)")+
  fav_theme

ggplot(StartDay29)+
  geom_boxplot(aes(x=Spores,y=LogTotDens,color=Nutrients),notch=TRUE)+
  geom_point(aes(x=Spores,y=LogTotDens,color=Nutrients),size=0.5,position = position_jitterdodge(),alpha=0.5)+ 
  scale_color_manual(values=NutColors)+
  ylab("Log (Total host density, #/m^2)")+
  labs(title="Day 29-46")+
  fav_theme

m29.1<-lme(LogTotDens~Nutrients*Mixing*Spores+Raft,random=~1|Bag,method="REML",data=StartDay29) 
summary(m29.1)
Anova(m29.1,type=3)

m29.2<-update(m29.1,~.-Nutrients:Mixing:Spores)
Anova(m29.2,type=3)

m29.3<-update(m29.2,~.-Mixing:Spores)
Anova(m29.3,type=3)

m29.4<-update(m29.3,~.-Nutrients:Spores)
Anova(m29.4,type=3)

m29.5<-update(m29.4,~.-Nutrients:Mixing)
Anova(m29.5,type=3)
summary(m29.5)

ggplot(StartDay29Spores)+
  geom_boxplot(aes(x=Mixing,y=LogInfDens,color=Nutrients),notch=TRUE)+
  geom_point(aes(x=Mixing,y=LogInfDens,color=Nutrients),size=0.5,position = position_jitterdodge(),alpha=0.5)+ 
  scale_color_manual(values=NutColors)+
  ylab("Log (Infected host density, #/m^2)")+
  labs(title="Day 29-46")+
  fav_theme

ggplot(StartDay29Spores)+
  geom_boxplot(aes(x=Nutrients,y=LogInfDens,color=Nutrients),notch=TRUE)+
  geom_point(aes(x=Nutrients,y=LogInfDens,color=Nutrients),size=0.5,position = position_jitter(),alpha=0.5)+ 
  scale_color_manual(values=NutColors)+
  ylab("Log (Infected host density, #/m^2)")+
  labs(title="Day 29-46")+
  fav_theme

m29.1<-lme(LogInfDens~Nutrients*Mixing+Raft,random=~1|Bag,method="REML",data=StartDay29Spores) 
summary(m29.1)
Anova(m29.1,type=3)

m29.2<-update(m29.1,~.-Nutrients:Mixing)
Anova(m29.2,type=3)
summary(m29.2)

##################################################################################################
## Subset data to begin at Day 32
StartDay32<-subset(DensData,Day>=32) # all bags
StartDay32Spores<-subset(DensDataSpores,Day>=32) # +spores only

SporesNames<-c("No"="-Spores","Spores"="+Spores")
ggplot(StartDay32)+
  geom_boxplot(aes(x=factor(Day),y=LogTotDens,color=Nutrients))+
  geom_point(aes(x=factor(Day),y=LogTotDens,color=Nutrients),size=0.5,position = position_jitterdodge(),alpha=0.5)+ 
  facet_wrap(~Spores,labeller = labeller(Spores = SporesNames))+
  scale_color_manual(values=NutColors)+
  ylab("Log (Total host density, #/m^2)")+
  fav_theme

ggplot(StartDay32)+
  geom_boxplot(aes(x=Spores,y=LogTotDens,color=Nutrients),notch=TRUE)+
  geom_point(aes(x=Spores,y=LogTotDens,color=Nutrients),size=0.5,position = position_jitterdodge(),alpha=0.5)+ 
  scale_color_manual(values=NutColors)+
  ylab("Log (Total host density, #/m^2)")+
  labs(title="Day 32-46")+
  fav_theme

m32.1<-lme(LogTotDens~Nutrients*Mixing*Spores+Raft,random=~1|Bag,method="REML",data=StartDay32) 
summary(m32.1)
Anova(m32.1,type=3)

m32.2<-update(m32.1,~.-Nutrients:Mixing:Spores)
Anova(m32.2,type=3)

m32.3<-update(m32.2,~.-Mixing:Spores)
Anova(m32.3,type=3)

m32.4<-update(m32.3,~.-Nutrients:Spores)
Anova(m32.4,type=3)

m32.5<-update(m32.4,~.-Nutrients:Mixing)
Anova(m32.5,type=3)
summary(m32.5)

# m32.6<-update(m32.5,~.-Mixing)
# Anova(m32.6,type=3)

ggplot(StartDay32Spores)+
  geom_boxplot(aes(x=Mixing,y=LogInfDens,color=Nutrients),notch=TRUE)+
  geom_point(aes(x=Mixing,y=LogInfDens,color=Nutrients),size=0.5,position = position_jitterdodge(),alpha=0.5)+ 
  scale_color_manual(values=NutColors)+
  ylim(c(0,12))+
  ylab("Log (Infected host density, #/m^2)")+
  labs(title="Day 32-46")+
  fav_theme

ggplot(StartDay32Spores)+
  geom_boxplot(aes(x=Nutrients,y=LogInfDens,color=Nutrients),notch=TRUE)+
  geom_point(aes(x=Nutrients,y=LogInfDens,color=Nutrients),size=0.5,position = position_jitter(),alpha=0.5)+ 
  scale_color_manual(values=NutColors)+
  ylim(c(0,12))+
  ylab("Log (Infected host density, #/m^2)")+
  labs(title="Day 32-46")+
  fav_theme

##################################################################################################
## Subset data to begin at Day 35
StartDay35<-subset(DensData,Day>=35) # all bags
StartDay35Spores<-subset(DensDataSpores,Day>=35) # +spores only

SporesNames<-c("No"="-Spores","Spores"="+Spores")
ggplot(StartDay35)+
  geom_boxplot(aes(x=factor(Day),y=LogTotDens,color=Nutrients))+
  facet_wrap(~Spores,labeller = labeller(Spores = SporesNames))+
  scale_color_manual(values=NutColors)+
  ylab("Log (Total host density, #/m^2)")+
  fav_theme

ggplot(StartDay35)+
  geom_boxplot(aes(x=Spores,y=LogTotDens,color=Nutrients),notch=TRUE)+
  geom_point(aes(x=Spores,y=LogTotDens,color=Nutrients),size=0.5,position = position_jitterdodge(),alpha=0.5)+ 
  scale_color_manual(values=NutColors)+
  ylab("Log (Total host density, #/m^2)")+
  labs(title="Day 35-46")+
  fav_theme

m35.1<-lme(LogTotDens~Nutrients*Mixing*Spores+Raft,random=~1|Bag,method="REML",data=StartDay35) 
summary(m35.1)
Anova(m35.1,type=3)

m35.2<-update(m35.1,~.-Nutrients:Mixing:Spores)
Anova(m35.2,type=3)

m35.3<-update(m35.2,~.-Mixing:Spores)
Anova(m35.3,type=3)

m35.4<-update(m35.3,~.-Nutrients:Mixing)
Anova(m35.4,type=3)

m35.5<-update(m35.4,~.-Nutrients:Spores)
Anova(m35.5,type=3)
summary(m35.5)

# m35.6<-update(m35.5,~.-Mixing)
# Anova(m35.6,type=3)

ggplot(StartDay35Spores)+
  geom_boxplot(aes(x=Nutrients,y=LogInfDens,color=Nutrients),notch=TRUE)+
  geom_point(aes(x=Nutrients,y=LogInfDens,color=Nutrients),size=0.5,position = position_jitter(),alpha=0.5)+ 
  scale_color_manual(values=NutColors)+
  ylim(c(0,12))+
  ylab("Log (Infected host density, #/m^2)")+
  labs(title="Day 35-46")+
  fav_theme

##################################################################################################
## Subset data to begin at Day 40
StartDay40<-subset(DensData,Day>=40) # all bags
StartDay40Spores<-subset(DensDataSpores,Day>=40) # +spores only

SporesNames<-c("No"="-Spores","Spores"="+Spores")
ggplot(StartDay40)+
  geom_boxplot(aes(x=factor(Day),y=LogTotDens,color=Nutrients))+
  facet_wrap(~Spores,labeller = labeller(Spores = SporesNames))+
  scale_color_manual(values=NutColors)+
  ylab("Log (Total host density, #/m^2)")+
  fav_theme

ggplot(StartDay40)+
  geom_boxplot(aes(x=Spores,y=LogTotDens,color=Nutrients),notch=TRUE)+
  geom_point(aes(x=Spores,y=LogTotDens,color=Nutrients),size=0.5,position = position_jitterdodge(),alpha=0.5)+ 
  scale_color_manual(values=NutColors)+
  ylab("Log (Total host density, #/m^2)")+
  labs(title="Day 40-46")+
  fav_theme

m40.1<-lme(LogTotDens~Nutrients*Mixing*Spores+Raft,random=~1|Bag,method="REML",data=StartDay40) 
summary(m40.1)
Anova(m40.1,type=3)

m40.2<-update(m40.1,~.-Nutrients:Mixing:Spores)
Anova(m40.2,type=3)

m40.3<-update(m40.2,~.-Mixing:Spores)
Anova(m40.3,type=3)

m40.4<-update(m40.3,~.-Nutrients:Mixing)
Anova(m40.4,type=3)

m40.5<-update(m40.4,~.-Nutrients:Spores)
Anova(m40.5,type=3)
summary(m40.5)

ggplot(StartDay40Spores)+
  geom_boxplot(aes(x=Mixing,y=LogInfDens,color=Nutrients),notch=TRUE)+
  geom_point(aes(x=Mixing,y=LogInfDens,color=Nutrients),size=0.5,position = position_jitterdodge(),alpha=0.5)+ 
  scale_color_manual(values=NutColors)+
  ylab("Log (Infected host density, #/m^2)")+
  labs(title="Day 40-46")+
  fav_theme

m40.1<-lme(LogInfDens~Nutrients*Mixing+Raft,random=~1|Bag,method="REML",data=StartDay40Spores) 
summary(m40.1)
Anova(m40.1,type=3)

m40.2<-update(m40.1,~.-Nutrients:Mixing)
Anova(m40.2,type=3)
summary(m40.2)

ggplot(StartDay40Spores)+
  geom_boxplot(aes(x=Nutrients,y=LogInfDens,color=Nutrients),notch=TRUE)+
  geom_point(aes(x=Nutrients,y=LogInfDens,color=Nutrients),size=0.5,position = position_jitterdodge(),alpha=0.5)+ 
  scale_color_manual(values=NutColors)+
  ylab("Log (Infected host density, #/m^2)")+
  labs(title="Day 40-46")+
  fav_theme

ggplot(StartDay40Spores)+
  geom_point(aes(x=LogEdChl,y=LogInfDens,color=Nutrients),size=0.5,alpha=0.5)+ 
  scale_color_manual(values=NutColors)+
  ylab("Log (Infected host density, #/m^2)")+
  labs(title="Day 40-46")+
  fav_theme

ggplot(StartDay40Spores)+
  geom_point(aes(x=LogTP,y=LogInfDens,color=Nutrients),size=0.5,alpha=0.5)+ 
  scale_color_manual(values=NutColors)+
  ylab("Log (Infected host density, #/m^2)")+
  labs(title="Day 40-46")+
  fav_theme

##################################################################################################
## Analyze time series of Log(Total host density), starting from Day 1
## There are 15 data points from Day 1 (first sampling day) to Day 46 (last sampling day)
## In this analysis, NDay = Days numbered 1 through 15 (i.e., pretending the time intervals are evenly spaced)
##################################################################################################

## Full model with no autocorrelation structure
TotDensMod.p0.1<-lme(LogTotDens~Raft+Nutrients*Mixing*Spores*factor(NDay),random=~1|Bag,method="ML",data=DensData) 
summary(TotDensMod.p0.1)

ACF(TotDensMod.p0.1)
plot(ACF(TotDensMod.p0.1),alpha=0.05) # suggests significant positive autocorrelation at time lags 1 and 2...but then increasingly negative autocorrelation at greater time lags? With a time series of only 15 points, obviously, we don't have enough data to actually allow for very large lags.

## Full model with AR(2) autocorrelation structure, followed by stepwise simplification
TotDensMod.p2.1<-lme(LogTotDens~Raft+Nutrients*Mixing*Spores*factor(NDay),random=~1|Bag,
                     correlation=corARMA(form=~NDay,p=2),method="ML",data=DensData) 

## In the next line, I'm just checking whether it matters if corARMA form is ~NDay or ~NDay|Bag in our analysis (it shouldn't matter, since all bags were sampled on the same days)
TotDensMod.p2.1.NDayWithinBag<-lme(LogTotDens~Raft+Nutrients*Mixing*Spores*factor(NDay),random=~1|Bag,
                                  correlation=corARMA(form=~NDay|Bag,p=2),method="ML",data=DensData) 

summary(TotDensMod.p2.1)
summary(TotDensMod.p2.1.NDayWithinBag) # all parameter estimates are the same as in TotDensMod.p2.1, which confirms that it is not necessary for us to code corARMA form as ~NDay|Bag
Anova(TotDensMod.p2.1,type=3)

plot(ACF(TotDensMod.p2.1,resType="normalized"),alpha=0.05) # Whoa, I don't understand why the autocorrelation at the first few time lags is WORSE in the model accounting for autocorrelation  

anova(TotDensMod.p2.1,TotDensMod.p0.1) # Is this a valid thing to do, to test whether the corARMA model is better than the non-corARMA model? If it is valid, it suggests that we should use the corARMA model.

TotDensMod.p2.2<-update(TotDensMod.p2.1,~.-Nutrients:Mixing:Spores:factor(NDay)) # 4-way interaction removed
anova(TotDensMod.p2.1,TotDensMod.p2.2) # 4-way interaction p > 0.05, so proceeding with removing it
Anova(TotDensMod.p2.2,type=3)

TotDensMod.p2.3<-update(TotDensMod.p2.2,~.-Nutrients:Mixing:factor(NDay)) # least significant 3-way interaction removed
anova(TotDensMod.p2.2,TotDensMod.p2.3) 
Anova(TotDensMod.p2.3,type=3) 

TotDensMod.p2.4<-update(TotDensMod.p2.3,~.-Nutrients:Mixing:Spores) # next least significant 3-way interaction removed
anova(TotDensMod.p2.3,TotDensMod.p2.4) 
Anova(TotDensMod.p2.4,type=3) 

TotDensMod.p2.5<-update(TotDensMod.p2.4,~.-Mixing:Spores:factor(NDay)) # next least significant 3-way interaction removed
anova(TotDensMod.p2.4,TotDensMod.p2.5) 
Anova(TotDensMod.p2.5,type=3) 

TotDensMod.p2.6<-update(TotDensMod.p2.5,~.-Nutrients:Spores:factor(NDay)) # next least significant 3-way interaction removed
anova(TotDensMod.p2.5,TotDensMod.p2.6) 
Anova(TotDensMod.p2.6,type=3) 

TotDensMod.p2.7<-update(TotDensMod.p2.6,~.-Mixing:Spores) # least significant 2-way interaction removed
anova(TotDensMod.p2.6,TotDensMod.p2.7) 
Anova(TotDensMod.p2.7,type=3)

TotDensMod.p2.8<-update(TotDensMod.p2.7,~.-Nutrients:Mixing) # next least significant 2-way interaction removed
anova(TotDensMod.p2.7,TotDensMod.p2.8) 
Anova(TotDensMod.p2.8,type=3)

TotDensMod.p2.9<-update(TotDensMod.p2.8,~.-Nutrients:Spores) # next least significant 2-way interaction removed
anova(TotDensMod.p2.8,TotDensMod.p2.9) 
Anova(TotDensMod.p2.9,type=3)

TotDensMod.p2.10<-update(TotDensMod.p2.9,~.-Nutrients:factor(NDay)) # next least significant 2-way interaction removed
anova(TotDensMod.p2.9,TotDensMod.p2.10) 
Anova(TotDensMod.p2.10,type=3)

TotDensMod.p2.11<-update(TotDensMod.p2.10,~.-Spores:factor(NDay)) # next least significant 2-way interaction removed
anova(TotDensMod.p2.10,TotDensMod.p2.11) 
Anova(TotDensMod.p2.11,type=3)

TotDensMod.p2.12<-update(TotDensMod.p2.11,~.-Mixing:factor(NDay)) # final 2-way interaction removed
anova(TotDensMod.p2.11,TotDensMod.p2.12) 
Anova(TotDensMod.p2.12,type=3)
summary(TotDensMod.p2.12)

plot(ACF(TotDensMod.p2.12,resType="normalized"),alpha=0.05) # Hmmm...I expected the ACF plot with resType="normalized" to show no remaining signif autocorrelation, but it appears there is still large autocorrelation, even beyond time lag 2.
pacf(residuals(TotDensMod.p2.12, retype="normalized")) # Or maybe the pacf is showing no correlation for lag values beyond 1?

##################################################################################################
## Analyze time series of Log(Total host density), starting from Day 1
## There are 15 data points from Day 1 (first sampling day) to Day 46 (last sampling day)
## In this analysis, intervals between values of Day are the actual intervals between sampling days.
##################################################################################################

## Note that I named the models in this section the same as in the previous section!!!
## The "p" in the model name signifies the time lag of the autocorrelation structure included in the model.
## If the model doesn't have corARMA specified, then the model has "p0" as part of its name.

## Full model with no autocorrelation structure
TotDensMod.p0.1<-lme(LogTotDens~Raft+Nutrients*Mixing*Spores*factor(Day),random=~1|Bag,method="ML",data=DensData) 
summary(TotDensMod.p0.1)
ACF(TotDensMod.p0.1)
plot(ACF(TotDensMod.p0.1),alpha=0.05) # suggests significant positive autocorrelation at time lags 1 and 2...but then increasingly negative autocorrelation at greater time lags? With a time series of only 15 points, obviously, we don't have enough data to actually allow for very large lags.

## Spoiler alert: stepwise simplification from full model with no autocorrelation structure specified yields this reduced model (after systematically deleting non-significant interaction terms):
TotDensMod.p0.11<-lme(LogTotDens~Raft+Nutrients+Mixing+Spores+factor(Day)+Mixing:factor(Day),random=~1|Bag,method="ML",data=DensData) 
Anova(TotDensMod.p0.11,type=3)
# Analysis of Deviance Table (Type III tests)
# 
# Response: LogTotDens
# Chisq Df Pr(>Chisq)    
# (Intercept)        1221.6689  1    < 2e-16 ***
# Raft                  7.8930  3    0.04828 *  
# Nutrients             4.7494  1    0.02931 *  
# Mixing                0.0206  1    0.88578    
# Spores                5.6913  1    0.01705 *  
# factor(Day)         547.2889 14    < 2e-16 ***
# Mixing:factor(Day)   25.3014 14    0.03171 *  
summary(TotDensMod.p0.11)
# Random effects:
#   Formula: ~1 | Bag
#         (Intercept) Residual
# StdDev:   0.3137198 0.732411
ACF(TotDensMod.p0.11)
plot(ACF(TotDensMod.p0.11),alpha=0.05) # suggests significant autocorrelation at time lags 1 and 2

## Full model with AR(2) autocorrelation structure, followed by stepwise simplification
TotDensMod.p2.1<-lme(LogTotDens~Raft+Nutrients*Mixing*Spores*factor(Day),random=~1|Bag,
                correlation=corARMA(form=~Day,p=2),method="ML",data=DensData) 

TotDensMod.p2.1.DayWithinBag<-lme(LogTotDens~Raft+Nutrients*Mixing*Spores*factor(Day),random=~1|Bag,
                     correlation=corARMA(form=~Day|Bag,p=2),method="ML",data=DensData) # just confirming that it doesn't matter whether we specify corARMA(form=~Day|Bag,p=2) vs. corARMA(form=~Day,p=2)
summary(TotDensMod.p2.1)
summary(TotDensMod.p2.1.DayWithinBag)
Anova(TotDensMod.p2.1,type=3)

TotDensMod.p2.2<-update(TotDensMod.p2.1,~.-Nutrients:Mixing:Spores:factor(Day)) # 4-way interaction removed
anova(TotDensMod.p2.1,TotDensMod.p2.2) # 4-way interaction p > 0.05 (though not by much!!), so proceeding with removing it
Anova(TotDensMod.p2.2,type=3)

TotDensMod.p2.3<-update(TotDensMod.p2.2,~.-Nutrients:Mixing:factor(Day)) # least significant 3-way interaction removed
anova(TotDensMod.p2.2,TotDensMod.p2.3) 
Anova(TotDensMod.p2.3,type=3) 

TotDensMod.p2.4<-update(TotDensMod.p2.3,~.-Nutrients:Mixing:Spores) # next least significant 3-way interaction removed
anova(TotDensMod.p2.3,TotDensMod.p2.4) 
Anova(TotDensMod.p2.4,type=3) 

TotDensMod.p2.5<-update(TotDensMod.p2.4,~.-Mixing:Spores:factor(Day)) # next least significant 3-way interaction removed
anova(TotDensMod.p2.4,TotDensMod.p2.5) 
Anova(TotDensMod.p2.5,type=3) 

TotDensMod.p2.6<-update(TotDensMod.p2.5,~.-Nutrients:Spores:factor(Day)) # next least significant 3-way interaction removed
anova(TotDensMod.p2.5,TotDensMod.p2.6) 
Anova(TotDensMod.p2.6,type=3) 

TotDensMod.p2.7<-update(TotDensMod.p2.6,~.-Mixing:Spores) # least significant 2-way interaction removed
anova(TotDensMod.p2.6,TotDensMod.p2.7) 
Anova(TotDensMod.p2.7,type=3)

TotDensMod.p2.8<-update(TotDensMod.p2.7,~.-Nutrients:Mixing) # next least significant 2-way interaction removed
anova(TotDensMod.p2.7,TotDensMod.p2.8) 
Anova(TotDensMod.p2.8,type=3)

TotDensMod.p2.9<-update(TotDensMod.p2.8,~.-Nutrients:Spores) # next least significant 2-way interaction removed
anova(TotDensMod.p2.8,TotDensMod.p2.9) 
Anova(TotDensMod.p2.9,type=3)

TotDensMod.p2.10<-update(TotDensMod.p2.9,~.-Spores:factor(Day)) # next least significant 2-way interaction removed
anova(TotDensMod.p2.9,TotDensMod.p2.10) 
Anova(TotDensMod.p2.10,type=3)

TotDensMod.p2.11<-update(TotDensMod.p2.10,~.-Nutrients:factor(Day)) # next least significant 2-way interaction removed
anova(TotDensMod.p2.10,TotDensMod.p2.11) 
Anova(TotDensMod.p2.11,type=3)

TotDensMod.p2.12<-update(TotDensMod.p2.11,~.-Mixing:factor(Day)) # final 2-way interaction removed
anova(TotDensMod.p2.11,TotDensMod.p2.12) 
Anova(TotDensMod.p2.12,type=3)
summary(TotDensMod.p2.12)
plot(ACF(TotDensMod.p2.12),alpha=0.05) 
plot(ACF(TotDensMod.p2.12,resType="normalized"),alpha=0.05) # Hmmm...I expected the ACF plot with resType="normalized" to show no remaining signif autocorr.
pacf(residuals(TotDensMod.p2.12, retype="normalized")) # Or maybe the pacf is showing no correlation for lag values beyond 1? 

##################################################################################################
## Analyze time series of Log(Total host density), starting from Day 12 (start of 1st epidemic)
## There are 11 data points from Day 12 (start of epidemic) to Day 46 (last sampling day)
##################################################################################################

## Full model with no autocorrelation structure
TotDensEpiMod.p0.1<-lme(LogTotDens~Raft+Nutrients*Mixing*Spores*factor(Day),random=~1|Bag,method="ML",data=DensDataStartDay12) 
summary(TotDensEpiMod.p0.1)

ACF(TotDensEpiMod.p0.1)
plot(ACF(TotDensEpiMod.p0.1),alpha=0.05) # suggests significant autocorrelation at time lags 1 and 2

## Full model with AR(2) autocorrelation structure
TotDensEpiMod.p2.1<-lme(LogTotDens~Raft+Nutrients*Mixing*Spores*factor(Day),random=~1|Bag,
                    correlation=corARMA(form=~Day|Bag,p=2),method="ML",data=DensDataStartDay12) 
summary(TotDensEpiMod.p2.1) # but estimate of phi = 0 for time lag 1 and 2
anova(TotDensEpiMod.p0.1,TotDensEpiMod.p2.1) # if this is a legit model comparison, it suggests there's no significant AR(2) temporal autocorrelation

## Stepwise simplification from full model with no autocorrelation structure specified
TotDensEpiMod.p0.2<-update(TotDensEpiMod.p0.1,~.-Nutrients:Mixing:Spores:factor(Day)) # 4-way interaction removed
anova(TotDensEpiMod.p0.1,TotDensEpiMod.p0.2) 
Anova(TotDensEpiMod.p0.2,type=3) 

TotDensEpiMod.p0.3<-update(TotDensEpiMod.p0.2,~.-Nutrients:Mixing:factor(Day)) # 3-way interaction removed
anova(TotDensEpiMod.p0.2,TotDensEpiMod.p0.3) 
Anova(TotDensEpiMod.p0.3,type=3) 

TotDensEpiMod.p0.4<-update(TotDensEpiMod.p0.3,~.-Mixing:Spores:factor(Day)) # 3-way interaction removed
anova(TotDensEpiMod.p0.3,TotDensEpiMod.p0.4) 
Anova(TotDensEpiMod.p0.4,type=3) 

TotDensEpiMod.p0.5<-update(TotDensEpiMod.p0.4,~.-Nutrients:Spores:factor(Day)) # 3-way interaction removed
anova(TotDensEpiMod.p0.4,TotDensEpiMod.p0.5) 
Anova(TotDensEpiMod.p0.5,type=3) 

TotDensEpiMod.p0.6<-update(TotDensEpiMod.p0.5,~.-Nutrients:Mixing:Spores) # 3-way interaction removed
anova(TotDensEpiMod.p0.5,TotDensEpiMod.p0.6) 
Anova(TotDensEpiMod.p0.6,type=3) 

TotDensEpiMod.p0.7<-update(TotDensEpiMod.p0.6,~.-Nutrients:factor(Day)) # 2-way interaction removed
anova(TotDensEpiMod.p0.6,TotDensEpiMod.p0.7) 
Anova(TotDensEpiMod.p0.7,type=3) 

TotDensEpiMod.p0.8<-update(TotDensEpiMod.p0.7,~.-Spores:factor(Day)) # 2-way interaction removed
anova(TotDensEpiMod.p0.7,TotDensEpiMod.p0.8) 
Anova(TotDensEpiMod.p0.8,type=3) 

TotDensEpiMod.p0.9<-update(TotDensEpiMod.p0.8,~.-Nutrients:Mixing) # 2-way interaction removed
anova(TotDensEpiMod.p0.8,TotDensEpiMod.p0.9) 
Anova(TotDensEpiMod.p0.9,type=3) 

TotDensEpiMod.p0.10<-update(TotDensEpiMod.p0.9,~.-Mixing:Spores) # 2-way interaction removed
anova(TotDensEpiMod.p0.9,TotDensEpiMod.p0.10) 
Anova(TotDensEpiMod.p0.10,type=3) 

TotDensEpiMod.p0.11<-update(TotDensEpiMod.p0.10,~.-Nutrients:Spores) # 2-way interaction removed
anova(TotDensEpiMod.p0.10,TotDensEpiMod.p0.11) 
Anova(TotDensEpiMod.p0.11,type=3) 
# Analysis of Deviance Table (Type III tests)
# 
# Response: LogTotDens
# Chisq Df Pr(>Chisq)    
# (Intercept)        1494.1627  1    < 2e-16 ***
# Raft                  7.5975  3    0.05511 .  
# Nutrients             6.5925  1    0.01024 *  
# Mixing                0.2251  1    0.63515    
# Spores                6.4868  1    0.01087 *  
# factor(Day)         175.6739 10    < 2e-16 ***
# Mixing:factor(Day)   20.9612 10    0.02137 *  
summary(TotDensEpiMod.p0.11)
# Random effects:
#   Formula: ~1 | Bag
#         (Intercept) Residual
# StdDev:   0.3941676 0.762639

##################################################################################################
## Analyze time series of Log(Total host density), starting from Day 26 (start of 2nd epidemic)
## There are 7 data points from Day 26 (start of 2nd epidemic) to Day 46 (last sampling day)
##################################################################################################

## Full model with no autocorrelation structure
TotDensEpi2Mod.p0.1<-lme(LogTotDens~Raft+Nutrients*Mixing*Spores*factor(Day),random=~1|Bag,method="ML",data=DensDataStartDay26) 
summary(TotDensEpi2Mod.p0.1)

ACF(TotDensEpi2Mod.p0.1)
plot(ACF(TotDensEpi2Mod.p0.1),alpha=0.05)

## Full model with AR(1) autocorrelation structure
TotDensEpi2Mod.p1.1<-lme(LogTotDens~Raft+Nutrients*Mixing*Spores*factor(Day),random=~1|Bag,
                        correlation=corARMA(form=~Day|Bag,p=1),method="ML",data=DensDataStartDay26) 
summary(TotDensEpi2Mod.p1.1) # but estimate of phi = 0 for time lag 1
anova(TotDensEpi2Mod.p0.1,TotDensEpi2Mod.p1.1) # if this is a legit model comparison, it suggests there's no significant AR(1) temporal autocorrelation

## Stepwise simplification from full model with no autocorrelation structure specified
TotDensEpi2Mod.p0.2<-update(TotDensEpi2Mod.p0.1,~.-Nutrients:Mixing:Spores:factor(Day)) # 4-way interaction removed
anova(TotDensEpi2Mod.p0.1,TotDensEpi2Mod.p0.2) 
Anova(TotDensEpi2Mod.p0.2,type=3) 

TotDensEpi2Mod.p0.3<-update(TotDensEpi2Mod.p0.2,~.-Mixing:Spores:factor(Day)) # 3-way interaction removed
anova(TotDensEpi2Mod.p0.2,TotDensEpi2Mod.p0.3) 
Anova(TotDensEpi2Mod.p0.3,type=3) 

TotDensEpi2Mod.p0.4<-update(TotDensEpi2Mod.p0.3,~.-Nutrients:Mixing:factor(Day)) # 3-way interaction removed
anova(TotDensEpi2Mod.p0.3,TotDensEpi2Mod.p0.4) 
Anova(TotDensEpi2Mod.p0.4,type=3) 

TotDensEpi2Mod.p0.5<-update(TotDensEpi2Mod.p0.4,~.-Nutrients:Spores:factor(Day)) # 3-way interaction removed
anova(TotDensEpi2Mod.p0.4,TotDensEpi2Mod.p0.5) 
Anova(TotDensEpi2Mod.p0.5,type=3) 

TotDensEpi2Mod.p0.6<-update(TotDensEpi2Mod.p0.5,~.-Nutrients:Mixing:Spores) # 3-way interaction removed
anova(TotDensEpi2Mod.p0.5,TotDensEpi2Mod.p0.6) 
Anova(TotDensEpi2Mod.p0.6,type=3) 

TotDensEpi2Mod.p0.7<-update(TotDensEpi2Mod.p0.6,~.-Mixing:Spores) # 2-way interaction removed
anova(TotDensEpi2Mod.p0.6,TotDensEpi2Mod.p0.7) 
Anova(TotDensEpi2Mod.p0.7,type=3) 

TotDensEpi2Mod.p0.8<-update(TotDensEpi2Mod.p0.7,~.-Nutrients:Spores) # 2-way interaction removed
anova(TotDensEpi2Mod.p0.7,TotDensEpi2Mod.p0.8) 
Anova(TotDensEpi2Mod.p0.8,type=3) 

TotDensEpi2Mod.p0.9<-update(TotDensEpi2Mod.p0.8,~.-Nutrients:Mixing) # 2-way interaction removed
anova(TotDensEpi2Mod.p0.8,TotDensEpi2Mod.p0.9) 
Anova(TotDensEpi2Mod.p0.9,type=3) 

TotDensEpi2Mod.p0.10<-update(TotDensEpi2Mod.p0.9,~.-Nutrients:factor(Day)) # 2-way interaction removed
anova(TotDensEpi2Mod.p0.9,TotDensEpi2Mod.p0.10) 
Anova(TotDensEpi2Mod.p0.10,type=3) 

TotDensEpi2Mod.p0.11<-update(TotDensEpi2Mod.p0.10,~.-Spores:factor(Day)) # 2-way interaction removed
anova(TotDensEpi2Mod.p0.10,TotDensEpi2Mod.p0.11) 
Anova(TotDensEpi2Mod.p0.11,type=3) 
# Analysis of Deviance Table (Type III tests)
# 
# Response: LogTotDens
# Chisq Df Pr(>Chisq)    
# (Intercept)        1368.4673  1  < 2.2e-16 ***
# Raft                  9.8886  3   0.019537 *  
# Nutrients             4.9170  1   0.026594 *  
# Mixing                1.6008  1   0.205791    
# Spores                5.7291  1   0.016686 *  
# factor(Day)         144.6853  6  < 2.2e-16 ***
# Mixing:factor(Day)   17.9518  6   0.006354 ** 
summary(TotDensEpi2Mod.p0.11)
# Random effects:
#   Formula: ~1 | Bag
#         (Intercept)  Residual
# StdDev:   0.5372503 0.7764517

##################################################################################################
## Analyze time series of Log(Infected host density), starting from Day 1
##################################################################################################

## Full model with no autocorrelation structure
InfDensMod.p0.1<-lme(LogInfDens~Raft+Nutrients*Mixing*factor(Day),random=~1|Bag,method="ML",data=DensDataSpores) 
summary(InfDensMod.p0.1)
ACF(InfDensMod.p0.1)
plot(ACF(InfDensMod.p0.1),alpha=0.05) # suggests significant autocorrelation at time lag 1

## Full model with AR(1) autocorrelation structure
InfDensMod.p1.1<-lme(LogInfDens~Raft+Nutrients*Mixing*factor(Day),random=~1|Bag,
                     correlation=corARMA(form=~Day,p=1),method="ML",data=DensDataSpores) 
summary(InfDensMod.p1.1) # but estimate for phi = 0

## Stepwise simplification of model with no autocorrelation structure
InfDensMod.p0.2<-update(InfDensMod.p0.1,~.-Nutrients:Mixing:factor(Day)) # 3-way interaction removed
anova(InfDensMod.p0.1,InfDensMod.p0.2)
Anova(InfDensMod.p0.2,type=3)

InfDensMod.p0.3<-update(InfDensMod.p0.2,~.-Nutrients:Mixing) # least significant 2-way interaction removed
anova(InfDensMod.p0.2,InfDensMod.p0.3)
Anova(InfDensMod.p0.3,type=3)

InfDensMod.p0.4<-update(InfDensMod.p0.3,~.-Mixing:factor(Day)) # next least significant 2-way interaction removed
anova(InfDensMod.p0.3,InfDensMod.p0.4)
Anova(InfDensMod.p0.4,type=3)
  # Analysis of Deviance Table (Type III tests)
  # 
  # Response: LogInfDens
  #                     Chisq Df Pr(>Chisq)    
  # (Intercept)     262.7569  1    < 2e-16 ***
  # Raft              4.1150  3    0.24931    
  # Nutrients         0.0280  1    0.86703    
  # Mixing            4.4214  1    0.03549 *  
  # factor(Day)    1358.9125 14    < 2e-16 ***
  # Nutrients:factor(Day)   26.5605 14    0.02195 * 
summary(InfDensMod.p0.4)
  # Random effects:
  # Formula: ~1 | Bag
  #         (Intercept)  Residual
  # StdDev:   0.1733052 0.6357495

##################################################################################################
## Analyze time series of Log(Infected host density), starting from Day 12 (start of 1st epidemic)
##################################################################################################

## Full model with no autocorrelation structure
InfDensEpiMod.p0.1<-lme(LogInfDens~Raft+Nutrients*Mixing*factor(Day),random=~1|Bag,method="ML",data=DensDataStartDay12Spores) 
summary(InfDensEpiMod.p0.1)
ACF(InfDensEpiMod.p0.1)
plot(ACF(InfDensEpiMod.p0.1),alpha=0.05) # suggests significant autocorrelation at time lag 1

## Full model with AR(1) autocorrelation structure
InfDensEpiMod.p1.1<-lme(LogInfDens~Raft+Nutrients*Mixing*factor(Day),random=~1|Bag,
                     correlation=corARMA(form=~Day|Bag,p=1),method="ML",data=DensDataStartDay12Spores) 
summary(InfDensEpiMod.p1.1) # but estimate for phi = 0

## Stepwise simplification of model with no autocorrelation structure
InfDensEpiMod.p0.2<-update(InfDensEpiMod.p0.1,~.-Nutrients:Mixing:factor(Day)) # 3-way interaction removed
anova(InfDensEpiMod.p0.1,InfDensEpiMod.p0.2)
Anova(InfDensEpiMod.p0.2,type=3)

InfDensEpiMod.p0.3<-update(InfDensEpiMod.p0.2,~.-Nutrients:Mixing) # least significant 2-way interaction removed
anova(InfDensEpiMod.p0.2,InfDensEpiMod.p0.3)
Anova(InfDensEpiMod.p0.3,type=3)

InfDensEpiMod.p0.4<-update(InfDensEpiMod.p0.3,~.-Mixing:factor(Day)) # next least significant 2-way interaction removed
anova(InfDensEpiMod.p0.3,InfDensEpiMod.p0.4)
Anova(InfDensEpiMod.p0.4,type=3)

InfDensEpiMod.p0.5<-update(InfDensEpiMod.p0.4,~.-Nutrients:factor(Day)) # next least significant 2-way interaction removed
anova(InfDensEpiMod.p0.4,InfDensEpiMod.p0.5)
Anova(InfDensEpiMod.p0.5,type=3)
  # Analysis of Deviance Table (Type III tests)
  # 
  # Response: LogInfDens
  #                 Chisq Df Pr(>Chisq)    
  # (Intercept) 311.3904  1  < 2.2e-16 ***
  # Raft          4.1150  3  0.2493136    
  # Nutrients    12.6941  1  0.0003668 ***
  # Mixing        4.4214  1  0.0354913 *  
  # factor(Day) 625.3889 10  < 2.2e-16 ***
summary(InfDensEpiMod.p0.5)
  # Random effects:
  #   Formula: ~1 | Bag
  #         (Intercept)  Residual
  # StdDev:   0.2326118 0.7551865
plot(ACF(InfDensEpiMod.p0.5),alpha=0.05)

##################################################################################################
## Analyze time series of Log(Infected host density), starting from Day 26 (start of 2nd epidemic)
##################################################################################################

## Full model with no autocorrelation structure
InfDensEpi2Mod.p0.1<-lme(LogInfDens~Raft+Nutrients*Mixing*factor(Day),random=~1|Bag,method="ML",data=DensDataStartDay26Spores) 
summary(InfDensEpi2Mod.p0.1)
ACF(InfDensEpi2Mod.p0.1)
plot(ACF(InfDensEpi2Mod.p0.1),alpha=0.05) # suggests significant autocorrelation at time lag 1

## Full model with AR(1) autocorrelation structure
InfDensEpi2Mod.p1.1<-lme(LogInfDens~Raft+Nutrients*Mixing*factor(Day),random=~1|Bag,
                        correlation=corARMA(form=~Day|Bag,p=1),method="ML",data=DensDataStartDay26Spores) 
summary(InfDensEpi2Mod.p1.1) # but estimate for phi = 0

## Stepwise simplification of model with no autocorrelation structure
InfDensEpi2Mod.p0.2<-update(InfDensEpi2Mod.p0.1,~.-Nutrients:Mixing:factor(Day)) # 3-way interaction removed
anova(InfDensEpi2Mod.p0.1,InfDensEpi2Mod.p0.2)
Anova(InfDensEpi2Mod.p0.2,type=3)

InfDensEpi2Mod.p0.3<-update(InfDensEpi2Mod.p0.2,~.-Mixing:factor(Day)) # least significant 2-way interaction removed
anova(InfDensEpi2Mod.p0.2,InfDensEpi2Mod.p0.3)
Anova(InfDensEpi2Mod.p0.3,type=3)

InfDensEpi2Mod.p0.4<-update(InfDensEpi2Mod.p0.3,~.-Nutrients:Mixing) # next least significant 2-way interaction removed
anova(InfDensEpi2Mod.p0.3,InfDensEpi2Mod.p0.4)
Anova(InfDensEpi2Mod.p0.4,type=3)

InfDensEpi2Mod.p0.5<-update(InfDensEpi2Mod.p0.4,~.-Nutrients:factor(Day)) # next least significant 2-way interaction removed
anova(InfDensEpi2Mod.p0.4,InfDensEpi2Mod.p0.5)
Anova(InfDensEpi2Mod.p0.5,type=3)
  # Analysis of Deviance Table (Type III tests)
  # 
  # Response: LogInfDens
  #                 Chisq Df Pr(>Chisq)    
  # (Intercept) 658.6463  1    < 2e-16 ***
  # Raft          7.2337  3    0.06481 .  
  # Nutrients     9.3496  1    0.00223 ** 
  # Mixing        5.2600  1    0.02182 *  
  # factor(Day) 221.2576  6    < 2e-16 ***
summary(InfDensEpi2Mod.p0.5)
  # Random effects:
  #   Formula: ~1 | Bag
  #         (Intercept)  Residual
  # StdDev:   0.3012317 0.7600628

##################################################################################################
## Prepare data for analyses of infection prevalence
##################################################################################################

PrevTotData<-Spores[,c("Raft","Nutrients","Mixing","Spores","Trt","Bag","Day","TotalCount","PrevTotal")]
sum(complete.cases(PrevTotData$PrevTotal))

PrevAdultData<-Spores[,c("Raft","Nutrients","Mixing","Spores","Trt","Bag","Day","AdultCount","PrevAdult")]
sum(complete.cases(PrevAdultData$PrevAdult))

## Omit missing data
PrevTotData<-na.omit(PrevTotData)
PrevAdultData<-na.omit(PrevAdultData)

## Subset data to begin at start of epidemics (Day 12 onward)
PrevTotDataEpi<-subset(PrevTotData,Day>=12) 
hist(PrevTotDataEpi$PrevTotal)

PrevAdultDataEpi<-subset(PrevAdultData,Day>=12) 
hist(PrevAdultDataEpi$PrevAdult)

## Subset data to begin at start of 2nd epidemics (Day 26 onward)
PrevTotDataEpi2<-subset(PrevTotData,Day>=26) 
hist(PrevTotDataEpi2$PrevTotal)

PrevAdultDataEpi2<-subset(PrevAdultData,Day>=26) 
hist(PrevAdultDataEpi2$PrevAdult)

PrevTotDataDay35<-subset(PrevTotData,Day==35) 
hist(PrevTotDataDay35$PrevTotal)

ggplot(PrevTotDataDay35)+
  geom_boxplot(aes(x=Mixing,y=PrevTotal,color=Nutrients))+
  scale_colour_manual(values=NutColors)+ 
  ylab("Total prevalence at Day 35")+
  fav_theme

m35.1<-glm(PrevTotal~Raft+Nutrients*Mixing,
      weights=TotalCount,
      data=PrevTotDataDay35,family=quasibinomial)
Anova(m35.1,type=3)
summary(m35.1)

m35.2<-glm(PrevTotal~Raft+Nutrients+Mixing,
           weights=TotalCount,
           data=PrevTotDataDay35,family=quasibinomial)
Anova(m35.2,type=3)
summary(m35.2)
##################################################################################################
## Plot and analyze time series of total infection prevalence, starting from Day 1
## (Ben Bolker helped with this analysis in Aug 2012)
##################################################################################################

(p7 <- ggplot(PrevTotData,aes(x=Day,y=PrevTotal,colour=Nutrients,shape=Mixing,
                              linetype=Mixing))+
   geom_point(alpha=0.7,aes(size=TotalCount))+ 
   scale_colour_manual(values=NutColors)+ 
   scale_y_continuous(limits=c(0,0.6),oob=squish)+
   ylab("Total prevalence")+
   xlab("Day of experiment")+ 
   fav_theme)

## Smooth lines by Mixing and Nutrients treatments
p7 + geom_smooth(alpha=0.2,method="loess",span=0.4,aes(weight=TotalCount))

## Smooth lines by bag
p7 + geom_smooth(alpha=0.2,method="loess",span=0.4,
                 aes(group=Bag,weight=TotalCount),se=FALSE)

## Lines (unsmoothed) by bag
p7 + geom_line(aes(group=Bag))

## Add splines 
p7 + geom_smooth(alpha=0.2,method="glm",family="binomial",
                 aes(weight=TotalCount),
                 formula=y~ns(x,3)) 
## Visually: spline with 3 df does a bad job fitting shape of epidemic humps

p7 + geom_smooth(alpha=0.2,method="glm",family="binomial",
                 aes(weight=TotalCount),
                 formula=y~ns(x,4)) 
## Visually: spline with 4 df does a good job fitting shape of epidemic humps

## Analyze using binomial glmm with 4 df spline
PrevTotMod.ns4.1 <- glmer(PrevTotal~Raft+Nutrients*Mixing*ns(Day,4)+(1|Bag),
                          weights=TotalCount,
                          data=PrevTotData,family=binomial)
# Warning message:
#   In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#     Model failed to converge with max|grad| = 0.00395817 (tol = 0.001, component 1)

relgrad.PrevTotMod.ns4.1 <- with(PrevTotMod.ns4.1@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad.PrevTotMod.ns4.1)) # 0.0028 (acceptable if < 0.001)

## Get the following error message if try to run glmer with 5 df spline
# Error in (function (fr, X, reTrms, family, nAGQ = 1L, verbose = 0L, maxit = 100L,  : 
#   (maxstephalfit) PIRLS step-halvings failed to reduce deviance in pwrssUpdate

summary(PrevTotMod.ns4.1)
Anova(PrevTotMod.ns4.1,type=3) # looks like 3-way interaction is significant...but can we trust model results if model failed to converge?

## Checking if model converges if we remove higher-order interactions
PrevTotMod.ns4.2<-update(PrevTotMod.ns4.1,~.-Nutrients:Mixing:ns(Day,4))
relgrad.PrevTotMod.ns4.2 <- with(PrevTotMod.ns4.2@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad.PrevTotMod.ns4.2)) # 0.0016 (acceptable if < 0.001)

anova(PrevTotMod.ns4.1,PrevTotMod.ns4.2) # significant 3-way interaction
Anova(PrevTotMod.ns4.2,type=3)

PrevTotMod.ns4.3<-update(PrevTotMod.ns4.2,~.-Nutrients:Mixing)
relgrad.PrevTotMod.ns4.3 <- with(PrevTotMod.ns4.3@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad.PrevTotMod.ns4.3)) # 0.0017
Anova(PrevTotMod.ns4.3,type=3)

pframe <- with(PrevTotData,expand.grid(Raft=levels(Raft),Nutrients=levels(Nutrients),
                                       Mixing=levels(Mixing),
                                       Bag=1,
                                       Day=seq(1,max(Day),length=max(Day))))

pframe$PrevTotal <- plogis(predict(PrevTotMod.ns4.1,newdata=pframe,re.form=NA))

## Model has Raft as a fixed effect. Here we plot model predictions for each bag separately.
p7 + geom_line(data=pframe,aes(group=interaction(Nutrients,Mixing,Raft)))

## Here we plot model predictions averaged by Raft
pframe2 <- ddply(pframe,c("Mixing","Nutrients","Day"),
                 function(x) {
                   x <- transform(x,PrevTotal_avg=mean(PrevTotal))
                 })

p7 + geom_line(data=pframe2,aes(y=PrevTotal_avg))

##################################################################################################
## Plot and analyze time series of total infection prevalence, starting from Day 12 (1st epidemic)
##################################################################################################

(p8 <- ggplot(PrevTotDataEpi,aes(x=Day,y=PrevTotal,colour=Nutrients,shape=Mixing,
                                 linetype=Mixing))+
   geom_point(alpha=0.7,aes(size=TotalCount))+ 
   scale_colour_manual(values=NutColors)+ 
   scale_y_continuous(limits=c(0,0.6),oob=squish)+
   ylab("Total prevalence")+
   xlab("Day of experiment")+ 
   fav_theme)

## Add splines with 3 df
p8 + geom_smooth(alpha=0.2,method="glm",family="binomial",
                 aes(weight=TotalCount),
                 formula=y~ns(x,3)) 

## Add splines with 4 df
p8 + geom_smooth(alpha=0.2,method="glm",family="binomial",
                 aes(weight=TotalCount),
                 formula=y~ns(x,4)) 

## Analyze using binomial glmm with 3 df spline
PrevTotModEpi.ns3.1 <- glmer(PrevTotal~Raft+Nutrients*Mixing*ns(Day,3)+(1|Bag),
                             weights=TotalCount,
                             data=PrevTotDataEpi,family=binomial)
Anova(PrevTotModEpi.ns3.1,Type=3) # significant 3-way interaction

##################################################################################################
## Plot and analyze time series of total infection prevalence, starting from Day 26 (2nd epidemic)
##################################################################################################

(p9 <- ggplot(PrevTotDataEpi2,aes(x=Day,y=PrevTotal,colour=Nutrients,shape=Mixing,
                                  linetype=Mixing))+
   geom_point(alpha=0.7,aes(size=TotalCount))+ 
   scale_colour_manual(values=NutColors)+ 
   scale_y_continuous(limits=c(0,0.6),oob=squish)+
   ylab("Total prevalence")+
   xlab("Day of experiment")+ 
   fav_theme)

## Add splines with 2 df
p9 + geom_smooth(alpha=0.2,method="glm",family="binomial",
                 aes(weight=TotalCount),
                 formula=y~ns(x,2)) 

## Analyze using binomial glmm with 2 df spline
PrevTotModEpi2.ns2.1 <- glmer(PrevTotal~Raft+Nutrients*Mixing*ns(Day,2)+(1|Bag),
                              weights=TotalCount,
                              data=PrevTotDataEpi2,family=binomial)
Anova(PrevTotModEpi2.ns2.1,Type=3) # significant 3-way interaction
summary(PrevTotModEpi2.ns2.1)

## Add splines with 3 df
p9 + geom_smooth(alpha=0.2,method="glm",family="binomial",
                 aes(weight=TotalCount),
                 formula=y~ns(x,3)) 

## Analyze using binomial glmm with 3 df spline
PrevTotModEpi2.ns3.1 <- glmer(PrevTotal~Raft+Nutrients*Mixing*ns(Day,3)+(1|Bag),
                              weights=TotalCount,
                              data=PrevTotDataEpi2,family=binomial)
Anova(PrevTotModEpi2.ns3.1,Type=3) # significant 3-way interaction

##################################################################################################
## Analyses of TotChl, EdChl, and InedChl at start days 29, 32, 35, 40
##################################################################################################

## Subset data 
ChlData<-mydata[,c("Raft","Nutrients","Mixing","Spores","Trt","Bag","Day",
                   "TotChl","EdChl","InedChl","PropInedChl","LogTotChl","LogEdChl","LogInedChl")]

## Omit missing data
ChlData<-na.omit(ChlData)

## Subset data for -spores (uninf) and +spores (inf) bags
ChlDataUninf<-subset(ChlData,Spores=="No")
ChlDataInf<-subset(ChlData,Spores=="Spores")

##################################################################################################
StartDay29Chl<-subset(ChlData,Day>=29) # all bags
StartDay29ChlSpores<-subset(ChlDataInf,Day>=29) # +spores only

SporesNames<-c("No"="-Spores","Spores"="+Spores")
ggplot(StartDay29Chl)+
  geom_boxplot(aes(x=factor(Day),y=LogTotChl,color=Nutrients))+
  facet_wrap(~Spores,labeller = labeller(Spores = SporesNames))+
  scale_color_manual(values=NutColors)+
  ylab("Log (Total chlorophyll a, ug/L)")+
  fav_theme

ggplot(StartDay29Chl)+
  geom_boxplot(aes(x=Spores,y=LogTotChl,color=Nutrients),notch=TRUE)+
  geom_point(aes(x=Spores,y=LogTotChl,color=Nutrients),position=position_jitterdodge(),alpha=0.8,size=0.8)+
  scale_color_manual(values=NutColors)+
  ylab("Log (Total chlorophyll a, ug/L)")+
  labs(title="Day 29-46")+
  fav_theme

mchl29.1<-lme(LogTotChl~Nutrients*Mixing*Spores+Raft,random=~1|Bag,method="REML",data=StartDay29Chl) 
summary(mchl29.1)
Anova(mchl29.1,type=3)

mchl29.2<-update(mchl29.1,~.-Nutrients:Mixing:Spores)
Anova(mchl29.2,type=3)

mchl29.3<-update(mchl29.2,~.-Nutrients:Spores)
Anova(mchl29.3,type=3)

mchl29.4<-update(mchl29.3,~.-Mixing:Spores)
Anova(mchl29.4,type=3)

mchl29.5<-update(mchl29.4,~.-Nutrients:Mixing)
Anova(mchl29.5,type=3)

##################################################################################################
StartDay32Chl<-subset(ChlData,Day>=32) # all bags
StartDay32ChlSpores<-subset(ChlDataInf,Day>=32) # +spores only

SporesNames<-c("No"="-Spores","Spores"="+Spores")
ggplot(StartDay32Chl)+
  geom_boxplot(aes(x=factor(Day),y=LogTotChl,color=Nutrients))+
  facet_wrap(~Spores,labeller = labeller(Spores = SporesNames))+
  scale_color_manual(values=NutColors)+
  ylab("Log (Total chlorophyll a, ug/L)")+
  fav_theme

ggplot(StartDay32Chl)+
  geom_boxplot(aes(x=Spores,y=LogTotChl,color=Nutrients),notch=TRUE)+
  geom_point(aes(x=Spores,y=LogTotChl,color=Nutrients),position=position_jitterdodge(),alpha=0.8,size=0.8)+
  scale_color_manual(values=NutColors)+
  ylab("Log (Total chlorophyll a, ug/L)")+
  labs(title="Day 32-46")+
  fav_theme

ggplot(StartDay32Chl)+
  geom_boxplot(aes(x=Spores,y=LogEdChl,color=Nutrients),notch=TRUE)+
  scale_color_manual(values=NutColors)+
  ylab("Log (Edible chlorophyll a, ug/L)")+
  labs(title="Day 32-46")+
  fav_theme

ggplot(StartDay32Chl)+
  geom_boxplot(aes(x=Spores,y=LogInedChl,color=Nutrients),notch=TRUE)+
  scale_color_manual(values=NutColors)+
  ylab("Log (Inedible chlorophyll a, ug/L)")+
  labs(title="Day 32-46")+
  fav_theme

mchl32.1<-lme(LogTotChl~Nutrients*Mixing*Spores+Raft,random=~1|Bag,method="REML",data=StartDay32Chl) 
summary(mchl32.1)
Anova(mchl32.1,type=3)

mchl32.2<-update(mchl32.1,~.-Nutrients:Mixing:Spores)
Anova(mchl32.2,type=3)

mchl32.3<-update(mchl32.2,~.-Nutrients:Spores)
Anova(mchl32.3,type=3)

mchl32.4<-update(mchl32.3,~.-Mixing:Spores)
Anova(mchl32.4,type=3)

mchl32.5<-update(mchl32.4,~.-Nutrients:Mixing)
Anova(mchl32.5,type=3)

##################################################################################################
StartDay35Chl<-subset(ChlData,Day>=35) # all bags
StartDay35ChlSpores<-subset(ChlDataInf,Day>=35) # +spores only

SporesNames<-c("No"="-Spores","Spores"="+Spores")
ggplot(StartDay35Chl)+
  geom_boxplot(aes(x=factor(Day),y=LogTotChl,color=Nutrients))+
  facet_wrap(~Spores,labeller = labeller(Spores = SporesNames))+
  scale_color_manual(values=NutColors)+
  ylab("Log (Total chlorophyll a, ug/L)")+
  fav_theme

ggplot(StartDay35Chl)+
  geom_boxplot(aes(x=Spores,y=LogTotChl,color=Nutrients),notch=TRUE)+
  geom_point(aes(x=Spores,y=LogTotChl,color=Nutrients),position=position_jitterdodge(),alpha=0.8,size=0.8)+
  scale_color_manual(values=NutColors)+
  ylab("Log (Total chlorophyll a, ug/L)")+
  labs(title="Day 35-46")+
  fav_theme

ggplot(StartDay35Chl)+
  geom_boxplot(aes(x=Nutrients,y=LogTotChl,color=Nutrients),notch=TRUE)+
  geom_point(aes(x=Nutrients,y=LogTotChl,color=Nutrients),position=position_jitter(),alpha=0.8,size=0.8)+
  scale_color_manual(values=NutColors)+
  ylab("Log (Total chlorophyll a, ug/L)")+
  labs(title="Day 35-46")+
  fav_theme

ggplot(StartDay35Chl)+
  geom_boxplot(aes(x=Spores,y=LogEdChl,color=Nutrients),notch=TRUE)+
  scale_color_manual(values=NutColors)+
  ylab("Log (Edible chlorophyll a, ug/L)")+
  labs(title="Day 35-46")+
  fav_theme

mchl35.1<-lme(LogTotChl~Nutrients*Mixing*Spores+Raft,random=~1|Bag,method="REML",data=StartDay35Chl) 
summary(mchl35.1)
Anova(mchl35.1,type=3)

mchl35.2<-update(mchl35.1,~.-Nutrients:Mixing:Spores)
Anova(mchl35.2,type=3)

mchl35.3<-update(mchl35.2,~.-Nutrients:Spores)
Anova(mchl35.3,type=3)

mchl35.4<-update(mchl35.3,~.-Mixing:Spores)
Anova(mchl35.4,type=3)

mchl35.5<-update(mchl35.4,~.-Nutrients:Mixing)
Anova(mchl35.5,type=3)

##################################################################################################
StartDay40Chl<-subset(ChlData,Day>=40) # all bags
StartDay40ChlSpores<-subset(ChlDataInf,Day>=40) # +spores only

SporesNames<-c("No"="-Spores","Spores"="+Spores")
ggplot(StartDay40Chl)+
  geom_boxplot(aes(x=factor(Day),y=LogTotChl,color=Nutrients))+
  facet_wrap(~Spores,labeller = labeller(Spores = SporesNames))+
  scale_color_manual(values=NutColors)+
  ylab("Log (Total chlorophyll a, ug/L)")+
  fav_theme

ggplot(StartDay40Chl)+
  geom_boxplot(aes(x=Spores,y=LogTotChl,color=Nutrients),notch=TRUE)+
  geom_point(aes(x=Spores,y=LogTotChl,color=Nutrients),position=position_jitterdodge(),alpha=0.8,size=0.8)+
  scale_color_manual(values=NutColors)+
  ylab("Log (Total chlorophyll a, ug/L)")+
  labs(title="Day 40-46")+
  fav_theme

ggplot(StartDay40Chl)+
  geom_boxplot(aes(x=Spores,y=LogEdChl,color=Nutrients),notch=TRUE)+
  scale_color_manual(values=NutColors)+
  ylab("Log (Edible chlorophyll a, ug/L)")+
  labs(title="Day 40-46")+
  fav_theme

mchl40.1<-lme(LogTotChl~Nutrients*Mixing*Spores+Raft,random=~1|Bag,method="REML",data=StartDay40Chl) 
summary(mchl40.1)
Anova(mchl40.1,type=3)

mchl40.2<-update(mchl40.1,~.-Nutrients:Mixing:Spores)
Anova(mchl40.2,type=3)

mchl40.3<-update(mchl40.2,~.-Nutrients:Spores)
Anova(mchl40.3,type=3)

mchl40.4<-update(mchl40.3,~.-Mixing:Spores)
Anova(mchl40.4,type=3)

mchl40.5<-update(mchl40.4,~.-Nutrients:Mixing)
Anova(mchl40.5,type=3)
summary(mchl40.5)

##################################################################################################
## Plot time series of Log(TotChl), starting from Day 1
##################################################################################################

## First, plot each bag as its own line

p18<-(ggplot(ChlDataUninf,aes(x=Day,y=LogTotChl,group=Bag,color=Nutrients))+
        scale_color_manual(values=NutColors)+
        geom_line(aes(linetype=Mixing))+
        ylim(c(0,7))+
        ylab("Log (Total chlorophyll a, ug/L)")+
        xlab("Day of experiment")+ 
        labs(title="- Spores")+
        guides(color="none",linetype="none")+ 
        fav_theme)

p19<-(ggplot(ChlDataInf,aes(x=Day,y=LogTotChl,group=Bag,color=Nutrients))+
        scale_color_manual(values=NutColors)+
        geom_line(aes(linetype=Mixing))+
        ylim(c(0,7))+
        ylab("Log (Total chlorophyll a, ug/L)")+
        xlab("Day of experiment")+ 
        labs(title="+ Spores")+
        guides(color="none",linetype="none")+  
        fav_theme)

## Then, plot treatment means +/- SE

UninfLogTotChlSummary<-summarySE(ChlDataUninf, measurevar="LogTotChl", groupvars=c("Nutrients","Mixing","Day"))

p20<-(ggplot(UninfLogTotChlSummary)+
        geom_point(aes(x=Day, y=LogTotChl, color=Nutrients)) +  
        geom_line(aes(x=Day, y=LogTotChl, color=Nutrients, linetype=Mixing))+
        geom_errorbar(aes(x=Day, ymin=LogTotChl-se, ymax=LogTotChl+se,color=Nutrients), width=0) + 
        scale_color_manual(values=NutColors)+
        fav_theme+
        ylim(c(0,7))+
        ylab("Log (Total chlorophyll a, ug/L)")+
        xlab("Day of experiment")+ 
        labs(title="- Spores"))

InfLogTotChlSummary<-summarySE(ChlDataInf, measurevar="LogTotChl", groupvars=c("Nutrients","Mixing","Day"))

p21<-(ggplot(InfLogTotChlSummary)+
        geom_point(aes(x=Day, y=LogTotChl, color=Nutrients)) +  
        geom_line(aes(x=Day, y=LogTotChl, color=Nutrients, linetype=Mixing))+
        geom_errorbar(aes(x=Day, ymin=LogTotChl-se, ymax=LogTotChl+se,color=Nutrients), width=0) + 
        scale_color_manual(values=NutColors)+
        fav_theme+
        ylim(c(0,7))+
        ylab("Log (Total chlorophyll a, ug/L)")+
        xlab("Day of experiment")+ 
        labs(title="+ Spores")) 

## Put panels together into one figure
grid.arrange(p18,p20,p19,p21, nrow=2,ncol=2)

##################################################################################################
## Plot time series of Log(EdChl), starting from Day 1
##################################################################################################

## First, plot each bag as its own line

p22<-(ggplot(ChlDataUninf,aes(x=Day,y=LogEdChl,group=Bag,color=Nutrients))+
        scale_color_manual(values=NutColors)+
        geom_line(aes(linetype=Mixing))+
        ylim(c(-0.5,5))+
        ylab("Log (Edible chlorophyll a, ug/L)")+
        xlab("Day of experiment")+ 
        labs(title="- Spores")+
        guides(color="none",linetype="none")+ 
        fav_theme)

p23<-(ggplot(ChlDataInf,aes(x=Day,y=LogEdChl,group=Bag,color=Nutrients))+
        scale_color_manual(values=NutColors)+
        geom_line(aes(linetype=Mixing))+
        ylim(c(-0.5,5))+
        ylab("Log (Edible chlorophyll a, ug/L)")+
        xlab("Day of experiment")+ 
        labs(title="+ Spores")+
        guides(color="none",linetype="none")+  
        fav_theme)

## Then, plot treatment means +/- SE

UninfLogEdChlSummary<-summarySE(ChlDataUninf, measurevar="LogEdChl", groupvars=c("Nutrients","Mixing","Day"))

p24<-(ggplot(UninfLogEdChlSummary)+
        geom_point(aes(x=Day, y=LogEdChl, color=Nutrients)) +  
        geom_line(aes(x=Day, y=LogEdChl, color=Nutrients, linetype=Mixing))+
        geom_errorbar(aes(x=Day, ymin=LogEdChl-se, ymax=LogEdChl+se,color=Nutrients), width=0) + 
        scale_color_manual(values=NutColors)+
        fav_theme+
        ylim(c(-0.5,5))+
        ylab("Log (Edible chlorophyll a, ug/L)")+
        xlab("Day of experiment")+ 
        labs(title="- Spores"))

InfLogEdChlSummary<-summarySE(ChlDataInf, measurevar="LogEdChl", groupvars=c("Nutrients","Mixing","Day"))

p25<-(ggplot(InfLogEdChlSummary)+
        geom_point(aes(x=Day, y=LogEdChl, color=Nutrients)) +  
        geom_line(aes(x=Day, y=LogEdChl, color=Nutrients, linetype=Mixing))+
        geom_errorbar(aes(x=Day, ymin=LogEdChl-se, ymax=LogEdChl+se,color=Nutrients), width=0) + 
        scale_color_manual(values=NutColors)+
        fav_theme+
        ylim(c(-0.5,5))+
        ylab("Log (Edible chlorophyll a, ug/L)")+
        xlab("Day of experiment")+ 
        labs(title="+ Spores")) 

## Put panels together into one figure
grid.arrange(p22,p24,p23,p25, nrow=2,ncol=2)

##################################################################################################
## Plot time series of Log(InedChl), starting from Day 1
##################################################################################################

## First, plot each bag as its own line

p26<-(ggplot(ChlDataUninf,aes(x=Day,y=LogInedChl,group=Bag,color=Nutrients))+
        scale_color_manual(values=NutColors)+
        geom_line(aes(linetype=Mixing))+
        ylim(c(-4,7))+
        ylab("Log (Inedible chlorophyll a, ug/L)")+
        xlab("Day of experiment")+ 
        labs(title="- Spores")+
        guides(color="none",linetype="none")+ 
        fav_theme)

p27<-(ggplot(ChlDataInf,aes(x=Day,y=LogInedChl,group=Bag,color=Nutrients))+
        scale_color_manual(values=NutColors)+
        geom_line(aes(linetype=Mixing))+
        ylim(c(-4,7))+
        ylab("Log (Inedible chlorophyll a, ug/L)")+
        xlab("Day of experiment")+ 
        labs(title="+ Spores")+
        guides(color="none",linetype="none")+  
        fav_theme)

## Then, plot treatment means +/- SE

UninfLogInedChlSummary<-summarySE(ChlDataUninf, measurevar="LogInedChl", groupvars=c("Nutrients","Mixing","Day"))

p28<-(ggplot(UninfLogInedChlSummary)+
        geom_point(aes(x=Day, y=LogInedChl, color=Nutrients)) +  
        geom_line(aes(x=Day, y=LogInedChl, color=Nutrients, linetype=Mixing))+
        geom_errorbar(aes(x=Day, ymin=LogInedChl-se, ymax=LogInedChl+se,color=Nutrients), width=0) + 
        scale_color_manual(values=NutColors)+
        fav_theme+
        ylim(c(-4,7))+
        ylab("Log (Inedible chlorophyll a, ug/L)")+
        xlab("Day of experiment")+ 
        labs(title="- Spores"))

InfLogInedChlSummary<-summarySE(ChlDataInf, measurevar="LogInedChl", groupvars=c("Nutrients","Mixing","Day"))

p29<-(ggplot(InfLogInedChlSummary)+
        geom_point(aes(x=Day, y=LogInedChl, color=Nutrients)) +  
        geom_line(aes(x=Day, y=LogInedChl, color=Nutrients, linetype=Mixing))+
        geom_errorbar(aes(x=Day, ymin=LogInedChl-se, ymax=LogInedChl+se,color=Nutrients), width=0) + 
        scale_color_manual(values=NutColors)+
        fav_theme+
        ylim(c(-4,7))+
        ylab("Log (Inedible chlorophyll a, ug/L)")+
        xlab("Day of experiment")+ 
        labs(title="+ Spores")) 

## Put panels together into one figure
grid.arrange(p26,p28,p27,p29, nrow=2,ncol=2)

##################################################################################################
## Plot time series of PropInedChl (proportion of chl a in inedible size fraction), starting from Day 1
##################################################################################################

## First, plot each bag as its own line

p30<-(ggplot(ChlDataUninf,aes(x=Day,y=PropInedChl,group=Bag,color=Nutrients))+
        scale_color_manual(values=NutColors)+
        geom_line(aes(linetype=Mixing))+
        ylim(c(0,1))+
        ylab("Proportion inedible chl a")+
        xlab("Day of experiment")+ 
        labs(title="- Spores")+
        guides(color="none",linetype="none")+ 
        fav_theme)

p31<-(ggplot(ChlDataInf,aes(x=Day,y=PropInedChl,group=Bag,color=Nutrients))+
        scale_color_manual(values=NutColors)+
        geom_line(aes(linetype=Mixing))+
        ylim(c(0,1))+
        ylab("Proportion inedible chl a")+
        xlab("Day of experiment")+ 
        labs(title="+ Spores")+
        guides(color="none",linetype="none")+  
        fav_theme)

## Then, plot treatment means +/- SE

UninfPropInedChlSummary<-summarySE(ChlDataUninf, measurevar="PropInedChl", groupvars=c("Nutrients","Mixing","Day"))

p32<-(ggplot(UninfPropInedChlSummary)+
        geom_point(aes(x=Day, y=PropInedChl, color=Nutrients)) +  
        geom_line(aes(x=Day, y=PropInedChl, color=Nutrients, linetype=Mixing))+
        geom_errorbar(aes(x=Day, ymin=PropInedChl-se, ymax=PropInedChl+se,color=Nutrients), width=0) + 
        scale_color_manual(values=NutColors)+
        fav_theme+
        ylim(c(0,1))+
        ylab("Proportion inedible chl a")+
        xlab("Day of experiment")+ 
        labs(title="- Spores"))

InfPropInedChlSummary<-summarySE(ChlDataInf, measurevar="PropInedChl", groupvars=c("Nutrients","Mixing","Day"))

p33<-(ggplot(InfPropInedChlSummary)+
        geom_point(aes(x=Day, y=PropInedChl, color=Nutrients)) +  
        geom_line(aes(x=Day, y=PropInedChl, color=Nutrients, linetype=Mixing))+
        geom_errorbar(aes(x=Day, ymin=PropInedChl-se, ymax=PropInedChl+se,color=Nutrients), width=0) + 
        scale_color_manual(values=NutColors)+
        fav_theme+
        ylim(c(0,1))+
        ylab("Proportion inedible chl a")+
        xlab("Day of experiment")+ 
        labs(title="+ Spores")) 

## Put panels together into one figure
grid.arrange(p30,p32,p31,p33, nrow=2,ncol=2)

(p34 <- ggplot(ChlDataUninf,aes(x=Day,y=PropInedChl,colour=Nutrients,shape=Mixing,
                                linetype=Mixing))+
    geom_point(alpha=0.7,aes(size=TotChl))+ 
    scale_colour_manual(values=NutColors)+ 
    scale_y_continuous(limits=c(0,1),oob=squish)+
    ylab("Proportion inedible chl a")+
    xlab("Day of experiment")+ 
    fav_theme)

## Add splines with 3 df
p34 + geom_smooth(alpha=0.2,method="glm",family="binomial",
                  aes(weight=TotChl),
                  formula=y~ns(x,3)) 

(p35 <- ggplot(ChlDataInf,aes(x=Day,y=PropInedChl,colour=Nutrients,shape=Mixing,
                              linetype=Mixing))+
    geom_point(alpha=0.7,aes(size=TotChl))+ 
    scale_colour_manual(values=NutColors)+ 
    scale_y_continuous(limits=c(0,1),oob=squish)+
    ylab("Proportion inedible chl a")+
    xlab("Day of experiment")+ 
    fav_theme)

## Add splines with 3 df
p35 + geom_smooth(alpha=0.2,method="glm",family="binomial",
                  aes(weight=TotChl),
                  formula=y~ns(x,3)) 

##################################################################################################
## Analyze time series of Log(Tot Chl), starting from Day 1
##################################################################################################

## Full model with no autocorrelation structure
TotChlMod.p0.1<-lme(LogTotChl~Raft+Nutrients*Mixing*Spores*factor(Day),random=~1|Bag,method="ML",data=ChlData) 
summary(TotChlMod.p0.1)
ACF(TotChlMod.p0.1)
plot(ACF(TotChlMod.p0.1),alpha=0.05) # suggests significant autocorrelation at time lags 1-3?

## Full model with AR(3) autocorrelation structure
TotChlMod.p3.1<-lme(LogTotChl~Raft+Nutrients*Mixing*Spores*factor(Day),random=~1|Bag,
                     correlation=corARMA(form=~Day,p=3),method="ML",data=ChlData) 
summary(TotChlMod.p3.1) 

## Stepwise simplification of model
TotChlMod.p3.2<-update(TotChlMod.p3.1,~.-Nutrients:Mixing:Spores:factor(Day)) # 4-way interaction removed
anova(TotChlMod.p3.2,TotChlMod.p3.1)
Anova(TotChlMod.p3.2,type=3)

TotChlMod.p3.3<-update(TotChlMod.p3.2,~.-Mixing:Spores:factor(Day))
anova(TotChlMod.p3.3,TotChlMod.p3.2)
Anova(TotChlMod.p3.3,type=3)

TotChlMod.p3.4<-update(TotChlMod.p3.3,~.-Nutrients:Mixing:Spores)
anova(TotChlMod.p3.4,TotChlMod.p3.3)
Anova(TotChlMod.p3.4,type=3)

TotChlMod.p3.5<-update(TotChlMod.p3.4,~.-Nutrients:Spores:factor(Day))
anova(TotChlMod.p3.5,TotChlMod.p3.4)
Anova(TotChlMod.p3.5,type=3)

TotChlMod.p3.6<-update(TotChlMod.p3.5,~.-Nutrients:Mixing:factor(Day))
anova(TotChlMod.p3.6,TotChlMod.p3.5)
Anova(TotChlMod.p3.6,type=3)

TotChlMod.p3.7<-update(TotChlMod.p3.6,~.-Nutrients:Spores)
anova(TotChlMod.p3.7,TotChlMod.p3.6)
Anova(TotChlMod.p3.7,type=3)

TotChlMod.p3.8<-update(TotChlMod.p3.7,~.-Mixing:Spores)
anova(TotChlMod.p3.8,TotChlMod.p3.7)
Anova(TotChlMod.p3.8,type=3)

TotChlMod.p3.9<-update(TotChlMod.p3.8,~.-Spores:factor(Day))
anova(TotChlMod.p3.9,TotChlMod.p3.8)
Anova(TotChlMod.p3.9,type=3)

TotChlMod.p3.10<-update(TotChlMod.p3.9,~.-Nutrients:Mixing)
anova(TotChlMod.p3.10,TotChlMod.p3.9)
Anova(TotChlMod.p3.10,type=3)

TotChlMod.p3.11<-update(TotChlMod.p3.10,~.-Mixing:factor(Day))
anova(TotChlMod.p3.11,TotChlMod.p3.10)
Anova(TotChlMod.p3.11,type=3)

##################################################################################################
## Regression of TP on Day 15 vs. each bag's mean infected host density from Days 15-46 
##################################################################################################

BagMeanEpiDay15LogInfDens<-ddply(DensDataStartDay15Spores,c("Bag"),summarize,
                             BagMeanEpiDay15LogInfDens=mean(LogInfDens)
)

BagTPDayData<-mydata[,c("Nutrients","Mixing","Spores","Bag","Day","TP")]
BagTPDay15Data<-subset(BagTPDayData,Day==15)
BagTPDay15DataInf<-subset(BagTPDay15Data,Spores=="Spores")

str(BagTPDay15DataInf)
str(BagMeanEpiDay15LogInfDens)
TPDay15BagMeanEpiDay15LogInfDens<-merge(BagMeanEpiDay15LogInfDens,BagTPDay15DataInf,by="Bag")
TPDay15BagMeanEpiDay15LogInfDens

ggplot(TPDay15BagMeanEpiDay15LogInfDens)+
  geom_point(aes(x=log(TP),y=BagMeanEpiDay15LogInfDens,color=Nutrients,shape=Mixing))+
  stat_smooth(aes(x=log(TP), y=BagMeanEpiDay15LogInfDens),method='lm',color='black')+
  scale_color_manual(values=NutColors)+
  fav_theme

summary(lm(BagMeanEpiDay15LogInfDens~log(TP),data=TPDay15BagMeanEpiDay15LogInfDens))

cor.test(log(TPDay15BagMeanEpiDay15LogInfDens$TP),TPDay15BagMeanEpiDay15LogInfDens$BagMeanEpiDay15LogInfDens)

##################################################################################################
## Regression of TN on Day 15 vs. each bag's mean infected host density from Days 15-46 
##################################################################################################

BagMeanEpiDay15LogInfDens<-ddply(DensDataStartDay15Spores,c("Bag"),summarize,
                                 BagMeanEpiDay15LogInfDens=mean(LogInfDens)
)

BagTNDayData<-mydata[,c("Nutrients","Mixing","Spores","Bag","Day","TN")]
BagTNDay15Data<-subset(BagTNDayData,Day==15)
BagTNDay15DataInf<-subset(BagTNDay15Data,Spores=="Spores")

str(BagTNDay15DataInf)
str(BagMeanEpiDay15LogInfDens)
TNDay15BagMeanEpiDay15LogInfDens<-merge(BagMeanEpiDay15LogInfDens,BagTNDay15DataInf,by="Bag")
TNDay15BagMeanEpiDay15LogInfDens

ggplot(TNDay15BagMeanEpiDay15LogInfDens)+
  geom_point(aes(x=log(TN),y=BagMeanEpiDay15LogInfDens,color=Nutrients,shape=Mixing))+
  stat_smooth(aes(x=log(TN), y=BagMeanEpiDay15LogInfDens),method='lm',color='black')+
  scale_color_manual(values=NutColors)+
  fav_theme

summary(lm(BagMeanEpiDay15LogInfDens~log(TN),data=TNDay15BagMeanEpiDay15LogInfDens))

cor.test(log(TNDay15BagMeanEpiDay15LogInfDens$TN),TNDay15BagMeanEpiDay15LogInfDens$BagMeanEpiDay15LogInfDens)

##################################################################################################
## Prepare data for analyses of TP and TN
##################################################################################################

## Subset data 
NutData<-mydata[,c("Raft","Nutrients","Mixing","Spores","Trt","Bag","Day","TP","TN","LogTP","LogTN")]
TPData<-mydata[,c("Raft","Nutrients","Mixing","Spores","Trt","Bag","Day","TP","LogTP")]
TNData<-mydata[,c("Raft","Nutrients","Mixing","Spores","Trt","Bag","Day","TN","LogTN")]

## Omit missing data
NutData<-na.omit(NutData)
TPData<-na.omit(TPData)
TNData<-na.omit(TNData)

## Subset data for -spores (uninf) and +spores (inf) bags
TPDataUninf<-subset(TPData,Spores=="No")
TPDataInf<-subset(TPData,Spores=="Spores")
TNDataUninf<-subset(TNData,Spores=="No")
TNDataInf<-subset(TNData,Spores=="Spores")

##################################################################################################
## Plot time series of Log(TP)
##################################################################################################

## First, plot each bag as its own line

p10<-(ggplot(TPDataUninf,aes(x=Day,y=LogTP,group=Bag,color=Nutrients))+
       scale_color_manual(values=NutColors)+
       geom_line(aes(linetype=Mixing))+
       ylim(c(1,6))+
       ylab("Log (TP, ug/L)")+
       xlab("Day of experiment")+ 
       labs(title="- Spores")+
       guides(color="none",linetype="none")+ 
       fav_theme)

p11<-(ggplot(TPDataInf,aes(x=Day,y=LogTP,group=Bag,color=Nutrients))+
       scale_color_manual(values=NutColors)+
       geom_line(aes(linetype=Mixing))+
       ylim(c(1,6))+
       ylab("Log (TP, ug/L)")+
       xlab("Day of experiment")+ 
       labs(title="+ Spores")+
       guides(color="none",linetype="none")+  
       fav_theme)

## Then, plot treatment means +/- SE

UninfLogTPSummary<-summarySE(TPDataUninf, measurevar="LogTP", groupvars=c("Nutrients","Mixing","Day"))

p12<-(ggplot(UninfLogTPSummary)+
       geom_point(aes(x=Day, y=LogTP, color=Nutrients)) +  
       geom_line(aes(x=Day, y=LogTP, color=Nutrients, linetype=Mixing))+
       geom_errorbar(aes(x=Day, ymin=LogTP-se, ymax=LogTP+se,color=Nutrients), width=0) + 
       scale_color_manual(values=NutColors)+
       fav_theme+
       ylim(c(1,6))+
       ylab("Log (TP, ug/L)")+
       xlab("Day of experiment")+ 
       labs(title="- Spores"))

InfLogTPSummary<-summarySE(TPDataInf, measurevar="LogTP", groupvars=c("Nutrients","Mixing","Day"))

p13<-(ggplot(InfLogTPSummary)+
       geom_point(aes(x=Day, y=LogTP, color=Nutrients)) +  
       geom_line(aes(x=Day, y=LogTP, color=Nutrients, linetype=Mixing))+
       geom_errorbar(aes(x=Day, ymin=LogTP-se, ymax=LogTP+se,color=Nutrients), width=0) + 
       scale_color_manual(values=NutColors)+
       fav_theme+
       ylim(c(1,6))+
       ylab("Log (TP, ug/L)")+
       xlab("Day of experiment")+ 
       labs(title="+ Spores")) 

## Put panels together into one figure
grid.arrange(p10,p12,p11,p13, nrow=2,ncol=2)

##################################################################################################
## Plot time series of Log(TN)
##################################################################################################

## First, plot each bag as its own line

p14<-(ggplot(TNDataUninf,aes(x=Day,y=LogTN,group=Bag,color=Nutrients))+
        scale_color_manual(values=NutColors)+
        geom_line(aes(linetype=Mixing))+
        ylim(c(5,9))+
        ylab("Log (TN, ug/L)")+
        xlab("Day of experiment")+ 
        labs(title="- Spores")+
        guides(color="none",linetype="none")+ 
        fav_theme)

p15<-(ggplot(TNDataInf,aes(x=Day,y=LogTN,group=Bag,color=Nutrients))+
        scale_color_manual(values=NutColors)+
        geom_line(aes(linetype=Mixing))+
        ylim(c(5,9))+
        ylab("Log (TN, ug/L)")+
        xlab("Day of experiment")+ 
        labs(title="+ Spores")+
        guides(color="none",linetype="none")+  
        fav_theme)

## Then, plot treatment means +/- SE

UninfLogTNSummary<-summarySE(TNDataUninf, measurevar="LogTN", groupvars=c("Nutrients","Mixing","Day"))

p16<-(ggplot(UninfLogTNSummary)+
        geom_point(aes(x=Day, y=LogTN, color=Nutrients)) +  
        geom_line(aes(x=Day, y=LogTN, color=Nutrients, linetype=Mixing))+
        geom_errorbar(aes(x=Day, ymin=LogTN-se, ymax=LogTN+se,color=Nutrients), width=0) + 
        scale_color_manual(values=NutColors)+
        fav_theme+
        ylim(c(5,9))+
        ylab("Log (TN, ug/L)")+
        xlab("Day of experiment")+ 
        labs(title="- Spores"))

InfLogTNSummary<-summarySE(TNDataInf, measurevar="LogTN", groupvars=c("Nutrients","Mixing","Day"))

p17<-(ggplot(InfLogTNSummary)+
        geom_point(aes(x=Day, y=LogTN, color=Nutrients)) +  
        geom_line(aes(x=Day, y=LogTN, color=Nutrients, linetype=Mixing))+
        geom_errorbar(aes(x=Day, ymin=LogTN-se, ymax=LogTN+se,color=Nutrients), width=0) + 
        scale_color_manual(values=NutColors)+
        fav_theme+
        ylim(c(5,9))+
        ylab("Log (TN, ug/L)")+
        xlab("Day of experiment")+ 
        labs(title="+ Spores")) 

## Put panels together into one figure
grid.arrange(p14,p16,p15,p17, nrow=2,ncol=2)

##################################################################################################
## Plot infection prevalence vs. proportion inedible chl a
##################################################################################################

PrevInedChlData<-mydata[,c("Raft","Nutrients","Mixing","Spores","Trt","Bag","Day","TotChl","EdChl","PrevTotal","PrevAdult")]
PrevInedChlData<-na.omit(PrevInedChlData)

## If I had defined InedChl at the very top of the script within mydata, I wouldn't have to do this again here:
for (i in 1:nrow(PrevInedChlData)){
  if(PrevInedChlData[i,"EdChl"]>PrevInedChlData[i,"TotChl"]){
    PrevInedChlData[i,"InedChl"]<-0
  } else {
    PrevInedChlData[i,"InedChl"]<-PrevInedChlData[i,"TotChl"]-PrevInedChlData[i,"EdChl"]
  }
}

NonzeroInedPrevInedChlData<-subset(PrevInedChlData,InedChl>0)
min(NonzeroInedPrevInedChlData$InedChl) 

for (i in 1:nrow(PrevInedChlData)){
  if(PrevInedChlData[i,"InedChl"]==0){
    PrevInedChlData[i,"InedChl"]<-min(NonzeroInedPrevInedChlData$InedChl)/2 
  } else {
    PrevInedChlData[i,"InedChl"]<-PrevInedChlData[i,"InedChl"]
  }
}

PrevInedChlData$PropInedChl<-PrevInedChlData$InedChl/PrevInedChlData$TotChl

## Subset data for -spores (uninf) and +spores (inf) bags
PrevInedChlDataUninf<-subset(PrevInedChlData,Spores=="No")
PrevInedChlDataInf<-subset(PrevInedChlData,Spores=="Spores")

PrevInedChlDataInfDisease<-subset(PrevInedChlDataInf,PrevTotal>0)

ggplot(PrevInedChlDataInfDisease,aes(x=PrevTotal,y=PropInedChl,color=Nutrients))+
  scale_color_manual(values=NutColors)+
  geom_point(aes(shape=Mixing))+
  ylim(c(0,1))+
  ylab("Proportion inedible chl a")+
  xlab("Infection prevalence")+ 
  labs(title="+ Spores")+
  fav_theme

PrevTotPropInedSummary<-ddply(PrevInedChlDataInf,c("Mixing","Spores","Nutrients","Bag"),summarize,
                            MaxPrevTotal=max(PrevTotal),
                            MaxPropInedChl=max(PropInedChl),
                            MeanPrevTotal=mean(PrevTotal),
                            MeanPropInedChl=mean(PropInedChl)
)

PrevTotPropInedSummary

ggplot(PrevTotPropInedSummary,aes(x=MaxPrevTotal,y=MaxPropInedChl,color=Nutrients))+
  scale_color_manual(values=NutColors)+
  geom_point(aes(shape=Mixing))+
  ylab("Max. proportion inedible chl a")+
  xlab("Max. infection prevalence")+ 
  labs(title="+ Spores")+
  fav_theme

ggplot(PrevTotPropInedSummary,aes(x=MeanPrevTotal,y=MeanPropInedChl,color=Nutrients))+
  scale_color_manual(values=NutColors)+
  geom_point(aes(shape=Mixing))+
  ylab("Mean proportion inedible chl a")+
  xlab("Mean infection prevalence")+ 
  labs(title="+ Spores")+
  fav_theme

##################################################################################################
## C, N, and P in edible seston
##################################################################################################

## Subset data 
CNPLightData<-mydata[,c("Raft","Nutrients","Mixing","Spores","Trt","Bag","Day","C","N","P","CN","CP","NP","Light")]

## Omit missing data
CNPLightData<-na.omit(CNPLightData)

## Subset data for -spores (uninf) and +spores (inf) bags
CNPLightDataUninf<-subset(CNPLightData,Spores=="No")
CNPLightDataInf<-subset(CNPLightData,Spores=="Spores")

## Subset data for Days 1, 31, and 46
CNPLightDataDay1<-subset(CNPLightData,Day==1)
CNPLightDataDay31<-subset(CNPLightData,Day==31)
CNPLightDataDay46<-subset(CNPLightData,Day==46)

##################################################################################################
## C:P molar ratio
##################################################################################################

Day1<-subset(mydata,Day==1)
Day31<-subset(mydata,Day==31)
Day46<-subset(mydata,Day==46)

hist(Day1$CP)
hist(Day31$CP)
hist(Day46$CP)

hist(Day1$CN)
hist(Day31$CN)
hist(Day46$CN)

hist(Day1$NP)
hist(Day31$NP)
hist(Day46$NP)

ggplot(Day1)+
  geom_boxplot(aes(x=Trt,y=CP,color=Nutrients))+
  scale_color_manual(values=NutColors)+
  fav_theme+
  ylim(0,600)

ggplot(Day31)+
  geom_boxplot(aes(x=Trt,y=CP,color=Nutrients))+
  scale_color_manual(values=NutColors)+
  fav_theme+
  ylim(0,600)

ggplot(Day46)+
  geom_boxplot(aes(x=Trt,y=CP,color=Nutrients))+
  scale_color_manual(values=NutColors)+
  fav_theme+
  ylim(0,600)

hist(Day1$CP)
m1Day1CP<-glm(CP~Nutrients*Mixing*Spores+Raft,data=Day1)
summary(m1Day1CP)
m2Day1CP<-glm(CP~Nutrients*Mixing+Nutrients*Spores+Mixing*Spores+Raft,data=Day1)
summary(m2Day1CP)
m3Day1CP<-glm(CP~Nutrients*Mixing+Mixing*Spores+Raft,data=Day1)
summary(m3Day1CP)
m4Day1CP<-glm(CP~Nutrients*Mixing+Spores+Raft,data=Day1)
summary(m4Day1CP)
m5Day1CP<-glm(CP~Nutrients+Mixing+Spores+Raft,data=Day1)
summary(m5Day1CP)
Anova(m5Day1CP)

hist(Day31$CP)
m1Day31CP<-glm(CP~Nutrients*Mixing*Spores+Raft,data=Day31)
summary(m1Day31CP)
m2Day31CP<-glm(CP~Nutrients*Mixing+Nutrients*Spores+Mixing*Spores+Raft,data=Day31)
summary(m2Day31CP)
m3Day31CP<-glm(CP~Nutrients*Mixing+Nutrients*Spores+Raft,data=Day31)
summary(m3Day31CP)
m4Day31CP<-glm(CP~Nutrients*Mixing+Spores+Raft,data=Day31)
summary(m4Day31CP)
m5Day31CP<-glm(CP~Nutrients+Mixing+Spores+Raft,data=Day31)
summary(m5Day31CP)
Anova(m5Day31CP)

ggplot(Day1)+
  geom_boxplot(aes(x=Trt,y=CN,color=Nutrients))+
  scale_color_manual(values=NutColors)+
  fav_theme+
  ylim(5,25)

ggplot(Day31)+
  geom_boxplot(aes(x=Trt,y=CN,color=Nutrients))+
  scale_color_manual(values=NutColors)+
  fav_theme+
  ylim(5,25)

ggplot(Day46)+
  geom_boxplot(aes(x=Trt,y=CN,color=Nutrients))+
  scale_color_manual(values=NutColors)+
  fav_theme+
  ylim(5,25)


ggplot(Day1)+
  geom_boxplot(aes(x=Trt,y=NP,color=Nutrients))+
  scale_color_manual(values=NutColors)+
  fav_theme+
  ylim(0,45)

ggplot(Day31)+
  geom_boxplot(aes(x=Trt,y=NP,color=Nutrients))+
  scale_color_manual(values=NutColors)+
  fav_theme+
  ylim(0,45)

ggplot(Day46)+
  geom_boxplot(aes(x=Trt,y=NP,color=Nutrients))+
  scale_color_manual(values=NutColors)+
  fav_theme+
  ylim(0,45)

UninfCPSummary<-summarySE(CNPLightDataUninf, measurevar="CP", groupvars=c("Nutrients","Mixing","Day"))

p1<-(ggplot(UninfCPSummary)+
       geom_point(aes(x=Day, y=CP, color=Nutrients)) +  
       geom_line(aes(x=Day, y=CP, color=Nutrients, linetype=Mixing))+
       geom_errorbar(aes(x=Day, ymin=CP-se, ymax=CP+se,color=Nutrients), width=0) + 
       scale_color_manual(values=NutColors)+
       fav_theme+
       ylab("C:P of edible seston (molar ratio)")+
       xlab("Day of experiment")+ 
       labs(title="- Spores")+
       ylim(c(50,500))
)

InfCPSummary<-summarySE(CNPLightDataInf, measurevar="CP", groupvars=c("Nutrients","Mixing","Day"))

p2<-(ggplot(InfCPSummary)+
       geom_point(aes(x=Day, y=CP, color=Nutrients)) +  
       geom_line(aes(x=Day, y=CP, color=Nutrients, linetype=Mixing))+
       geom_errorbar(aes(x=Day, ymin=CP-se, ymax=CP+se,color=Nutrients), width=0) + 
       scale_color_manual(values=NutColors)+
       fav_theme+
       ylab("C:P of edible seston (molar ratio)")+
       xlab("Day of experiment")+ 
       labs(title="+ Spores")+
       ylim(c(50,500))
)

## Put panels together into one figure
grid.arrange(p1,p2, nrow=1,ncol=2)

##################################################################################################
## C:N molar ratio
##################################################################################################

## Plot treatment means +/- SE

UninfCNSummary<-summarySE(CNPLightDataUninf, measurevar="CN", groupvars=c("Nutrients","Mixing","Day"))

p3<-(ggplot(UninfCNSummary)+
       geom_point(aes(x=Day, y=CN, color=Nutrients)) +  
       geom_line(aes(x=Day, y=CN, color=Nutrients, linetype=Mixing))+
       geom_errorbar(aes(x=Day, ymin=CN-se, ymax=CN+se,color=Nutrients), width=0) + 
       scale_color_manual(values=NutColors)+
       fav_theme+
       ylab("C:N of edible seston (molar ratio)")+
       xlab("Day of experiment")+ 
       labs(title="- Spores")+
       ylim(c(5,20))
)

InfCNSummary<-summarySE(CNPLightDataInf, measurevar="CN", groupvars=c("Nutrients","Mixing","Day"))

p4<-(ggplot(InfCNSummary)+
       geom_point(aes(x=Day, y=CN, color=Nutrients)) +  
       geom_line(aes(x=Day, y=CN, color=Nutrients, linetype=Mixing))+
       geom_errorbar(aes(x=Day, ymin=CN-se, ymax=CN+se,color=Nutrients), width=0) + 
       scale_color_manual(values=NutColors)+
       fav_theme+
       ylab("C:N of edible seston (molar ratio)")+
       xlab("Day of experiment")+ 
       labs(title="+ Spores")+
       ylim(c(5,20))
)

## Put panels together into one figure
grid.arrange(p3,p4, nrow=1,ncol=2)

##################################################################################################
## N:P molar ratio
##################################################################################################

## Plot treatment means +/- SE

UninfNPSummary<-summarySE(CNPLightDataUninf, measurevar="NP", groupvars=c("Nutrients","Mixing","Day"))

p5<-(ggplot(UninfNPSummary)+
       geom_point(aes(x=Day, y=NP, color=Nutrients)) +  
       geom_line(aes(x=Day, y=NP, color=Nutrients, linetype=Mixing))+
       geom_errorbar(aes(x=Day, ymin=NP-se, ymax=NP+se,color=Nutrients), width=0) + 
       scale_color_manual(values=NutColors)+
       fav_theme+
       ylab("N:P of edible seston (molar ratio)")+
       xlab("Day of experiment")+ 
       labs(title="- Spores")+
       ylim(c(10,40))
)

InfNPSummary<-summarySE(CNPLightDataInf, measurevar="NP", groupvars=c("Nutrients","Mixing","Day"))

p6<-(ggplot(InfNPSummary)+
       geom_point(aes(x=Day, y=NP, color=Nutrients)) +  
       geom_line(aes(x=Day, y=NP, color=Nutrients, linetype=Mixing))+
       geom_errorbar(aes(x=Day, ymin=NP-se, ymax=NP+se,color=Nutrients), width=0) + 
       scale_color_manual(values=NutColors)+
       fav_theme+
       ylab("N:P of edible seston (molar ratio)")+
       xlab("Day of experiment")+ 
       labs(title="+ Spores")+
       ylim(c(10,40))
)

## Put panels together into one figure
grid.arrange(p5,p6, nrow=1,ncol=2)

## Put panels together into one figure
grid.arrange(p1,p2,p3,p4,p5,p6, nrow=3,ncol=2)

##################################################################################################
## Log(C) vs. light extinction, Day 31
##################################################################################################

ggplot(CNPLightDataDay31)+
  geom_point(aes(x=LogC, y=Light, color=Nutrients,shape=Mixing,group=Spores)) +  
  facet_wrap(~Spores)+
  stat_smooth(aes(x=LogC, y=Light),method='lm',color='black')+
  scale_color_manual(values=NutColors)+
  fav_theme+
  ylab("Light extinction coefficient, k (umol quanta cm-2 s-1 m-1)")+
  xlab("Log (C, ug/L) in edible seston")

ggplot(CNPLightDataDay31)+
  geom_point(aes(x=LogC, y=Light, color=Nutrients,shape=Mixing)) +  
  stat_smooth(aes(x=LogC, y=Light),method='lm',color='black')+
  scale_color_manual(values=NutColors)+
  fav_theme+
  ylab("Light extinction coefficient, k (umol quanta cm-2 s-1 m-1)")+
  xlab("Log (C, ug/L) in edible seston")

##################################################################################################
## Log(C) vs. C:P, Day 31
##################################################################################################

ggplot(CNPLightDataDay31)+
  geom_point(aes(x=LogC, y=CP, color=Nutrients,shape=Mixing)) +  
  stat_smooth(aes(x=LogC, y=CP),method='lm',color='black')+
  scale_color_manual(values=NutColors)+
  fav_theme+
  ylab("C:P of edible seston (molar ratio)")+
  xlab("Log (C, ug/L) in edible seston")

##################################################################################################
## Light vs. C:P, Day 31
##################################################################################################

ggplot(CNPLightDataDay31)+
  geom_point(aes(x=Light, y=CP, color=Nutrients,shape=Mixing)) +  
  stat_smooth(aes(x=Light, y=CP),method='lm',color='black')+
  scale_color_manual(values=NutColors)+
  fav_theme+
  ylab("C:P of edible seston (molar ratio)")+
  xlab("Light")

##################################################################################################
## To do still:  analyses of infected and total host density vs. the various size fractions of chlorophyll a
##################################################################################################

ggplot(StartDay26,aes(x=LogTotDens,y=LogTotChl,color=Nutrients))+
  geom_point()+
  geom_smooth(method=lm)+
  scale_color_manual(values=NutColors)+
  facet_wrap(~Spores)+
  fav_theme

ggplot(StartDay12,aes(x=LogTotDens,y=LogTotChl,color=Nutrients))+
  geom_point()+
  geom_smooth(method=lm)+
  scale_color_manual(values=NutColors)+
  facet_wrap(~Spores)+
  fav_theme

## Correlate hosts and algae with in each bag, and compare correlation coefficients among treatments???

##################################################################################################
## End of script
##################################################################################################