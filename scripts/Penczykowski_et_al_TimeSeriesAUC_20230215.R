## Code to produce the time series and AUC figures in manuscript titled:
## "Pathways linking nutrient enrichment, habitat structure, and parasitism to host–resource interactions"
## Penczykowski et al. submitting to Oecologia on February 15, 2023

## code written by: Rachel M. Penczykowski (modified AUC calculations written by Michelle L. Fearon)
## last updated: February 15, 2023

# load packages
library(plyr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(DescTools)
library(lme4)
library(car)
library(multcomp)
library(cowplot)
library(here)

# set the path to the script relative to the project root directory
here::i_am("scripts/Penczykowski_et_al_TimeSeriesAUC_20230215.R")

# load data
mydata <- read.csv(here("data/Penczykowski_et_al_BagExptData.csv"), stringsAsFactors = F, header = T)
str(mydata)

# create dummy variables for categorical variables 
mydata <- mutate(mydata,  
                 NutrientTrt = ifelse(Nutrients == "Low", 0, 1),
                 MixingTrt = ifelse(Mixing == "No", 0, 1),
                 SporeTrt = ifelse(Spores == "No", 0, 1))

mydata$Bag<-as.factor(mydata$Bag)
levels(mydata$Bag) # checking that the excluded bags were already removed (for reasons explained in methods section of manuscript)

# remove NAs from focal variables
mydata2 <- filter(mydata, !is.na(PrevTotal))
mydata2 <- filter(mydata2, !is.na(InfDensity))
mydata2 <- filter(mydata2, !is.na(LogInfDens))
mydata2 <- filter(mydata2, !is.na(TotalDensity))
mydata2 <- filter(mydata2, !is.na(LogTotDens))
mydata2 <- filter(mydata2, !is.na(EdChl))
mydata2 <- filter(mydata2, !is.na(LogEdChl))
mydata2 <- filter(mydata2, !is.na(TP))
mydata2 <- filter(mydata2, !is.na(LogTP))
mydata2 <- filter(mydata2, !is.na(TN))
mydata2 <- filter(mydata2, !is.na(LogTN))

##################################################################################################
## Calculate AUC (Area Under Curve) for all response variables during 2nd epidemic (days 22-43)
##################################################################################################

loop_data <- mydata2[,c("Nutrients","Mixing","Spores","NutrientTrt","MixingTrt","SporeTrt","Raft","Bag","Day",
                        "PrevTotal","LogInfDens","LogTotDens","LogEdChl","LogTP","LogTN")]

# arrange data 
loop_data <- loop_data %>%
  group_by(NutrientTrt,MixingTrt,SporeTrt) %>%
  dplyr::arrange( .by_group = T)
str(loop_data)

# set treatments, Bag IDs, and days of experiment as correct class
loop_data$NutrientTrt <- as.factor(loop_data$NutrientTrt)
loop_data$MixingTrt <- as.factor(loop_data$MixingTrt)
loop_data$SporeTrt <- as.factor(loop_data$SporeTrt)
loop_data$Raft <- as.factor(loop_data$Raft)
loop_data$Bag <- as.factor(loop_data$Bag)
loop_data$Day <- as.integer(loop_data$Day)

# order everything by day of experiment
loop_data <- loop_data[order(loop_data$Bag), ] #make sure days are in order

# make sure it's in data frame format
loop_data <- as.data.frame(loop_data)

# loop through each nutrient, mixing, spore treatment, and bag ID to calculate AUCs
OUT <- NULL
NUTRIENTS<-unique(loop_data$NutrientTrt)
for (i in NUTRIENTS) {
  thisnutrients <- loop_data[loop_data$NutrientTrt == i, ] 
  MIXING <- unique(thisnutrients$MixingTrt)
  for (j in MIXING) {
    thismixing <- thisnutrients[thisnutrients$MixingTrt == j, ] 
    SPORES <- unique(thismixing$SporeTrt)
    for (h in SPORES) {
      thisspores <- thismixing[thismixing$SporeTrt == h, ] 
      RAFT <- unique(thisspores$Raft)
        for (r in RAFT) {
        thisraft <- thisspores[thisspores$Raft == r, ] 
        BAG <- unique(thisraft$Bag)
    
              for (p in BAG) {
              thisbag <- thisraft[thisraft$Bag == p, ] 
              
              TOTS <- NULL #make an empty matrix to include sequential values. This will be used to calculate integrated areas.
              for (m in 1:length(thisbag$PrevTotal)) {
                total <- thisbag[m, "PrevTotal"] 
                day <- thisbag[m, "Day"] 
                tots <- c(m, day, total) #new row in dataframe with m, day, and total
                TOTS <- rbind(TOTS, tots) #bind all rows together (different one for each bag)
                colnames(TOTS) <- c("m", "day", "total")
              }
              rownames(TOTS) <- NULL
              TOTS <- as.data.frame(TOTS)
              AREA_Day22Thru43_PrevTotal <- AUC(TOTS[ , "day"], TOTS[ , "total"], from = min(22, na.rm = TRUE), to = max(43, na.rm = TRUE), method = "trapezoid")
             
              TOTS <- NULL #make an empty matrix to include sequential values. This will be used to calculate integrated areas.
              for (m in 1:length(thisbag$LogInfDens)) {
                total <- thisbag[m, "LogInfDens"] 
                day <- thisbag[m, "Day"] 
                tots <- c(m, day, total) #new row in dataframe with m, day, and total
                TOTS <- rbind(TOTS, tots) #bind all rows together (different one for each bag)
                colnames(TOTS) <- c("m", "day", "total")
              }
              rownames(TOTS) <- NULL
              TOTS <- as.data.frame(TOTS)
              AREA_Day22Thru43_LogInfDens <- AUC(TOTS[ , "day"], TOTS[ , "total"], from = min(22, na.rm = TRUE), to = max(43, na.rm = TRUE), method = "trapezoid")
              
              TOTS <- NULL #make an empty matrix to include sequential values. This will be used to calculate integrated areas.
              for (m in 1:length(thisbag$LogTotDens)) {
                total <- thisbag[m, "LogTotDens"] 
                day <- thisbag[m, "Day"] 
                tots <- c(m, day, total) #new row in dataframe with m, day, and total
                TOTS <- rbind(TOTS, tots) #bind all rows together (different one for each bag)
                colnames(TOTS) <- c("m", "day", "total")
              }
              rownames(TOTS) <- NULL
              TOTS <- as.data.frame(TOTS)
              AREA_Day22Thru43_LogTotDens <- AUC(TOTS[ , "day"], TOTS[ , "total"], from = min(22, na.rm = TRUE), to = max(43, na.rm = TRUE), method = "trapezoid")
              
              TOTS <- NULL #make an empty matrix to include sequential values. This will be used to calculate integrated areas.
              for (m in 1:length(thisbag$LogEdChl)) {
                total <- thisbag[m, "LogEdChl"] 
                day <- thisbag[m, "Day"] 
                tots <- c(m, day, total) #new row in dataframe with m, day, and total
                TOTS <- rbind(TOTS, tots) #bind all rows together (different one for each bag)
                colnames(TOTS) <- c("m", "day", "total")
              }
              rownames(TOTS) <- NULL
              TOTS <- as.data.frame(TOTS)
              AREA_Day22Thru43_LogEdChl <- AUC(TOTS[ , "day"], TOTS[ , "total"], from = min(22, na.rm = TRUE), to = max(43, na.rm = TRUE), method = "trapezoid")
              
              TOTS <- NULL #make an empty matrix to include sequential values. This will be used to calculate integrated areas.
              for (m in 1:length(thisbag$LogTP)) {
                total <- thisbag[m, "LogTP"] 
                day <- thisbag[m, "Day"] 
                tots <- c(m, day, total) #new row in dataframe with m, day, and total
                TOTS <- rbind(TOTS, tots) #bind all rows together (different one for each bag)
                colnames(TOTS) <- c("m", "day", "total")
              }
              rownames(TOTS) <- NULL
              TOTS <- as.data.frame(TOTS)
              AREA_Day22Thru43_LogTP <- AUC(TOTS[ , "day"], TOTS[ , "total"], from = min(22, na.rm = TRUE), to = max(43, na.rm = TRUE), method = "trapezoid")
              
              TOTS <- NULL #make an empty matrix to include sequential values. This will be used to calculate integrated areas.
              for (m in 1:length(thisbag$LogTN)) {
                total <- thisbag[m, "LogTN"] 
                day <- thisbag[m, "Day"] 
                tots <- c(m, day, total) #new row in dataframe with m, day, and total
                TOTS <- rbind(TOTS, tots) #bind all rows together (different one for each bag)
                colnames(TOTS) <- c("m", "day", "total")
              }
              rownames(TOTS) <- NULL
              TOTS <- as.data.frame(TOTS)
              AREA_Day22Thru43_LogTN <- AUC(TOTS[ , "day"], TOTS[ , "total"], from = min(22, na.rm = TRUE), to = max(43, na.rm = TRUE), method = "trapezoid")
              
              output <- c(i, j, h, r, p, AREA_Day22Thru43_PrevTotal, AREA_Day22Thru43_LogInfDens, AREA_Day22Thru43_LogTotDens,
                          AREA_Day22Thru43_LogEdChl, AREA_Day22Thru43_LogTP, AREA_Day22Thru43_LogTN)
              OUT <- rbind(OUT, output)
        }
      }
    }
  }
}

colnames(OUT) <- c("NutrientTrt", "MixingTrt", "SporeTrt", "Raft","Bag", 
                   "AUC.Days22Thru43.PrevTotal","AUC.Days22Thru43.LogInfDens","AUC.Days22Thru43.LogTotDens",
                   "AUC.Days22Thru43.LogEdChl","AUC.Days22Thru43.LogTP","AUC.Days22Thru43.LogTN")

auc_data <- as.data.frame(OUT)
row.names(auc_data) <- c()

sapply(auc_data, class)
c <- c("NutrientTrt", "MixingTrt", "SporeTrt", "Raft","Bag")
auc_data[c] <- lapply(auc_data[c], as.factor)
d <- c("AUC.Days22Thru43.PrevTotal","AUC.Days22Thru43.LogInfDens","AUC.Days22Thru43.LogTotDens",
       "AUC.Days22Thru43.LogEdChl","AUC.Days22Thru43.LogTP","AUC.Days22Thru43.LogTN")
auc_data[d] <- lapply(auc_data[d], as.numeric)
str(auc_data)

# put original treatment names/levels back in
auc_data2 <- mutate(auc_data,  
                 Nutrients = ifelse(NutrientTrt == 0, "Low","High"),
                 Mixing = ifelse(MixingTrt == 0, "No","Mix"),
                 Spores = ifelse(SporeTrt == 0,"No","Spores"))

##################################################################################################
## Set standard formatting for figures 
##################################################################################################

fav_theme<-theme_bw()+
  theme(panel.border = element_rect(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        strip.text.x=element_text(size=12),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=10), 
        axis.title=element_text(size=12),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        axis.title.y=element_text(vjust=1.2),
        axis.title.x=element_text(vjust=0.4))

SporesNames<-c(No="No parasite",Spores="Parasite")
NutColors<-c("darkgreen","#56B4E9") # green for high nutrients (color of algae bloom); blue for low nutrients. These are from a colorblind-friendly palette.

##################################################################################################
## Plots and statistics for Fig. 1 panels and Table 1

##################################################################################################
## Infection prevalence
##################################################################################################
InfPrevSummary<-ddply(mydata2, c("Nutrients", "Mixing","Spores","Day"), summarise,
                         N    = length(PrevTotal),
                         mean = mean(PrevTotal),
                         sd   = sd(PrevTotal),
                         se   = sd / sqrt(N)
)

InfPrev_Open<-subset(InfPrevSummary, Day<22|Day>43)
InfPrev_Closed<-subset(InfPrevSummary,Day>=22&Day<=43)
InfPrev_OpenMixed<-subset(InfPrev_Open,Mixing=="Mix")
InfPrev_OpenNoMixed<-subset(InfPrev_Open,Mixing=="No")
InfPrev_ClosedMixed<-subset(InfPrev_Closed,Mixing=="Mix")
InfPrev_ClosedNoMixed<-subset(InfPrev_Closed,Mixing=="No")

p1<-(ggplot()+
  geom_point(data=InfPrev_OpenMixed,aes(x=Day+0.4,y=mean,color=Nutrients,shape=Mixing),size=2)+
  geom_point(data=InfPrev_ClosedMixed,aes(x=Day+0.4,y=mean,color=Nutrients,fill=Nutrients,shape=Mixing),size=2)+
  geom_errorbar(data=InfPrev_OpenMixed,aes(x=Day+0.4, ymin=mean-se, ymax=mean+se,color=Nutrients,linetype=Mixing), width=1) + 
  geom_errorbar(data=InfPrev_ClosedMixed,aes(x=Day+0.4, ymin=mean-se, ymax=mean+se,color=Nutrients,linetype=Mixing), width=1) + 
  geom_point(data=InfPrev_OpenNoMixed,aes(x=Day-0.4,y=mean,color=Nutrients,shape=Mixing),size=2)+
  geom_point(data=InfPrev_ClosedNoMixed,aes(x=Day-0.4,y=mean,color=Nutrients,fill=Nutrients,shape=Mixing),size=2)+
  geom_errorbar(data=InfPrev_OpenNoMixed,aes(x=Day-0.4, ymin=mean-se, ymax=mean+se,color=Nutrients,linetype=Mixing), width=1) + 
  geom_errorbar(data=InfPrev_ClosedNoMixed,aes(x=Day-0.4, ymin=mean-se, ymax=mean+se,color=Nutrients,linetype=Mixing), width=1) + 
  scale_shape_manual(values=c(21,24))+
  scale_color_manual(values=NutColors)+
  scale_fill_manual(values=NutColors)+
  facet_wrap(~Spores,labeller=labeller(Spores=SporesNames))+
  ylab("Infection prevalence")+
  xlab("Day of experiment")+
  fav_theme+
  annotate("rect", xmin = 22-0.8, xmax = 43+0.8, ymin = 0, ymax = 0.43,fill="gray90",alpha=0.3)
)

p2<-(ggplot(auc_data2)+
  geom_boxplot(aes(x=Mixing,y=AUC.Days22Thru43.PrevTotal,color=Nutrients))+
  geom_point(aes(x=Mixing,y=AUC.Days22Thru43.PrevTotal,color=Nutrients,fill=Nutrients,shape=Mixing),position=position_dodge(width=0.75),size=2)+
  scale_shape_manual(values=c(21,24))+
  scale_color_manual(values=NutColors)+
  scale_fill_manual(values=NutColors)+
  facet_wrap(~Spores,labeller=labeller(Spores=SporesNames))+
  ylab("AUC, infection prevalence")+
  fav_theme
)

auc_SporesOnly<-subset(auc_data2,Spores=="Spores")

mPrevTotal.glmm1<-lmer(AUC.Days22Thru43.PrevTotal~Nutrients*Mixing+(1|Raft),data=auc_SporesOnly)
Anova(mPrevTotal.glmm1)
mPrevTotal.glmm2<-lmer(AUC.Days22Thru43.PrevTotal~Nutrients+Mixing+(1|Raft),data=auc_SporesOnly)
summary(mPrevTotal.glmm2)
Anova(mPrevTotal.glmm2)

##################################################################################################
## Log(Infected host density)
##################################################################################################
LogInfDensSummary<-ddply(mydata2, c("Nutrients", "Mixing","Spores","Day"), summarise,
                         N    = length(LogInfDens),
                         mean = mean(LogInfDens),
                         sd   = sd(LogInfDens),
                         se   = sd / sqrt(N)
)

LogInfDens_Open<-subset(LogInfDensSummary, Day<22|Day>43)
LogInfDens_Closed<-subset(LogInfDensSummary,Day>=22&Day<=43)
LogInfDens_OpenMixed<-subset(LogInfDens_Open,Mixing=="Mix")
LogInfDens_OpenNoMixed<-subset(LogInfDens_Open,Mixing=="No")
LogInfDens_ClosedMixed<-subset(LogInfDens_Closed,Mixing=="Mix")
LogInfDens_ClosedNoMixed<-subset(LogInfDens_Closed,Mixing=="No")

p3<-(ggplot()+
  geom_point(data=LogInfDens_OpenMixed,aes(x=Day+0.4,y=mean,color=Nutrients,shape=Mixing),size=2)+
  geom_point(data=LogInfDens_ClosedMixed,aes(x=Day+0.4,y=mean,color=Nutrients,fill=Nutrients,shape=Mixing),size=2)+
  geom_errorbar(data=LogInfDens_OpenMixed,aes(x=Day+0.4, ymin=mean-se, ymax=mean+se,color=Nutrients,linetype=Mixing), width=1) + 
  geom_errorbar(data=LogInfDens_ClosedMixed,aes(x=Day+0.4, ymin=mean-se, ymax=mean+se,color=Nutrients,linetype=Mixing), width=1) + 
  geom_point(data=LogInfDens_OpenNoMixed,aes(x=Day-0.4,y=mean,color=Nutrients,shape=Mixing),size=2)+
  geom_point(data=LogInfDens_ClosedNoMixed,aes(x=Day-0.4,y=mean,color=Nutrients,fill=Nutrients,shape=Mixing),size=2)+
  geom_errorbar(data=LogInfDens_OpenNoMixed,aes(x=Day-0.4, ymin=mean-se, ymax=mean+se,color=Nutrients,linetype=Mixing), width=1) + 
  geom_errorbar(data=LogInfDens_ClosedNoMixed,aes(x=Day-0.4, ymin=mean-se, ymax=mean+se,color=Nutrients,linetype=Mixing), width=1) + 
  scale_shape_manual(values=c(21,24))+
  scale_color_manual(values=NutColors)+
  scale_fill_manual(values=NutColors)+
  facet_wrap(~Spores,labeller=labeller(Spores=SporesNames))+
  ylab("Log(infected host density)")+
  xlab("Day of experiment")+
  fav_theme+
  annotate("rect", xmin = 22-0.8, xmax = 43+0.8, ymin = 4.25, ymax = 11.5,fill="gray90",alpha=0.3)
)

p4<-(ggplot(auc_data2)+
  geom_boxplot(aes(x=Mixing,y=AUC.Days22Thru43.LogInfDens,color=Nutrients))+
  geom_point(aes(x=Mixing,y=AUC.Days22Thru43.LogInfDens,color=Nutrients,fill=Nutrients,shape=Mixing),position=position_dodge(width=0.75),size=2)+
  scale_shape_manual(values=c(21,24))+
  scale_color_manual(values=NutColors)+
  scale_fill_manual(values=NutColors)+
  facet_wrap(~Spores,labeller=labeller(Spores=SporesNames))+
  ylab("AUC, log(infected host density)")+
  fav_theme
)

auc_SporesOnly<-subset(auc_data2,Spores=="Spores")

mLogInfDens.glmm1<-lmer(AUC.Days22Thru43.LogInfDens~Nutrients*Mixing+(1|Raft),data=auc_SporesOnly)
summary(mLogInfDens.glmm1)
Anova(mLogInfDens.glmm1)
mLogInfDens.glmm2<-lmer(AUC.Days22Thru43.LogInfDens~Nutrients+Mixing+(1|Raft),data=auc_SporesOnly)
summary(mLogInfDens.glmm2)
Anova(mLogInfDens.glmm2)

##################################################################################################
## Log(Total host density)
##################################################################################################
LogTotDensSummary<-ddply(mydata2, c("Nutrients", "Mixing","Spores","Day"), summarise,
                         N    = length(LogTotDens),
                         mean = mean(LogTotDens),
                         sd   = sd(LogTotDens),
                         se   = sd / sqrt(N)
)

LogTotDens_Open<-subset(LogTotDensSummary, Day<22|Day>43)
LogTotDens_Closed<-subset(LogTotDensSummary,Day>=22&Day<=43)
LogTotDens_OpenMixed<-subset(LogTotDens_Open,Mixing=="Mix")
LogTotDens_OpenNoMixed<-subset(LogTotDens_Open,Mixing=="No")
LogTotDens_ClosedMixed<-subset(LogTotDens_Closed,Mixing=="Mix")
LogTotDens_ClosedNoMixed<-subset(LogTotDens_Closed,Mixing=="No")

p5<-(ggplot()+
  geom_point(data=LogTotDens_OpenMixed,aes(x=Day+0.4,y=mean,color=Nutrients,shape=Mixing),size=2)+
  geom_point(data=LogTotDens_ClosedMixed,aes(x=Day+0.4,y=mean,color=Nutrients,fill=Nutrients,shape=Mixing),size=2)+
  geom_errorbar(data=LogTotDens_OpenMixed,aes(x=Day+0.4, ymin=mean-se, ymax=mean+se,color=Nutrients,linetype=Mixing), width=1) + 
  geom_errorbar(data=LogTotDens_ClosedMixed,aes(x=Day+0.4, ymin=mean-se, ymax=mean+se,color=Nutrients,linetype=Mixing), width=1) + 
  geom_point(data=LogTotDens_OpenNoMixed,aes(x=Day-0.4,y=mean,color=Nutrients,shape=Mixing),size=2)+
  geom_point(data=LogTotDens_ClosedNoMixed,aes(x=Day-0.4,y=mean,color=Nutrients,fill=Nutrients,shape=Mixing),size=2)+
  geom_errorbar(data=LogTotDens_OpenNoMixed,aes(x=Day-0.4, ymin=mean-se, ymax=mean+se,color=Nutrients,linetype=Mixing), width=1) + 
  geom_errorbar(data=LogTotDens_ClosedNoMixed,aes(x=Day-0.4, ymin=mean-se, ymax=mean+se,color=Nutrients,linetype=Mixing), width=1) + 
  scale_shape_manual(values=c(21,24))+
  scale_color_manual(values=NutColors)+
  scale_fill_manual(values=NutColors)+
  facet_wrap(~Spores,labeller=labeller(Spores=SporesNames))+
  ylab("Log(total host density)")+
  xlab("Day of experiment")+
  fav_theme+
  annotate("rect", xmin = 22-0.8, xmax = 43+0.8, ymin = 6.8, ymax = 13,fill="gray90",alpha=0.3)
)

p6<-(ggplot(auc_data2)+
  geom_boxplot(aes(x=Mixing,y=AUC.Days22Thru43.LogTotDens,color=Nutrients))+
  geom_point(aes(x=Mixing,y=AUC.Days22Thru43.LogTotDens,color=Nutrients,fill=Nutrients,shape=Mixing),position=position_dodge(width=0.75),size=2)+
  scale_shape_manual(values=c(21,24))+
  scale_color_manual(values=NutColors)+
  scale_fill_manual(values=NutColors)+
  facet_wrap(~Spores,labeller=labeller(Spores=SporesNames))+
  ylab("AUC, log(total host density)")+
  fav_theme
)

mLogTotDens.glmm1<-lmer(AUC.Days22Thru43.LogTotDens~Nutrients*Mixing*Spores+(1|Raft),data=auc_data2)
summary(mLogTotDens.glmm1)
Anova(mLogTotDens.glmm1)

## performed stepwise deletion of non-significant interaction terms (not shown); all interactions were non-significant
mLogTotDens.glmm2<-lmer(AUC.Days22Thru43.LogTotDens~Nutrients+Mixing+Spores+(1|Raft),data=auc_data2)
summary(mLogTotDens.glmm2)
Anova(mLogTotDens.glmm2)

##################################################################################################
## Complete Fig. 1 
##################################################################################################
Fig1<-plot_grid(p1+theme(legend.position="none"),
          p2+theme(legend.position="none"),
          p3+theme(legend.position="none"),
          p4+theme(legend.position="none"),
          p5+theme(legend.position="none"),
          p6+theme(legend.position="none"), 
          labels = c('a','b','c','d','e','f'),
          ncol=2,rel_widths = c(1,0.5,1,0.5,1,0.5))

print(Fig1)
ggsave(here("figures/Fig1.tiff"), plot = Fig1, dpi = 300, width = 8, height = 8.5, units = "in", compression="lzw")

library(ggpubr)
# Extract the legend. Returns a gtable
leg1 <- get_legend(p1)
leg2 <- get_legend(p2)

# Convert to a ggplot and print
print(as_ggplot(leg1))
print(as_ggplot(leg2))

##################################################################################################
## Plots and statistics for Fig. S1 panels and Table S1

##################################################################################################
## Log(Edible chlorophyll a)
##################################################################################################
LogEdChlSummary<-ddply(mydata2, c("Nutrients", "Mixing","Spores","Day"), summarise,
                         N    = length(LogEdChl),
                         mean = mean(LogEdChl),
                         sd   = sd(LogEdChl),
                         se   = sd / sqrt(N)
)

LogEdChl_Open<-subset(LogEdChlSummary, Day<22|Day>43)
LogEdChl_Closed<-subset(LogEdChlSummary,Day>=22&Day<=43)
LogEdChl_OpenMixed<-subset(LogEdChl_Open,Mixing=="Mix")
LogEdChl_OpenNoMixed<-subset(LogEdChl_Open,Mixing=="No")
LogEdChl_ClosedMixed<-subset(LogEdChl_Closed,Mixing=="Mix")
LogEdChl_ClosedNoMixed<-subset(LogEdChl_Closed,Mixing=="No")

PS1<-(ggplot()+
  geom_point(data=LogEdChl_OpenMixed,aes(x=Day+0.4,y=mean,color=Nutrients,shape=Mixing),size=2)+
  geom_point(data=LogEdChl_ClosedMixed,aes(x=Day+0.4,y=mean,color=Nutrients,fill=Nutrients,shape=Mixing),size=2)+
  geom_errorbar(data=LogEdChl_OpenMixed,aes(x=Day+0.4, ymin=mean-se, ymax=mean+se,color=Nutrients,linetype=Mixing), width=1) + 
  geom_errorbar(data=LogEdChl_ClosedMixed,aes(x=Day+0.4, ymin=mean-se, ymax=mean+se,color=Nutrients,linetype=Mixing), width=1) + 
  geom_point(data=LogEdChl_OpenNoMixed,aes(x=Day-0.4,y=mean,color=Nutrients,shape=Mixing),size=2)+
  geom_point(data=LogEdChl_ClosedNoMixed,aes(x=Day-0.4,y=mean,color=Nutrients,fill=Nutrients,shape=Mixing),size=2)+
  geom_errorbar(data=LogEdChl_OpenNoMixed,aes(x=Day-0.4, ymin=mean-se, ymax=mean+se,color=Nutrients,linetype=Mixing), width=1) + 
  geom_errorbar(data=LogEdChl_ClosedNoMixed,aes(x=Day-0.4, ymin=mean-se, ymax=mean+se,color=Nutrients,linetype=Mixing), width=1) + 
  scale_shape_manual(values=c(21,24))+
  scale_color_manual(values=NutColors)+
  scale_fill_manual(values=NutColors)+
  facet_wrap(~Spores,labeller=labeller(Spores=SporesNames))+
  ylab(expression("Log(edible chlorophyll "~italic(a)~")"))+
  xlab("Day of experiment")+
  fav_theme+
  annotate("rect", xmin = 22-0.8, xmax = 43+0.8, ymin = 0, ymax = 4.2,fill="gray90",alpha=0.3)
)

PS2<-(ggplot(auc_data2)+
  geom_boxplot(aes(x=Mixing,y=AUC.Days22Thru43.LogEdChl,color=Nutrients))+
  geom_point(aes(x=Mixing,y=AUC.Days22Thru43.LogEdChl,color=Nutrients,fill=Nutrients,shape=Mixing),position=position_dodge(width=0.75),size=2)+
  scale_shape_manual(values=c(21,24))+
  scale_color_manual(values=NutColors)+
  scale_fill_manual(values=NutColors)+
  facet_wrap(~Spores,labeller=labeller(Spores=SporesNames))+
  ylab(expression("AUC, log(edible chlorophyll "~italic(a)~")"))+
  fav_theme)

mLogEdChl.glmm1<-lmer(AUC.Days22Thru43.LogEdChl~Nutrients*Mixing*Spores+(1|Raft),data=auc_data2)
summary(mLogEdChl.glmm1)
Anova(mLogEdChl.glmm1)

## performed stepwise deletion of non-significant interaction terms (not shown); all interactions were non-significant
mLogEdChl.glmm2<-lmer(AUC.Days22Thru43.LogEdChl~Nutrients+Mixing+Spores+(1|Raft),data=auc_data2)
summary(mLogEdChl.glmm2)
Anova(mLogEdChl.glmm2)

##################################################################################################
## Log(Total phosphorus)
##################################################################################################
LogTPSummary<-ddply(mydata2, c("Nutrients", "Mixing","Spores","Day"), summarise,
                       N    = length(LogTP),
                       mean = mean(LogTP),
                       sd   = sd(LogTP),
                       se   = sd / sqrt(N)
)

LogTP_Open<-subset(LogTPSummary, Day<22|Day>43)
LogTP_Closed<-subset(LogTPSummary,Day>=22&Day<=43)
LogTP_OpenMixed<-subset(LogTP_Open,Mixing=="Mix")
LogTP_OpenNoMixed<-subset(LogTP_Open,Mixing=="No")
LogTP_ClosedMixed<-subset(LogTP_Closed,Mixing=="Mix")
LogTP_ClosedNoMixed<-subset(LogTP_Closed,Mixing=="No")

PS3<-(ggplot()+
  geom_point(data=LogTP_OpenMixed,aes(x=Day+0.4,y=mean,color=Nutrients,shape=Mixing),size=2)+
  geom_point(data=LogTP_ClosedMixed,aes(x=Day+0.4,y=mean,color=Nutrients,fill=Nutrients,shape=Mixing),size=2)+
  geom_errorbar(data=LogTP_OpenMixed,aes(x=Day+0.4, ymin=mean-se, ymax=mean+se,color=Nutrients,linetype=Mixing), width=1) + 
  geom_errorbar(data=LogTP_ClosedMixed,aes(x=Day+0.4, ymin=mean-se, ymax=mean+se,color=Nutrients,linetype=Mixing), width=1) + 
  geom_point(data=LogTP_OpenNoMixed,aes(x=Day-0.4,y=mean,color=Nutrients,shape=Mixing),size=2)+
  geom_point(data=LogTP_ClosedNoMixed,aes(x=Day-0.4,y=mean,color=Nutrients,fill=Nutrients,shape=Mixing),size=2)+
  geom_errorbar(data=LogTP_OpenNoMixed,aes(x=Day-0.4, ymin=mean-se, ymax=mean+se,color=Nutrients,linetype=Mixing), width=1) + 
  geom_errorbar(data=LogTP_ClosedNoMixed,aes(x=Day-0.4, ymin=mean-se, ymax=mean+se,color=Nutrients,linetype=Mixing), width=1) + 
  scale_shape_manual(values=c(21,24))+
  scale_color_manual(values=NutColors)+
  scale_fill_manual(values=NutColors)+
  facet_wrap(~Spores,labeller=labeller(Spores=SporesNames))+
  ylab("Log(total phosphorus)")+
  xlab("Day of experiment")+
  fav_theme+
  annotate("rect", xmin = 22-0.8, xmax = 43+0.8, ymin = 1.7, ymax = 4.95,fill="gray90",alpha=0.3)
)

PS4<-(ggplot(auc_data2)+
  geom_boxplot(aes(x=Mixing,y=AUC.Days22Thru43.LogTP,color=Nutrients))+
  geom_point(aes(x=Mixing,y=AUC.Days22Thru43.LogTP,color=Nutrients,fill=Nutrients,shape=Mixing),position=position_dodge(width=0.75),size=2)+
  scale_shape_manual(values=c(21,24))+
  scale_color_manual(values=NutColors)+
  scale_fill_manual(values=NutColors)+
  facet_wrap(~Spores,labeller=labeller(Spores=SporesNames))+
  ylab("AUC, log(total phosphorus)")+
  fav_theme)

mLogTP.glmm1<-lmer(AUC.Days22Thru43.LogTP~Nutrients*Mixing*Spores+(1|Raft),data=auc_data2)
summary(mLogTP.glmm1)
Anova(mLogTP.glmm1)

mLogTP.glmm2<-lmer(AUC.Days22Thru43.LogTP~Nutrients*Mixing+Nutrients*Spores+Mixing*Spores+(1|Raft),data=auc_data2)
Anova(mLogTP.glmm2)

mLogTP.glmm3<-lmer(AUC.Days22Thru43.LogTP~Nutrients*Mixing+Mixing*Spores+(1|Raft),data=auc_data2)
Anova(mLogTP.glmm3)

mLogTP.glmm4<-lmer(AUC.Days22Thru43.LogTP~Nutrients+Mixing*Spores+(1|Raft),data=auc_data2)
Anova(mLogTP.glmm4)

mLogTP.glmm5<-lmer(AUC.Days22Thru43.LogTP~Nutrients+Mixing+Spores+(1|Raft),data=auc_data2)
Anova(mLogTP.glmm5)

##################################################################################################
## Log(Total nitrogen)
##################################################################################################
LogTNSummary<-ddply(mydata2, c("Nutrients", "Mixing","Spores","Day"), summarise,
                    N    = length(LogTN),
                    mean = mean(LogTN),
                    sd   = sd(LogTN),
                    se   = sd / sqrt(N)
)

LogTN_Open<-subset(LogTNSummary, Day<22|Day>43)
LogTN_Closed<-subset(LogTNSummary,Day>=22&Day<=43)
LogTN_OpenMixed<-subset(LogTN_Open,Mixing=="Mix")
LogTN_OpenNoMixed<-subset(LogTN_Open,Mixing=="No")
LogTN_ClosedMixed<-subset(LogTN_Closed,Mixing=="Mix")
LogTN_ClosedNoMixed<-subset(LogTN_Closed,Mixing=="No")

PS5<-(ggplot()+
  geom_point(data=LogTN_OpenMixed,aes(x=Day+0.4,y=mean,color=Nutrients,shape=Mixing),size=2)+
  geom_point(data=LogTN_ClosedMixed,aes(x=Day+0.4,y=mean,color=Nutrients,fill=Nutrients,shape=Mixing),size=2)+
  geom_errorbar(data=LogTN_OpenMixed,aes(x=Day+0.4, ymin=mean-se, ymax=mean+se,color=Nutrients,linetype=Mixing), width=1) + 
  geom_errorbar(data=LogTN_ClosedMixed,aes(x=Day+0.4, ymin=mean-se, ymax=mean+se,color=Nutrients,linetype=Mixing), width=1) + 
  geom_point(data=LogTN_OpenNoMixed,aes(x=Day-0.4,y=mean,color=Nutrients,shape=Mixing),size=2)+
  geom_point(data=LogTN_ClosedNoMixed,aes(x=Day-0.4,y=mean,color=Nutrients,fill=Nutrients,shape=Mixing),size=2)+
  geom_errorbar(data=LogTN_OpenNoMixed,aes(x=Day-0.4, ymin=mean-se, ymax=mean+se,color=Nutrients,linetype=Mixing), width=1) + 
  geom_errorbar(data=LogTN_ClosedNoMixed,aes(x=Day-0.4, ymin=mean-se, ymax=mean+se,color=Nutrients,linetype=Mixing), width=1) + 
  scale_shape_manual(values=c(21,24))+
  scale_color_manual(values=NutColors)+
  scale_fill_manual(values=NutColors)+
  facet_wrap(~Spores,labeller=labeller(Spores=SporesNames))+
  ylab("Log(total nitrogen)")+
  xlab("Day of experiment")+
  fav_theme+
  annotate("rect", xmin = 22-0.8, xmax = 43+0.8, ymin = 5.85, ymax = 8,fill="gray90",alpha=0.3)
)

PS6<-(ggplot(auc_data2)+
  geom_boxplot(aes(x=Mixing,y=AUC.Days22Thru43.LogTN,color=Nutrients))+
  geom_point(aes(x=Mixing,y=AUC.Days22Thru43.LogTN,color=Nutrients,fill=Nutrients,shape=Mixing),position=position_dodge(width=0.75),size=2)+
  scale_shape_manual(values=c(21,24))+
  scale_color_manual(values=NutColors)+
  scale_fill_manual(values=NutColors)+
  facet_wrap(~Spores,labeller=labeller(Spores=SporesNames))+
  ylab("AUC, log(total nitrogen)")+
  fav_theme)

mLogTN.glmm1<-lmer(AUC.Days22Thru43.LogTN~Nutrients*Mixing*Spores+(1|Raft),data=auc_data2)
summary(mLogTN.glmm1)
Anova(mLogTN.glmm1)

mLogTN.glmm2<-lmer(AUC.Days22Thru43.LogTN~Nutrients*Mixing+Nutrients*Spores+Mixing*Spores+(1|Raft),data=auc_data2)
Anova(mLogTN.glmm2)

mLogTN.glmm3<-lmer(AUC.Days22Thru43.LogTN~Nutrients*Mixing+Nutrients*Spores+(1|Raft),data=auc_data2)
Anova(mLogTN.glmm3)

mLogTN.glmm4<-lmer(AUC.Days22Thru43.LogTN~Mixing+Nutrients*Spores+(1|Raft),data=auc_data2)
Anova(mLogTN.glmm4)

mLogTN.glmm5<-lmer(AUC.Days22Thru43.LogTN~Nutrients+Mixing+Spores+(1|Raft),data=auc_data2)
Anova(mLogTN.glmm5)
summary(mLogTN.glmm5)

mLogTN.glm1<-glm(AUC.Days22Thru43.LogTN~Nutrients+Mixing+Spores,data=auc_data2)
Anova(mLogTN.glm1)
summary(mLogTN.glm1)

mLogTN.glm2<-glm(AUC.Days22Thru43.LogTN~Nutrients+Mixing+Spores,family="quasipoisson",data=auc_data2)
Anova(mLogTN.glm2)
summary(mLogTN.glm2)

##################################################################################################
## Complete Fig. S1 
##################################################################################################
FigS1<-plot_grid(PS1+theme(legend.position="none"),
                PS2+theme(legend.position="none"),
                PS3+theme(legend.position="none"),
                PS4+theme(legend.position="none"),
                PS5+theme(legend.position="none"),
                PS6+theme(legend.position="none"), 
                labels = c('a','b','c','d','e','f'),
                ncol=2,rel_widths = c(1,0.5,1,0.5,1,0.5))

print(FigS1)
ggsave(here("figures/FigS1.tiff"), plot = FigS1, dpi = 300, width = 8, height = 8.5, units = "in", compression="lzw")

##################################################################################################
## END OF SCRIPT 
##################################################################################################