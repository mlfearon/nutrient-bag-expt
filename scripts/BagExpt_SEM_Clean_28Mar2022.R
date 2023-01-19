## Code to produce the path models included in "Nutrient enrichment, habitat structure, and disease in the plankton"
# Penczykowski et al. submitted to Oecologia


# code written by: Michelle L Fearon
# last updated: Jan 19, 2023


# load packages
library(piecewiseSEM)
library(ggplot2)
library(ggeffects)
library(nlme)
library(lme4)
library(glmmTMB)
library(lmerTest)
library(dplyr)
library(tidyr)
library(car)
library(here)


# set the path to the script relative to the project root directory
here::i_am("scripts/BagExpt_SEM_Clean_28Mar2022.R")

# load data
mydata <- read.csv(here("data/BagExpt_clean.csv"), stringsAsFactors = F, header = T)
str(mydata)
#View(mydata)

# scale and center key variables
# mydata$LogInfDens_z <- as.numeric(scale(mydata$LogInfDens))
# mydata$LogTotDens_z <- as.numeric(scale(mydata$LogTotDens))
# mydata$LogTotChl_z <- as.numeric(scale(mydata$LogTotChl))
# mydata$LogEdChl_z <- as.numeric(scale(mydata$LogEdChl))

# create dummy variables for categorical variables
mydata <- mutate(mydata, SporeTrt = ifelse(Spores == "No", 0, 1), MixingTrt = ifelse(Mixing == "No", 0, 1), 
                 NutrientTrt = ifelse(Nutrients == "Low", 0, 1))

mydata$NDay_fac <- as.factor(mydata$NDay)

# remove NAs
mydata2 <- filter(mydata, !is.na(UA))
mydata2 <- filter(mydata2, !is.na(TotChl))
mydata2 <- filter(mydata2, !is.na(LogTotDens))
is.na(mydata2$TotChl)


# filter data to only include day 22 to 43 of experiment (second epidemic)
mydata22_43 <- mydata2 %>%
  filter(Day >= 22 & Day <= 43)


# filter data for +Spores and -Spores treatments
mydata22_43_spores <- mydata22_43 %>%
  filter(Spores == "Spores")
mydata22_43_NOspores <- mydata22_43 %>%
  filter(Spores == "No")



# log and scale TP and TN, remove NA in nutrient data
dim(mydata22_43)
mydata22_43_nutrients <- mydata22_43 %>%
  filter(!is.na(TP)) %>%
  filter(!is.na(TN))
dim(mydata22_43_nutrients)
hist(log(mydata22_43_nutrients$TP))
hist(log(mydata22_43_nutrients$TN))
mydata22_43_nutrients$LogTP <- as.numeric(log(mydata22_43_nutrients$TP))
mydata22_43_nutrients$LogTN <- as.numeric(log(mydata22_43_nutrients$TN))
mydata22_43_nutrients$LogTP_z <- as.numeric(scale(log(mydata22_43_nutrients$TP)))
mydata22_43_nutrients$LogTN_z <- as.numeric(scale(log(mydata22_43_nutrients$TN)))
str(mydata22_43_nutrients)




# filter data for +Spores and -Spores treatments
mydata22_43_nutrients_spores <- mydata22_43_nutrients %>%
  filter(Spores == "Spores") %>%
  mutate(OLRE = 1:79)
mydata22_43_nutrients_NOspores <- mydata22_43_nutrients %>%
  filter(Spores == "No") %>%
  mutate(OLRE = 1:88)

#View(mydata22_43_nutrients_spores)

# overdispersion function by ben bolker
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}



### Simple SEM with data from day 22 to 43 (second epidemic)
# Run separate SEMs for +Spore and -Spore data because the infection densities are all zero in the -Spore treatments

###### +Spore model

# update columns in data set with re-scaled variables
mydata22_43_nutrients_spores$TotalDensity2 <- round(mydata22_43_nutrients_spores$TotalDensity/10000)
mydata22_43_nutrients_spores$InfDensity2 <- round(mydata22_43_nutrients_spores$InfDensity/1000)
mydata22_43_nutrients_spores$EdChl2 <- round(mydata22_43_nutrients_spores$EdChl)
mydata22_43_nutrients_spores$TN2 <- (mydata22_43_nutrients_spores$TN/100)
mydata22_43_nutrients_spores$TP2 <- (mydata22_43_nutrients_spores$TP/10)
head(mydata22_43_nutrients_spores)


# Initial model for + spores
exp.sem.fit_spores <- psem(
  glmer(InfDensity2~TotalDensity2+MixingTrt+EdChl2+(1|Bag)+(1|NDay_fac)+(1|OLRE),control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),family="poisson",data=mydata22_43_nutrients_spores),
  glmer(TotalDensity2~EdChl2+MixingTrt+(1|Bag)+(1|NDay_fac)+(1|OLRE), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),family="poisson",data=mydata22_43_nutrients_spores),
  glmer(EdChl2~TN2+TP2+MixingTrt+(1|Bag)+(1|NDay_fac)+(1|OLRE),family="poisson",data=mydata22_43_nutrients_spores),
  TotalDensity2 %~~% TP2,
  TN2%~~%TP2,
  mydata22_43_nutrients_spores
)
summary(exp.sem.fit_spores) # seems like there is a collinearity issue with including both TN and TP, causing effects in opposite directions but highly correlated.


cor.test(mydata22_43_nutrients_spores$TN, mydata22_43_nutrients_spores$TP) # Yes, Pearson's correlation is 0.809


# Run separate SEMs testing TP and TN independently


### BEST + SPORE and TP only model
# Need to add add a correlated error between TP and total density to fix the collinearity problem!! 
exp.sem.fit_spores_TP <- psem(
  glmer(InfDensity2~TotalDensity2+MixingTrt+EdChl2+(1|Bag)+(1|NDay_fac)+(1|OLRE),control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),family="poisson",data=mydata22_43_nutrients_spores),
  glmer(TotalDensity2~EdChl2+MixingTrt+(1|Bag)+(1|NDay_fac)+(1|OLRE), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),family="poisson",data=mydata22_43_nutrients_spores),
  glmer(EdChl2~TP2+MixingTrt+(1|Bag)+(1|NDay_fac)+(1|OLRE),family="poisson",data=mydata22_43_nutrients_spores),
  TotalDensity2 %~~% TP2,
  mydata22_43_nutrients_spores
)
summary(exp.sem.fit_spores_TP)

# Table S1 in Appendix
# coefficient table for +Spore path model with only TP (coefficients, st error, critical value, p-value and standard coefficients)
std_scale_coefs_sporesTP <- stdCoefs(exp.sem.fit_spores_TP, data=mydata22_43_nutrients_spores, 
                                         standardize="scale",
                                         standardize.type = "latent.linear",
                                         intercepts = F)
write.csv(std_scale_coefs_sporesTP, here("tables/sem_std_scale_coefs_SporesTP_model.csv"), row.names=FALSE)




### BEST +SPORE and TN only model
exp.sem.fit_spores_TN <- psem(
  glmer(InfDensity2~TotalDensity2+MixingTrt+EdChl2+(1|Bag)+(1|NDay_fac)+(1|OLRE),control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),family="poisson",data=mydata22_43_nutrients_spores),
  glmer(TotalDensity2~EdChl2+MixingTrt+(1|Bag)+(1|NDay_fac)+(1|OLRE), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),family="poisson",data=mydata22_43_nutrients_spores),
  glmer(EdChl2~TN2+MixingTrt+(1|Bag)+(1|NDay_fac)+(1|OLRE),family="poisson",data=mydata22_43_nutrients_spores),
  mydata22_43_nutrients_spores
)
summary(exp.sem.fit_spores_TN) #model converges and fits the data, does not have a link between TN and TotalDensity


# Table S3 in Appendix
# coefficient table for +Spore path model with only TP (coefficients, st error, critical value, p-value and standard coefficients)
std_scale_coefs_sporesTN <- stdCoefs(exp.sem.fit_spores_TN, data=mydata22_43_nutrients_spores, 
                                     standardize="scale",
                                     standardize.type = "latent.linear",
                                     intercepts = F)
write.csv(std_scale_coefs_sporesTN, here("tables/sem_std_scale_coefs_SporesTN_model.csv"), row.names=FALSE)





#Test of +Spores with TP model, removing Infection density (epidemic size) to compare power with the best model above
exp.sem.fit_spores_TP_NoInf <- psem(
  glmer(TotalDensity2~EdChl2+MixingTrt+(1|Bag)+(1|NDay_fac)+(1|OLRE), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),family="poisson",data=mydata22_43_nutrients_spores),
  glmer(EdChl2~TP2+MixingTrt+(1|Bag)+(1|NDay_fac)+(1|OLRE),family="poisson",data=mydata22_43_nutrients_spores),
  TotalDensity2 %~~% TP2,
  mydata22_43_nutrients_spores
)
summary(exp.sem.fit_spores_TP_NoInf)

# huh, the model without infection density has a much lower AIC, I'm guessing because it is simpler and has fewer relationships...
# but I think including epidemic is important to explain the dynamics happening, because we know adding spores will change population dynamics and the model can't know that without us telling it
AIC(exp.sem.fit_spores_TP, exp.sem.fit_spores_TP_NoInf) 



#Test of +Spores with TN model, removing Infection density (epidemic size) to compare power with the best model above
exp.sem.fit_spores_TN_NoInf <- psem(
  glmer(TotalDensity2~EdChl2+MixingTrt+(1|Bag)+(1|NDay_fac)+(1|OLRE), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),family="poisson",data=mydata22_43_nutrients_spores),
  glmer(EdChl2~TN2+MixingTrt+(1|Bag)+(1|NDay_fac)+(1|OLRE),family="poisson",data=mydata22_43_nutrients_spores),
  mydata22_43_nutrients_spores
)
summary(exp.sem.fit_spores_TN_NoInf) 

AIC(exp.sem.fit_spores_TN, exp.sem.fit_spores_TN_NoInf)



#Test of +Spores with TP model, replacing Infection density with infection prevalence to evaluate epidemic size and reduce collinearity between infection density and total host density
exp.sem.fit_spores_TP_Prev <- psem(
  glmer(PrevTotal~TotalDensity2+MixingTrt+EdChl2+(1|Bag)+(1|NDay_fac)+(1|OLRE), weights = TotalCount, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),family="binomial", data=mydata22_43_nutrients_spores),
  glmer(TotalDensity2~EdChl2+MixingTrt+(1|Bag)+(1|NDay_fac)+(1|OLRE), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),family="poisson",data=mydata22_43_nutrients_spores),
  glmer(EdChl2~TP2+MixingTrt+(1|Bag)+(1|NDay_fac)+(1|OLRE),family="poisson",data=mydata22_43_nutrients_spores),
  PrevTotal %~~% TP2,
  TotalDensity2 %~~% TP2,
  mydata22_43_nutrients_spores
)
summary(exp.sem.fit_spores_TP_Prev)

cor.test(mydata22_43_nutrients_spores$PrevTotal, mydata22_43_nutrients_spores$InfDensity)

# comparison of different model versions and tests of assumptions
mod_test <- glmer(PrevTotal~TotalDensity2+MixingTrt+EdChl2+(1|Bag)+(1|NDay_fac)+(1|OLRE), weights = TotalCount, family="binomial", data=mydata22_43_nutrients_spores)
mod_test <- glmer(cbind(InfCount,UninfCount)~TotalDensity2+MixingTrt+EdChl2+(1|Bag)+(1|NDay_fac)+(1|OLRE), family="binomial", data=mydata22_43_nutrients_spores)
summary(mod_test)
overdisp_fun(mod_test)
plot(mod_test)
qqnorm(resid(mod_test))
qqline(resid(mod_test))

library(ggeffects)
predict_test <- ggpredict(mod_test, c("TotalDensity2 [all]"))
plot(predict_test, add.data = T)

AIC(exp.sem.fit_spores_TP, exp.sem.fit_spores_TP_Prev)


#Test of +Spores with TN model, replacing Infection density with infection prevalence to evaluate epidemic size and reduce collinearity between infection density and total host density
exp.sem.fit_spores_TN_Prev <- psem(
  glmer(PrevTotal~TotalDensity2+MixingTrt+EdChl2+(1|Bag)+(1|NDay_fac)+(1|OLRE),weights = TotalCount,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),family="binomial",data=mydata22_43_nutrients_spores),
  glmer(TotalDensity2~EdChl2+MixingTrt+(1|Bag)+(1|NDay_fac)+(1|OLRE), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),family="poisson",data=mydata22_43_nutrients_spores),
  glmer(EdChl2~TN2+MixingTrt+(1|Bag)+(1|NDay_fac)+(1|OLRE),family="poisson",data=mydata22_43_nutrients_spores),
  mydata22_43_nutrients_spores
)
summary(exp.sem.fit_spores_TN_Prev) 


mod_test2 <- glmer(PrevTotal~TotalDensity2+MixingTrt+EdChl2+(1|Bag)+(1|NDay_fac)+(1|OLRE),weights = TotalCount,family="binomial",data=mydata22_43_nutrients_spores)
summary(mod_test2)
overdisp_fun(mod_test2)
plot(mod_test2)
qqnorm(resid(mod_test2))
qqline(resid(mod_test2))

AIC(exp.sem.fit_spores_TN, exp.sem.fit_spores_TN_Prev)





###### -Spore models

#update columns with calculations
mydata22_43_nutrients_NOspores$TotalDensity2 <- round(mydata22_43_nutrients_NOspores$TotalDensity/10000)
mydata22_43_nutrients_NOspores$InfDensity2 <- round(mydata22_43_nutrients_NOspores$InfDensity/1000)
mydata22_43_nutrients_NOspores$EdChl2 <- round(mydata22_43_nutrients_NOspores$EdChl)
mydata22_43_nutrients_NOspores$TN2 <- (mydata22_43_nutrients_NOspores$TN/100)
mydata22_43_nutrients_NOspores$TP2 <- (mydata22_43_nutrients_NOspores$TP/10)
head(mydata22_43_nutrients_NOspores)

# similar as above, we will run separate No spore models for TP and TN because they are too highly correlated to keep within the same model


### BEST NO SPORE TP MODEL
# -spore model, only TP
exp.sem.fit_NOspores_TP <- psem(
  glmer(TotalDensity2~EdChl2+MixingTrt+(1|Bag)+(1|NDay_fac)+(1|OLRE), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),family="poisson",data=mydata22_43_nutrients_NOspores),
  glmer(EdChl2~TP2+MixingTrt+(1|Bag)+(1|NDay_fac)+(1|OLRE),family="poisson",data=mydata22_43_nutrients_NOspores),
  mydata22_43_nutrients_NOspores
)
summary(exp.sem.fit_NOspores_TP)

# Table S2 in Appendix
# coefficient table for -Spore path model with only TP (coefficients, st error, critical value, p-value and standard coefficients)
std_scale_coefs_NOsporesTP <- stdCoefs(exp.sem.fit_NOspores_TP, data=mydata22_43_nutrients_NOspores, 
                                     standardize="scale",
                                     standardize.type = "latent.linear",
                                     intercepts = F)
write.csv(std_scale_coefs_NOsporesTP, here("tables/sem_std_scale_coefs_NOSporesTP_model.csv"), row.names=FALSE)


### BEST NO SPORE TN MODEL
# needed to add a correlated error between TN and total density
exp.sem.fit_NOspores_TN <- psem(
  glmer(TotalDensity2~EdChl2+MixingTrt+(1|Bag)+(1|NDay_fac)+(1|OLRE), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),family="poisson",data=mydata22_43_nutrients_NOspores),
  glmer(EdChl2~TN2+MixingTrt+(1|Bag)+(1|NDay_fac)+(1|OLRE),family="poisson",data=mydata22_43_nutrients_NOspores),
  TotalDensity2 %~~% TN2,
  mydata22_43_nutrients_NOspores
)
summary(exp.sem.fit_NOspores_TN)



# Table S4 in Appendix
# coefficient table for -Spore path model with only TN (coefficients, st error, critical value, p-value and standard coefficients)
std_scale_coefs_NOsporesTN <- stdCoefs(exp.sem.fit_NOspores_TN, data=mydata22_43_nutrients_NOspores, 
                                       standardize="scale",
                                       standardize.type = "latent.linear",
                                       intercepts = F)
write.csv(std_scale_coefs_NOsporesTN, here("tables/sem_std_scale_coefs_NOSporesTN_model.csv"), row.names=FALSE)


# Testing the NO Spore TN model for goodness of fit by removing a non-sig pathway from the model so that it is not fully saturdated
exp.sem.fit_NOspores_TN_test <- psem(
  glmer(TotalDensity2~EdChl2+(1|Bag)+(1|NDay_fac)+(1|OLRE), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),family="poisson",data=mydata22_43_nutrients_NOspores),
  glmer(EdChl2~TN2+MixingTrt+(1|Bag)+(1|NDay_fac)+(1|OLRE),family="poisson",data=mydata22_43_nutrients_NOspores),
  TotalDensity2 %~~% TN2,
  mydata22_43_nutrients_NOspores
)
summary(exp.sem.fit_NOspores_TN_test)




# correlations between key variables, some are quite high!
cor(mydata22_43_nutrients_NOspores$TN2, mydata22_43_nutrients_NOspores$TP2)
cor(mydata22_43_nutrients_NOspores$EdChl2, mydata22_43_nutrients_NOspores$TP2)
cor(mydata22_43_nutrients_NOspores$EdChl2, mydata22_43_nutrients_NOspores$TN2)
cor(mydata22_43_nutrients_NOspores$EdChl2, mydata22_43_nutrients_NOspores$TotalDensity2)
cor(mydata22_43_nutrients_NOspores$TP2, mydata22_43_nutrients_NOspores$TotalDensity2)
cor(mydata22_43_nutrients_NOspores$TN2, mydata22_43_nutrients_NOspores$TotalDensity2)






####################################################
# Check of component models and assumptions
####################################################
# Raft is often singular, removed from  all models

### Infection density models

infdens <- glmer(InfDensity2~TotalDensity2+MixingTrt+EdChl2+(1|Bag)+(1|NDay_fac)+(1|OLRE),control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),family="poisson",data=mydata22_43_nutrients_spores)
summary(infdens)
plot(infdens)
qqnorm(resid(infdens)) 
qqline(resid(infdens))
overdisp_fun(infdens) ## added Observation level random effect because original model was over-dispersed. No longer over-dispersed with OLRE included.




### Host density models

hostdens_spores <- glmer(TotalDensity2~EdChl2+MixingTrt+(1|Bag)+(1|NDay_fac)+(1|OLRE), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),family="poisson",data=mydata22_43_nutrients_spores)
summary(hostdens_spores)
plot(hostdens_spores)
qqnorm(resid(hostdens_spores)) 
qqline(resid(hostdens_spores))
overdisp_fun(hostdens_spores) ### added Observation level random effect because original model was over-dispersed. No longer over-dispersed with OLRE included..

hostdens_NOspores <- glmer(TotalDensity2~EdChl2+MixingTrt+(1|Bag)+(1|NDay_fac)+(1|OLRE), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),family="poisson",data=mydata22_43_nutrients_NOspores)
summary(hostdens_NOspores)
plot(hostdens_NOspores)
qqnorm(resid(hostdens_NOspores)) 
qqline(resid(hostdens_NOspores))
overdisp_fun(hostdens_NOspores)  ## added Observation level random effect because original model was over-dispersed. No longer over-dispersed with OLRE included.



### Edible Chl algae models

chl_sporesTP <- glmer(EdChl2~TP2+MixingTrt+(1|Bag)+(1|NDay_fac)+(1|OLRE),family="poisson",data=mydata22_43_nutrients_spores)
summary(chl_sporesTP)
plot(chl_sporesTP)
qqnorm(resid(chl_sporesTP)) 
qqline(resid(chl_sporesTP))
overdisp_fun(chl_sporesTP)  # no significant overdispersion, but added OLRE to the model to be consistent with the other versions of the model


chl_sporesTN <- glmer(EdChl2~TN2+MixingTrt+(1|Bag)+(1|NDay_fac)+(1|OLRE),family="poisson",data=mydata22_43_nutrients_spores)
summary(chl_sporesTN)
plot(chl_sporesTN)
qqnorm(resid(chl_sporesTN)) 
qqline(resid(chl_sporesTN))
overdisp_fun(chl_sporesTN) # no significant overdispersion, but added OLRE to the model to be consistent with the other versions of the model


chl_NOsporesTP <- glmer(EdChl2~TP2+MixingTrt+(1|Bag)+(1|NDay_fac)+(1|OLRE),family="poisson",data=mydata22_43_nutrients_NOspores)
summary(chl_NOsporesTP)
plot(chl_NOsporesTP)
qqnorm(resid(chl_NOsporesTP)) 
qqline(resid(chl_NOsporesTP))
overdisp_fun(chl_NOsporesTP)  ## added Observation level random effect because original model was over-dispersed. No longer over-dispersed with OLRE included.

chl_NOsporesTN <- glmer(EdChl2~TN2+MixingTrt+(1|Bag)+(1|NDay_fac)+(1|OLRE),family="poisson",data=mydata22_43_nutrients_NOspores)
summary(chl_NOsporesTN)
plot(chl_NOsporesTN)
qqnorm(resid(chl_NOsporesTN)) 
qqline(resid(chl_NOsporesTN))
overdisp_fun(chl_NOsporesTN)  ## added Observation level random effect because original model was over-dispersed. No longer over-dispersed with OLRE included.





# APPENDIX S1 Table S5
# correlation table among factors in the SEM for +Spore data
Spore_subset <- mydata22_43_nutrients_spores %>%
  select(InfDensity2, TotalDensity2, EdChl2, TP2, TN2)

Spore_correlations <- cor(Spore_subset)

# add p-values for each correlation above the diagonal
Spore_correlations["InfDensity2", "TotalDensity2"] <- unlist(cor.test(Spore_subset$InfDensity2, Spore_subset$TotalDensity2)[3])
Spore_correlations["InfDensity2", "EdChl2"] <- unlist(cor.test(Spore_subset$InfDensity2, Spore_subset$EdChl2)[3])
Spore_correlations["InfDensity2", "TP2"] <- unlist(cor.test(Spore_subset$InfDensity2, Spore_subset$TP2)[3])
Spore_correlations["InfDensity2", "TN2"] <- unlist(cor.test(Spore_subset$InfDensity2, Spore_subset$TN2)[3])
Spore_correlations["TotalDensity2", "EdChl2"] <- unlist(cor.test(Spore_subset$TotalDensity2, Spore_subset$EdChl2)[3])
Spore_correlations["TotalDensity2", "TP2"] <- unlist(cor.test(Spore_subset$TotalDensity2, Spore_subset$TP2)[3])
Spore_correlations["TotalDensity2", "TN2"] <- unlist(cor.test(Spore_subset$TotalDensity2, Spore_subset$TN2)[3])
Spore_correlations["EdChl2", "TP2"] <- unlist(cor.test(Spore_subset$EdChl2, Spore_subset$TP2)[3])
Spore_correlations["EdChl2", "TN2"] <- unlist(cor.test(Spore_subset$EdChl2, Spore_subset$TN2)[3])
Spore_correlations["TP2", "TN2"] <- unlist(cor.test(Spore_subset$TP2, Spore_subset$TN2)[3])


Spore_correlations <- format(round(Spore_correlations, 3))
# set diagonal to NA
diag(Spore_correlations) <- "--"
rownames(Spore_correlations) <- c("Infection Density", "Total Host Density", "Edible Chlorophyll", "Total Phosphorus", "Total Nitrogen")
colnames(Spore_correlations) <- c("Infection Density", "Total Host Density", "Edible Chlorophyll", "Total Phosphorus", "Total Nitrogen")
write.csv(Spore_correlations, "Spore_SEM_factor_correlations.csv", quote = F, row.names = T)



# APPENDIX S1 Table S6
# correlation table among factors in the SEM for -Spore data
NoSpore_subset <- mydata22_43_nutrients_NOspores %>%
  select(TotalDensity2, EdChl2, TP2, TN2)

NoSpore_correlations <- cor(NoSpore_subset)

# add p-values for each correlation above the diagonal
NoSpore_correlations["TotalDensity2", "EdChl2"] <- unlist(cor.test(NoSpore_subset$TotalDensity2, NoSpore_subset$EdChl2)[3])
NoSpore_correlations["TotalDensity2", "TP2"] <- unlist(cor.test(NoSpore_subset$TotalDensity2, NoSpore_subset$TP2)[3])
NoSpore_correlations["TotalDensity2", "TN2"] <- unlist(cor.test(NoSpore_subset$TotalDensity2, NoSpore_subset$TN2)[3])
NoSpore_correlations["EdChl2", "TP2"] <- unlist(cor.test(NoSpore_subset$EdChl2, NoSpore_subset$TP2)[3])
NoSpore_correlations["EdChl2", "TN2"] <- unlist(cor.test(NoSpore_subset$EdChl2, NoSpore_subset$TN2)[3])
NoSpore_correlations["TP2", "TN2"] <- unlist(cor.test(NoSpore_subset$TP2, NoSpore_subset$TN2)[3])


NoSpore_correlations <- format(round(NoSpore_correlations, 3))
# set diagonal to NA
diag(NoSpore_correlations) <- "--"
rownames(NoSpore_correlations) <- c("Total Host Density", "Edible Chlorophyll", "Total Phosphorus", "Total Nitrogen")
colnames(NoSpore_correlations) <- c("Total Host Density", "Edible Chlorophyll", "Total Phosphorus", "Total Nitrogen")
write.csv(NoSpore_correlations, "NoSpore_SEM_factor_correlations.csv", quote = F, row.names = T)

