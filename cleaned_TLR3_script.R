# Charli Davies

# R code and analysis for: Contemporary evolution of the viral-sensing TLR3 gene in an isolated vertebrate population

# last modified 04/12/2020

# using R version 3.6.1

# Need to first put all files from dryad into the same folder and set this as you're working directory

# clear R of all objects
rm(list=ls())

# load libraries  ####
library(tidyverse) 
library(dplyr) 
library(plyr) 
library(grid) 
library(gridExtra) 
library(ggplot2) 

library(survival) 
library(survminer) 
library(rms) 
library(coxme) 

library(glmmTMB) 
library(MuMIn) 
library(DHARMa) 
library(sjPlot) 
library(arm) 


######### Q1: Temporal patterns of TLR3 variation on Cousin############################################################################################

# Load data####
GENO<-read_csv("Genotype_frequencies.csv")
head(GENO)
View(GENO)
str(GENO)

#Tidy data and convert variables to factors etc
GENO$Genotype <- as.factor(GENO$Genotype)
GENO$Allele <- as.factor(GENO$Allele)
GENO$Year <- as.numeric(GENO$year)
GENO$Adult_frequency <- as.numeric(GENO$Adult_frequency )
GENO$YoB_frequency <- as.numeric(GENO$YoB_frequency )
GENO$Adult_population <- as.numeric(GENO$Adult_population )
GENO$Adult_allele_frequency <- as.numeric(GENO$Adult_allele_frequency )
GENO$YoB_allele_frequency <- as.numeric(GENO$YoB_allele_frequency )

# Do allele frequencies change over time in the adult population?####

Allele_C<- GENO %>% filter(Allele == "C")%>% droplevels()%>%filter(Year != "1992") %>% droplevels()
mod_Allele_C <-lm(Adult_allele_frequency~Year, data=Allele_C)
summary(mod_Allele_C)

Allele_A<- GENO %>% filter(Allele == "A")%>% droplevels()%>%filter(Year != "1992") %>% droplevels()
mod_Allele_A <-lm(Adult_allele_frequency~Year, data=Allele_A)
summary(mod_Allele_A)

# Do allele frequencies change over time in the juvenile population?####

Hatch_C<- GENO %>% filter(Allele == "C")%>% droplevels()%>%filter(Year != "1992") %>% droplevels()%>%filter(YoB_allele_frequency != "NA") %>% droplevels()
view(Hatch_C)
mod_HATCH_C <-lm(YoB_allele_frequency~Year, data=Hatch_C)
summary(mod_HATCH_C)

Hatch_A<- GENO %>% filter(Allele == "A")%>% droplevels()%>%filter(Year != "1992") %>% droplevels()%>%filter(YoB_allele_frequency != "NA") %>% droplevels()
view(Hatch_A)
mod_HATCH_A <-lm(YoB_allele_frequency~Year, data=Hatch_A)
summary(mod_HATCH_A)

# Plot Figure 1 ####

allele_graph<-
  GENO %>%   filter(Year != "1992") %>% droplevels()%>%filter(Allele != "na") %>% droplevels()%>%
  ggplot(aes(x = Year, y = Adult_allele_frequency, group=Allele, colour=Allele))+
  geom_point() + 
  geom_smooth(method='lm', se= FALSE, size=1.5)+
  scale_color_manual(values = c("#034a1b", "#ffca47"))+
  scale_linetype_manual(values = c("longdash","solid", "solid", "solid"))+
  scale_x_continuous(breaks=c(1993, 1998, 2003, 2008, 2013, 2018))+
  labs(y = "Allele frequency", x = "Year") +
  theme_classic()+theme(axis.title.x=element_text(size=20),axis.text.x =element_text(size=16))+theme(axis.title.y=element_text(size=20),axis.text.y=element_text(size=16), legend.text = element_text(size=16), legend.title = element_text(size=16))+ 
  theme(legend.position=c(0.2, 0.2))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)))+
  geom_area(aes( x=Year,y = Adult_population),fill= "grey", colour="grey", alpha=0.1)+
  geom_line(aes(x = Year, y = YoB_allele_frequency, group=Allele,colour=Allele), linetype = "dashed", size=1)+ 
  scale_y_continuous(sec.axis = dup_axis(~ . * 100,name = "Adult population sampled (%)"))+
  theme(axis.title.y.right = element_text(margin = margin(t = 0, r = 0, b = 0, l = 5)))+
  theme(axis.title.y.right = element_text(angle = 90))

allele_graph

combo_plot<-grid.arrange(allele_graph, ncol = 1)
ggsave("Figure_1.png", combo_plot, width = 10, height=6, dpi = 300)
ggsave("Figure_1.pdf", combo_plot, width = 10, height=6, dpi = 300)


######### Q2:Testing for contemporary selection on TLR3 variation on Cousin  ##########################################################################


# import data file####
ALLTLR3<-read_csv("TLR3_CLEANED.csv")
head(ALLTLR3)
View(ALLTLR3)
str(ALLTLR3)

#Tidy data and convert variables to factors etc
ALLTLR3$SexEstimate <- as.factor(ALLTLR3$SexEstimate)
ALLTLR3$HatchYearFactor <- as.factor(ALLTLR3$HatchYear)
ALLTLR3$CohortYear <- as.factor(ALLTLR3$CohortYear)
ALLTLR3$SeasonBorn <- as.factor(ALLTLR3$SeasonBorn)
ALLTLR3$AdultAtFirstCapture <- as.factor(ALLTLR3$AdultAtFirstCapture)
ALLTLR3$Translocated <- as.factor(ALLTLR3$Translocated)
ALLTLR3$MinorBreedingSeason <- as.factor(ALLTLR3$MinorBreedingSeason)
ALLTLR3$AgeAtDeathCousinBiAnnualSeason <- as.numeric(ALLTLR3$AgeAtDeathCousinBiAnnualSeason)
ALLTLR3$SurviveToFirstYear_BiAnnual <- as.factor(ALLTLR3$SurviveToFirstYear_BiAnnual)
ALLTLR3$Survival_Status <- as.numeric(ALLTLR3$Survival_Status)
ALLTLR3$Survival_Status_WS2018 <- as.numeric(ALLTLR3$Survival_Status_WS2018)
ALLTLR3$CountOfSurvivingOffspringTo1Year <- as.numeric(ALLTLR3$CountOfSurvivingOffspringTo1Year)
ALLTLR3$CountOfSurvivingOffspringToOFL <- as.numeric(ALLTLR3$CountOfSurvivingOffspringToOFL)
ALLTLR3$TLR3 <- as.factor(ALLTLR3$TLR3)
ALLTLR3$A_allele_presence <- as.factor(ALLTLR3$A_allele_presence)
ALLTLR3$C_allele_presence <- as.factor(ALLTLR3$C_allele_presence)
ALLTLR3$Hs_obs = as.numeric(ALLTLR3$Hs_obs)
ALLTLR3$Maternal_Hs_obs <- as.numeric(ALLTLR3$Maternal_Hs_obs)
ALLTLR3$MHCDiversity <- as.numeric(ALLTLR3$MHCDiversity)
ALLTLR3$Aseua4 <- as.factor(ALLTLR3$Aseua4)
str(ALLTLR3)

# to reorder levels  
ALLTLR3$CohortYear <- factor(ALLTLR3$CohortYear, levels(ALLTLR3$CohortYear)[c(2,1)])
View(ALLTLR3)
ALLTLR3$TLR3b <- factor(ALLTLR3$TLR3, levels(ALLTLR3$TLR3)[c(3, 1, 2)])
ALLTLR3$TLR3c <- factor(ALLTLR3$TLR3, levels(ALLTLR3$TLR3)[c(2, 1, 3)])
View(ALLTLR3)



# Q2A: Are there any survival differences####

#      Survival model using COXME - TLR3 reference level set as AA ####
survival1 <- coxme(Surv(AgeAtDeathCousinBiAnnualSeason,Survival_Status) ~ SexEstimate+TLR3 + Aseua4+MHCDiversity+Hs_obs+Maternal_Hs_obs+SeasonBorn+(1|HatchYearFactor),data = ALLTLR3)
survival1

# check no difference in results if exclude birds missing a minor season
ALLTLR3_nominor<- ALLTLR3 %>% filter(MinorBreedingSeason != "NOMINOR")%>% droplevels()
survival1_nominor <- coxme(Surv(AgeAtDeathCousinBiAnnualSeason,Survival_Status) ~ SexEstimate+TLR3 + Aseua4+MHCDiversity+Hs_obs+Maternal_Hs_obs+SeasonBorn+(1|HatchYearFactor),data = ALLTLR3_nominor)
survival1_nominor # results do not meaningfully change - use dataset with all birds included

# Same model but with TLR3 reference level set as CC
survival1b <- coxme(Surv(AgeAtDeathCousinBiAnnualSeason,Survival_Status) ~ SexEstimate+TLR3b + Aseua4+MHCDiversity+Hs_obs+Maternal_Hs_obs+SeasonBorn+(1|HatchYearFactor),data = ALLTLR3)
survival1b

# Same model but with TLR3 reference level set as AC
survival1c <- coxme(Surv(AgeAtDeathCousinBiAnnualSeason,Survival_Status) ~ SexEstimate+TLR3c + Aseua4+MHCDiversity+Hs_obs+Maternal_Hs_obs+SeasonBorn+(1|HatchYearFactor),data = ALLTLR3)
survival1c

#      Check assumptions - note this is just using coxph model####

# first run model not including random factor
surv_object <- Surv(time = ALLTLR3$AgeAtDeathCousinBiAnnualSeason, event = ALLTLR3$Survival_Status)
survival_coxph <- coxph(surv_object ~ TLR3 + SexEstimate+Aseua4+MHCDiversity+Maternal_Hs_obs+Hs_obs+SeasonBorn, data = ALLTLR3)
ggforest(survival_coxph, data = ALLTLR3)

# check assumptions
test.ph <- cox.zph(survival_coxph)
test.ph # test is not statistically significant for each of the covariates, and the global test is also not statistically significant. Therefore, we can assume the proportional hazards.
ggcoxzph(test.ph) # all looks fine
ggcoxdiagnostics(survival_coxph, type = "dfbeta", linear.predictions = FALSE, ggtheme = theme_bw())  #also all look fine

#test autoccrelation
cvif <- vif(survival_coxph)
cvif #looks fine

#      Plot Figure 2 ####
surv_object_plot <- Surv(time = ALLTLR3$AgeAtDeathCousinBiAnnualSeason, event = ALLTLR3$Survival_Status)

surv_plot_data <- survfit(surv_object_plot ~ TLR3 , data = ALLTLR3)

survival_graph <- ggsurvplot(surv_plot_data, data = ALLTLR3, pval = FALSE, pval.coord = c(1, 0.9), conf.int = TRUE, conf.int.alpha=0.1, break.time.by = 2, surv.median.line="hv",   risk.table = TRUE, tables.height = 0.2,tables.theme = theme_cleantable(),legend.title = "TLR3 genotype",legend.labs = c("AA", "AC", "CC"), legend = c(0.85, 0.8), palette = c("#034a1b", "#56b82e", "#ffca47"),linetype=c("solid", "dotted", "longdash"), xlab = "Years", ylab = "Survival probability", ggtheme = theme_classic()+theme(axis.title.x=element_text(size=20),axis.text.x =element_text(size=16))+theme(legend.key.width = grid::unit(1.4, "cm"))+theme(axis.title.y=element_text(size=20),axis.text.y=element_text(size=16),legend.text = element_text(size=12), legend.title = element_text(size=16)))

survival_graph$plot <- survival_graph$plot+ ggplot2::annotate("text",x = 0.7, y = 0.05, label = "1", size = 5, colour="#ffca47")
survival_graph$plot <- survival_graph$plot+ ggplot2::annotate("text",x = 1.7, y = 0.05, label = "2", size = 5, colour="#56b82e")
survival_graph$plot <- survival_graph$plot+ ggplot2::annotate("text",x = 3, y = 0.05, label = "2.5", size = 5, colour="#034a1b")
survival_graph

ggsave("Figure_2.png", survival_graph, width = 8.5, height=5.75, dpi = 300)
ggsave("Figure_2.pdf", survival_graph, width = 8.5, height=5.75, dpi = 300)


# Q2B: Are there any reproductive success differences?  ####

#      Create reduced dataset####

# seperate by sex
# exclude translocated birds
# exclude birds still alive after winter 2018 season
# for reproductive rate also want to exclude individuals who died before 1 year of age

# dataset for male LRS
ALLTLR3_male_LRS<- ALLTLR3 %>% filter(Survival_Status_WS2018 != "0")%>% droplevels()%>% filter(SexEstimate != "0")%>% droplevels()%>% filter(Translocated != "Yes")%>% droplevels()

# dataset for female LRS
ALLTLR3_female_LRS<- ALLTLR3 %>% filter(Survival_Status_WS2018 != "0")%>% droplevels()%>% filter(SexEstimate != "1")%>% droplevels()%>% filter(Translocated != "Yes")%>% droplevels()

# dataset for male reproductive rate
ALLTLR3_male_rate<- ALLTLR3 %>% filter(Survival_Status_WS2018 != "0")%>% droplevels()%>% filter(SexEstimate != "0")%>% droplevels()%>% filter(Translocated != "Yes")%>% droplevels() %>% filter(SurviveToFirstYear_BiAnnual != "N")%>% droplevels()

# dataset for female reproductive rate
ALLTLR3_female_rate<- ALLTLR3 %>% filter(Survival_Status_WS2018 != "0")%>% droplevels()%>% filter(SexEstimate != "1")%>% droplevels()%>% filter(Translocated != "Yes")%>% droplevels() %>% filter(SurviveToFirstYear_BiAnnual != "N")%>% droplevels()

#      Scale and center variables ####

ALLTLR3_male_rate$ScaledMHCDiversity<-rescale(ALLTLR3_male_rate$MHCDiversity, binary.inputs="full")
ALLTLR3_male_rate$ScaledHs_obs<-rescale(ALLTLR3_male_rate$Hs_obs, binary.inputs="full")
ALLTLR3_male_rate$ScaledAge<-rescale(ALLTLR3_male_rate$AgeAtDeathCousinBiAnnualSeason, binary.inputs="full")
ALLTLR3_female_rate$ScaledMHCDiversity<-rescale(ALLTLR3_female_rate$MHCDiversity, binary.inputs="full")
ALLTLR3_female_rate$ScaledHs_obs<-rescale(ALLTLR3_female_rate$Hs_obs, binary.inputs="full")
ALLTLR3_female_rate$ScaledAge<-rescale(ALLTLR3_female_rate$AgeAtDeathCousinBiAnnualSeason, binary.inputs="full")
ALLTLR3_male_LRS$ScaledMHCDiversity<-rescale(ALLTLR3_male_LRS$MHCDiversity, binary.inputs="full")
ALLTLR3_male_LRS$ScaledHs_obs<-rescale(ALLTLR3_male_LRS$Hs_obs, binary.inputs="full")
ALLTLR3_male_LRS$ScaledAge<-rescale(ALLTLR3_male_LRS$AgeAtDeathCousinBiAnnualSeason, binary.inputs="full")
ALLTLR3_female_LRS$ScaledMHCDiversity<-rescale(ALLTLR3_female_LRS$MHCDiversity, binary.inputs="full")
ALLTLR3_female_LRS$ScaledHs_obs<-rescale(ALLTLR3_female_LRS$Hs_obs, binary.inputs="full")
ALLTLR3_female_LRS$ScaledAge<-rescale(ALLTLR3_female_LRS$AgeAtDeathCousinBiAnnualSeason, binary.inputs="full")

# add in quadratic age as a variable
ALLTLR3_female_LRS$ScaledAge2 <- as.numeric((ALLTLR3_female_LRS$ScaledAge)^2)
ALLTLR3_male_LRS$ScaledAge2 <- as.numeric((ALLTLR3_male_LRS$ScaledAge)^2)
ALLTLR3_female_rate$ScaledAge2 <- as.numeric((ALLTLR3_female_rate$ScaledAge)^2)
ALLTLR3_male_rate$ScaledAge2 <- as.numeric((ALLTLR3_male_rate$ScaledAge)^2)

#    Lifetime reproductive success ####
#      LRS male####

male_LRS <- glmmTMB(CountOfSurvivingOffspringToOFL ~ TLR3+ScaledHs_obs+ScaledMHCDiversity+Aseua4 + (1 | HatchYearFactor), zi=~1, data = ALLTLR3_male_LRS, family=poisson)

# check zero inflation
hist(ALLTLR3_male_LRS$CountOfSurvivingOffspringToOFL,breaks=30)
male_LRS_notzero <- glmmTMB(CountOfSurvivingOffspringToOFL ~ TLR3+ScaledHs_obs+ScaledMHCDiversity+Aseua4 + (1 | HatchYearFactor), data = ALLTLR3_male_LRS, family=poisson)
AIC(male_LRS, male_LRS_notzero) # zero inflated model much better - check other assumptions

plot_model(male_LRS)
plot_model(male_LRS,type = "re")

# check model assumptions
resp = simulateResiduals(male_LRS)
plot(resp, rank = T)
simulationOutput <- simulateResiduals(fittedModel = male_LRS, n = 1000)
testDispersion(simulationOutput = simulationOutput, alternative ="less")
testDispersion(simulationOutput = simulationOutput, alternative ="greater")
testUniformity(simulationOutput = simulationOutput)

# check variance inflation - all under 2
vif(glm(CountOfSurvivingOffspringToOFL ~ TLR3+ScaledHs_obs+ScaledMHCDiversity+Aseua4  ,data = ALLTLR3_male_LRS))

#use MUMIN::DREDGE to test all models options and find best one and average best ones
male_LRS_dredge<-dredge(male_LRS)
male_LRS_dredge
model.avg(male_LRS_dredge)
plot(male_LRS_dredge)
averaged_male_LRS<-summary(model.avg(get.models(male_LRS_dredge, subset = delta < 7)))
averaged_male_LRS
sw(averaged_male_LRS)

# How do AC males differ from CC males?
male_LRSc <- glmmTMB(CountOfSurvivingOffspringToOFL ~ TLR3c+ScaledHs_obs+ScaledMHCDiversity+Aseua4 + (1 |HatchYearFactor), zi=~1, data = ALLTLR3_male_LRS, family=poisson)
#use MUMIN::DREDGE to test all models options and find best one and average best ones
male_LRSc_dredge<-dredge(male_LRSc)
averaged_male_LRSc<-summary(model.avg(get.models(male_LRSc_dredge, subset = delta < 7)))
averaged_male_LRSc

#      LRS female####

female_LRS <- glmmTMB(CountOfSurvivingOffspringToOFL ~ TLR3+ScaledHs_obs+ScaledMHCDiversity+Aseua4 + (1 | HatchYearFactor), zi=~1, data = ALLTLR3_female_LRS, family=poisson)

# check zero inflation
hist(ALLTLR3_female_LRS$CountOfSurvivingOffspringToOFL,breaks=30)
female_LRS_notzero <- glmmTMB(CountOfSurvivingOffspringToOFL ~ TLR3+ScaledHs_obs+ScaledMHCDiversity+Aseua4 + (1 | HatchYearFactor), data = ALLTLR3_female_LRS, family=poisson)
AIC(female_LRS, female_LRS_notzero) # zero inflated model much better - check other assumptions

plot_model(female_LRS)
plot_model(female_LRS,type = "re")

# check model assumptions
resp = simulateResiduals(female_LRS)
plot(resp, rank = T)
simulationOutput <- simulateResiduals(fittedModel = female_LRS, n = 1000)
testDispersion(simulationOutput = simulationOutput, alternative ="less")
testDispersion(simulationOutput = simulationOutput, alternative ="greater")
testUniformity(simulationOutput = simulationOutput)

# check variance inflation - all under 2
vif(glm(CountOfSurvivingOffspringToOFL ~ TLR3+ScaledHs_obs+ScaledMHCDiversity+Aseua4  ,data = ALLTLR3_female_LRS))

#use MUMIN::DREDGE to test all models options and find best one and average best ones
female_LRS_dredge<-dredge(female_LRS)
female_LRS_dredge
model.avg(female_LRS_dredge)
plot(female_LRS_dredge)
averaged_female_LRS<-summary(model.avg(get.models(female_LRS_dredge, subset = delta < 7)))
averaged_female_LRS
sw(averaged_female_LRS)

#    Reproductive rate ####
#      Reproductive rate male####

male_rate <- glmmTMB(CountOfSurvivingOffspringToOFL ~ TLR3+ScaledHs_obs+ScaledMHCDiversity+Aseua4 + ScaledAge+ScaledAge2+(1 | HatchYearFactor), zi=~1, data = ALLTLR3_male_rate, family=poisson)

# check zero inflation
hist(ALLTLR3_male_rate$CountOfSurvivingOffspringToOFL,breaks=30)
male_rate_notzero <- glmmTMB(CountOfSurvivingOffspringToOFL ~ TLR3+ScaledHs_obs+ScaledMHCDiversity+Aseua4 + ScaledAge + ScaledAge2 + (1 | HatchYearFactor), data = ALLTLR3_male_rate, family=poisson)
AIC(male_rate, male_rate_notzero) # zero inflated model slightly better - check other assumptions - none zero inflatted model is slightly overdispersed - proceed with zero inflatted model

plot_model(male_rate)
plot_model(male_rate,type = "re") # hatch year does not account for much variance

# check model assumptions - also check for none zero inflated model
resp = simulateResiduals(male_rate)
plot(resp, rank = T)
simulationOutput <- simulateResiduals(fittedModel = male_rate, n = 1000)
testDispersion(simulationOutput = simulationOutput, alternative ="less")
testDispersion(simulationOutput = simulationOutput, alternative ="greater")
testUniformity(simulationOutput = simulationOutput)

# check variance inflation - all under 2
vif(glm(CountOfSurvivingOffspringToOFL ~ TLR3+ScaledHs_obs+ScaledMHCDiversity+Aseua4+ScaledAge+ScaledAge2  ,data = ALLTLR3_male_rate))

#use MUMIN::DREDGE to test all models options and find best one and average best ones
male_rate_dredge<-dredge(male_rate)
male_rate_dredge
model.avg(male_rate_dredge)
plot(male_rate_dredge)
averaged_male_rate<-summary(model.avg(get.models(male_rate_dredge, subset = delta < 7)))
averaged_male_rate
sw(averaged_male_rate)

# How do AC males differ from CC males?
male_rateb <- glmmTMB(CountOfSurvivingOffspringToOFL ~ TLR3b+ScaledHs_obs+ScaledMHCDiversity+Aseua4 + ScaledAge+ScaledAge2+(1 | HatchYearFactor), zi=~1, data = ALLTLR3_male_rate, family=poisson)
#use MUMIN::DREDGE to test all models options and find best one and average best ones
male_rateb_dredge<-dredge(male_rateb)
averaged_male_rateb<-summary(model.avg(get.models(male_rateb_dredge, subset = delta < 7)))
averaged_male_rateb


#      Reproductive rate female####

female_rate <- glmmTMB(CountOfSurvivingOffspringToOFL ~ TLR3+ScaledHs_obs+ScaledMHCDiversity+Aseua4+ScaledAge+ScaledAge2 + (1 | HatchYearFactor), zi=~1, data = ALLTLR3_female_rate, family=poisson)

# check zero inflation
hist(ALLTLR3_female_rate$CountOfSurvivingOffspringToOFL,breaks=30)
female_rate_notzero <- glmmTMB(CountOfSurvivingOffspringToOFL ~ TLR3+ScaledHs_obs+ScaledMHCDiversity+Aseua4 + ScaledAge + ScaledAge2 + (1 | HatchYearFactor), data = ALLTLR3_female_rate, family=poisson)
AIC(female_rate, female_rate_notzero) # none zero inflated model slightly better - check other assumptions - none zero inflatted model is slightly overdispersed - and results do not differ between models- proceed with zero inflatted model

plot_model(female_rate)
plot_model(female_rate,type = "re") # hatch year does not account for much variance

# check model assumptions
resp = simulateResiduals(female_rate)
plot(resp, rank = T)
simulationOutput <- simulateResiduals(fittedModel = female_rate, n = 1000)
testDispersion(simulationOutput = simulationOutput, alternative ="less")
testDispersion(simulationOutput = simulationOutput, alternative ="greater")
testUniformity(simulationOutput = simulationOutput)

# check variance inflation - all under 2
vif(glm(CountOfSurvivingOffspringToOFL ~ TLR3+ScaledHs_obs+ScaledMHCDiversity+Aseua4 +ScaledAge+ScaledAge2 ,data = ALLTLR3_female_rate))

#use MUMIN::DREDGE to test all models options and find best one and average best ones
female_rate_dredge<-dredge(female_rate)
female_rate_dredge
model.avg(female_rate_dredge)
plot(female_rate_dredge)
averaged_female_rate<-summary(model.avg(get.models(female_rate_dredge, subset = delta < 7)))
averaged_female_rate
sw(averaged_female_rate)


#    Plot Figure 3 ####

# dataset for LRS - both sexes
ALLTLR3_LRS<- ALLTLR3 %>% filter(Survival_Status_WS2018 != "0")%>% droplevels()%>% droplevels()%>% filter(Translocated != "Yes")%>% droplevels()

# dataset for reproductive rate - both sexes
ALLTLR3_rate<- ALLTLR3 %>% filter(Survival_Status_WS2018 != "0")%>% droplevels()%>% droplevels()%>% filter(Translocated != "Yes")%>% droplevels() %>% filter(SurviveToFirstYear_BiAnnual != "N")%>% droplevels()

# Figure 3a - Lifetime reproductive success 
# summary statistics
LRS <- ddply(ALLTLR3_LRS, .(TLR3, SexEstimate), summarise, callmean=mean(CountOfSurvivingOffspringToOFL),  N    = length(CountOfSurvivingOffspringToOFL),sd   = sd(CountOfSurvivingOffspringToOFL),se   = sd / sqrt(N))

pd <- position_dodge(0.4)

LRSplot<-ggplot(LRS, aes(x=TLR3, y=callmean,  group=SexEstimate, color= SexEstimate))+geom_errorbar(aes(ymin=callmean-se, ymax=callmean+se), width=.2, size=1, position = pd)+geom_point(position = pd, size=3)+
  scale_color_manual(values=c('grey60','black'), labels = c("Female", "Male") ,name="Sex")+
  theme_classic()+
  labs(x="Genotype", y="Lifetime reproductive success (mean ± SE)") + 
  theme(axis.title.x=element_text(size=20, color="black"),axis.text.x =element_text(size=16, color="black"), legend.text = element_text(size=16), legend.title = element_text(size=20))+  theme(legend.position = 'none')+
  theme(axis.title.y=element_text(color="black", size=20, margin = margin(-5,10,0,0)),axis.text.y=element_text(size=16, color="black"), axis.line = element_line(colour = "black",size = 1, linetype = "solid")) +ylim(0, 2)+
  scale_x_discrete(labels=c("AA","AC", "CC"))+labs(fill = "Sex")

LRSplot

#  Figure 3b - Rate of reproductive success

# create annual reproductive success
ALLTLR3_rate$AnnualOffspringOFL <- as.numeric(ALLTLR3_rate$CountOfSurvivingOffspringToOFL/ALLTLR3_rate$AgeAtDeathCousinBiAnnualSeason)

# summary statistics
RRS  <- ddply(ALLTLR3_rate, .(TLR3, SexEstimate), summarise, callmean=mean(AnnualOffspringOFL),  N    = length(AnnualOffspringOFL),sd   = sd(AnnualOffspringOFL),se   = sd / sqrt(N))

pd <- position_dodge(0.4)

RRSplot<-ggplot(RRS, aes(x=TLR3, y=callmean,  group=SexEstimate, color= SexEstimate))+geom_errorbar(aes(ymin=callmean-se, ymax=callmean+se), width=.2, size=1, position = pd)+geom_point(position = pd, size=3)+
  scale_color_manual(values=c('grey60','black'), labels = c("Female", "Male") ,name="Sex")+
  theme_classic()+
  labs(x="Genotype", y="Rate of reproduction (mean ± SE)") + 
  theme(axis.title.x=element_text(size=20, color="black"),axis.text.x =element_text(size=16, color="black"), legend.text = element_text(size=16), legend.title = element_text(size=20))+  
  theme(axis.title.y=element_text(color="black", size=20, margin = margin(-5,10,0,0)),axis.text.y=element_text(size=16, color="black"), axis.line = element_line(colour = "black",size = 1, linetype = "solid")) +ylim(0, 0.4)+theme(legend.position=c(0.8, 0.9))+
  scale_x_discrete(labels=c("AA","AC", "CC"))+labs(fill = "Sex")

RRSplot

# add in labels
RRSplot <- arrangeGrob(RRSplot, top = textGrob("B", x = unit(0, "npc"), y   = unit(1, "npc"), just=c("left","top"),gp=gpar(col="black", type ="bold", fontsize=24, fontfamily="Times Roman")))

LRSplot <- arrangeGrob(LRSplot, top = textGrob("A", x = unit(0, "npc"), y   = unit(1, "npc"), just=c("left","top"),gp=gpar(col="black", type ="bold", fontsize=24, fontfamily="Times Roman")))

# create combined figure 3:
combo_plot<-grid.arrange( LRSplot,RRSplot, ncol = 2)
ggsave("Figure_3_RS.png", combo_plot, width = 14, height = 7, dpi = 300)





# Q2c: Selection coefficients ####

# Calculate selection coefficents over a course estimate of 3 overlapping generations
Selection_LRS <- ddply(ALLTLR3_LRS, .(TLR3), summarise, callmean=mean(CountOfSurvivingOffspringTo1Year))

Selection_LRS

# TLR3: AA set as 1

TLR3_AC<-1-(0.6230366/0.9363636)
TLR3_AC

TLR3_CC<-1-(0.5/0.9363636)
TLR3_CC
                       
# Calculate mean LRS for each genotype







######### Q3: Spatial patterns of TLR3 variation accross islands############################################################################################

# Import dataset - note all analysis for this section conducted using https://genepop.curtin.edu.au/

ISLAND<-read_csv("Genotype_frequencies_islands.csv")
head(ISLAND)
View(ISLAND)
str(ISLAND)

#Tidy data and convert variables to factors etc
ISLAND$Island <- as.factor(ISLAND$Island)
ISLAND$A_allele_frequency <- as.numeric(ISLAND$A_allele_frequency)
ISLAND$Year <- as.numeric(ISLAND$Year)
ISLAND$C_allele_frequency = as.numeric(ISLAND$C_allele_frequency )
ISLAND$Island <- factor(ISLAND$Island, levels(ISLAND$Island)[c(2,1,3,4,5)])

# Plot Figure 4 ####

Figure_4_Island<-ISLAND%>%
  ggplot(aes(x = Year, y = C_allele_frequency, group=Island, colour=Island))+
  geom_point(size = 3) +
  geom_line(aes(linetype=Island), size = 1.5)+
  scale_color_manual(values = c("#ffca47","grey80" , "grey60", "gray35", "gray1"))+
  scale_linetype_manual(values = c("solid","longdash" , "dotdash","dotted", "dashed"), name ="Island")+
  scale_x_continuous(breaks=c(1993, 1998, 2003, 2008, 2013, 2018))+
  labs(y = "C allele frequency", x = "Year") +ylim(0.20, 0.45)+ 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)))+
  theme_classic()+theme(axis.title.x=element_text(size=20),axis.text.x =element_text(size=16),legend.text = element_text(size=14), legend.title = element_text(size=16))+  theme(legend.position=c(0.2, 0.3))+theme(axis.title.y=element_text(size=20),axis.text.y=element_text(size=16))+ theme(legend.key.width = grid::unit(1.5, "cm"))+guides(colour = guide_legend(override.aes = list(shape = NA)))+annotate("text", x = 2016.5, y = 0.22, label = "-0.014", size=5, color="gray35")+annotate("text", x = 2019.5, y = 0.266, label = "-0.006", size=5, color="#ffca47")+annotate("text", x = 2019.5, y = 0.287, label = "-0.011", size=5, color="gray1")+annotate("text", x = 2015.5, y = 0.39, label = "-0.002", size=5, color="gray80")+annotate("text", x = 2020.5, y = 0.34, label = "-0.003", size=5, color="gray60")

combo_plot<-grid.arrange(Figure_4_Island, ncol = 1)
ggsave("Figure_4_Island.png", combo_plot, width = 10, height=6, dpi = 300)



######### Q4: HWE in fledglings sampled on Cousin###############################################################################################################
# Import dataset - note all analysis for this section conducted using https://genepop.curtin.edu.au/

HWE<-read_csv("HWE_graph_data.csv")
head(HWE)
View(HWE)
str(HWE)


######### Supplementary Material #######################################################

# Table S1: Survival model using COXME - presence/absence of alleles ####
survival1 <- coxme(Surv(AgeAtDeathCousinBiAnnualSeason,Survival_Status) ~ SexEstimate+A_allele_presence+C_allele_presence + Aseua4+MHCDiversity+Hs_obs+Maternal_Hs_obs+SeasonBorn+(1|HatchYearFactor),data = ALLTLR3)
survival1

# check assumptions - note this is just using coxph model

# first run model not including random factor
surv_object <- Surv(time = ALLTLR3$AgeAtDeathCousinBiAnnualSeason, event = ALLTLR3$Survival_Status)
survival_coxph_allele <- coxph(surv_object ~ A_allele_presence+C_allele_presence + SexEstimate+Aseua4+MHCDiversity+Maternal_Hs_obs+Hs_obs+SeasonBorn, data = ALLTLR3)
ggforest(survival_coxph, data = ALLTLR3)

# check assumptions
test.ph <- cox.zph(survival_coxph_allele)
test.ph # test is not statistically significant for each of the covariates, and the global test is also not statistically significant. Therefore, we can assume the proportional hazards.
ggcoxzph(test.ph) # all looks fine
ggcoxdiagnostics(survival_coxph_allele, type = "dfbeta", linear.predictions = FALSE, ggtheme = theme_bw())  #also all look fine

#test autoccrelation
cvif <- vif(survival_coxph_allele)
cvif #looks fine



# Table S2: Reproductive success - presence/absence of alleles ####
#           LRS male - presence/absence of allele ####

male_LRS_allele <- glmmTMB(CountOfSurvivingOffspringToOFL ~ A_allele_presence+C_allele_presence+ScaledHs_obs+ScaledMHCDiversity+Aseua4 + (1 | HatchYearFactor), zi=~1, data = ALLTLR3_male_LRS, family=poisson)

plot_model(male_LRS_allele)
plot_model(male_LRS_allele,type = "re")

# check model assumptions
resp = simulateResiduals(male_LRS_allele)
plot(resp, rank = T)
simulationOutput <- simulateResiduals(fittedModel = male_LRS_allele, n = 1000)
testDispersion(simulationOutput = simulationOutput, alternative ="less")
testDispersion(simulationOutput = simulationOutput, alternative ="greater")
testUniformity(simulationOutput = simulationOutput)

# check variance inflation - all under 2
vif(glm(CountOfSurvivingOffspringToOFL ~ A_allele_presence+C_allele_presence+ScaledHs_obs+ScaledMHCDiversity+Aseua4  ,data = ALLTLR3_male_LRS))

#use MUMIN::DREDGE to test all models options and find best one and average best ones
male_LRS_allele_dredge<-dredge(male_LRS_allele)
male_LRS_allele_dredge
model.avg(male_LRS_allele_dredge)
plot(male_LRS_allele_dredge)
averaged_male_LRS_allele<-summary(model.avg(get.models(male_LRS_allele_dredge, subset = delta < 7)))
averaged_male_LRS_allele
sw(averaged_male_LRS)

#           LRS female - presence/absence of allele ####

female_LRS_allele <- glmmTMB(CountOfSurvivingOffspringToOFL ~ A_allele_presence+C_allele_presence+ScaledHs_obs+ScaledMHCDiversity+Aseua4 + (1 | HatchYearFactor), zi=~1, data = ALLTLR3_female_LRS, family=poisson)

plot_model(female_LRS_allele)
plot_model(female_LRS_allele,type = "re")

# check model assumptions
resp = simulateResiduals(female_LRS_allele)
plot(resp, rank = T)
simulationOutput <- simulateResiduals(fittedModel = female_LRS_allele, n = 1000)
testDispersion(simulationOutput = simulationOutput, alternative ="less")
testDispersion(simulationOutput = simulationOutput, alternative ="greater")
testUniformity(simulationOutput = simulationOutput)

# check variance inflation - all under 2
vif(glm(CountOfSurvivingOffspringToOFL ~ A_allele_presence+C_allele_presence+ScaledHs_obs+ScaledMHCDiversity+Aseua4  ,data = ALLTLR3_female_LRS))

#use MUMIN::DREDGE to test all models options and find best one and average best ones
female_LRS_allele_dredge<-dredge(female_LRS_allele)
female_LRS_allele_dredge
model.avg(female_LRS_allele_dredge)
plot(female_LRS_allele_dredge)
averaged_female_LRS_allele<-summary(model.avg(get.models(female_LRS_allele_dredge, subset = delta < 7)))
averaged_female_LRS_allele
sw(averaged_female_LRS_allele)

#           Reproductive rate male - presence/absence of allele ####

male_rate_allele <- glmmTMB(CountOfSurvivingOffspringToOFL ~ A_allele_presence+C_allele_presence+ScaledHs_obs+ScaledMHCDiversity+Aseua4 + ScaledAge+ScaledAge2+(1 | HatchYearFactor), zi=~1, data = ALLTLR3_male_rate, family=poisson)

plot_model(male_rate_allele)
plot_model(male_rate_allele,type = "re")

# check model assumptions
resp = simulateResiduals(male_rate_allele)
plot(resp, rank = T)
simulationOutput <- simulateResiduals(fittedModel = male_rate_allele, n = 1000)
testDispersion(simulationOutput = simulationOutput, alternative ="less")
testDispersion(simulationOutput = simulationOutput, alternative ="greater")
testUniformity(simulationOutput = simulationOutput)

# check variance inflation - all under 2
vif(glm(CountOfSurvivingOffspringToOFL ~ A_allele_presence+C_allele_presence+ScaledHs_obs+ScaledMHCDiversity+Aseua4+ScaledAge+ScaledAge2  ,data = ALLTLR3_male_rate))

#use MUMIN::DREDGE to test all models options and find best one and average best ones
male_rate_allele_dredge<-dredge(male_rate_allele)
male_rate_allele_dredge
model.avg(male_rate_allele_dredge)
plot(male_rate_allele_dredge)
averaged_male_rate_allele<-summary(model.avg(get.models(male_rate_allele_dredge, subset = delta < 7)))
averaged_male_rate_allele
sw(averaged_male_rate_allele)

#           Reproductive rate female - presence/absence of allele ####

female_rate_allele <- glmmTMB(CountOfSurvivingOffspringToOFL ~ A_allele_presence+C_allele_presence+ScaledHs_obs+ScaledMHCDiversity+Aseua4+ScaledAge+ScaledAge2 + (1 | HatchYearFactor), zi=~1, data = ALLTLR3_female_rate, family=poisson)

plot_model(female_rate_allele)
plot_model(female_rate_allele,type = "re")# Hatch year does not account for much variance

# check model assumptions
resp = simulateResiduals(female_rate_allele)
plot(resp, rank = T)
simulationOutput <- simulateResiduals(fittedModel = female_rate_allele, n = 1000)
testDispersion(simulationOutput = simulationOutput, alternative ="less")
testDispersion(simulationOutput = simulationOutput, alternative ="greater")
testUniformity(simulationOutput = simulationOutput)

# check variance inflation - all under 2
vif(glm(CountOfSurvivingOffspringToOFL ~ A_allele_presence+C_allele_presence+ScaledHs_obs+ScaledMHCDiversity+Aseua4 +ScaledAge+ScaledAge2 ,data = ALLTLR3_female_rate))

#use MUMIN::DREDGE to test all models options and find best one and average best ones
female_rate_allele_dredge<-dredge(female_rate_allele) # models do not converge

# re-run analysis - but not including hatch year as a random factor
female_rate_allele_norandom <- glmmTMB(CountOfSurvivingOffspringToOFL ~ A_allele_presence+C_allele_presence+ScaledHs_obs+ScaledMHCDiversity+Aseua4+ScaledAge+ScaledAge2,zi=~1, data = ALLTLR3_female_rate, family=poisson)

plot_model(female_rate_allele_norandom)

summary(female_rate_allele)
summary(female_rate_allele_norandom)
# model outputs seem very similar - check assumptions to confirm which model a better fit

# check model assumptions
resp = simulateResiduals(female_rate_allele_norandom)
plot(resp, rank = T)
simulationOutput <- simulateResiduals(fittedModel = female_rate_allele_norandom, n = 1000)
testDispersion(simulationOutput = simulationOutput, alternative ="less")
testDispersion(simulationOutput = simulationOutput, alternative ="greater")
testUniformity(simulationOutput = simulationOutput)

# check variance inflation - all under 2
vif(glm(CountOfSurvivingOffspringToOFL ~ A_allele_presence+C_allele_presence+ScaledHs_obs+ScaledMHCDiversity+Aseua4 +ScaledAge+ScaledAge2 ,data = ALLTLR3_female_rate))

#use MUMIN::DREDGE to test all models options and find best one and average best ones
female_rate_allele_norandom_dredge<-dredge(female_rate_allele_norandom)
female_rate_allele_norandom_dredge
model.avg(female_rate_allele_norandom_dredge)
plot(female_rate_allele_norandom_dredge)
averaged_female_rate_allele_norandom<-summary(model.avg(get.models(female_rate_allele_norandom_dredge, subset = delta < 7)))
averaged_female_rate_allele_norandom
sw(averaged_female_rate_allele_norandom)

# Table S3-S4 analysed using genepop as in Q3 above####

# Plot Supplementary figure S1####
HWE_plot<-ggplot(HWE, aes(x=Genotype, y=value, fill=factor(Obs_Exp)))+geom_col(position = "dodge", colour="black")+
  scale_fill_manual(values=c('black','lightgrey'), labels = c("Observed", "Expected"))+
  theme_classic()+
  labs(x="Genotype", y="Individuals") + 
  theme(axis.title.x=element_text(size=20, color="black"),axis.text.x =element_text(size=16, color="black"), legend.text = element_text(size=16),legend.title = element_blank())+ theme(legend.position =c(0.9, 0.85))+
  theme(axis.title.y=element_text(color="black", size=20, margin = margin(-5,10,0,0)),axis.text.y=element_text(size=16, color="black"), axis.line = element_line(colour = "black",size = 1, linetype = "solid")) +
  scale_x_discrete(labels=c("AA","AC", "CC"))+facet_wrap(~Category) +theme(strip.text = element_text(size = 16))

HWE_plot

ggsave("Figure_S1.png", HWE_plot, width = 14, height = 7, dpi = 300)

# View session info
sessionInfo()
