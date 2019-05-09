#Clear Environment#
rm(list=ls())

#Load Required Packages#
library(ggplot2)
library(visreg)
library(lmerTest)
library(plyr)
library(gridExtra)
library(reshape2)
library(nplr)
library(drc)
library(sandwich)
library(lmtest)
library(multcomp)

#Loadfiles#
setwd("~/Box Sync/Lynch Lab/Code Respositories/Figure 1_20190301/")

#Read in Average Concentrations#
BronchialData <- read.csv("INPUT-BronchialInfiltration_20180817.csv")
BronchialData$Treatment <- factor(BronchialData$Treatment, c("DiHOME","PBS","CRA"))
BronchialData

#Summarize Count Data for Bronchials#
BronchialDataSummary <- ddply(BronchialData, c("Assay","Treatment", "Group"), summarise,
                           N    = length(Avg_Neutrophil.Score),
                           MeanScore= mean(Avg_Neutrophil.Score)
)
BronchialDataSummary


#Run Test#
Input<- BronchialDataSummary
Input$Treatment <- factor(Input$Treatment, c("PBS","CRA","DiHOME"))

#Use LME to Correct for Run Variability#
mod.MeanScore<- lmer(MeanScore~Treatment+(1|Assay),data=Input)
y.MeanScore<-summary(mod.MeanScore)
y.MeanScore
coef(mod.MeanScore)
mod.MeanScore$res
vMeanScore<-visreg(fit=mod.MeanScore,xvar="Treatment",type="conditional")

#summarize each individual point from the model#
resMeanScore<-vMeanScore$res
resMeanScore["Don"] <- Input$Assay
resMeanScore
resMeanScore$Repeat <- ifelse(resMeanScore$Don=="P1A2",15,ifelse(resMeanScore$Don=="P1A3", 16, ifelse(resMeanScore$Don=="P1A4",17, ifelse(resMeanScore$Don=="P1A5", 18, ifelse(resMeanScore$Don=="P1A6",3,0)))))

#Create a table of the points taht will be included in the graph#
write.table(resMeanScore,file="lme_graphpoints_MeanScore.txt",sep="\t",quote=F,row.names=F)

#Read in Table#
lmepoint<-read.table(file="lme_graphpoints_MeanScore.txt",header=T, sep="\t")
#Factor the table by Treatment#
lmepoint$Treatment = factor(lmepoint$Treatment,c("PBS","CRA", "DiHOME"))
lmepoint

#plot the table in a dot plot#
#Generate plot space with lmepoint data#
MeanScore <- ggplot(data=lmepoint, aes(x=Treatment, y=visregRes)) + 
  geom_boxplot()+
  geom_jitter(data=lmepoint, aes(x=Treatment, y=visregRes), position=position_jitter(0.2), size=5, shape=lmepoint$Repeat)+
  theme_classic() + scale_fill_brewer(palette = "Greys", direction=1)+
  guides(fill = guide_legend(reverse=FALSE))
MeanScore
ggsave(filename="GRAPH-MeanScore_MouseBronchial.pdf", plot=MeanScore,  useDingbats=FALSE, width=8.5, height=11, units="in")

#Model Statistics#
#Compared to PBS#
Input$Treatment <- factor(Input$Treatment, c("PBS","CRA","DiHOME"))
model <- lmer(MeanScore ~ Treatment + (1|Assay), data = Input)
summary(model)

#Compared to DiHOME#
Input$Treatment <- factor(Input$Treatment, c("DiHOME","PBS","CRA"))
model <- lmer(MeanScore ~ Treatment + (1|Assay), data = Input)
summary(model)

#Remove Outlier and Repeat Statistics#
NoOutlier <- droplevels(subset(Input, Input$Group!="P1A6_12"))

#Compared to PBS#
NoOutlier$Treatment <- factor(NoOutlier$Treatment, c("PBS","CRA","DiHOME"))
model <- lmer(MeanScore ~ Treatment + (1|Assay), data = NoOutlier)
summary(model)

#Compared to DiHOME#
NoOutlier$Treatment <- factor(NoOutlier$Treatment, c("DiHOME","PBS","CRA"))
model <- lmer(MeanScore ~ Treatment + (1|Assay), data = NoOutlier)
summary(model)

