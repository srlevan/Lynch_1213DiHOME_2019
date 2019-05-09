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
            N    = length(Count),
            MeanCount= mean(Count)
)
BronchialDataSummary


#Use a LME to correct for Run and Assay variability and compare to DiHOME-Treatment
x <- matrix(data = NA, nrow =  3, ncol = 3, dimnames=list(1:3, c("Estimate", "Standard Error", "P")))
Input <- BronchialDataSummary
Input

for (i in 5){
  model <- lmer(Input[[i]] ~ Treatment + (1|Assay), data = Input)
  y <- summary(model)
  x[,1] <- y$coefficients[,1]
  x[,2] <- y$coefficients[,2]
  x[,3] <- y$coefficients[,5]
  print(colnames(Input[i]))
  print(x)
}

#Run Test#
Input<- BronchialDataSummary
Input$Treatment <- factor(Input$Treatment, c("PBS","CRA","DiHOME"))

#Use LME to Correct for Run Variability#
mod.MeanCount<- lmer(MeanCount~Treatment+(1|Assay),data=Input)
y.MeanCount<-summary(mod.MeanCount)
y.MeanCount
coef(mod.MeanCount)
mod.MeanCount$res
vMeanCount<-visreg(fit=mod.MeanCount,xvar="Treatment",type="conditional")

#summarize each individual point from the model#
resMeanCount<-vMeanCount$res
resMeanCount["Don"] <- Input$Assay
resMeanCount
levels(resMeanCount$Don)
resMeanCount$Repeat <- ifelse(resMeanCount$Don=="P1A2",15,ifelse(resMeanCount$Don=="P1A3", 16, ifelse(resMeanCount$Don=="P1A4",17, ifelse(resMeanCount$Don=="P1A5", 18, ifelse(resMeanCount$Don=="P1A6",3,0)))))

#Create a table of the points taht will be included in the graph#
write.table(resMeanCount,file="lme_graphpoints_MeanCount.txt",sep="\t",quote=F,row.names=F)

#Read in Table#
lmepoint<-read.table(file="lme_graphpoints_MeanCount.txt",header=T, sep="\t")
#Factor the table by Treatment#
lmepoint$Treatment = factor(lmepoint$Treatment,c("PBS","CRA", "DiHOME"))
lmepoint$Repeat

#plot the table in a dot plot#
#Generate plot space with lmepoint data#
MeanCount <- ggplot(lmepoint, aes(x=Treatment, y=visregRes)) +
  geom_boxplot()+
  geom_jitter(data=lmepoint, aes(x=Treatment, y=visregRes), position=position_jitter(0.2), size=5, shape=lmepoint$Repeat)+
  theme_classic() + scale_fill_brewer(palette = "Greys", direction=1)+
  guides(fill = guide_legend(reverse=FALSE))
MeanCount
ggsave(filename="GRAPH-MeanCount_MouseBronchial.pdf", plot=MeanCount,  useDingbats=FALSE, width=8.5, height=11, units="in")

#Model Statistics#
#Compared to PBS#
Input$Treatment <- factor(Input$Treatment, c("PBS","CRA","DiHOME"))
model <- lmer(MeanCount ~ Treatment + (1|Assay), data = Input)
summary(model)

#Compared to DiHOME#
Input$Treatment <- factor(Input$Treatment, c("DiHOME","PBS","CRA"))
model <- lmer(MeanCount ~ Treatment + (1|Assay), data = Input)
summary(model)

