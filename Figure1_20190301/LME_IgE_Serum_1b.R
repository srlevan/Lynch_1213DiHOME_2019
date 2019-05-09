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
SerumData <- read.csv("INPUT-MouseIgE_20180724.csv")
SerumData$Group <- factor(SerumData$Group, c("DiHOME","PBS","CRA"))
SerumData

#Run Test#
Input<- SerumData
Input$Group <- factor(Input$Group, c("PBS","CRA","DiHOME"))

#Serum.igE#
#Use LME to Correct for Run Variability#
mod.Serum.igE<- lmer(Serum.igE~Group+(1|Assay),data=Input)
y.Serum.igE<-summary(mod.Serum.igE)
y.Serum.igE
coef(mod.Serum.igE)
mod.Serum.igE$res
vSerum.igE<-visreg(fit=mod.Serum.igE,xvar="Group",type="conditional")

#summarize each individual point from the model#
resSerum.igE<-vSerum.igE$res
resSerum.igE["Don"] <- Input$Assay
resSerum.igE
resSerum.igE$Repeat <- ifelse(resSerum.igE$Don=="P1A2",15,ifelse(resSerum.igE$Don=="P1A3", 16, ifelse(resSerum.igE$Don=="P1A4",17, ifelse(resSerum.igE$Don=="P1A5", 18, ifelse(resSerum.igE$Don=="P1A6",3,0)))))

#Create a table of the points taht will be included in the graph#
write.table(resSerum.igE,file="lme_graphpoints_Serum.igE.txt",sep="\t",quote=F,row.names=F)

#Read in Table#
lmepoint<-read.table(file="lme_graphpoints_Serum.igE.txt",header=T, sep="\t")
#Factor the table by Group#
lmepoint$Group = factor(lmepoint$Group,c("PBS","CRA", "DiHOME"))
lmepoint

#plot the table in a dot plot#
#Generate plot space with lmepoint data#
Serum.igE <- ggplot() + 
  geom_boxplot(data=lmepoint, aes(x=Group, y=visregRes))+
  geom_jitter(data=lmepoint, aes(x=Group, y=visregRes), position=position_jitter(0.2), size=5, shape=lmepoint$Repeat)+
  theme_classic() + scale_fill_brewer(palette = "Greys", direction=1)+
  guides(fill = guide_legend(reverse=FALSE))+
  ylim(0,125)
Serum.igE
ggsave(filename="GRAPH-Serum.IgE_MouseLung.pdf", plot=Serum.igE,  useDingbats=FALSE, width=8.5, height=11, units="in")


#Model Statistics#
#Compared to PBS#
Input$Group <- factor(Input$Group, c("PBS","CRA","DiHOME"))
model <- lmer(Serum.igE ~ Group + (1|Assay), data = Input)
summary(model)

#Compared to DiHOME#
Input$Group <- factor(Input$Group, c("DiHOME","PBS","CRA"))
model <- lmer(Serum.igE ~ Group + (1|Assay), data = Input)
summary(model)



#Remove Outlier and Repeat Statistics#
NoOutlier <- droplevels(subset(Input, Input$Sample!="P1A2-2" & Input$Sample!="P1A3-4" & Input$Sample!="P1A5-C1" & Input$Sample!="P1A5-C3"))

#Compared to PBS#
NoOutlier$Group <- factor(NoOutlier$Group, c("PBS","CRA","DiHOME"))
model <- lmer(Serum.igE ~ Group + (1|Assay), data = NoOutlier)
summary(model)

#Compared to DiHOME#
NoOutlier$Group <- factor(NoOutlier$Group, c("DiHOME","PBS","CRA"))
model <- lmer(Serum.igE ~ Group + (1|Assay), data = NoOutlier)
summary(model)
