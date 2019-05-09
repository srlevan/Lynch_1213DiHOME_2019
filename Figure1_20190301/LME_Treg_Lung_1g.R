#Clear Environment#
rm(list=ls())

#load required packages#
library(lmerTest)
library(visreg)
library(ggplot2)
library(plyr)

#Load File#

setwd("~/Box Sync/Lynch Lab/Code Respositories/Figure 1_20190301/")

dihome=read.csv("INPUT_Mouse_Treg_Lung.csv")
dihome


#Use LME to correct for assay variability#
Input <- dihome
Input$Treatment<-factor(Input$Treatment,c("PBS","CRA","DiHOME"))

#Model Treg Frequency#
mod.Treg<- lmer(Treg~Treatment+(1|Assay),data=Input)
y.Treg<-summary(mod.Treg)
y.Treg
coef(mod.Treg)
mod.Treg$res
vTreg<-visreg(fit=mod.Treg,xvar="Treatment",type="conditional")

#summarize each individual point from the model#
resTreg<-vTreg$res
resTreg["Don"] <- Input$Assay
resTreg
resTreg$Repeat <- ifelse(resTreg$Don=="P1A2",15,ifelse(resTreg$Don=="P1A3", 16, ifelse(resTreg$Don=="P1A4",17, ifelse(resTreg$Don=="P1A5", 18, ifelse(resTreg$Don=="P1A6",3,0)))))


#Create a table of the points taht will be included in the graph#
write.table(resTreg,file="lme_graphpoints_Treg.txt",sep="\t",quote=F,row.names=F)

#read the table of graphing points#
lmepoint<-read.table("lme_graphpoints_Treg.txt",header=T,sep="\t")

#Factor the table by treatment#
lmepoint$Treatment = factor(lmepoint$Treatment,c("PBS","CRA", "DiHOME"))
lmepoint
lmepoint$Sample <- Input$Sample.

#Generate a Bar Graph of the donor-corrected Treg Data#
Treg <- ggplot(data=lmepoint, aes(x=Treatment, y=visregRes)) + 
  geom_boxplot()+
  geom_jitter(data=lmepoint, aes(x=Treatment, y=visregRes), position=position_jitter(0.2), size=5, shape=lmepoint$Repeat)+
  theme_classic() + scale_fill_brewer(palette = "Greys", direction=1)+
  guides(fill = guide_legend(reverse=FALSE))
Treg
ggsave(filename="GRAPH-MouseLungTreg.pdf", plot=Treg,  useDingbats=FALSE, width=8.5, height=11, units="in")

#Compared to PBS#
Input$Treatment <- factor(Input$Treatment, c("PBS","CRA","DiHOME"))
model <- lmer(Treg ~ Treatment + (1|Assay), data = Input)
summary(model)

#Compared to DiHOME#
Input$Treatment <- factor(Input$Treatment, c("DiHOME","PBS","CRA"))
model <- lmer(Treg ~ Treatment + (1|Assay), data = Input)
summary(model)

#Remove Outliers and Repeat Statistics#
lmepoint
NoOutlier <- droplevels(subset(Input, Input$Sample.!="Lung_dc_005_039.fcs" & Input$Sample.!="Lung_dc_003_037.fcs" & Input$Sample.!="Lung-15_dc_012_026.fcs"))

#Compared to PBS#
NoOutlier$Treatment <- factor(NoOutlier$Treatment, c("PBS","CRA","DiHOME"))
model <- lmer(Treg ~ Treatment + (1|Assay), data = NoOutlier)
summary(model)

#Compared to DiHOME#
NoOutlier$Treatment <- factor(NoOutlier$Treatment, c("DiHOME","PBS","CRA"))
model <- lmer(Treg ~ Treatment + (1|Assay), data = NoOutlier)
summary(model)


