#Clear the Environment#
rm(list=ls())

#Load Required Packages#
library(lmerTest)
library(visreg)
library(ggplot2)
library(plyr)
library(gridExtra)
library(reshape2)
library(dplyr)

#set wd#
setwd("~/Box Sync/Lynch Lab/Code Respositories/Figure 1_20190301/")

#read in LCMS data for Mouse Lung and Plasma#
ConcP=read.csv("INPUT-LCMS-MousePlasma.csv")
ConcL=read.csv("INPUT-LCMS-MouseLung.csv")


#DiHOMEng.uL.2#
#Use LME to Correct for Run Variability#
Input <- ConcP
mod.DiHOMEng.uL.2<- lmer(DiHOMEng.uL.2~Sample.Type+(1|Assay),data=Input)
y.DiHOMEng.uL.2<-summary(mod.DiHOMEng.uL.2)
y.DiHOMEng.uL.2
coef(mod.DiHOMEng.uL.2)
mod.DiHOMEng.uL.2$res
vDiHOMEng.uL.2<-visreg(fit=mod.DiHOMEng.uL.2,xvar="Sample.Type",type="conditional")

#summarize each individual point from the model#
resDiHOMEng.uL.2<-vDiHOMEng.uL.2$res
resDiHOMEng.uL.2["Don"] <- Input$Assay
resDiHOMEng.uL.2
resDiHOMEng.uL.2$Repeat <- ifelse(resDiHOMEng.uL.2$Don=="P2A1",4,ifelse(resDiHOMEng.uL.2$Don=="P2A2", 8, 0))

#Create a table of the points taht will be included in the graph#
write.table(resDiHOMEng.uL.2,file="lme_graphpoints_DiHOMEng.uL.2.txt",sep="\t",quote=F,row.names=F)

#Read in Table#
lmepoint<-read.table(file="lme_graphpoints_DiHOMEng.uL.2.txt",header=T, sep="\t")
#Factor the table by Sample.Type#
lmepoint$Sample.Type = factor(lmepoint$Sample.Type,c("Control", "DiHOME"))
lmepoint

#plot the table in a dot plot#
#Generate plot space with lmepoint data#
DiHOMEng.uL.2 <- ggplot(data=lmepoint, aes(x=Sample.Type, y=visregRes)) + 
  geom_boxplot()+
   geom_jitter(data=lmepoint, aes(x=Sample.Type, y=visregRes), position=position_jitter(0.2), size=5, shape=lmepoint$Repeat)+
  theme_classic() + scale_fill_brewer(palette = "Greys", direction=1)+
  guides(fill = guide_legend(reverse=FALSE))
DiHOMEng.uL.2
ggsave(filename="GRAPH-DiHOMEng.uL.2_MouseLung.pdf", plot=DiHOMEng.uL.2,  useDingbats=FALSE, width=8.5, height=11, units="in")

#Compared to PBS#
Input$Sample.Type <- factor(Input$Sample.Type, c("Control","DiHOME"))
model <- lmer(DiHOMEng.uL.2 ~ Sample.Type + (1|Assay), data = Input)
summary(model)

#Remove Outliers and Repeat Statistics#
lmepoint
NoOutlier <- droplevels(subset(Input, Input$Sample.Name!="P D6"))

#Compared to PBS#
NoOutlier$Sample.Type <- factor(NoOutlier$Sample.Type, c("Control","DiHOME"))
model <- lmer(DiHOMEng.uL.2 ~ Sample.Type + (1|Assay), data = NoOutlier)
summary(model)


#DiHOMEng.mg.1#
#Use LME to Correct for Run Variability#
Input <- ConcL
mod.DiHOMEng.mg.1<- lmer(DiHOMEng.mg.1~Sample.Type+(1|Assay),data=Input)
y.DiHOMEng.mg.1<-summary(mod.DiHOMEng.mg.1)
y.DiHOMEng.mg.1
coef(mod.DiHOMEng.mg.1)
mod.DiHOMEng.mg.1$res
vDiHOMEng.mg.1<-visreg(fit=mod.DiHOMEng.mg.1,xvar="Sample.Type",type="conditional")

#summarize each individual point from the model#
resDiHOMEng.mg.1<-vDiHOMEng.mg.1$res
resDiHOMEng.mg.1["Don"] <- Input$Assay
resDiHOMEng.mg.1
resDiHOMEng.mg.1$Repeat <- ifelse(resDiHOMEng.mg.1$Don=="P2A1",4,ifelse(resDiHOMEng.mg.1$Don=="P2A2", 8, 0))

#Create a table of the points taht will be included in the graph#
write.table(resDiHOMEng.mg.1,file="lme_graphpoints_DiHOMEng.mg.1.txt",sep="\t",quote=F,row.names=F)

#Read in Table#
lmepoint<-read.table(file="lme_graphpoints_DiHOMEng.mg.1.txt",header=T, sep="\t")
#Factor the table by Sample.Type#
lmepoint$Sample.Type = factor(lmepoint$Sample.Type,c("Control", "DiHOME"))
lmepoint

#plot the table in a dot plot#
#Generate plot space with lmepoint data#
DiHOMEng.mg.1 <- ggplot(data=lmepoint, aes(x=Sample.Type, y=visregRes)) + 
  geom_boxplot()+
    geom_jitter(data=lmepoint, aes(x=Sample.Type, y=visregRes), position=position_jitter(0.2), size=5, shape=lmepoint$Repeat)+
  theme_classic() + scale_fill_brewer(palette = "Greys", direction=1)+
  guides(fill = guide_legend(reverse=FALSE))
DiHOMEng.mg.1
ggsave(filename="GRAPH-DiHOMEng.mg.1_MouseLung.pdf", plot=DiHOMEng.mg.1,  useDingbats=FALSE, width=8.5, height=11, units="in")

#Compared to PBS#
Input$Sample.Type <- factor(Input$Sample.Type, c("Control","DiHOME"))
model <- lmer(DiHOMEng.mg.1 ~ Sample.Type + (1|Assay), data = Input)
summary(model)


