#Clear Environment#
rm(list=ls())

#load required packages#
library(lmerTest)
library(visreg)
library(ggplot2)

#Load File#
setwd("~/Figure2_20190301/")

dihome=read.csv("INPUT-IL10CBA.csv")
dihome

#Subset raw data to include treatments and tcell subtypes of interest#
S1 = dihome
S1$Treatment = factor(S1$Treatment, c("DMSO", "DiHOME-75", "DiHOME-130", "DiHOME-200"))
S1

#Include only Treatment, Donor and Final Concentration#
Input<-S1[,c("Treatment","Donor","Sample.Name","Final.CC")]

#Plot Linear Mixed Effects Model as a Dot plot#
#model IL10 Final.CC#
mod.Final.CC<- lmer(Final.CC~Treatment+(1|Donor),data=Input)
y.Final.CC<-summary(mod.Final.CC)
y.Final.CC
coef(mod.Final.CC)
mod.Final.CC$res
vFinal.CC<-visreg(fit=mod.Final.CC,xvar="Treatment",type="conditional")
summary(y.Final.CC)

#summarize each individual point from the model#
resFinal.CC<-vFinal.CC$res
resFinal.CC["Don"] <- Input$Donor
resFinal.CC
levels(resFinal.CC$Don) <- c("16","17","18")

#Create a table of the points taht will be included in the graph#
write.table(resFinal.CC,file="lme_graphpoints_IL10_Final.CC.txt",sep="\t",quote=F,row.names=F)

#read the table of graphing points#
lmepoint<-read.table("lme_graphpoints_IL10_Final.CC.txt",header=T,sep="\t")
#Factor the table by treatment#
lmepoint$Treatment = factor(lmepoint$Treatment,c("DMSO","DiHOME-75", "DiHOME-130", "DiHOME-200"))

lmeeffect

#plot the table in a dot plot#
#Generate plot space with lmepoint data#
Final.CC <- ggplot(data=lmepoint, aes(x=Treatment, y=visregRes)) + 
  geom_boxplot()+
  geom_jitter(data=lmepoint, aes(x=Treatment, y=visregRes), shape=lmepoint$Don, position=position_jitter(0.3), size=8)+
   theme(axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.border = element_blank(),
         panel.background = element_blank()) +
   labs(title="Final.CC_SubsetSamples", x="Treatment", y="Final.CC Frequency")
Final.CC
ggsave(filename="GRAPH-IL10Concentration.pdf", plot=Final.CC,  useDingbats=FALSE,width=11, height=8.5, units="in")

#Compared to DMSO#
Input$Group <- factor(Input$Treatment,c("DMSO","DiHOME-75", "DiHOME-130", "DiHOME-200"))
model <- lmer(Final.CC ~ Group + (1|Donor), data = Input)
summary(model)

#Compared to DiHOME-200#
Input$Group <- factor(Input$Treatment,c("DiHOME-200","DMSO","DiHOME-75", "DiHOME-130"))
model <- lmer(Final.CC ~ Group + (1|Donor), data = Input)
summary(model)

#Compared to DiHOME-75#
Input$Group <- factor(Input$Treatment,c("DiHOME-75","DMSO","DiHOME-130", "DiHOME-200"))
model <- lmer(Final.CC ~ Group + (1|Donor), data = Input)
summary(model)


#Remove Outliers and Repeat Statistics#
lmepoint
NoOutlier <- droplevels(subset(Input, Input$Sample.Name!="DL_C04"))

#Compared to DMSO#
NoOutlier$Group <- factor(NoOutlier$Treatment,c("DMSO","DiHOME-75", "DiHOME-130", "DiHOME-200"))
model <- lmer(Final.CC ~ Group + (1|Donor), data = NoOutlier)
summary(model)

#Compared to DiHOME-200#
NoOutlier$Group <- factor(NoOutlier$Treatment,c("DiHOME-200","DMSO","DiHOME-75", "DiHOME-130"))
model <- lmer(Final.CC ~ Group + (1|Donor), data = NoOutlier)
summary(model)