#Clear Environment#
rm(list=ls())

#load required packages#
library(lmerTest)
library(visreg)
library(ggplot2)
library(plyr)

#Load File#

setwd("~/Box Sync/Lynch Lab/Code Respositories/Figure 1_20190301/")

dihome=read.csv("INPUT-MouseCD3RLung.csv")
dihome


#Use LME to correct for assay variability#
Input <- dihome
Input$Treatment<-factor(Input$Treatment,c("PBS","CRA","DiHOME"))

#Model CD3R Frequency#
mod.CD3R<- lmer(CD3R~Treatment+(1|Assay),data=Input)
y.CD3R<-summary(mod.CD3R)
y.CD3R
coef(mod.CD3R)
mod.CD3R$res
vCD3R<-visreg(fit=mod.CD3R,xvar="Treatment",type="conditional")

#summarize each individual point from the model#
resCD3R<-vCD3R$res
resCD3R["Don"] <- Input$Assay
resCD3R
resCD3R$Repeat <- ifelse(resCD3R$Don=="P1A2",15,ifelse(resCD3R$Don=="P1A3", 16, ifelse(resCD3R$Don=="P1A4",17, ifelse(resCD3R$Don=="P1A5", 18, ifelse(resCD3R$Don=="P1A6",3,0)))))


#Create a table of the points taht will be included in the graph#
write.table(resCD3R,file="lme_graphpoints_CD3R.txt",sep="\t",quote=F,row.names=F)

#read the table of graphing points#
lmepoint<-read.table("lme_graphpoints_CD3R.txt",header=T,sep="\t")

#Factor the table by treatment#
lmepoint$Treatment = factor(lmepoint$Treatment,c("PBS","CRA", "DiHOME"))
lmepoint
lmepoint$Sample <- Input$Sample.

#Generate a Bar Graph of the donor-corrected CD3R Data#
CD3R <- ggplot(data=lmepoint, aes(x=Treatment, y=visregRes)) + 
  geom_boxplot()+
  geom_jitter(data=lmepoint, aes(x=Treatment, y=visregRes), position=position_jitter(0.2), size=5, shape=lmepoint$Repeat)+
  theme_classic() + scale_fill_brewer(palette = "Greys", direction=1)+
  guides(fill = guide_legend(reverse=FALSE))
CD3R
ggsave(filename="GRAPH-MouseLungCD3R.pdf", plot=CD3R,  useDingbats=FALSE, width=8.5, height=11, units="in")

#Compared to PBS#
Input$Treatment <- factor(Input$Treatment, c("PBS","CRA","DiHOME"))
model <- lmer(CD3R ~ Treatment + (1|Assay), data = Input)
summary(model)

#Compared to DiHOME#
Input$Treatment <- factor(Input$Treatment, c("DiHOME","PBS","CRA"))
model <- lmer(CD3R ~ Treatment + (1|Assay), data = Input)
summary(model)

