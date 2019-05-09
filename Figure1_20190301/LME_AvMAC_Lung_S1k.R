#Clear the Environment#
rm(list=ls())

#load required packages#
library(lmerTest)
library(visreg)
library(ggplot2)
library(plyr)

#Set Working Directory#
setwd("~/Box Sync/Lynch Lab/Code Respositories/Figure 1_20190301/")

#Read in AvMac Data#
dihome=read.csv("INPUT-Mouse-AvMAC.csv")
dihome

#Graph LME of AvMACs in Bar Graph#
dihome$Treatment <- factor(dihome$Treatment, c("PBS","CRA","DiHOME"))
Input<- dihome


#FreqLive#
#model Competition + Titration_FreqLive#
mod.FreqLive<- lmer(FreqLive~Treatment+(1|Assay),data=Input)
y.FreqLive<-summary(mod.FreqLive)
y.FreqLive
coef(mod.FreqLive)
mod.FreqLive$res
vFreqLive<-visreg(fit=mod.FreqLive,xvar="Treatment",type="conditional")

#summarize each individual point from the model#
resFreqLive<-vFreqLive$res
resFreqLive["Don"] <- Input$Assay
resFreqLive
resFreqLive$Repeat <- ifelse(resFreqLive$Don=="P1A2",15,ifelse(resFreqLive$Don=="P1A3", 16, ifelse(resFreqLive$Don=="P1A4",17, ifelse(resFreqLive$Don=="P1A5", 18, ifelse(resFreqLive$Don=="P1A6",3,0)))))

#Create a table of the points that will be included in the graph#
write.table(resFreqLive,file="lme_graphpoints_AvMAC_FreqLive.txt",sep="\t",quote=F,row.names=F)

#read the table of graphing points#
lmepoint<-read.table("lme_graphpoints_AvMAC_FreqLive.txt",header=T,sep="\t")

#Factor the table by treatment#
lmepoint$Treatment = factor(lmepoint$Treatment,c("PBS","CRA", "DiHOME"))
lmepoint

#plot the table in a dot plot#
#Generate plot space with lmepoint data#
FreqLive <- ggplot(data=lmepoint, aes(x=Treatment, y=visregRes)) + 
  geom_boxplot()+
  geom_jitter(data=lmepoint, aes(x=Treatment, y=visregRes), position=position_jitter(0.2), size=5, shape=lmepoint$Repeat)+
  theme_classic() + scale_fill_brewer(palette = "Greys", direction=1)+
  guides(fill = guide_legend(reverse=FALSE))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(title="FreqLive_SubsetSamples", x="Treatment", y="FreqLive Frequency")
FreqLive
ggsave(filename="GRAPH-AvMAC_MouseLung.pdf", plot=FreqLive,  useDingbats=FALSE, width=8.5, height=11, units="in")

#Compared to PBS#
Input$Treatment <- factor(Input$Treatment, c("PBS","CRA","DiHOME"))
model <- lmer(FreqLive ~ Treatment + (1|Assay), data = Input)
summary(model)

#Compared to DiHOME#
Input$Treatment <- factor(Input$Treatment, c("DiHOME","PBS","CRA"))
model <- lmer(FreqLive ~ Treatment + (1|Assay), data = Input)
summary(model)
