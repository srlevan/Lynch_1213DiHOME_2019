#Clear Envorinment#
rm(list=ls())

#load required packages#
library(lmerTest)
library(visreg)
library(ggplot2)
library(plyr)
library(gridExtra)
library(reshape2)

#Set the Working Directory
setwd("~/Figure2_20190301/")

#Import csv with effects of 12,13 DiHOME-treated DCs on Tcell Subsets#
dihome=read.csv("INPUT-TcellSubsets.csv")
dihome

#Include only Frequency Data#
S3 <- dihome[,c("Donor","Treatment","Live","CD4","Th1","Th2","Th17","Treg")]
S4 <- dihome[,c("Donor","Treatment","CD4Count")]

#Input for Linear Mixed Effects Model#
Input<-S3
Input

#Use LME to model theT cell subsets#
#Treg#
#LME Model of Treg + DONOR#
mod.Treg<- lmer(Treg~Treatment+(1|Donor),data=Input)
y.Treg<-summary(mod.Treg)
y.Treg
coef(mod.Treg)
mod.Treg$res
vTreg<-visreg(fit=mod.Treg,xvar="Treatment",type="conditional")

#summarize each individual point from the model#
resTreg<-vTreg$res
resTreg["Tcell"] <- "Treg"
colnames(resTreg)[3]<-"Freq"
resTreg$DonorShape <- dihome$Donor
resTreg

#Th1#
#LME Model of Th1 + DONOR#
mod.Th1<- lmer(Th1~Treatment+(1|Donor),data=Input)
y.Th1<-summary(mod.Th1)
y.Th1
coef(mod.Th1)
mod.Th1$res
vTh1<-visreg(fit=mod.Th1,xvar="Treatment",type="conditional")

#summarize each individual point from the model#
resTh1<-vTh1$res
resTh1["Tcell"] <- "Th1"
colnames(resTh1)[3]<-"Freq"
resTh1$DonorShape <- na.omit(dihome)$Donor
resTh1

#Th2#
#LME Model of Th2 + DONOR#
mod.Th2<- lmer(Th2~Treatment+(1|Donor),data=Input)
y.Th2<-summary(mod.Th2)
y.Th2
coef(mod.Th2)
mod.Th2$res
vTh2<-visreg(fit=mod.Th2,xvar="Treatment",type="conditional")

#summarize each individual point from the model#
resTh2<-vTh2$res
resTh2["Tcell"] <- "Th2"
colnames(resTh2)[3]<-"Freq"
resTh2$DonorShape <- na.omit(dihome)$Donor
resTh2

#Th17#
#LME Model of Th17 + DONOR#
mod.Th17<- lmer(Th17~Treatment+(1|Donor),data=Input)
y.Th17<-summary(mod.Th17)
y.Th17
coef(mod.Th17)
mod.Th17$res
vTh17<-visreg(fit=mod.Th17,xvar="Treatment",type="conditional")

#summarize each individual point from the model#
resTh17<-vTh17$res
resTh17["Tcell"] <- "Th17"
colnames(resTh17)[3]<-"Freq"
resTh17$DonorShape <- na.omit(dihome)$Donor
resTh17

#Combine Tcell Frequencies (Adjusted for Donor) to Make a Single Dataframe#
Bind1to2 <- rbind(resTh1,resTh2)
Bind12to17 <- rbind(Bind1to2,resTh17)
BindFreq <- rbind(Bind12to17,resTreg)
BindFreq
BindFreq$Tcell = factor(BindFreq$Tcell, c("Th1", "Th2","Th17","Treg"))
levels(BindFreq$DonorShape) <- c(18,16,17,15)

#Save Data Tables#
write.table(BindFreq,file="lme_graphpoints_BindFreq.txt",sep="\t",quote=F,row.names=F)

#read the table of graphing points#
BindFreq<-read.table("lme_graphpoints_BindFreq.txt",header=T,sep="\t")

#Seperate Th1 from Th2,Th17, and Treg#
Th1 <- subset(BindFreq, Tcell=="Th1")
OtherT <- droplevels(subset(BindFreq, Tcell!="Th1"))
OtherT$Tcell = factor(OtherT$Tcell, c("Th2","Th17","Treg"))

#Plot Th1 with boxplot and jitter#
Th1Graph <- ggplot() + 
  geom_boxplot(data=Th1, aes(x=Treatment, y=visregRes, fill=Treatment))+
  geom_jitter(data=Th1, aes(x=Treatment, y=visregRes), shape=Th1$DonorShape, position=position_jitter(0.3), size=8)+
  theme_classic()+scale_fill_brewer(palette = "Greys", direction=1) + 
  theme(legend.position="none")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(title="Th1_SubsetSamples", x="Treatment", y="Th1 Frequency")
Th1Graph

ggsave(filename="GRAPH-DistributedTh1Frequency.pdf", plot=Th1Graph,  useDingbats=FALSE,width=2.75, height=8.5, units="in")


#Plot as OtherMean T with Boxplot and Jitterr#
OtherT

#Plot as OtherT with 95% CI and Scatter#
OtherTGraph <- ggplot() + 
  geom_boxplot(data=OtherT, aes(x=Tcell, y=visregRes, fill=Treatment))+
  geom_jitter(data=OtherT, aes(x=Tcell, y=visregRes, fill=Treatment, shape=as.factor(OtherT$DonorShape)), position=position_dodge(width=0.8), size=8)+
  theme_classic()+scale_fill_brewer(palette = "Greys", direction=1) + 
  theme(legend.position="none")+
  scale_shape_manual(values = c(18,16,15,17))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(title="OtherT_SubsetSamples", x="Treatment", y="OtherT Frequency")
OtherTGraph

ggsave(filename="GRAPH-DistributedOtherTFrequency.pdf", plot=OtherTGraph,  useDingbats=FALSE,width=8.25, height=8.5, units="in")


#LME Model of Live + DONOR#
mod.Live<- lmer(Live~Treatment+(1|Donor),data=Input)
y.Live<-summary(mod.Live)
y.Live
coef(mod.Live)
mod.Live$res
vLive<-visreg(fit=mod.Live,xvar="Treatment",type="conditional")

#summarize each individual point from the model#
resLive<-vLive$res
resLive["Tcell"] <- "Live"
colnames(resLive)[3]<-"Freq"
resLive$DonorShape <- dihome$Donor
resLive

#Graph Live with Box Plot and Jitter#
LivePlot<-ggplot()+
  geom_boxplot(data=resLive, aes(x=Treatment, y=visregRes))+
geom_jitter(data=resLive,aes(x=Treatment, y=visregRes, shape=resLive$DonorShape),position=position_jitter(0.2),size=5)+
  theme_classic()+scale_fill_brewer(palette = "Greys", direction=1) + 
  scale_shape_manual(values = c(18,16,15,17))+
  scale_y_continuous(limits=c(0,100))+
  theme(legend.position="none")
LivePlot
ggsave(filename="GRAPH-LiveCD3Freq.pdf", plot=LivePlot,  useDingbats=FALSE,width=4, height=8.5, units="in")

#Use LME to model the CD4 T cells#
#LME Model of CD4 + DONOR#
mod.CD4<- lmer(CD4~Treatment+(1|Donor),data=Input)
y.CD4<-summary(mod.CD4)
y.CD4
coef(mod.CD4)
mod.CD4$res
vCD4<-visreg(fit=mod.CD4,xvar="Treatment",type="conditional")

#summarize each individual point from the model#
resCD4<-vCD4$res
resCD4["Tcell"] <- "CD4"
colnames(resCD4)[3]<-"Freq"
resCD4$DonorShape <- dihome$Donor
resCD4

#Graph CD4 Freq with Boxplot and Jitter#
limits <- aes(ymax = mean + sd, ymin=resp - sd)
CD4Plot=ggplot()+
  geom_boxplot(data=resCD4, aes(x=Treatment, y=visregRes))+
  geom_jitter(data=resCD4,aes(x=Treatment, y=visregRes, shape=resCD4$DonorShape),position=position_jitter(0.2),size=5)+
  theme_classic()+scale_fill_brewer(palette = "Greys", direction=1) + 
  scale_shape_manual(values = c(18,16,15,17))+
  scale_y_continuous(limits=c(50,80))+
  theme(legend.position="none")
CD4Plot
ggsave(filename="GRAPH-CD4Freq.pdf", plot=CD4Plot,  useDingbats=FALSE,width=4, height=8.5, units="in")

#Use LME to model the CD4 T cell Counts#
#LME Model of CD4 + DONOR#
Input <- S4
mod.CD4Count<- lmer(CD4Count~Treatment+(1|Donor),data=Input)
y.CD4Count<-summary(mod.CD4Count)
y.CD4Count
coef(mod.CD4Count)
mod.CD4Count$res
vCD4Count<-visreg(fit=mod.CD4Count,xvar="Treatment",type="conditional")

#summarize each individual point from the model#
resCD4Count<-vCD4Count$res
resCD4Count["Tcell"] <- "CD4Count"
colnames(resCD4Count)[3]<-"Count"
resCD4Count$DonorShape <- dihome$Donor
resCD4Count

#Graph Boxplot and Jitter#
CD4CountPlot=ggplot()+
  geom_boxplot(data=resCD4Count, aes(x=Treatment, y=visregRes))+
  geom_jitter(data=resCD4Count,aes(x=Treatment, y=visregRes, shape=resCD4Count$DonorShape),position=position_jitter(0.2),size=5)+
  theme_classic()+scale_fill_brewer(palette = "Greys", direction=1) + 
  scale_shape_manual(values = c(18,16,15,17))+
  theme(legend.position="none")
CD4CountPlot
ggsave(filename="GRAPH-CD4Count.pdf", plot=CD4CountPlot,  useDingbats=FALSE,width=4, height=8.5, units="in")

