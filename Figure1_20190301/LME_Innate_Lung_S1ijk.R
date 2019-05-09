#Clear Environment#
rm(list=ls())

#load required packages#
library(lmerTest)
library(visreg)
library(ggplot2)
library(plyr)

#Set Working Directory#
setwd("~/Box Sync/Lynch Lab/Code Respositories/Figure 1_20190301/")

#Load Mouse Lung Innate Cells#
dihome=read.csv("INPUT-Mouse-Innate.csv")

#Run Test#
Input<- dihome
Input$Treatment <- factor(Input$Treatment, c("PBS","CRA","DiHOME"))

#Resident#
#Use LME to Correct for Assay Variability#
mod.Resident<- lmer(Resident~Treatment+(1|Assay),data=Input)
y.Resident<-summary(mod.Resident)
y.Resident
coef(mod.Resident)
mod.Resident$res
vResident<-visreg(fit=mod.Resident,xvar="Treatment",type="conditional")

#summarize each individual point from the model#
resResident<-vResident$res
resResident["Don"] <- Input$Assay
resResident
resResident$Repeat <- ifelse(resResident$Don=="P1A2",15,ifelse(resResident$Don=="P1A3", 16, ifelse(resResident$Don=="P1A4",17, ifelse(resResident$Don=="P1A5", 18, ifelse(resResident$Don=="P1A6",3,0)))))

#Create a table of the points taht will be included in the graph#
write.table(resResident,file="lme_graphpoints_Tcells_Resident.txt",sep="\t",quote=F,row.names=F)
lmepoint<-read.table(file="lme_graphpoints_Tcells_Resident.txt",header=T, sep="\t")

#Factor the table by treatment#
lmepoint$Treatment = factor(lmepoint$Treatment,c("PBS","CRA", "DiHOME"))
lmepoint

#plot the table in a dot plot#
#Generate plot space with lmepoint data#
Resident <- ggplot(data=lmepoint, aes(x=Treatment, y=visregRes)) +
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2), size=5, shape=lmepoint$Repeat)+
  theme_classic() + scale_fill_brewer(palette = "Greys", direction=1)+
  guides(fill = guide_legend(reverse=FALSE))
Resident
ggsave(filename="GRAPH-Resident_MouseLung.pdf", plot=Resident,  useDingbats=FALSE, width=8.5, height=11, units="in")

#Model Statistics#
#Compared to PBS#
Input$Treatment <- factor(Input$Treatment, c("PBS","CRA","DiHOME"))
model <- lmer(Resident ~ Treatment + (1|Assay), data = Input)
summary(model)

#Compared to DiHOME#
Input$Treatment <- factor(Input$Treatment, c("DiHOME","PBS","CRA"))
model <- lmer(Resident ~ Treatment + (1|Assay), data = Input)
summary(model)



#Remove Outlier and Repeat Statistics#
NoOutlier <- droplevels(subset(Input, Input$Sample.!="Lung-15_dc_009_021.fcs" & Input$Sample.!="Lung-15_dc_010_022.fcs"))

#Compared to PBS#
NoOutlier$Treatment <- factor(NoOutlier$Treatment, c("PBS","CRA","DiHOME"))
model <- lmer(Resident ~ Treatment + (1|Assay), data = NoOutlier)
summary(model)

#Compared to DiHOME#
NoOutlier$Treatment <- factor(NoOutlier$Treatment, c("DiHOME","PBS","CRA"))
model <- lmer(Resident ~ Treatment + (1|Assay), data = NoOutlier)
summary(model)



#Live#
#Use LME to Correct for Assay Variability#
mod.Live<- lmer(Live~Treatment+(1|Assay),data=Input)
y.Live<-summary(mod.Live)
y.Live
coef(mod.Live)
mod.Live$res
vLive<-visreg(fit=mod.Live,xvar="Treatment",type="conditional")

#summarize each individual point from the model#
resLive<-vLive$res
resLive["Don"] <- Input$Assay
resLive
resLive$Repeat <- ifelse(resLive$Don=="P1A2",15,ifelse(resLive$Don=="P1A3", 16, ifelse(resLive$Don=="P1A4",17, ifelse(resLive$Don=="P1A5", 18, ifelse(resLive$Don=="P1A6",3,0)))))

#Create a table of the points taht will be included in the graph#
write.table(resLive,file="lme_graphpoints_Tcells_Live.txt",sep="\t",quote=F,row.names=F)

#read the table of graphing points#
#Graph Live#
lmepoint<-read.table("lme_graphpoints_Tcells_Live.txt",header=T,sep="\t")

#Factor the table by treatment#
lmepoint$Treatment = factor(lmepoint$Treatment,c("PBS","CRA", "DiHOME"))
lmepoint

#plot the table in a dot plot#
#Generate plot space with lmepoint data#
Live <- ggplot(data=lmepoint, aes(x=Treatment, y=visregRes)) + 
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2), size=5, shape=lmepoint$Repeat)+
  theme_classic() + scale_fill_brewer(palette = "Greys", direction=1)+
  guides(fill = guide_legend(reverse=FALSE))
Live
ggsave(filename="GRAPH-Live_MouseLung.pdf", plot=Live,  useDingbats=FALSE, width=8.5, height=11, units="in")

#Compared to PBS#
Input$Treatment <- factor(Input$Treatment, c("PBS","CRA","DiHOME"))
model <- lmer(Live ~ Treatment + (1|Assay), data = Input)
summary(model)

#Compared to DiHOME#
Input$Treatment <- factor(Input$Treatment, c("DiHOME","PBS","CRA"))
model <- lmer(Live ~ Treatment + (1|Assay), data = Input)
summary(model)



#Remove Outlier and Repeat Statistics#
NoOutlier <- droplevels(subset(Input, Input$Sample.!="Lung-15_dc_009_021.fcs" & Input$Sample.!="Lung-15_dc_010_022.fcs"))

#Compared to PBS#
NoOutlier$Treatment <- factor(NoOutlier$Treatment, c("PBS","CRA","DiHOME"))
model <- lmer(Live ~ Treatment + (1|Assay), data = NoOutlier)
summary(model)

#Compared to DiHOME#
NoOutlier$Treatment <- factor(NoOutlier$Treatment, c("DiHOME","PBS","CRA"))
model <- lmer(Live ~ Treatment + (1|Assay), data = NoOutlier)
summary(model)


#Neutro#
#Use LME to Correct for Assay Variability#
#Drop Outliers from P1A4#
Input$Treatment <- factor(Input$Treatment, c("PBS","CRA","DiHOME"))

mod.Neutro<- lmer(Neutro~Treatment+(1|Assay),data=Input)
y.Neutro<-summary(mod.Neutro)
y.Neutro
coef(mod.Neutro)
mod.Neutro$res
vNeutro<-visreg(fit=mod.Neutro,xvar="Treatment",type="conditional")

#summarize each individual point from the model#
resNeutro<-vNeutro$res
resNeutro["Don"] <- Input$Assay
resNeutro
resNeutro$Repeat <- ifelse(resNeutro$Don=="P1A2",15,ifelse(resNeutro$Don=="P1A3", 16, ifelse(resNeutro$Don=="P1A4",17, ifelse(resNeutro$Don=="P1A5", 18, ifelse(resNeutro$Don=="P1A6",3,0)))))


#Create a table of the points taht will be included in the graph#
write.table(resNeutro,file="lme_graphpoints_Tcells_Neutro.txt",sep="\t",quote=F,row.names=F)

#read the table of graphing points#
#Graph Neutro#
lmepoint<-read.table("lme_graphpoints_Tcells_Neutro.txt",header=T,sep="\t")
#Factor the table by treatment#
lmepoint$Treatment = factor(lmepoint$Treatment,c("PBS","CRA", "DiHOME"))
lmepoint


#plot the table in a dot plot#
#Generate plot space with lmepoint data#
Neutro <- ggplot(data=lmepoint, aes(x=Treatment, y=visregRes)) +
  geom_boxplot()+
  geom_jitter(data=lmepoint, aes(x=Treatment, y=visregRes), position=position_jitter(0.2), size=5, shape=lmepoint$Repeat)+
  theme_classic() + scale_fill_brewer(palette = "Greys", direction=1)+
  guides(fill = guide_legend(reverse=FALSE))
Neutro
ggsave(filename="GRAPH-Neutro-MouseLung.pdf", plot=Neutro,  useDingbats=FALSE, width=8.5, height=11, units="in")

#Compared to PBS#
Input$Treatment <- factor(Input$Treatment, c("PBS","CRA","DiHOME"))
model <- lmer(Neutro ~ Treatment + (1|Assay), data = Input)
summary(model)

#Compared to DiHOME#
Input$Treatment <- factor(Input$Treatment, c("DiHOME","PBS","CRA"))
model <- lmer(Neutro ~ Treatment + (1|Assay), data = Input)
summary(model)



#Remove Outlier and Repeat Statistics#
NoOutlier <- droplevels(subset(Input, Input$Sample.!="Lung-15_dc_012_024.fcs"))

#Compared to PBS#
NoOutlier$Treatment <- factor(NoOutlier$Treatment, c("PBS","CRA","DiHOME"))
model <- lmer(Neutro ~ Treatment + (1|Assay), data = NoOutlier)
summary(model)

#Compared to DiHOME#
NoOutlier$Treatment <- factor(NoOutlier$Treatment, c("DiHOME","PBS","CRA"))
model <- lmer(Neutro ~ Treatment + (1|Assay), data = NoOutlier)
summary(model)


#Eos#
#Use LME to Correct for Assay Variability#
Input <- dihome
Input$Treatment <- factor(Input$Treatment, c("PBS","CRA","DiHOME"))

mod.Eos<- lmer(Eos~Treatment+(1|Assay),data=Input)
y.Eos<-summary(mod.Eos)
y.Eos
coef(mod.Eos)
mod.Eos$res
vEos<-visreg(fit=mod.Eos,xvar="Treatment",type="conditional")

#summarize each individual point from the model#
resEos<-vEos$res
resEos["Don"] <- Input$Assay
resEos
resEos$Repeat <- ifelse(resEos$Don=="P1A2",15,ifelse(resEos$Don=="P1A3", 16, ifelse(resEos$Don=="P1A4",17, ifelse(resEos$Don=="P1A5", 18, ifelse(resEos$Don=="P1A6",3,0)))))


#Create a table of the points taht will be included in the graph#
write.table(resEos,file="lme_graphpoints_Tcells_Eos.txt",sep="\t",quote=F,row.names=F)

#Graph Eos#
lmepoint<-read.table("lme_graphpoints_Tcells_Eos.txt",header=T,sep="\t")

#Factor the table by treatment#
lmepoint$Treatment = factor(lmepoint$Treatment,c("PBS","CRA", "DiHOME"))
lmepoint

#plot the table in a dot plot#
#Generate plot space with lmepoint data#
Eos <- ggplot(data=lmepoint, aes(x=Treatment, y=visregRes)) + 
  geom_boxplot()+
  geom_jitter(data=lmepoint, aes(x=Treatment, y=visregRes), position=position_jitter(0.2), size=5, shape=lmepoint$Repeat)+
  theme_classic() + scale_fill_brewer(palette = "Greys", direction=1)+
  guides(fill = guide_legend(reverse=FALSE))
Eos
ggsave(filename="GRAPH-Eos_MouseLung.pdf", plot=Eos,  useDingbats=FALSE, width=8.5, height=11, units="in")


#Compared to PBS#
Input$Treatment <- factor(Input$Treatment, c("PBS","CRA","DiHOME"))
model <- lmer(Eos ~ Treatment + (1|Assay), data = Input)
summary(model)

#Compared to DiHOME#
Input$Treatment <- factor(Input$Treatment, c("DiHOME","PBS","CRA"))
model <- lmer(Eos ~ Treatment + (1|Assay), data = Input)
summary(model)


#Mono#
#Use LME to Correct for Assay Variability#
Input <- dihome
Input$Treatment <- factor(Input$Treatment, c("PBS","CRA","DiHOME"))

Input$Treatment <- factor(Input$Treatment, c("PBS","CRA","DiHOME"))
mod.Mono<- lmer(Mono~Treatment+(1|Assay),data=Input)
y.Mono<-summary(mod.Mono)
y.Mono
coef(mod.Mono)
mod.Mono$res
vMono<-visreg(fit=mod.Mono,xvar="Treatment",type="conditional")

#summarize each individual point from the model#
resMono<-vMono$res
resMono["Don"] <- Input$Assay
resMono
resMono$Repeat <- ifelse(resMono$Don=="P1A2",15,ifelse(resMono$Don=="P1A3", 16, ifelse(resMono$Don=="P1A4",17, ifelse(resMono$Don=="P1A5", 18, ifelse(resMono$Don=="P1A6",3,0)))))


#Create a table of the points taht will be included in the graph#
write.table(resMono,file="lme_graphpoints_Tcells_Mono.txt",sep="\t",quote=F,row.names=F)

#Graph Mono#
lmepoint<-read.table("lme_graphpoints_Tcells_Mono.txt",header=T,sep="\t")
#Factor the table by treatment#
lmepoint$Treatment = factor(lmepoint$Treatment,c("PBS","CRA", "DiHOME"))
lmepoint

#plot the table in a dot plot#
#Generate plot space with lmepoint data#
Mono <- ggplot(data=lmepoint, aes(x=Treatment, y=visregRes)) + 
  geom_boxplot()+
  geom_jitter(data=lmepoint, aes(x=Treatment, y=visregRes), position=position_jitter(0.2), size=5, shape=lmepoint$Repeat)+
  theme_classic() + scale_fill_brewer(palette = "Greys", direction=1)+
  guides(fill = guide_legend(reverse=FALSE))
Mono
ggsave(filename="GRAPH-Mono_MouseLung.pdf", plot=Mono,  useDingbats=FALSE, width=8.5, height=11, units="in")

#Compared to PBS#
Input$Treatment <- factor(Input$Treatment, c("PBS","CRA","DiHOME"))
model <- lmer(Mono ~ Treatment + (1|Assay), data = Input)
summary(model)

#Compared to DiHOME#
Input$Treatment <- factor(Input$Treatment, c("DiHOME","PBS","CRA"))
model <- lmer(Mono ~ Treatment + (1|Assay), data = Input)
summary(model)



#Remove Outlier and Repeat Statistics#
NoOutlier <- droplevels(subset(Input, Input$Sample.!="Lung-15_dc_012_024.fcs"))

#Compared to PBS#
NoOutlier$Treatment <- factor(NoOutlier$Treatment, c("PBS","CRA","DiHOME"))
model <- lmer(Mono ~ Treatment + (1|Assay), data = NoOutlier)
summary(model)

#Compared to DiHOME#
NoOutlier$Treatment <- factor(NoOutlier$Treatment, c("DiHOME","PBS","CRA"))
model <- lmer(Mono ~ Treatment + (1|Assay), data = NoOutlier)
summary(model)