#Clear the Environment#
rm(list=ls())

#Load Required Packages#
library(ggplot2)
library(visreg)
library(plyr)
library(gridExtra)
library(reshape2)
library(lmerTest)

#Set WD#
setwd("~/Figure2_20190301/")

#Read in qPCR Data#
qPCR=read.csv("INPUT-qPCRDC_Updated20190103.csv")
head(qPCR)
levels(qPCR$Target.Name)

#Subset by Gene#
CD36 <- subset(qPCR, Target.Name =="CD36")
CD36$Treatment <- factor(CD36$Treatment, c("Vehicle", "DiHOME"))

CD1a <- droplevels(subset(qPCR, Target.Name == "CD1a"))
CD1a$Treatment <- factor(CD1a$Treatment, c("Vehicle", "DiHOME"))

HADH <- droplevels(subset(qPCR, Target.Name == "HADH"))
HADH$Treatment <- factor(HADH$Treatment, c("Vehicle", "DiHOME"))

FABP4 <- droplevels(subset(qPCR, Target.Name == "FABP4"))
FABP4$Treatment <- factor(FABP4$Treatment, c("Vehicle", "DiHOME"))

#Select Gene for Analysis#
Input<-FABP4
Input

#Model and Graph Dotplot of ddCT for the Gene of interest#
mod.ddCT<- lmer(ddCT~Treatment+(1|Donor),data=Input)
y.ddCT<-summary(mod.ddCT)
y.ddCT
coef(mod.ddCT)
mod.ddCT$res
vddCT<-visreg(fit=mod.ddCT,xvar="Treatment",type="conditional")

#summarize each individual point from the model#
resddCT<-vddCT$res
resddCT["Don"] <- as.factor(Input$Assay)
resddCT$Treatment
levels(resddCT$Don) <- c(18,16,17,15)
resddCT
#Add a column with Fold Change#
resddCT$FoldChange <- 2^(-resddCT$visregRes)
#Simplify Table#
SimpleddCT <- resddCT[,c(1,4,6,7)]
names(SimpleddCT)[2]<-paste("ddCT")
names(SimpleddCT)[3]<-paste("Donor")

#Create a table of the points taht will be included in the graph#
write.table(SimpleddCT,file="lme_graphpoints_ddCT_FABP4.txt",sep="\t",quote=F,row.names=F)

#read the table of graphing points#
lmepoint<-read.table("lme_graphpoints_ddCT_FABP4.txt",header=T,sep="\t")
#Factor the table by treatment#
lmepoint$Treatment = factor(lmepoint$Treatment,c("Vehicle","DiHOME"))
lmepoint
#plot the table in a dot plot#
FoldChange <- ggplot(lmepoint, aes(x=Treatment, y=FoldChange)) + 
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2), size=8, shape=lmepoint$Donor) +
  stat_summary(data=lmepoint, aes(ymax=..y..,ymin=..y.., line = factor(Treatment)), fun.y = "median", geom = "crossbar", colour = "black", width = 0.5)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(title="FABP4 Fold Change", x="Treatment", y="FoldChange")
FoldChange
ggsave(filename="GRAPH-FABP4_FoldChange.pdf", plot=ddCT,  useDingbats=FALSE,width=8.5, height=11, units="in")

#Select Gene for Analysis#
Input<-HADH
Input

#Model and Graph Dotplot of ddCT for the Gene of interest#
mod.ddCT<- lmer(ddCT~Treatment+(1|Donor),data=Input)
y.ddCT<-summary(mod.ddCT)
y.ddCT
coef(mod.ddCT)
mod.ddCT$res
vddCT<-visreg(fit=mod.ddCT,xvar="Treatment",type="conditional")

#summarize each individual point from the model#
resddCT<-vddCT$res
resddCT["Don"] <- as.factor(Input$Assay)
resddCT$Treatment
levels(resddCT$Don) <- c(18,16,17,15)
resddCT
#Add a column with Fold Change#
resddCT$FoldChange <- 2^(-resddCT$visregRes)
#Simplify Table#
SimpleddCT <- resddCT[,c(1,4,6,7)]
names(SimpleddCT)[2]<-paste("ddCT")
names(SimpleddCT)[3]<-paste("Donor")

#Create a table of the points taht will be included in the graph#
write.table(SimpleddCT,file="lme_graphpoints_ddCT_HADH.txt",sep="\t",quote=F,row.names=F)

#read the table of graphing points#
lmepoint<-read.table("lme_graphpoints_ddCT_HADH.txt",header=T,sep="\t")
#Factor the table by treatment#
lmepoint$Treatment = factor(lmepoint$Treatment,c("Vehicle","DiHOME"))
lmepoint
#plot the table in a dot plot#
FoldChange <- ggplot(lmepoint, aes(x=Treatment, y=FoldChange)) + 
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2), size=8, shape=lmepoint$Donor) +
  stat_summary(data=lmepoint, aes(ymax=..y..,ymin=..y.., line = factor(Treatment)), fun.y = "median", geom = "crossbar", colour = "black", width = 0.5)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(title="HADH Fold Change", x="Treatment", y="FoldChange")
FoldChange
ggsave(filename="GRAPH-HADH_FoldChange.pdf", plot=ddCT,  useDingbats=FALSE,width=8.5, height=11, units="in")


#Select Gene for Analysis#
Input<-CD36
Input

#Model and Graph Dotplot of ddCT for the Gene of interest#
mod.ddCT<- lmer(ddCT~Treatment+(1|Donor),data=Input)
y.ddCT<-summary(mod.ddCT)
y.ddCT
coef(mod.ddCT)
mod.ddCT$res
vddCT<-visreg(fit=mod.ddCT,xvar="Treatment",type="conditional")
levels(Input$Assay)
#summarize each individual point from the model#
resddCT<-vddCT$res
resddCT["Don"] <- as.factor(Input$Donor)
resddCT$Treatment
levels(resddCT$Don) <- c(18,16,17,15)
resddCT
#Add a column with Fold Change#
resddCT$FoldChange <- 2^(-resddCT$visregRes)
#Simplify Table#
SimpleddCT <- resddCT[,c(1,4,6,7)]
names(SimpleddCT)[2]<-paste("ddCT")
names(SimpleddCT)[3]<-paste("Donor")

#Create a table of the points taht will be included in the graph#
write.table(SimpleddCT,file="lme_graphpoints_ddCT_CD36.txt",sep="\t",quote=F,row.names=F)

#read the table of graphing points#
lmepoint<-read.table("lme_graphpoints_ddCT_CD36.txt",header=T,sep="\t")
#Factor the table by treatment#
lmepoint$Treatment = factor(lmepoint$Treatment,c("Vehicle","DiHOME"))
lmepoint
#plot the table in a dot plot#
FoldChange <- ggplot(lmepoint, aes(x=Treatment, y=FoldChange)) + 
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2), size=8, shape=lmepoint$Donor) +
  stat_summary(data=lmepoint, aes(ymax=..y..,ymin=..y.., line = factor(Treatment)), fun.y = "median", geom = "crossbar", colour = "black", width = 0.5)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(title="CD36 Fold Change", x="Treatment", y="FoldChange")
FoldChange
ggsave(filename="GRAPH-CD36_FoldChange.pdf", plot=ddCT,  useDingbats=FALSE,width=8.5, height=11, units="in")

#Select Gene for Analysis#
Input<-CD1a
Input

#Model and Graph Dotplot of ddCT for the Gene of interest#
mod.ddCT<- lmer(ddCT~Treatment+(1|Donor),data=Input)
y.ddCT<-summary(mod.ddCT)
y.ddCT
coef(mod.ddCT)
mod.ddCT$res
vddCT<-visreg(fit=mod.ddCT,xvar="Treatment",type="conditional")

#summarize each individual point from the model#
resddCT<-vddCT$res
resddCT["Don"] <- as.factor(Input$Assay)
resddCT$Treatment
levels(resddCT$Don) <- c(18,16,17,15)
resddCT
#Add a column with Fold Change#
resddCT$FoldChange <- 2^(-resddCT$visregRes)
#Simplify Table#
SimpleddCT <- resddCT[,c(1,4,6,7)]
names(SimpleddCT)[2]<-paste("ddCT")
names(SimpleddCT)[3]<-paste("Donor")

#Create a table of the points taht will be included in the graph#
write.table(SimpleddCT,file="lme_graphpoints_ddCT_CD1a.txt",sep="\t",quote=F,row.names=F)

#read the table of graphing points#
lmepoint<-read.table("lme_graphpoints_ddCT_CD1a.txt",header=T,sep="\t")
#Factor the table by treatment#
lmepoint$Treatment = factor(lmepoint$Treatment,c("Vehicle","DiHOME"))
lmepoint
#plot the table in a dot plot#
FoldChange <- ggplot(lmepoint, aes(x=Treatment, y=FoldChange)) + 
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2), size=8, shape=lmepoint$Donor) +
  stat_summary(data=lmepoint, aes(ymax=..y..,ymin=..y.., line = factor(Treatment)), fun.y = "median", geom = "crossbar", colour = "black", width = 0.5)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(title="CD1a Fold Change", x="Treatment", y="FoldChange")
FoldChange
ggsave(filename="GRAPH-CD1a_FoldChange.pdf", plot=ddCT,  useDingbats=FALSE,width=8.5, height=11, units="in")


#Read in Modeled Fold Changes#
CD1aFC=read.table("lme_graphpoints_ddCT_CD1a.txt", sep="\t", header=TRUE)
CD1aFC$Target.Name <- "CD1a"
CD36FC=read.table("lme_graphpoints_ddCT_CD36.txt", sep="\t", header=TRUE)
CD36FC$Target.Name <- "CD36"
HADHFC=read.table("lme_graphpoints_ddCT_HADH.txt", sep="\t", header=TRUE)
HADHFC$Target.Name <- "HADH"
FABP4FC=read.table("lme_graphpoints_ddCT_FABP4.txt", sep="\t", header=TRUE)
FABP4FC$Target.Name <- "FABP4"

#Make a Single Data Frame with All Modelled Fold Changes#
FC <- rbind(CD1aFC, CD36FC, FABP4FC, HADHFC)
FC
FC1 <- FC[,c(1,3,4,5)]
FC1


#Convert to Long Form#
FoldChange <- melt(FC1, id.vars=c("Treatment","Target.Name", "Donor"))
FoldChange
FoldChange$Treatment <-factor(FoldChange$Treatment, c("Vehicle","DiHOME"))


#Unable to reproduce CD36 data in additional donors#

#Plot bar Graph that includes CD1a, FABP4, HADH#
Points <- droplevels(subset(FoldChange, Target.Name!="CD36"))
write.table(Points,file="lme_graphpoints_FC.txt",sep="\t",quote=F,row.names=F)
#read the table of graphing points#
Points<-read.table("lme_graphpoints_FC.txt",header=T,sep="\t")
Points$Donor<-as.factor(Points$Donor)
Points$Treatment <-factor(Points$Treatment, c("Vehicle","DiHOME"))


Input <-Plot
qPCR=ggplot()+
  geom_boxplot(data=Points, aes(x=Target.Name, y=value, fill=Treatment))+
  geom_jitter(data=Points, aes(x=Target.Name, y=value, fill=Treatment, shape=Points$Donor),position=position_dodge(width=0.8),size=5)+
  scale_shape_manual(values = c(18,16,15,17))+
  theme_classic()+scale_fill_brewer(palette = "Greys", direction=1) + 
  theme(legend.position="none")+
  scale_y_log10(limits=c(0.0015,30),breaks=c(0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,20,30))
qPCR

ggsave(filename="GRAPH-FC.pdf", plot=qPCR,  useDingbats=FALSE,width=10, height=11, units="in")
