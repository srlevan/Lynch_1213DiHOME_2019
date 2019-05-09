#test for expression changes in Mouse Lung# 
rm(list=ls())

#Load Required Packages#
library(ggplot2)
library(visreg)
library(plyr)
library(lmerTest)


#Set WD#
setwd("~/Box Sync/Lynch Lab/Code Respositories/Figure 1_20190301/")

#Read in qPCR Data#
qPCR=read.csv("INPUT_MouseLungQPCR.csv", na.strings="NA")
qPCR

Target.Name <- as.factor(levels(qPCR$Target.Name))
Target.Name

#Subset by Gene#
IL1a <- subset(qPCR, Target.Name == "IL1a")
IL1b <- subset(qPCR, Target.Name == "IL1b")
TNF <- droplevels(subset(qPCR, Target.Name == "TNF"))
CCL11 <- droplevels(subset(qPCR, Target.Name == "CCL11"))
Muc5ac <- droplevels(subset(qPCR, Target.Name == "Muc5ac"))

#Select Gene for Analysis#
#Analysis of IL1a#
Input<-IL1a
Input$Group <- factor(Input$Group, c("PBS","CRA","DiHOME"))
Input

#Model and Graph Dotplot of ddCT for the Gene of interest#
mod.ddCT<- lmer(ddCT~Group+(1|Assay),data=Input)
y.ddCT<-summary(mod.ddCT)
y.ddCT
coef(mod.ddCT)
mod.ddCT$res
vddCT<-visreg(fit=mod.ddCT,xvar="Group",type="conditional")

#summarize each individual point from the model#
resddCT<-vddCT$res
resddCT["Don"] <- as.factor(Input$Assay)
resddCT$Group
resddCT$Repeat <- ifelse(resddCT$Don=="P1A2",15,ifelse(resddCT$Don=="P1A3", 16, ifelse(resddCT$Don=="P1A4",17, ifelse(resddCT$Don=="P1A5", 18, ifelse(resddCT$Don=="P1A6",3,0)))))
resddCT

#Add a column with Fold Change#
resddCT$FoldChange <- 2^(-resddCT$visregRes)

#Add a column with Sample ID#
resddCT$Sample  <- Input$Sample.Name
#Simplify Table#
SimpleFC <- resddCT[,c("Group","Assay","Repeat","Sample","FoldChange")]

#Create a table of the points taht will be included in the graph#
write.table(SimpleFC,file="lme_graphpoints_ddCT_IL1a.txt",sep="\t",quote=F,row.names=F)

#read the table of graphing points#
lmepoint<-read.table("lme_graphpoints_ddCT_IL1a.txt",header=T,sep="\t")

#Factor the table by Group#
lmepoint$Group = factor(lmepoint$Group,c("PBS","CRA","DiHOME"))
lmepoint

#plot the table in a dot plot#
FoldChange <- ggplot(lmepoint, aes(x=Group, y=FoldChange)) + 
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2), size=8, shape=lmepoint$Repeat) +
  stat_summary(data=lmepoint, aes(ymax=..y..,ymin=..y.., line = factor(Group)), fun.y = "median", geom = "crossbar", colour = "black", width = 0.5)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(title="IL1a Fold Change", x="Group", y="FoldChange")
FoldChange
ggsave(filename="GRAPH-IL1a_FoldChange.pdf", plot=FoldChange,  useDingbats=FALSE,width=8.5, height=11, units="in")

#Compared to PBS#
Input$Group <- factor(Input$Group, c("PBS","CRA","DiHOME"))
model <- lmer(ddCT ~ Group + (1|Assay), data = Input)
summary(model)

#Compared to DiHOME#
Input$Group <- factor(Input$Group, c("DiHOME","PBS","CRA"))
model <- lmer(ddCT ~ Group + (1|Assay), data = Input)
summary(model)

#Remove Outlier and Repeat Statistics#
lmepoint
NoOutlier <- droplevels(subset(Input, Input$Sample.Name!="4-D1"))

#Compared to PBS#
NoOutlier$Group <- factor(NoOutlier$Group, c("PBS","CRA","DiHOME"))
model <- lmer(ddCT ~ Group + (1|Assay), data = NoOutlier)
summary(model)

#Compared to DiHOME#
NoOutlier$Group <- factor(NoOutlier$Group, c("DiHOME","PBS","CRA"))
model <- lmer(ddCT ~ Group + (1|Assay), data = NoOutlier)
summary(model)

#Analysis of IL1b#
#Select Gene for Analysis#
Input<-IL1b
Input$Group <- factor(Input$Group, c("PBS","CRA","DiHOME"))
Input

#Model and Graph Dotplot of ddCT for the Gene of interest#
mod.ddCT<- lmer(ddCT~Group+(1|Assay),data=Input)
y.ddCT<-summary(mod.ddCT)
y.ddCT
coef(mod.ddCT)
mod.ddCT$res
vddCT<-visreg(fit=mod.ddCT,xvar="Group",type="conditional")

#summarize each individual point from the model#
resddCT<-vddCT$res
resddCT["Don"] <- as.factor(Input$Assay)
resddCT$Group
resddCT$Repeat <- ifelse(resddCT$Don=="P1A2",15,ifelse(resddCT$Don=="P1A3", 16, ifelse(resddCT$Don=="P1A4",17, ifelse(resddCT$Don=="P1A5", 18, ifelse(resddCT$Don=="P1A6",3,0)))))
resddCT

#Add a column with Fold Change#
resddCT$FoldChange <- 2^(-resddCT$visregRes)

#Add a column with Sample ID#
resddCT$Sample  <- Input$Sample.Name
#Simplify Table#
SimpleFC <- resddCT[,c("Group","Assay","Repeat","Sample","FoldChange")]

#Create a table of the points taht will be included in the graph#
write.table(SimpleFC,file="lme_graphpoints_ddCT_IL1b.txt",sep="\t",quote=F,row.names=F)

#read the table of graphing points#
lmepoint<-read.table("lme_graphpoints_ddCT_IL1b.txt",header=T,sep="\t")

#Factor the table by Group#
lmepoint$Group = factor(lmepoint$Group,c("PBS","CRA","DiHOME"))
lmepoint

#plot the table in a dot plot#
FoldChange <- ggplot(lmepoint, aes(x=Group, y=FoldChange)) + 
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2), size=8, shape=lmepoint$Repeat) +
  stat_summary(data=lmepoint, aes(ymax=..y..,ymin=..y.., line = factor(Group)), fun.y = "median", geom = "crossbar", colour = "black", width = 0.5)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(title="IL1b Fold Change", x="Group", y="FoldChange")
FoldChange
ggsave(filename="GRAPH-IL1b_FoldChange.pdf", plot=FoldChange,  useDingbats=FALSE,width=8.5, height=11, units="in")

#Compared to PBS#
Input$Group <- factor(Input$Group, c("PBS","CRA","DiHOME"))
model <- lmer(ddCT ~ Group + (1|Assay), data = Input)
summary(model)

#Compared to DiHOME#
Input$Group <- factor(Input$Group, c("DiHOME","PBS","CRA"))
model <- lmer(ddCT ~ Group + (1|Assay), data = Input)
summary(model)


#Analysis of TNF#
#Select Gene for Analysis#
Input<-droplevels(subset(TNF, TNF$Sample.Name!="5-D3"))

Input$Group <- factor(Input$Group, c("PBS","CRA","DiHOME"))
Input

#Model and Graph Dotplot of ddCT for the Gene of interest#
mod.ddCT<- lmer(ddCT~Group+(1|Assay),data=Input)
y.ddCT<-summary(mod.ddCT)
y.ddCT
coef(mod.ddCT)
mod.ddCT$res
vddCT<-visreg(fit=mod.ddCT,xvar="Group",type="conditional")

#summarize each individual point from the model#
resddCT<-vddCT$res
resddCT["Don"] <- as.factor(Input$Assay)
resddCT$Group
resddCT$Repeat <- ifelse(resddCT$Don=="P1A2",15,ifelse(resddCT$Don=="P1A3", 16, ifelse(resddCT$Don=="P1A4",17, ifelse(resddCT$Don=="P1A5", 18, ifelse(resddCT$Don=="P1A6",3,0)))))
resddCT

#Add a column with Fold Change#
resddCT$FoldChange <- 2^(-resddCT$visregRes)

#Add a column with Sample ID#
resddCT$Sample  <- Input$Sample.Name
#Simplify Table#
SimpleFC <- resddCT[,c("Group","Assay","Repeat","Sample","FoldChange")]

#Create a table of the points taht will be included in the graph#
write.table(SimpleFC,file="lme_graphpoints_ddCT_TNF.txt",sep="\t",quote=F,row.names=F)

#read the table of graphing points#
lmepoint<-read.table("lme_graphpoints_ddCT_TNF.txt",header=T,sep="\t")

#Factor the table by Group#
lmepoint$Group = factor(lmepoint$Group,c("PBS","CRA","DiHOME"))
lmepoint

#plot the table in a dot plot#
FoldChange <- ggplot(lmepoint, aes(x=Group, y=FoldChange)) + 
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2), size=8, shape=lmepoint$Repeat) +
  stat_summary(data=lmepoint, aes(ymax=..y..,ymin=..y.., line = factor(Group)), fun.y = "median", geom = "crossbar", colour = "black", width = 0.5)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(title="TNF Fold Change", x="Group", y="FoldChange")
FoldChange
ggsave(filename="GRAPH-TNF_FoldChange.pdf", plot=FoldChange,  useDingbats=FALSE,width=8.5, height=11, units="in")

#Compared to PBS#
Input$Group <- factor(Input$Group, c("PBS","CRA","DiHOME"))
model <- lmer(ddCT ~ Group + (1|Assay), data = Input)
summary(model)

#Compared to DiHOME#
Input$Group <- factor(Input$Group, c("DiHOME","PBS","CRA"))
model <- lmer(ddCT ~ Group + (1|Assay), data = Input)
summary(model)

#Remove Outlier and Repeat Statistics#
lmepoint
NoOutlier <- droplevels(subset(Input, Input$Sample.Name!="4-D2" & Input$Sample.Name!="4-D3" & Input$Sample.Name!="4-P3"))

#Compared to PBS#
NoOutlier$Group <- factor(NoOutlier$Group, c("PBS","CRA","DiHOME"))
model <- lmer(ddCT ~ Group + (1|Assay), data = NoOutlier)
summary(model)

#Compared to DiHOME#
NoOutlier$Group <- factor(NoOutlier$Group, c("DiHOME","PBS","CRA"))
model <- lmer(ddCT ~ Group + (1|Assay), data = NoOutlier)
summary(model)



#Select Gene for Analysis#
#Analysis of Muc5ac#
Input<-Muc5ac
Input$Group <- factor(Input$Group, c("PBS","CRA","DiHOME"))
Input

#Model and Graph Dotplot of ddCT for the Gene of interest#
mod.ddCT<- lmer(ddCT~Group+(1|Assay),data=Input)
y.ddCT<-summary(mod.ddCT)
y.ddCT
coef(mod.ddCT)
mod.ddCT$res
vddCT<-visreg(fit=mod.ddCT,xvar="Group",type="conditional")

#summarize each individual point from the model#
resddCT<-vddCT$res
resddCT["Don"] <- as.factor(Input$Assay)
resddCT$Group
resddCT$Repeat <- ifelse(resddCT$Don=="P1A2",15,ifelse(resddCT$Don=="P1A3", 16, ifelse(resddCT$Don=="P1A4",17, ifelse(resddCT$Don=="P1A5", 18, ifelse(resddCT$Don=="P1A6",3,0)))))
resddCT

#Add a column with Fold Change#
resddCT$FoldChange <- 2^(-resddCT$visregRes)

#Add a column with Sample ID#
resddCT$Sample  <- Input$Sample.Name
#Simplify Table#
SimpleFC <- resddCT[,c("Group","Assay","Repeat","Sample","FoldChange")]

#Create a table of the points taht will be included in the graph#
write.table(SimpleFC,file="lme_graphpoints_ddCT_Muc5ac.txt",sep="\t",quote=F,row.names=F)

#read the table of graphing points#
lmepoint<-read.table("lme_graphpoints_ddCT_Muc5ac.txt",header=T,sep="\t")

#Factor the table by Group#
lmepoint$Group = factor(lmepoint$Group,c("PBS","CRA","DiHOME"))
lmepoint

#plot the table in a dot plot#
FoldChange <- ggplot(lmepoint, aes(x=Group, y=FoldChange)) + 
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2), size=8, shape=lmepoint$Repeat) +
  stat_summary(data=lmepoint, aes(ymax=..y..,ymin=..y.., line = factor(Group)), fun.y = "median", geom = "crossbar", colour = "black", width = 0.5)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(title="Muc5ac Fold Change", x="Group", y="FoldChange")
FoldChange
ggsave(filename="GRAPH-Muc5ac_FoldChange.pdf", plot=FoldChange,  useDingbats=FALSE,width=8.5, height=11, units="in")

#Compared to PBS#
Input$Group <- factor(Input$Group, c("PBS","CRA","DiHOME"))
model <- lmer(ddCT ~ Group + (1|Assay), data = Input)
summary(model)

#Compared to DiHOME#
Input$Group <- factor(Input$Group, c("DiHOME","PBS","CRA"))
model <- lmer(ddCT ~ Group + (1|Assay), data = Input)
summary(model)


#Select Gene for Analysis#
#Analysis of CCL11#
Input<-CCL11
Input$Group <- factor(Input$Group, c("PBS","CRA","DiHOME"))
Input

#Model and Graph Dotplot of ddCT for the Gene of interest#
mod.ddCT<- lmer(ddCT~Group+(1|Assay),data=Input)
y.ddCT<-summary(mod.ddCT)
y.ddCT
coef(mod.ddCT)
mod.ddCT$res
vddCT<-visreg(fit=mod.ddCT,xvar="Group",type="conditional")

#summarize each individual point from the model#
resddCT<-vddCT$res
resddCT["Don"] <- as.factor(Input$Assay)
resddCT$Group
resddCT$Repeat <- ifelse(resddCT$Don=="P1A2",15,ifelse(resddCT$Don=="P1A3", 16, ifelse(resddCT$Don=="P1A4",17, ifelse(resddCT$Don=="P1A5", 18, ifelse(resddCT$Don=="P1A6",3,0)))))
resddCT

#Add a column with Fold Change#
resddCT$FoldChange <- 2^(-resddCT$visregRes)

#Add a column with Sample ID#
resddCT$Sample  <- Input$Sample.Name
#Simplify Table#
SimpleFC <- resddCT[,c("Group","Assay","Repeat","Sample","FoldChange")]

#Create a table of the points taht will be included in the graph#
write.table(SimpleFC,file="lme_graphpoints_ddCT_CCL11.txt",sep="\t",quote=F,row.names=F)

#read the table of graphing points#
lmepoint<-read.table("lme_graphpoints_ddCT_CCL11.txt",header=T,sep="\t")

#Factor the table by Group#
lmepoint$Group = factor(lmepoint$Group,c("PBS","CRA","DiHOME"))
lmepoint

#plot the table in a dot plot#
FoldChange <- ggplot(lmepoint, aes(x=Group, y=FoldChange)) + 
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2), size=8, shape=lmepoint$Repeat) +
  stat_summary(data=lmepoint, aes(ymax=..y..,ymin=..y.., line = factor(Group)), fun.y = "median", geom = "crossbar", colour = "black", width = 0.5)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(title="CCL11 Fold Change", x="Group", y="FoldChange")
FoldChange
ggsave(filename="GRAPH-CCL11_FoldChange.pdf", plot=FoldChange,  useDingbats=FALSE,width=8.5, height=11, units="in")

#Compared to PBS#
Input$Group <- factor(Input$Group, c("PBS","CRA","DiHOME"))
model <- lmer(ddCT ~ Group + (1|Assay), data = Input)
summary(model)

#Compared to DiHOME#
Input$Group <- factor(Input$Group, c("DiHOME","PBS","CRA"))
model <- lmer(ddCT ~ Group + (1|Assay), data = Input)
summary(model)

#Remove Outliers and Repeat Statistics#
lmepoint
NoOutlier <- droplevels(subset(Input, Input$Sample.Name!="5-D2"))

#Compared to PBS#
NoOutlier$Group <- factor(NoOutlier$Group, c("PBS","CRA","DiHOME"))
model <- lmer(ddCT ~ Group + (1|Assay), data = NoOutlier)
summary(model)

#Compared to DiHOME#
NoOutlier$Group <- factor(NoOutlier$Group, c("DiHOME","PBS","CRA"))
model <- lmer(ddCT ~ Group + (1|Assay), data = NoOutlier)
summary(model)

