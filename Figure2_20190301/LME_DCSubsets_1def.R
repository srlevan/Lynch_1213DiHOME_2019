#Clear the Environment#
rm(list=ls())

#load required packages#
library(lmerTest)
library(visreg)
library(ggplot2)

#Set Working Directory#
setwd("~/Figure2_20190301/")

#Read in Flow Data with DC Subsets#
dihome=read.csv("INPUT-DCSubsets.csv")
dihome

#Subset raw data to include treatments and cell subtypes of interest#
dihome$Treatment = factor(dihome$Treatment, c("DMSO", "DiHOME-75", "DiHOME-130", "DiHOME-200"))


#Select Frequencies#
S1<-dihome[,c("Treatment","Donor","Live","mDC","CCR7","CD36","CD80","CD1a")]
S1

#Compare to baseline for each donor#
D15 <- subset(S1, Donor =="D15")
D15$RelativeCCR7 <- D15$CCR7/mean(subset(D15, Treatment=="DMSO")$CCR7)
D15$RelativeCD36 <- D15$CD36/mean(subset(D15, Treatment=="DMSO")$CD36)
D15$RelativeCD80 <- D15$CD80/mean(subset(D15, Treatment=="DMSO")$CD80)
D15$RelativeCD1a <- D15$CD1a/mean(subset(D15, Treatment=="DMSO")$CD1a)

D17 <- subset(S1, Donor =="D17")
D17$RelativeCCR7 <- D17$CCR7/mean(subset(D17, Treatment=="DMSO")$CCR7)
D17$RelativeCD36 <- D17$CD36/mean(subset(D17, Treatment=="DMSO")$CD36)
D17$RelativeCD80 <- D17$CD80/mean(subset(D17, Treatment=="DMSO")$CD80)
D17$RelativeCD1a <- D17$CD1a/mean(subset(D17, Treatment=="DMSO")$CD1a)

D24 <- subset(S1, Donor =="D24")
D24$RelativeCCR7 <- D24$CCR7/mean(subset(D24, Treatment=="DMSO")$CCR7)
D24$RelativeCD36 <- D24$CD36/mean(subset(D24, Treatment=="DMSO")$CD36)
D24$RelativeCD80 <- D24$CD80/mean(subset(D24, Treatment=="DMSO")$CD80)
D24$RelativeCD1a <- D24$CD1a/mean(subset(D24, Treatment=="DMSO")$CD1a)

D26 <- subset(S1, Donor =="D26")
D26$RelativeCCR7 <- D26$CCR7/mean(subset(D26, Treatment=="DMSO")$CCR7)
D26$RelativeCD36 <- D26$CD36/mean(subset(D26, Treatment=="DMSO")$CD36)
D26$RelativeCD80 <- D26$CD80/mean(subset(D26, Treatment=="DMSO")$CD80)
D26$RelativeCD1a <- D26$CD1a/mean(subset(D26, Treatment=="DMSO")$CD1a)
D26
#Combine Relative Data#
S2 <- rbind(D15,D17,D24,D26)
Input <- S2[,c("Treatment","Donor","Live","mDC","RelativeCCR7","RelativeCD36","RelativeCD80","RelativeCD1a")]

#Run a Linear Mixed Effects Model correcting for Donor Variability on the DC Subsets#
Input$Treatment = factor(Input$Treatment, c("DMSO", "DiHOME-75", "DiHOME-130", "DiHOME-200"))
x <- matrix(data = NA, nrow =  4, ncol = 3, dimnames=list(1:4, c("Estimate", "Standard Error", "P")))
for (i in 3:8){
  model <- lmer(Input[[i]] ~ Treatment + (1|Donor), data = Input)
  y <- summary(model)
  x[,1] <- y$coefficients[,1]
  x[,2] <- y$coefficients[,2]
  x[,3] <- y$coefficients[,5]
  print(colnames(Input[i]))
  print(x)
}


#Plot LME Model of Live Cell Frequency as a Dotplot#
#model Live + Donorn#
mod.Live<- lmer(Live~Treatment+(1|Donor),data=Input)
y.Live<-summary(mod.Live)
y.Live
coef(mod.Live)
mod.Live$res

vLive<-visreg(fit=mod.Live,xvar="Treatment",type="conditional")

#summarize each individual point from the model#
resLive<-vLive$res
resLive$Donor <- Input$Donor
resLive$Don <- resLive$Donor
levels(resLive$Don) <- c("18","16","17","15")
resLive

#Save a table of the points that will be included in the graph#
write.table(resLive,file="lme_graphpoints_Live.txt",sep="\t",quote=F,row.names=F)

#read the table of graphing points#
lmepoint <- read.table("lme_graphpoints_Live.txt", header=T, sep="\t")

#Factor the table by treatment#
lmepoint$Treatment = factor(lmepoint$Treatment,c("DMSO","DiHOME-75", "DiHOME-130", "DiHOME-200"))
lmepoint

#Plot as a DotPlot + Boxplot#
Live <- ggplot(data=lmepoint, aes(x=Treatment, y=visregRes)) + 
  geom_boxplot()+
  geom_jitter(data=lmepoint, aes(x=Treatment, y=visregRes), shape=lmepoint$Don, position=position_jitter(0.3), size=8)+
  scale_shape_manual(values = c(18,16,15,17))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(title="Live_SubsetSamples", x="Treatment", y="Live Frequency")
Live

#Save the Plot of Live Cells#
ggsave(filename="GRAPH-LiveDCFrequency.pdf", plot=Live, useDingbats=FALSE,width=9, height=8.5, units="in")

#Plot LME Model of mDC Cell Frequency as a Dotplot#
#model mDC + Donorn#
mod.mDC<- lmer(mDC~Treatment+(1|Donor),data=Input)
y.mDC<-summary(mod.mDC)
y.mDC
coef(mod.mDC)
mod.mDC$res

vmDC<-visreg(fit=mod.mDC,xvar="Treatment",type="conditional")

#summarize each individual point from the model#
resmDC<-vmDC$res
resmDC$Donor <- Input$Donor
resmDC$Don <- resmDC$Donor
levels(resmDC$Don) <- c("18","16","17","15")
resmDC

#Save a table of the points that will be included in the graph#
write.table(resmDC,file="lme_graphpoints_mDC.txt",sep="\t",quote=F,row.names=F)


#read the table of graphing points#
lmepoint<-read.table("lme_graphpoints_mDC.txt",header=T,sep="\t")
#Factor the table by treatment#
lmepoint$Treatment = factor(lmepoint$Treatment,c("DMSO","DiHOME-75", "DiHOME-130", "DiHOME-200"))
lmepoint


#Plot as a DotPlot#
mDC <- ggplot(data=lmepoint, aes(x=Treatment, y=visregRes)) + 
  geom_boxplot()+
  geom_jitter(data=lmepoint, aes(x=Treatment, y=visregRes), shape=lmepoint$Don, position=position_jitter(0.3), size=8)+
  scale_shape_manual(values = c(18,16,15,17))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(title="mDC_SubsetSamples", x="Treatment", y="mDC Frequency")
mDC

#Save the Plot of mDC Cells#
ggsave(filename="GRAPH-mDCDCFrequency.pdf", plot=mDC, useDingbats=FALSE,width=9, height=8.5, units="in")


#Plot LME Model of RelativeCCR7 Cell Frequency as a Dotplot#
#model RelativeCCR7 + Donor#
mod.RelativeCCR7<- lmer(RelativeCCR7~Treatment+(1|Donor),data=Input)
y.RelativeCCR7<-summary(mod.RelativeCCR7)
y.RelativeCCR7
coef(mod.RelativeCCR7)
mod.RelativeCCR7$res

vRelativeCCR7<-visreg(fit=mod.RelativeCCR7,xvar="Treatment",type="conditional")

#summarize each individual point from the model#
resRelativeCCR7<-vRelativeCCR7$res
resRelativeCCR7$Donor <- Input$Donor
resRelativeCCR7$Don <- resRelativeCCR7$Donor
levels(resRelativeCCR7$Don) <- c("18","16","17","15")
resRelativeCCR7

#Save a table of the points that will be included in the graph#
write.table(resRelativeCCR7,file="lme_graphpoints_RelativeCCR7.txt",sep="\t",quote=F,row.names=F)

#read the table of graphing points#
lmepoint<-read.table("lme_graphpoints_RelativeCCR7.txt",header=T,sep="\t")
#Factor the table by treatment#
lmepoint$Treatment = factor(lmepoint$Treatment,c("DMSO","DiHOME-75", "DiHOME-130", "DiHOME-200"))
lmepoint

#Plot as a DotPlot + Boxplot#
RelativeCCR7 <- ggplot(data=lmepoint, aes(x=Treatment, y=visregRes)) + 
  geom_boxplot()+
  geom_jitter(data=lmepoint, aes(x=Treatment, y=visregRes), shape=lmepoint$Don, position=position_jitter(0.3), size=8)+
  scale_shape_manual(values = c(18,16,15,17))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(title="RelativeCCR7_SubsetSamples", x="Treatment", y="RelativeCCR7 Frequency")
RelativeCCR7

#Save the Plot of RelativeCCR7 Cells#
ggsave(filename="GRAPH-RelativeCCR7DCFrequency.pdf", plot=RelativeCCR7,useDingbats=FALSE,width=9, height=8.5, units="in")

#Compared to DMSO#
Input$Group <- factor(Input$Treatment,c("DMSO","DiHOME-75", "DiHOME-130", "DiHOME-200"))
model <- lmer(RelativeCCR7 ~ Group + (1|Donor), data = Input)
summary(model)

#Compared to DiHOME-200#
Input$Group <- factor(Input$Treatment,c("DiHOME-200","DMSO","DiHOME-75", "DiHOME-130"))
model <- lmer(RelativeCCR7 ~ Group + (1|Donor), data = Input)
summary(model)

#Compared to DiHOME-75#
Input$Group <- factor(Input$Treatment,c("DiHOME-75","DMSO", "DiHOME-130","DiHOME-200"))
model <- lmer(RelativeCCR7 ~ Group + (1|Donor), data = Input)
summary(model)

#Remove Outliers and Repeat Statistics#
Input
NoOutlier <- Input[c(1:6,8:26,28:43,45:48),]

#Compared to DMSO#
NoOutlier$Group <- factor(NoOutlier$Treatment,c("DMSO","DiHOME-75", "DiHOME-130", "DiHOME-200"))
model <- lmer(RelativeCCR7 ~ Group + (1|Donor), data = NoOutlier)
summary(model)

#Compared to DiHOME-200#
NoOutlier$Group <- factor(NoOutlier$Treatment,c("DiHOME-200","DMSO","DiHOME-75", "DiHOME-130"))
model <- lmer(RelativeCCR7 ~ Group + (1|Donor), data = NoOutlier)
summary(model)

#Compared to DiHOME-75#
NoOutlier$Group <- factor(NoOutlier$Treatment,c("DiHOME-75","DMSO","DiHOME-130", "DiHOME-200"))
model <- lmer(RelativeCCR7 ~ Group + (1|Donor), data = NoOutlier)
summary(model)

#Plot LME Model of RelativeCD80 Cell Frequency as a Dotplot#
#model RelativeCD80 + Donorn#
mod.RelativeCD80<- lmer(RelativeCD80~Treatment+(1|Donor),data=Input)
y.RelativeCD80<-summary(mod.RelativeCD80)
y.RelativeCD80
coef(mod.RelativeCD80)
mod.RelativeCD80$res

vRelativeCD80<-visreg(fit=mod.RelativeCD80,xvar="Treatment",type="conditional")

#summarize each individual point from the model#
resRelativeCD80<-vRelativeCD80$res
resRelativeCD80$Donor <- Input$Donor
resRelativeCD80$Don <- resRelativeCD80$Donor
levels(resRelativeCD80$Don) <- c("18","16","17","15")
resRelativeCD80

#Save a table of the points that will be included in the graph#
write.table(resRelativeCD80,file="lme_graphpoints_RelativeCD80.txt",sep="\t",quote=F,row.names=F)

#read the table of graphing points#
lmepoint<-read.table("lme_graphpoints_RelativeCD80.txt",header=T,sep="\t")
#Factor the table by treatment#
lmepoint$Treatment = factor(lmepoint$Treatment,c("DMSO","DiHOME-75", "DiHOME-130", "DiHOME-200"))
lmepoint


#Plot as a DotPlot#
RelativeCD80 <- ggplot(data=lmepoint, aes(x=Treatment, y=visregRes)) + 
  geom_boxplot()+
  geom_jitter(data=lmepoint, aes(x=Treatment, y=visregRes), shape=lmepoint$Don, position=position_jitter(0.3), size=8)+
  scale_shape_manual(values = c(18,16,15,17))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(title="RelativeCD80_SubsetSamples", x="Treatment", y="RelativeCD80 Frequency")
RelativeCD80

#Save the Plot of RelativeCD80 Cells#
ggsave(filename="GRAPH-RelativeCD80DCFrequency.pdf", plot=RelativeCD80,useDingbats=FALSE,width=9, height=8.5, units="in")

#Compared to DMSO#
Input$Group <- factor(Input$Treatment,c("DMSO","DiHOME-75", "DiHOME-130", "DiHOME-200"))
model <- lmer(RelativeCD80 ~ Group + (1|Donor), data = Input)
summary(model)

#Compared to DiHOME-200#
Input$Group <- factor(Input$Treatment,c("DiHOME-200","DMSO","DiHOME-75", "DiHOME-130"))
model <- lmer(RelativeCD80 ~ Group + (1|Donor), data = Input)
summary(model)

#Compared to DiHOME-75#
Input$Group <- factor(Input$Treatment,c("DiHOME-75","DMSO","DiHOME-130", "DiHOME-200"))
model <- lmer(RelativeCD80 ~ Group + (1|Donor), data = Input)
summary(model)


#Remove Outliers and Repeat Statistics#
#Outliers 1, 4,7#
Input
NoOutlier <- Input[c(2:3,5:6,8:48),]

#Compared to DMSO#
NoOutlier$Group <- factor(NoOutlier$Treatment,c("DMSO","DiHOME-75", "DiHOME-130", "DiHOME-200"))
model <- lmer(RelativeCD80 ~ Group + (1|Donor), data = NoOutlier)
summary(model)

#Compared to DiHOME-200#
NoOutlier$Group <- factor(NoOutlier$Treatment,c("DiHOME-200","DMSO","DiHOME-75", "DiHOME-130"))
model <- lmer(RelativeCD80 ~ Group + (1|Donor), data = NoOutlier)
summary(model)

#Compared to DiHOME-75#
NoOutlier$Group <- factor(NoOutlier$Treatment,c("DiHOME-75","DMSO","DiHOME-130", "DiHOME-200"))
model <- lmer(RelativeCD80 ~ Group + (1|Donor), data = NoOutlier)
summary(model)


#Plot LME Model of RelativeCD1a Cell Frequency as a Dotplot#
#model RelativeCD1a + Donorn#
mod.RelativeCD1a<- lmer(RelativeCD1a~Treatment+(1|Donor),data=Input)
y.RelativeCD1a<-summary(mod.RelativeCD1a)
y.RelativeCD1a
coef(mod.RelativeCD1a)
mod.RelativeCD1a$res

vRelativeCD1a<-visreg(fit=mod.RelativeCD1a,xvar="Treatment",type="conditional")

#summarize each individual point from the model#
resRelativeCD1a<-vRelativeCD1a$res
resRelativeCD1a$Donor <- Input$Donor
resRelativeCD1a$Don <- resRelativeCD1a$Donor
levels(resRelativeCD1a$Don) <- c("18","16","17","15")
resRelativeCD1a

#Save a table of the points that will be included in the graph#
write.table(resRelativeCD1a,file="lme_graphpoints_RelativeCD1a.txt",sep="\t",quote=F,row.names=F)

#read the table of graphing points#
lmepoint<-read.table("lme_graphpoints_RelativeCD1a.txt",header=T,sep="\t")
#Factor the table by treatment#
lmepoint$Treatment = factor(lmepoint$Treatment,c("DMSO","DiHOME-75", "DiHOME-130", "DiHOME-200"))
lmepoint


#Plot as a DotPlot#
RelativeCD1a <- ggplot(data=lmepoint, aes(x=Treatment, y=visregRes)) + 
  geom_boxplot()+
  geom_jitter(data=lmepoint, aes(x=Treatment, y=visregRes), shape=lmepoint$Don, position=position_jitter(0.3), size=8)+
  scale_shape_manual(values = c(18,16,15,17))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(title="RelativeCD1a_SubsetSamples", x="Treatment", y="RelativeCD1a Frequency")
RelativeCD1a

#Save the Plot of RelativeCD1a Cells#
ggsave(filename="GRAPH-RelativeCD1aDCFrequency.pdf", plot=RelativeCD1a,useDingbats=FALSE,width=9, height=8.5, units="in")

#Compared to DMSO#
Input$Group <- factor(Input$Treatment,c("DMSO","DiHOME-75", "DiHOME-130", "DiHOME-200"))
model <- lmer(RelativeCD1a ~ Group + (1|Donor), data = Input)
summary(model)

#Compared to DiHOME-200#
Input$Group <- factor(Input$Treatment,c("DiHOME-200","DMSO","DiHOME-75", "DiHOME-130"))
model <- lmer(RelativeCD1a ~ Group + (1|Donor), data = Input)
summary(model)

#Compared to DiHOME-75#
Input$Group <- factor(Input$Treatment,c("DiHOME-75","DMSO","DiHOME-130", "DiHOME-200"))
model <- lmer(RelativeCD1a ~ Group + (1|Donor), data = Input)
summary(model)

#Remove Outliers and Repeat Statistics#
#13, 41, 19
Input
NoOutlier <- Input[c(1:12, 13:18, 20:40, 42:48),]

#Compared to DMSO#
NoOutlier$Group <- factor(NoOutlier$Treatment,c("DMSO","DiHOME-75", "DiHOME-130", "DiHOME-200"))
model <- lmer(RelativeCD1a ~ Group + (1|Donor), data = NoOutlier)
summary(model)

#Compared to DiHOME-200#
NoOutlier$Group <- factor(NoOutlier$Treatment,c("DiHOME-200","DMSO","DiHOME-75", "DiHOME-130"))
model <- lmer(RelativeCD1a ~ Group + (1|Donor), data = NoOutlier)
summary(model)

#Compared to DiHOME-75#
NoOutlier$Group <- factor(NoOutlier$Treatment,c("DiHOME-75","DMSO","DiHOME-130", "DiHOME-200"))
model <- lmer(RelativeCD1a ~ Group + (1|Donor), data = NoOutlier)
summary(model)


#Plot LME Model of RelativeCD36 Cell Frequency as a Dotplot#
#model RelativeCD36 + Donorn#
mod.RelativeCD36<- lmer(RelativeCD36~Treatment+(1|Donor),data=Input)
y.RelativeCD36<-summary(mod.RelativeCD36)
y.RelativeCD36
coef(mod.RelativeCD36)
mod.RelativeCD36$res

vRelativeCD36<-visreg(fit=mod.RelativeCD36,xvar="Treatment",type="conditional")

#summarize each individual point from the model#
resRelativeCD36<-vRelativeCD36$res
resRelativeCD36$Donor <- Input$Donor
resRelativeCD36$Don <- resRelativeCD36$Donor
levels(resRelativeCD36$Don) <- c("18","16","17","15")
resRelativeCD36

#Save a table of the points that will be included in the graph#
write.table(resRelativeCD36,file="lme_graphpoints_RelativeCD36.txt",sep="\t",quote=F,row.names=F)

#read the table of graphing points#
lmepoint<-read.table("lme_graphpoints_RelativeCD36.txt",header=T,sep="\t")
#Factor the table by treatment#
lmepoint$Treatment = factor(lmepoint$Treatment,c("DMSO","DiHOME-75", "DiHOME-130", "DiHOME-200"))
lmepoint


#Plot as a DotPlot#
RelativeCD36 <- ggplot(data=lmepoint, aes(x=Treatment, y=visregRes)) + 
  geom_boxplot()+
  geom_jitter(data=lmepoint, aes(x=Treatment, y=visregRes), shape=lmepoint$Don, position=position_jitter(0.3), size=8)+
  scale_shape_manual(values = c(18,16,15,17))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(title="RelativeCD36_SubsetSamples", x="Treatment", y="RelativeCD36 Frequency")
RelativeCD36

#Save the Plot of RelativeCD36 Cells#
ggsave(filename="GRAPH-RelativeCD36DCFrequency.pdf", plot=RelativeCD36,useDingbats=FALSE,width=9, height=8.5, units="in")

