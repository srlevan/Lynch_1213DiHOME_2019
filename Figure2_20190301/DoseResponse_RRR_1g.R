rm(list=ls())

#Set working Directory#
setwd("~/Figure2_20190301/")

#Read in csv with Input Data#
INPUT <- read.csv("INPUT-Luciferase_20180729.csv")

#Load Required Packages#
library(plyr)
library(gridExtra)
library(reshape2)
library(nplr)
library(drc)
library(sandwich)
library(lmtest)
library(ggplot2)
library(multcomp)


#Subset each treatment#
Control<- subset(INPUT, Group=="Control")
DiHOME <- droplevels(subset(INPUT, Chem!="NoTF" & Chem!="NoPPAR" & Chem!="Rosi" & Chem!="GW1929" & Chem!="GW9662"))
GW1929 <- droplevels(subset(INPUT, Chem!="NoTF" & Chem!="NoPPAR" & Chem!="DiHOME" & Chem!="Rosi" & Chem!="GW9662"))
Rosi <- droplevels(subset(INPUT, Chem!="NoTF" & Chem!="NoPPAR" & Chem!="GW1929" & Chem!="DiHOME" & Chem!="GW9662"))

#Create SEM Function#
sem <- function(x) {
  sem<-sd(x)/sqrt(length(x))
  return(sem)
} 

meanDiHOME <- ddply(DiHOME, c("Dose"), summarise,
                  mean=mean(RRR),lower=mean(RRR)-sem(RRR),upper=mean(RRR)+sem(RRR))
meanDiHOME

meanGW1929 <- ddply(GW1929, c("Dose"), summarise,
                    mean=mean(RRR),lower=mean(RRR)-sem(RRR),upper=mean(RRR)+sem(RRR))
meanGW1929

meanRosi <- ddply(Rosi, c("Dose"), summarise,
                    mean=mean(RRR),lower=mean(RRR)-sem(RRR),upper=mean(RRR)+sem(RRR))
meanRosi


#Model GW1929#
model.GW1929 <- drm(RRR ~ Dose, data=GW1929, logDose=10,fct=LL2.3())
summary(model.GW1929)
coeftest(model.GW1929, vcov = sandwich)

# new dose levels as support for the line
GW1929newdata <- expand.grid(conc=exp(seq(log(0.000001), log(200), length=100))) # predictions and confidence intervals
pm <- predict(model.GW1929, newdata=GW1929newdata, interval="confidence") # new data with predictions
GW1929newdata$p <- pm[,1]
GW1929newdata$pmin <- pm[,2]
GW1929newdata$pmax <- pm[,3]

# need to shift conc == 0 a bit up, otherwise there are problems with coord_trans ryegrass$conc0 <- ryegrass$conc
GW1929$Dose[GW1929$Dose == 0] <- 0.000001
GW1929
# plotting the curve
GW1929.plot <- ggplot(GW1929, aes(x = Dose, y = RRR)) +
  geom_point() +
  geom_ribbon(data=GW1929newdata, aes(x=conc, y=p, ymin=pmin, ymax=pmax), alpha=0.2) +
  geom_line(data=GW1929newdata, aes(x=conc, y=p)) +
  coord_trans(x="log") +
  xlab("Dose (uM)") + ylab("Corrected Luminescence")+
  theme(axis.line = element_line(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank())
GW1929.plot  
ggsave(filename="GW1929.FitPlot.pdf", plot=GW1929.plot, useDingbats=FALSE)
summary(model.GW1929)

#Model Rosi#
Rosi
model.Rosi <- drm(RRR ~ Dose, data=Rosi, logDose=10,fct=LL2.3())
summary(model.Rosi)
coeftest(model.Rosi, vcov = sandwich)


# new dose levels as support for the line
Rosinewdata <- expand.grid(conc=exp(seq(log(0.000001), log(200), length=100))) # predictions and confidence intervals
pm <- predict(model.Rosi, newdata=Rosinewdata, interval="confidence") # new data with predictions
Rosinewdata$p <- pm[,1]
Rosinewdata$pmin <- pm[,2]
Rosinewdata$pmax <- pm[,3]

# need to shift conc == 0 a bit up, otherwise there are problems with coord_trans ryegrass$conc0 <- ryegrass$conc
Rosi$Dose[Rosi$Dose == 0] <- 0.000001
Rosi
# plotting the curve
Rosi.plot <- ggplot(Rosi, aes(x = Dose, y = RRR)) +
  geom_point() +
  geom_ribbon(data=Rosinewdata, aes(x=conc, y=p, ymin=pmin, ymax=pmax), alpha=0.2) +
  geom_line(data=Rosinewdata, aes(x=conc, y=p)) +
  coord_trans(x="log") +
  xlim(0.000001,100)+
  xlab("Dose (uM)") + ylab("Corrected Luminescence")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
Rosi.plot  
ggsave(filename="Rosi.FitPlot.pdf", plot=Rosi.plot, useDingbats=FALSE)
summary(model.Rosi)

#Model DiHOME#
model.DiHOME <- drm(RRR ~ Dose, data=DiHOME, logDose=10,fct=LL2.3())
summary(model.DiHOME)
coeftest(model.DiHOME, vcov = sandwich)


# new dose levels as support for the line
DiHOMEnewdata <- expand.grid(conc=exp(seq(log(0.000001), log(200), length=100))) # predictions and confidence intervals
pm <- predict(model.DiHOME, newdata=DiHOMEnewdata, interval="confidence") # new data with predictions
DiHOMEnewdata$p <- pm[,1]
DiHOMEnewdata$pmin <- pm[,2]
DiHOMEnewdata$pmax <- pm[,3]

# need to shift conc == 0 a bit up, otherwise there are problems with coord_trans ryegrass$conc0 <- ryegrass$conc
DiHOME$Dose[DiHOME$Dose == 0] <- 0.000001
DiHOME
# plotting the curve
DiHOME.plot <- ggplot(DiHOME, aes(x = Dose, y = RRR)) +
  geom_point() +
  geom_ribbon(data=DiHOMEnewdata, aes(x=conc, y=p, ymin=pmin, ymax=pmax), alpha=0.2) +
  geom_line(data=DiHOMEnewdata, aes(x=conc, y=p)) +
  coord_trans(x="log") +
  xlim(0.000001,300)+
  xlab("Dose (uM)") + ylab("Corrected Luminescence")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
DiHOME.plot  
ggsave(filename="DiHOME.FitPlot.pdf", plot=DiHOME.plot, useDingbats=FALSE)
summary(model.DiHOME)


#Limit X Range#
Combo.plot <- ggplot() +
  geom_point(data=meanGW1929, colour="red", aes(x = Dose, y = mean)) +
  geom_point(data=meanRosi, colour="blue", aes(x = Dose, y = mean))+
  geom_point(data=meanDiHOME, aes(x = Dose, y = mean))+
  geom_errorbar(data=meanGW1929, aes(x=Dose, ymin=lower, ymax=upper), colour="red", width=.01, position=position_dodge(0.8))+
  geom_errorbar(data=meanRosi, aes(x=Dose, ymin=lower, ymax=upper), colour="blue", width=.01, position=position_dodge(0.8))+
  geom_errorbar(data=meanDiHOME, aes(x=Dose, ymin=lower, ymax=upper), width=.01, position=position_dodge(0.8))+
  geom_line(data=GW1929newdata, aes(x=conc, y=p), colour="red") +
  geom_line(data=Rosinewdata, aes(x=conc, y=p), colour="blue") +
  geom_line(data=DiHOMEnewdata, aes(x=conc, y=p)) +
  geom_ribbon(data=GW1929newdata, aes(x=conc, y=p, ymin=pmin, ymax=pmax), fill="red",alpha=0.2) +
  geom_ribbon(data=Rosinewdata, aes(x=conc, y=p, ymin=pmin, ymax=pmax), fill="blue", alpha=0.2) +
  geom_ribbon(data=DiHOMEnewdata, aes(x=conc, y=p, ymin=pmin, ymax=pmax), fill="black", alpha=0.2) +
  coord_trans(x="log") +
  xlab("Dose (uM)") + ylab("Corrected Luminescence")+
  scale_x_continuous(breaks=c(0.01,0.1,1,10,100,200), trans="log1p", limits=c(0.075,200))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
Combo.plot
ggsave(filename="GRAPH-DoseResponseCurves.pdf", plot=Combo.plot, useDingbats=FALSE,width=11, height=8.5, units="in")


#Negative Controls#
NoPPAR <- subset(INPUT, Chem=="NoPPAR")
NoTF <- subset(INPUT, Chem=="NoTF")
UnTx <- subset(INPUT, Chem=="UnTx")

#Calculate mean and sem for UnTx#
UnTx
meanUnTx <- ddply(UnTx, c("Dose"), summarise,
                    mean=mean(RRR),lower=mean(RRR)-sem(RRR),upper=mean(RRR)+sem(RRR))
meanUnTx

#Calculate mean and sem for NoPPAR.DiHOME#
NoPPAR
NoPPAR.DiHOME <- subset(NoPPAR, Treatment=="200 uM  DiHOME")
NoPPAR.DiHOME$Dose[NoPPAR.DiHOME$Dose == 0] <- 200

meanNoPPAR.DiHOME <- ddply(NoPPAR.DiHOME, c("Dose"), summarise,
                  mean=mean(RRR),lower=mean(RRR)-sem(RRR),upper=mean(RRR)+sem(RRR))
meanNoPPAR.DiHOME

#Calculate mean and sem for No TF#
NoTF
meanNoTF <- ddply(NoTF, c("Dose"), summarise,
                           mean=mean(RRR),lower=mean(RRR)-sem(RRR),upper=mean(RRR)+sem(RRR))
meanNoTF

#Subset Data#
meanNoTF$Chem <- "NoTF"
meanNoPPAR.DiHOME$Chem <- "NoPPAR"
meanDiHOME.200 <- subset(meanDiHOME, Dose==200) 
meanDiHOME.200$Chem <- "DiHOME"
meanUnTx$Chem <- "UnTx"
NegativeControlMeans<-rbind(meanUnTx,meanNoTF,meanNoPPAR.DiHOME,meanDiHOME.200)
NegativeControlMeans$Chem <- factor(NegativeControlMeans$Chem, c("NoTF","NoPPAR","UnTx","DiHOME"))

DiHOME.200 <- subset(DiHOME, Dose==200)
NoTF
NoPPAR.DiHOME
UnTx
NegativeControlScatter <- rbind(UnTx,NoTF,NoPPAR.DiHOME,DiHOME.200)
NegativeControlScatter$Chem <- factor(NegativeControlScatter$Chem, c("NoTF","NoPPAR","UnTx","DiHOME"))

#Limit X Range#
NegativeControls.plot <- ggplot() +
  geom_bar(data=NegativeControlMeans, aes(x=Chem, y=mean,fill=Chem),stat='identity',position="dodge",colour="black")+
   geom_errorbar(data=NegativeControlMeans, aes(x=Chem,y=mean, ymin=lower, ymax=upper),position=position_dodge(0.8), colour="black", width=.1)+
  geom_jitter(data=NegativeControlScatter, aes(x=Chem, y=RRR), position=position_jitter(0.2), size=5)+
  xlab("Dose (uM)") + ylab("Corrected Luminescence")+
  theme_classic() + scale_fill_brewer(palette = "Greys", direction=1)+
  guides(fill = guide_legend(reverse=FALSE))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
NegativeControls.plot
ggsave(filename="GRAPH-NegativeControlsBarGraph.pdf", plot=NegativeControls.plot, useDingbats=FALSE,width=8.5, height=11, units="in")

