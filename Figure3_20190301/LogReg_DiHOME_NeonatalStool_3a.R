#Clear Environment#
rm(list=ls())

#Load Required Packages#
library(lmerTest)
library(visreg)
library(ggplot2)
library(plyr)
library(gridExtra)
library(reshape2)

#set wd#
setwd("~/Box Sync/Lynch Lab/Code Respositories/Figure3_20190301/")

#read csv#
HumanDB=read.csv("INPUT-HumanSamples_MasterDatabase.csv")
HumanDB$Test.1213.ngmg <- HumanDB$Test.1213.pgmg/1000
HumanDB$Test.910.ngmg <- HumanDB$Test.910.pgmg/1000


#Select only WHEALS DiHOME Data#
WHEALS <- subset(HumanDB, HumanDB$Study=="WHEALS")
colnames(WHEALS)

#Add Health Status and Condition#
WHEALS$Healthy <- ifelse(WHEALS$Diagnosis.AorA==1,"AorA","Healthy")
WHEALS$Condition <- ifelse(WHEALS$Diagnosis.Atopy+WHEALS$Diagnosis.Asthma>1,"Atopic Asthma",ifelse(WHEALS$Diagnosis.Atopy>0,"Atopic",ifelse(WHEALS$Diagnosis.Asthma>0,"Asthmatic","Healthy")))
WHEALS$Healthy<- factor(WHEALS$Healthy, c("Healthy","AorA"))
WHEALS$Condition <- factor(WHEALS$Condition, c("Healthy","Atopic","Asthmatic","Atopic Asthma"))

#Reduce to only Columns of Interest#
colnames(WHEALS)
WHEALS_META <- WHEALS[,c(24,25,3,22,23,14:21)]
WHEALS_DiHOME <- WHEALS_META[,c(1:5)]
head(WHEALS_DiHOME)

#Factor by condition and health status#
WHEALS_DiHOME$Healthy<- factor(WHEALS_DiHOME$Healthy, c("Healthy","AorA"))
WHEALS_DiHOME$Condition <- factor(WHEALS_DiHOME$Condition, c("Healthy","Atopic","Asthmatic","Atopic Asthma"))

#Convert to Long Form#
WHEALS_DiHOME_Long <- melt(WHEALS_DiHOME, id=c("Healthy","Condition","Lynch_ID"))
WHEALS_DiHOME_Long <- WHEALS_DiHOME_Long[,c(1:2,4:5)]

WHEALS_1213_Long <- subset(WHEALS_DiHOME_Long, WHEALS_DiHOME_Long$variable=="Test.1213.ngmg")
WHEALS_910_Long <- subset(WHEALS_DiHOME_Long, WHEALS_DiHOME_Long$variable=="Test.910.ngmg")


#Plot Box Plot of 12,13 DiHOME by Healthy Status#
BoxPlot <- ggplot() + 
  geom_boxplot(data=WHEALS_1213_Long, aes(x=Healthy, y=value))+
  geom_jitter(data=WHEALS_1213_Long, aes(x=Healthy, y=value,shape=WHEALS_1213_Long$Condition),position=position_jitter(0.2), size=5)+
  theme_classic() +scale_fill_brewer(palette = "Greys", direction=1)+
  scale_shape_manual(values = c(18,17,15,16))+
  guides(fill = guide_legend(reverse=FALSE))+
  labs(title="WHEALS_1213_Summary", y="ng.mg")
BoxPlot

ggsave(filename="GRAPH-1213DiHOME_NeonatalStool.pdf", plot=BoxPlot, useDingbats=FALSE,width=8.5, height=11, units="in")


#Plot Bar Graph of 9,10 DiHOME by Health Status#
BoxPlot <- ggplot() + 
  geom_boxplot(data=WHEALS_910_Long, aes(x=Healthy, y=value))+
  geom_jitter(data=WHEALS_910_Long, aes(x=Healthy, y=value,shape=WHEALS_910_Long$Condition),position=position_jitter(0.2),size=5)+
  theme_classic() +scale_fill_brewer(palette = "Greys", direction=1)+
  scale_shape_manual(values = c(18,17,15,16))+
  guides(fill = guide_legend(reverse=FALSE))+
  labs(title="WHEALS_910_Summary", y="ng.mg")
BoxPlot

ggsave(filename="GRAPH-910DiHOME_NeonatalStool.pdf", plot=BoxPlot, useDingbats=FALSE,width=8.5, height=11, units="in")
 
#Test the significance of the relationship between 12,13 DiHOME and Health Status using logistical regression#
Input<-WHEALS_META
head(Input)
model = glm(WHEALS_META$Healthy ~ WHEALS_META$Test.1213.ngmg, 
            data = Input, 
            family = binomial(link="logit"))
summary<-summary(model)
coefficients <- summary$coefficients
Estimate <- coefficients[c(2),c(1)]
Error <- coefficients[c(2),c(2)]
Pvalue <- coefficients[c(2),c(4)]
OR <- exp(Estimate)
lower <- exp(Estimate-1.96*Error)
upper <- exp(Estimate+1.96*Error)
AIC <- model$aic
Output <- matrix(data = NA, nrow =  7, ncol = 1, dimnames=list(c("Estimate","Error","pvalue","OR","lower CI","upper CI","AIC"), c(1)))
Output[1,1] <- Estimate
Output[2,1] <- Error
Output[3,1] <- Pvalue
Output[4,1] <- OR
Output[5,1] <- lower
Output[6,1] <- upper
Output[7,1] <- AIC
t(Output)
OR_Univariate_LogReg <- OR

#Test for Confounding using ChiSquared Test#
Input <- WHEALS_META

ConfounderTest <- function(Diagnosis,Test,SignificantInput){
  parent.model=glm(Diagnosis ~ Test,
                   data=Input,
                   family=binomial(link="logit"))
  extended.model = glm(Diagnosis ~ Test + SignificantInput, 
                       data = Input, 
                       family = binomial(link="logit"))
  testfordeviance <- anova(parent.model,extended.model, test="Chi")
  pvalue <- testfordeviance[2,5]
  print(pvalue)
}

colnames(WHEALS_META)

PotentialConfounders <- names(Input[,c(6:13)])

Output<- matrix(data = NA, nrow =  8, ncol = 3, dimnames=list(c(1:8),c("PotentialConfounder", "pvalue","Confounding?")))

for (i in 1:8){
  var <- PotentialConfounders[i]
  Stats <- ConfounderTest(Input$Healthy, Input$Test.1213.ngmg, Input[,var])
  Output[(i),2] <- Stats
  Output[(i),1] <- paste0(var,".AorA")
  Output[(i),3] <- ifelse(Stats<0.05,"YES","NO")
}
Confounder_WHEALS_1213_AorA <- as.data.frame(Output)
write.table(Confounder_WHEALS_1213_AorA,file="Confounders_WHEALS_1213_AorA.txt",sep="\t",quote=F,row.names=F)

#NO CONFOUNDERS IDENTIFIED for 12,13 DiHOME and A or A by variation of deviance#

#Test for Confounding by Change in OR#
ConfounderORChange <- function(Diagnosis,SignificantInput,Test){
  model = glm(Diagnosis ~ SignificantInput + Test, 
              data = Input, 
              family = binomial(link="logit"))
  summary<-summary(model)
  coefficients <- summary$coefficients
  Estimate <- coefficients[c(2),c(1)]
  Error <- coefficients[c(2),c(2)]
  Pvalue <- coefficients[c(2),c(4)]
  OR <- exp(Estimate)
  lower <- exp(Estimate-1.96*Error)
  upper <- exp(Estimate+1.96*Error)
  AIC <- model$aic
  ORChange <- 100*((OR/OR_SignificantInput_Diagnosis)-1)
  Confounder <- ifelse(ORChange>10,1,0)
  Output <- matrix(data = NA, nrow =  9, ncol = 1, dimnames=list(c("Estimate","Error","pvalue","OR","lower CI","upper CI","AIC","OR%Change","Confounder"), c(1)))
  Output[1,1] <- Estimate
  Output[2,1] <- Error
  Output[3,1] <- Pvalue
  Output[4,1] <- OR
  Output[5,1] <- lower
  Output[6,1] <- upper
  Output[7,1] <- AIC
  Output[8,1] <- ORChange
  Output[9,1] <- Confounder
  return(Output)
}

#Must Specify the OR for the Significant Test#
OR_Univariate_LogReg

#Select Potential Confounders of Test.1213.ngmg#
colnames(Input)
RiskFactors <- names(Input[,c(6:13)])

Output<- matrix(data = NA, nrow =  9, ncol = 8, dimnames=list(c("Estimate","Error","pvalue","OR","lower CI","upper CI","AIC","OR%Change","Confounder"),c(1:8)))

#Specify the OR of the Significant Test#
OR_SignificantInput_Diagnosis <- OR_Univariate_LogReg

for (i in 1:8){
  var <- RiskFactors[i]
  Stats <- ConfounderORChange(Input$Healthy, Input$Test.1213.ngmg, Input[,var])
  Output[,(i)] <- Stats
  colnames(Output) <- paste0(RiskFactors,".AorA")
}
Confounder_Test.1213.ngmg_Healthy <- as.data.frame(t(Output))
write.table(Confounder_Test.1213.ngmg_Healthy,file="WHEALS_ConfounderORChange_Test.1213.ngmg_AorA.txt",sep="\t",quote=F,row.names=T)

#Risk.Black Identified as a potential confounder#

#Count Data#
head(Input)

RiskFactors <- names(Input[,c(6:12)])

Summary <- matrix(data=NA, nrow=2, ncol=7,dimnames=list(c("N","Percent"),c(1:7)))
for (i in 1:7){
  var <- RiskFactors[i]
  Summary[1,(i)] <- sum(Input[,var])
  Summary[2,(i)] <- 100*sum(Input[,var])/length(Input[,var])  
  colnames(Summary) <- paste0(RiskFactors,".Summary")
}   
Summary 
median(Input$Age.Stool)


#Test the significance of the relationship between 9,10 DiHOME and Health Status using logistical regression#
Input<-WHEALS_META
head(Input)
model = glm(WHEALS_META$Healthy ~ WHEALS_META$Test.910.ngmg, 
            data = Input, 
            family = binomial(link="logit"))
summary<-summary(model)

#Determine pvalue#
coefficients <- summary$coefficients
Estimate <- coefficients[c(2),c(1)]
Error <- coefficients[c(2),c(2)]
Pvalue <- coefficients[c(2),c(4)]
OR <- exp(Estimate)
lower <- exp(Estimate-1.96*Error)
upper <- exp(Estimate+1.96*Error)
AIC <- model$aic
Output <- matrix(data = NA, nrow =  7, ncol = 1, dimnames=list(c("Estimate","Error","pvalue","OR","lower CI","upper CI","AIC"), c(1)))
Output[1,1] <- Estimate
Output[2,1] <- Error
Output[3,1] <- Pvalue
Output[4,1] <- OR
Output[5,1] <- lower
Output[6,1] <- upper
Output[7,1] <- AIC
t(Output)
OR_Univariate_LogReg <- OR

#Test for Confounding#
Input <- WHEALS_META

ConfounderTest <- function(Diagnosis,Test,SignificantInput){
  parent.model=glm(Diagnosis ~ Test,
                   data=Input,
                   family=binomial(link="logit"))
  extended.model = glm(Diagnosis ~ Test + SignificantInput, 
                       data = Input, 
                       family = binomial(link="logit"))
  testfordeviance <- anova(parent.model,extended.model, test="Chi")
  pvalue <- testfordeviance[2,5]
  print(pvalue)
}

colnames(WHEALS_META)

PotentialConfounders <- names(Input[,c(6:13)])

Output<- matrix(data = NA, nrow =  8, ncol = 3, dimnames=list(c(1:8),c("PotentialConfounder", "pvalue","Confounding?")))

for (i in 1:8){
  var <- PotentialConfounders[i]
  Stats <- ConfounderTest(Input$Healthy, Input$Test.910.ngmg, Input[,var])
  Output[(i),2] <- Stats
  Output[(i),1] <- paste0(var,".AorA")
  Output[(i),3] <- ifelse(Stats<0.05,"YES","NO")
}
Confounder_WHEALS_910_AorA <- as.data.frame(Output)
write.table(Confounder_WHEALS_910_AorA,file="Confounders_WHEALS_910_AorA.txt",sep="\t",quote=F,row.names=F)

#NO CONFOUNDERS IDENTIFIED for 9,10 DiHOME and A or A#

#Test for Confounding by Change in OR#
ConfounderORChange <- function(Diagnosis,SignificantInput,Test){
  model = glm(Diagnosis ~ SignificantInput + Test, 
              data = Input, 
              family = binomial(link="logit"))
  summary<-summary(model)
  coefficients <- summary$coefficients
  Estimate <- coefficients[c(2),c(1)]
  Error <- coefficients[c(2),c(2)]
  Pvalue <- coefficients[c(2),c(4)]
  OR <- exp(Estimate)
  lower <- exp(Estimate-1.96*Error)
  upper <- exp(Estimate+1.96*Error)
  AIC <- model$aic
  ORChange <- 100*((OR/OR_SignificantInput_Diagnosis)-1)
  Confounder <- ifelse(ORChange>10,1,0)
  Output <- matrix(data = NA, nrow =  9, ncol = 1, dimnames=list(c("Estimate","Error","pvalue","OR","lower CI","upper CI","AIC","OR%Change","Confounder"), c(1)))
  Output[1,1] <- Estimate
  Output[2,1] <- Error
  Output[3,1] <- Pvalue
  Output[4,1] <- OR
  Output[5,1] <- lower
  Output[6,1] <- upper
  Output[7,1] <- AIC
  Output[8,1] <- ORChange
  Output[9,1] <- Confounder
  return(Output)
}

#Must Specify the OR for the Significant Test#
OR_Univariate_LogReg

#Select Potential Confounders of Test.1213.ngmg#
colnames(Input)
RiskFactors <- names(Input[,c(6:13)])

Output<- matrix(data = NA, nrow =  9, ncol = 8, dimnames=list(c("Estimate","Error","pvalue","OR","lower CI","upper CI","AIC","OR%Change","Confounder"),c(1:8)))

#Specify the OR of the Significant Test#
OR_SignificantInput_Diagnosis <- OR_Univariate_LogReg

for (i in 1:8){
  var <- RiskFactors[i]
  Stats <- ConfounderORChange(Input$Healthy, Input$Test.910.ngmg, Input[,var])
  Output[,(i)] <- Stats
  colnames(Output) <- paste0(RiskFactors,".AorA")
}
Confounder_Test.910.ngmg_Healthy <- as.data.frame(t(Output))
write.table(Confounder_Test.910.ngmg_Healthy,file="WHEALS_ConfounderORChange_Test.910.ngmg_AorA.txt",sep="\t",quote=F,row.names=T)

#No confounders identified for 9,10 DiHOME#
