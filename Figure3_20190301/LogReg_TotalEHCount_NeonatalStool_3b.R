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
GeneCount=read.csv("INPUT-3EHCount_MetaData.csv")
head(GeneCount)
GeneCount$log2.TotalEHCount <- log2(GeneCount$TotalEHCount)

head(GeneCount)
colnames(GeneCount)

#Factor GeneCount by Condition and Health Status#
GeneCount$Healthy <- ifelse(GeneCount$Diagnosis.AorA==1,"AorA","Healthy")
GeneCount$Condition <- ifelse(GeneCount$Diagnosis.Atopy+GeneCount$Diagnosis.Asthma>1,"Atopic Asthmatic",ifelse(GeneCount$Diagnosis.Atopy>0,"Atopic",ifelse(GeneCount$Diagnosis.Asthma>0,"Asthmatic","Healthy")))
GeneCount$Healthy<- factor(GeneCount$Healthy, c("Healthy","AorA"))
GeneCount$Condition <- factor(GeneCount$Condition, c("Healthy","Atopic","Asthmatic","Atopic Asthmatic"))

#Select only Gene Count and Health Status#
head(GeneCount)
EHCountDF <- GeneCount[,c(4,21,22)]

#Factor EHCountDF by Condition and Health Status#
EHCountDF$Healthy<- factor(EHCountDF$Healthy, c("Healthy","AorA"))
EHCountDF$Condition <- factor(EHCountDF$Condition, c("Healthy","Atopic","Asthmatic","Atopic Asthmatic"))


#Factor Summary of EH Gene Count#
Summary_EHCountDF$Healthy<-factor(Summary_EHCountDF$Healthy,c("Healthy","AorA"))
EHCountDF

#Graph EH Gene Count as a Boxplot Graph#
Boxplot <- ggplot() + 
  geom_boxplot(data=EHCountDF, aes(x=Healthy, y=EHCountDF$TotalEHCount))+
  geom_jitter(data=EHCountDF, aes(x=Healthy, y=EHCountDF$TotalEHCount, shape=EHCountDF$Condition),position=position_jitter(0.2),size=5)+
  theme_classic() +scale_fill_brewer(palette = "Greys", direction=1)+
  scale_shape_manual(values = c(18,17,15,16))+
  guides(fill = guide_legend(reverse=FALSE))+
  labs(title="Epoxide Hydrolase Abundance", y="Copies/ng")
Boxplot
ggsave(filename="GRAPH-EHCount_NeonatalStool.pdf", plot=Boxplot,  useDingbats=FALSE,width=8.5, height=11, units="in")

#Test Relationship Between 3EH Gene Count and A or A with Simple Logistical Regression#
Input<-GeneCount
head(GeneCount)

model = glm(Input$Healthy ~ Input$log2.TotalEHCount, 
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
#Test for Confounding by Chi-squared test of deviance#
Input <- GeneCount

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

#Select potential confounders#
colnames(GeneCount)
PotentialConfounders <- names(Input[,c(12:19)])

Output<- matrix(data = NA, nrow =  8, ncol = 3, dimnames=list(c(1:8),c("PotentialConfounder", "pvalue","Confounding?")))
for (i in 1:8){
  var <- PotentialConfounders[i]
  Stats <- ConfounderTest(Input$Diagnosis.AorA, Input$log2.TotalEHCount, Input[,var])
  Output[(i),2] <- Stats
  Output[(i),1] <- paste0(var,".AorA")
  Output[(i),3] <- ifelse(Stats<0.05,"YES","NO")
}
Confounder_EHCount_AorA <- as.data.frame(Output)
write.table(Confounder_EHCount_AorA,file="Confounders_EHCount_AorA.txt",sep="\t",quote=F,row.names=F)

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

#Select Potential Confounders of Test.Log2.3EH.ngmg#
colnames(Input)
RiskFactors <- names(Input[,c(12:19)])

Output<- matrix(data = NA, nrow =  9, ncol = 8, dimnames=list(c("Estimate","Error","pvalue","OR","lower CI","upper CI","AIC","OR%Change","Confounder"),c(1:8)))

#Specify the OR of the Significant Test#
OR_SignificantInput_Diagnosis <- OR_Univariate_LogReg

for (i in 1:8){
  var <- RiskFactors[i]
  Stats <- ConfounderORChange(Input$Healthy, Input$log2.TotalEHCount, Input[,var])
  Output[,(i)] <- Stats
  colnames(Output) <- paste0(RiskFactors,".AorA")
}
Confounder_TotalEHCount_Healthy <- as.data.frame(t(Output))
write.table(Confounder_TotalEHCount_Healthy,file="WHEALS_ConfounderORChange_TotalEHCount_AorA.txt",sep="\t",quote=F,row.names=T)

#Risk.ExclBfed Identified as a potential confounder#

#Count Data#
head(Input)

AorA <- subset(GeneCount, Healthy =="AorA")
Healthy <- subset(GeneCount, Healthy=="Healthy")

Input <- GeneCount
RiskFactors <- names(Input[,c(12:18)])

Summary <- matrix(data=NA, nrow=2, ncol=7,dimnames=list(c("N","Percent"),c(1:7)))
for (i in 1:7){
  var <- RiskFactors[i]
  Summary[1,(i)] <- sum(Input[,var])
  Summary[2,(i)] <- 100*sum(Input[,var])/length(Input[,var])  
  colnames(Summary) <- paste0(RiskFactors,".Summary")
}   
Summary 

median(Input$Age.Stool)
