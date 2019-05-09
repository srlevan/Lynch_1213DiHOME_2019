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
#ADD LOG BASE 2 TO DATAFRAME#
HumanDB$Test.log2.3EH <- log2(HumanDB$Test.3EH.ng)
HumanDB$Test.1213.ngmg <- HumanDB$Test.1213.pgmg/1000

#Select only WHEALS Data#
WHEALS <- subset(HumanDB, HumanDB$Study=="WHEALS")

#Select on the TIPS Data#
TIPS <- subset(HumanDB, HumanDB$Study=="TIPS")

#Begin by Characterizing the WHEALS Dataset#

#Seperate Catagorical and Continuous Cohort Characteristics#
colnames(WHEALS)
CatWHEALS <- WHEALS[,c(1:3,10:13,14:20)]
ContWHEALS <- WHEALS[,c(1:3,10:13,21)]

#Summarize Catagorical Data for Whole Cohort#
Input<-CatWHEALS

RiskFactors <- names(Input[,c(8:14)])

Summary <- matrix(data=NA, nrow=2, ncol=7,dimnames=list(c("N","Percent"),c(1:7)))
for (i in 1:7){
  var <- RiskFactors[i]
  Summary[1,(i)] <- sum(Input[,var])
  Summary[2,(i)] <- 100*sum(Input[,var])/length(Input[,var])  
  colnames(Summary) <- paste0(RiskFactors,".Summary")
}   
Summary 


#Format of Table 
#     Disease+ Disease-  Total
# Test+ TP        FP      a+c
# Test- FN        TN      b+d
#   a+b       c+d     a+b+c+d
head(CatWHEALS)
#Make a Function to Calculate Table Statistics#
Table <- function(Diagnosis,Test){
  Input$TP <- ifelse(Test+Diagnosis=="2",1,0)
  TP <- sum(Input$TP)
  Input$TN <- ifelse(Test+Diagnosis=="0",1,0)
  TN <- sum(Input$TN)
  Input$FP <- ifelse(Test-Diagnosis=="1",1,0)
  FP <- sum(Input$FP)
  Input$FN <- ifelse(Test-Diagnosis<0,1,0)
  FN <- sum(Input$FN)
  Output <- matrix(data = NA, nrow =  4, ncol = 1, dimnames=list(c("TP","FP","FN","TN"), c(1)))
  Output[c(1),] <- c(TP)
  Output[c(2),] <- c(FP)
  Output[c(3),] <- c(FN)
  Output[c(4),] <- c(TN)
  return(Output)
}

TableStats <- function(Diagnosis,Test){
  Input$TP <- ifelse(Test+Diagnosis=="2",1,0)
  TP <- sum(Input$TP)
  Input$TN <- ifelse(Test+Diagnosis=="0",1,0)
  TN <- sum(Input$TN)
  Input$FP <- ifelse(Test-Diagnosis=="1",1,0)
  FP <- sum(Input$FP)
  Input$FN <- ifelse(Test-Diagnosis<0,1,0)
  FN <- sum(Input$FN)
  Total <- length(Diagnosis)
  Fishertest<-fisher.test(matrix(c(TP,FP,FN,TN),nrow=2),or=1)
  Output <- matrix(data = NA, nrow =  10, ncol = 1, dimnames=list(c("Disease+","Disease+Cat+","Percent","Disease-","Disease-Cat+", "Percent","Fisher OR estimate","lower CI", "upper CI","p-value"), c(1)))
  Output[c(1),] <- c(TP)+c(FN)
  Output[c(2),] <- c(TP)
  Output[c(3),] <- c(TP)/(c(TP)+c(FN))
  Output[c(4),] <- c(FP)+c(TN)
  Output[c(5),] <- c(FP)
  Output[c(6),] <- c(FP)/(c(FP)+c(TN))
  Output[c(7),] <- as.numeric(Fishertest$estimate)
  Output[c(8),] <- Fishertest$conf.int[1]
  Output[c(9),] <- Fishertest$conf.int[2]
  Output[c(10),] <- Fishertest$p.value
  return(Output)
}

#Calculate Table Statistics and p value using Fisher Exact Test for all Catagories for AorA#
Input<-CatWHEALS
head(Input)

RiskFactors <- names(Input[,c(8:14)])

Output<- matrix(data = NA, nrow =  10, ncol =7, dimnames=list(c("Disease+","Disease+Cat+","Percent","Disease-","Disease-Cat+", "Percent","Fisher OR estimate","lower CI", "upper CI","p-value"),c(1:7)))

for (i in 1:7){
  var <- RiskFactors[i]
  Stats <- TableStats(Input$Diagnosis.AorA, Input[,var])
  Output[,(i)] <- Stats
  colnames(Output) <- paste0(RiskFactors,".AorA")
}                                                    

AorACatChar <-t(Output)
write.table(AorACatChar,file="WHEALS_AorACatagoricalCharacteristics.txt",sep="\t",quote=F,row.names=T)

#Calculate statistics and Fisher estimates and pvalues for Atopy#
for (i in 1:7){
  var <- RiskFactors[i]
  Stats <- TableStats(Input$Diagnosis.Atopy, Input[,var])
  Output[,(i)] <- Stats
  colnames(Output) <- paste0(RiskFactors,".Atopy")
}                                                    

AtopyCatChar <-t(Output)
write.table(AtopyCatChar,file="WHEALS_AtopyCatagoricalCharacteristics.txt",sep="\t",quote=F,row.names=T)

#Calculate statistics and Fisher estimates and pvalues for Asthma#
for (i in 1:7){
  var <- RiskFactors[i]
  Stats <- TableStats(Input$Diagnosis.Asthma, Input[,var])
  Output[,(i)] <- Stats
  colnames(Output) <- paste0(RiskFactors,".Asthma")
}                                                    

AsthmaCatChar <-t(Output)
write.table(AsthmaCatChar,file="WHEALS_AsthmaCatagoricalCharacteristics.txt",sep="\t",quote=F,row.names=T)

#Calculate statistics and Fisher estimates and pvalues for Atopy#
for (i in 1:7){
  var <- RiskFactors[i]
  Stats <- TableStats(Input$Diagnosis.Atopy, Input[,var])
  Output[,(i)] <- Stats
  colnames(Output) <- paste0(RiskFactors,".Atopy")
}                                                    

AtopyCatChar <-t(Output)
write.table(AtopyCatChar,file="WHEALS_AtopyCatagoricalCharacteristics.txt",sep="\t",quote=F,row.names=T)

#Assess Age for Subset#
colnames(ContWHEALS)

#Calculate Median and SD for all samples in subset#
median(ContWHEALS$Age.Stool) 
sd(ContWHEALS$Age.Stool)


#Calculate Median, SD, and Student T for AorA by Age#
AorA_Age <- as.data.frame(ddply(ContWHEALS, c("Diagnosis.AorA"), summarise,
                                N=length(Diagnosis.AorA), median=median(Age.Stool), stdev=sd(Age.Stool)))
AorA_Age_p <- pairwise.t.test(ContWHEALS$Age.Stool, ContWHEALS$Diagnosis.AorA, alternative = "two.sided", pool.sd = T, paired = FALSE)
AorA_Age$pvalue <- as.numeric(AorA_Age_p[c(3)])

#Calculate Median, SD, and Student T for Atopy by Age#
Atopy_Age <- as.data.frame(ddply(ContWHEALS, c("Diagnosis.Atopy"), summarise,
                                N=length(Diagnosis.Atopy), median=median(Age.Stool), stdev=sd(Age.Stool)))
Atopy_Age_p <- pairwise.t.test(ContWHEALS$Age.Stool, ContWHEALS$Diagnosis.Atopy, alternative = "two.sided", pool.sd = T, paired = FALSE)
Atopy_Age$pvalue <- as.numeric(Atopy_Age_p[c(3)])
#Calculate Median, SD, and Student T for Asthma by Age#
Asthma_Age <- as.data.frame(ddply(ContWHEALS, c("Diagnosis.Asthma"), summarise,
                                 N=length(Diagnosis.Asthma), median=median(Age.Stool), stdev=sd(Age.Stool)))
Asthma_Age_p <- pairwise.t.test(ContWHEALS$Age.Stool, ContWHEALS$Diagnosis.Asthma, alternative = "two.sided", pool.sd = T, paired = FALSE)
Asthma_Age$pvalue <- as.numeric(Asthma_Age_p[c(3)])

#Make a Dataframe with All Age Statistics
AgeData <- cbind(Asthma_Age, Atopy_Age)
AgeData <- as.data.frame(cbind(AgeData, AorA_Age))
write.table(AgeData,file="WHEALS_AgeData.csv",sep=",",quote=F,row.names=F)




#Simple Logistical Regression#
SimpleLogReg <- function(Diagnosis,Test){
  model = glm(Diagnosis ~ Test, 
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
  Output <- matrix(data = NA, nrow =  6, ncol = 1, dimnames=list(c("Estimate","Error","OR","OR lower", "OR upper","pvalue"), c(1)))
  Output[1,1] <- Estimate
  Output[2,1] <- Error
  Output[3,1] <- OR
  Output[4,1] <- lower
  Output[5,1] <- upper
  Output[6,1] <- Pvalue
  return(Output)
}

#Run Simple Logistical Regression#
#Analyze Effect of Individual Variables on AorA Outcomes#
#USE Test.1213.ngmg and Test.Log2.3EH#
Input<-WHEALS
colnames(WHEALS)
RiskFactors <- names(Input[,c(22,23,14:21)])

Output<- matrix(data = NA, nrow =  6, ncol = 10, dimnames=list(c("Estimate","Error","OR","OR lower", "OR upper","pvalue"),c(1:10)))

for (i in 1:10){
  var <- RiskFactors[i]
  Stats <- SimpleLogReg(Input$Diagnosis.AorA, Input[,var])
  Output[,(i)] <- Stats
  colnames(Output) <- paste0(RiskFactors,".AorA")
}
SimpLogReg_AorA <- as.data.frame(t(Output))
write.table(SimpLogReg_AorA,file="WHEALS_SimpLogReg_AorA.txt",sep="\t",quote=F,row.names=T)

#Identify Potentially Significant Associations#
Significant <- subset(SimpLogReg_AorA, pvalue<0.1)

#Test for Confounding of Test.1213.ngmg#
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

colnames(Input)

PotentialConfounders <- names(Input[,c(14:21)])

Output<- matrix(data = NA, nrow =  8, ncol = 3, dimnames=list(c(1:8),c("PotentialConfounder", "pvalue","Confounding?")))

for (i in 1:8){
  var <- PotentialConfounders[i]
  Stats <- ConfounderTest(Input$Diagnosis.AorA, Input$Test.1213.ngmg, Input[,var])
  Output[(i),2] <- Stats
  Output[(i),1] <- paste0(var,".AorA")
  Output[(i),3] <- ifelse(Stats<0.05,"YES","NO")
}
Confounder_Test.1213.ngmg_AorA <- as.data.frame(Output)
write.table(Confounder_Test.1213.ngmg_AorA,file="WHEALS_ConfounderTest_log1213_AorA.txt",sep="\t",quote=F,row.names=F)

#Test for Confounding of Test.Log2.3EH#
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

colnames(Input)

PotentialConfounders <- names(Input[,c(14:21)])

Output<- matrix(data = NA, nrow =  8, ncol = 3, dimnames=list(c(1:8),c("PotentialConfounder", "pvalue","Confounding?")))

for (i in 1:8){
  var <- PotentialConfounders[i]
  Stats <- ConfounderTest(Input$Diagnosis.AorA, Input$Test.log2.3EH, Input[,var])
  Output[(i),2] <- Stats
  Output[(i),1] <- paste0(var,".AorA")
  Output[(i),3] <- ifelse(Stats<0.05,"YES","NO")
}
Confounder_Test.Log2.3EH_AorA <- as.data.frame(Output)
write.table(Confounder_Test.Log2.3EH_AorA,file="WHEALS_Confounders_Log2.3EH_AorA.txt",sep="\t",quote=F,row.names=F)

#Verify that potential confounders do not change the Odds Ratio#
#Want to add potential confounders one at a time and check that OR doesn't cance by more than 10%#

#Test for Confounding Logistical Regression#
#Potential Confounders are Indicated by a 1#
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
Significant
OR_Test.Log2.3EH <- Significant[c(1),c(3)]
OR_Test.1213.ngmg <- Significant[c(2),c(3)]

#Select Potential Confounders of Test.1213.ngmg#
colnames(Input)
RiskFactors <- names(Input[,c(14:21)])

Output<- matrix(data = NA, nrow =  9, ncol = 8, dimnames=list(c("Estimate","Error","pvalue","OR","lower CI","upper CI","AIC","OR%Change","Confounder"),c(1:8)))

#Specify the OR of the Significant Test#
OR_SignificantInput_Diagnosis <- OR_Test.1213.ngmg

for (i in 1:8){
  var <- RiskFactors[i]
  Stats <- ConfounderORChange(Input$Diagnosis.AorA, Input$Test.1213.ngmg, Input[,var])
  Output[,(i)] <- Stats
  colnames(Output) <- paste0(RiskFactors,".AorA")
}
Confounder_Test.1213.ngmg_AorA <- as.data.frame(t(Output))
write.table(Confounder_Test.1213.ngmg_AorA,file="WHEALS_ConfounderORChange_Test.1213.ngmg_AorA.txt",sep="\t",quote=F,row.names=T)

#Include Confounder in Model#
#Test Relationship Between Test.1213.ngmg, Risk.Black, and A or A with Multiple Logistical Regression#
model = glm(Input$Diagnosis.AorA ~ Input$Test.1213.ngmg + Input$Risk.Black, 
            data = Input, 
            family = binomial(link="logit"))
summary<-summary(model)

#Determine pvalue#
coefficients <- summary$coefficients
Pvalue <- coefficients[c(2),c(4)]
Pvalue

#Select Potential Confounders of Test.Log2.3EH#
colnames(Input)
RiskFactors <- names(Input[,c(14:21)])

Output<- matrix(data = NA, nrow =  9, ncol = 8, dimnames=list(c("Estimate","Error","pvalue","OR","lower CI","upper CI","AIC","OR%Change","Confounder"),c(1:8)))

#Specify the OR of the Significant Test#
OR_SignificantInput_Diagnosis <- OR_Test.Log2.3EH
Significant
for (i in 1:8){
  var <- RiskFactors[i]
  Stats <- ConfounderORChange(Input$Diagnosis.AorA, Input$Test.log2.3EH, Input[,var])
  Output[,(i)] <- Stats
  colnames(Output) <- paste0(RiskFactors,".AorA")
}
Confounder_Test.Log2.3EH_AorA <- as.data.frame(t(Output))
write.table(Confounder_Test.Log2.3EH_AorA,file="WHEALS_ConfounderORChange_Test.Log2.3EH_AorA.txt",sep="\t",quote=F,row.names=T)



#Check that the Analysis is not Driven by Outliers#
#Catagorize Biological Markers into Low and High Groups#
#Plot Distribution of 1213 DiHOME and the Q3 Marker#
head(WHEALS)
Input<-WHEALS
Quantile_1213 <- quantile(Input$Test.1213.ngmg)
Plot_1213 <- ggplot(data=Input, aes(x=Test.1213.ngmg))+
  geom_histogram(aes(y=..density..),binwidth = 0.1, fill=NA, color="black")+
  geom_density(alpha=0.2, fill="#FF6666")+
  geom_vline(aes(xintercept=Quantile_1213[4]),color="red", linetype="dashed", size=1)+ #Add Q3
  theme_classic() +scale_fill_brewer(palette = "Greys", direction=1)
Plot_1213
ggsave(filename="WHEALS_Test.1213.ngmg_Histogram.pdf", plot=Plot_1213, useDingbats=FALSE,width=11, height=8.5, units="in")


#Plot Distribution of EH DiHOME and the Q3 Marker#
Quantile_EH <- quantile(Input$Test.3EH.ng)
Plot_EH <- ggplot(data=Input, aes(x=Test.3EH.ng))+
  geom_histogram(aes(y=..density..),binwidth = 10000, fill=NA, color="black")+
  geom_density(alpha=0.2, fill="#FF6666")+
  geom_vline(aes(xintercept=Quantile_EH[4]),color="red", linetype="dashed", size=1)+ #Add Q3
  theme_classic() +scale_fill_brewer(palette = "Greys", direction=1)
Plot_EH
ggsave(filename="WHEALS_Test.3EH.ng_Histogram.pdf", plot=Plot_EH, useDingbats=FALSE,width=11, height=8.5, units="in")



#Generate Forrest Plots for AorA#
AorA_Graph<-SimpLogReg_AorA

#Add Adjusted OR for Test.1213.ngmg#
Confounder_Test.1213.ngmg_AorA
ADJUSTED_Test.1213.ngmg<-subset(Confounder_Test.1213.ngmg_AorA, Confounder>0)
ADJUSTED_Test.1213.ngmg<-ADJUSTED_Test.1213.ngmg[,c(1,2,4,5,6,3)]
names(ADJUSTED_Test.1213.ngmg)[4]<-"OR lower"
names(ADJUSTED_Test.1213.ngmg)[5]<-"OR upper"
rownames(ADJUSTED_Test.1213.ngmg) <- "AdjustedOR_Test.1213.ngmg"

AorA_Graph <- rbind(ADJUSTED_Test.1213.ngmg,SimpLogReg_AorA)
AorA_Graph

AorA_Graph$Test<-rownames(AorA_Graph)
AorA_Graph <- AorA_Graph[order(AorA_Graph$OR),]
rev(AorA_Graph$Test)

#Use Reversed List to Factor by Test#
AorA_Graph$Test <- factor(AorA_Graph$Test, c("AdjustedOR_Test.1213.ngmg","Test.1213.ngmg.AorA","Test.log2.3EH.AorA","Risk.MomAsthma.AorA","Risk.NoPets.AorA","Risk.Csect.AorA","Risk.MomSmoke.AorA","Age.Stool.AorA","Risk.Male.AorA","Risk.Black.AorA","Risk.ExclBfed.AorA" ))
AorA_Graph <- droplevels(subset(AorA_Graph, Test!="Age.Stool.AorA" & Test!="Test.1213.ngmg.AorA"))
Graph <- ggplot(data=AorA_Graph)+
  geom_point(aes(x=AorA_Graph$Test, y=AorA_Graph$OR))+
  geom_errorbar(aes(x=AorA_Graph$Test, ymin=AorA_Graph$'OR lower', ymax=AorA_Graph$'OR upper'))+
  geom_hline(yintercept = 1, linetype="dotted")+
  theme_classic() +scale_fill_brewer(palette = "Greys", direction=1)+
  scale_y_log10(limits=c(0.03,50),breaks=c(0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,20,30,40,50))+
  coord_flip()
Graph
ggsave(filename="WHEALS_AorA_ORs.pdf", plot=Graph, useDingbats=FALSE,width=11, height=8.5, units="in")


