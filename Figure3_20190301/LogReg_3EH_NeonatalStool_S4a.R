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
WHEALS_META <- WHEALS[,c(22,23,3,7,14:21)]
WHEALS_3EH <- WHEALS_META[,c(1:4)]
head(WHEALS_3EH)

#Factor by condition and health status#
WHEALS_3EH$Healthy<- factor(WHEALS_3EH$Healthy, c("Healthy","AorA"))
WHEALS_3EH$Condition <- factor(WHEALS_3EH$Condition, c("Healthy","Atopic","Asthmatic","Atopic Asthma"))

#Convert to Long Form#
WHEALS_3EH_Long <- melt(WHEALS_3EH, id=c("Healthy","Condition","Lynch_ID"))
WHEALS_3EH_Long <- WHEALS_3EH_Long[,c(1:2,4:5)]


#Make Lists of Healthy and AorA values#
Raw_Healthy <- subset(WHEALS_3EH_Long, Healthy=="Healthy")
Raw_AorA <- subset(WHEALS_3EH_Long, Healthy=="AorA")

Raw_Healthy.List <- Raw_Healthy[,c(4)]
Raw_AorA.List <- Raw_AorA[,c(4)]

#Create a matrix with Healthy and AorA values#
x <- matrix(data = NA, nrow =  22, ncol = 2, dimnames=list(1:22, c("Healthy", "AorA")))
x[,c(1)]<-Raw_Healthy.List
x[c(1:19),c(2)]<-Raw_AorA.List

#Plot a Boxplot with these values on the logscale using the baseplots boxplot function#
#The ggplot boxplot function gives errors when data was plotted on the log scale#
pdf('3EH_boxplot.pdf')
boxplot(x,log="y")
dev.off()

#Summarize Oxylipins by Health Status#
#Calculate SEM#
sem <- function(x) {
  sem<-sd(x)/sqrt(length(x))
  return(sem)
}
WHEALS_3EH_Summary <- as.data.frame(ddply(WHEALS_3EH_Long, c("variable","Healthy"), summarise,
                                                  mean=mean(value), median=median(value), lower=mean(value)-sem(value), upper=mean(value)+sem(value)))


#Make a Scatterplot with logscale data using ggplot#
Scatterplot <- ggplot() + 
  geom_crossbar(data=WHEALS_3EH_Summary, aes(x=Healthy,y=median, ymin=median, ymax=median),stat = 'identity',position="dodge",colour="black", width=0.8)+
  geom_jitter(data=WHEALS_3EH_Long, aes(x=Healthy, y=value,shape=WHEALS_3EH_Long$Condition),position=position_jitter(0.2),size=5)+
  scale_shape_manual(values = c(18,17,15,16))+
  theme_classic() +scale_fill_brewer(palette = "Greys", direction=1)+ 
  labs(title="Epoxide Hydrolase Abundance", y="Copies/ng")+
  theme(legend.position="none")+
  scale_y_log10(limits=c(1000,3000000),breaks=c(1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,20000,30000,40000,50000,60000,70000,80000,90000,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000,2000000,3000000))
Scatterplot
ggsave(filename="GRAPH-3EH_NeonatalStool_scatter_narrow_Log.pdf", plot=Scatterplot, useDingbats=FALSE,width=6, height=8, units="in")

#In adobe illustrator layer the scatterplot and the boxplot#
#Scales will align when width of the Scatter is adjusted to 5.55 and the height and witdth of the boxplot figure is adjusted to 6.05 inches 5.20 inches#
#Widths will align when adjusted to 5.55 inches#

#Add Log3EH to DF#
WHEALS_META$Log2.3EH <- log2(WHEALS_META$Test.3EH.ng)

#Test the significance of the relationship between 3EH and Health Status using logistical regression#
Input<-WHEALS_META
head(Input)

model = glm(Input$Healthy ~ Input$Log2.3EH, 
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

#Test for Confounding using Chi-Squared Test of Deviance#
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

PotentialConfounders <- names(Input[,c(5:12)])

Output<- matrix(data = NA, nrow =  8, ncol = 3, dimnames=list(c(1:8),c("PotentialConfounder", "pvalue","Confounding?")))

for (i in 1:8){
  var <- PotentialConfounders[i]
  Stats <- ConfounderTest(Input$Healthy, Input$Log2.3EH, Input[,var])
  Output[(i),2] <- Stats
  Output[(i),1] <- paste0(var,".AorA")
  Output[(i),3] <- ifelse(Stats<0.05,"YES","NO")
}
Confounder_WHEALS_Log2.3EH_AorA <- as.data.frame(Output)
write.table(Confounder_WHEALS_Log2.3EH_AorA,file="Confounders_WHEALS_Log2.3EH_AorA.txt",sep="\t",quote=F,row.names=F)

#Chi Squared Test of Deviance -> NO CONFOUNDERS IDENTIFIED for 3EH and A or A#

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
RiskFactors <- names(Input[,c(5:12)])

Output<- matrix(data = NA, nrow =  9, ncol = 8, dimnames=list(c("Estimate","Error","pvalue","OR","lower CI","upper CI","AIC","OR%Change","Confounder"),c(1:8)))

#Specify the OR of the Significant Test#
OR_SignificantInput_Diagnosis <- OR_Univariate_LogReg

for (i in 1:8){
  var <- RiskFactors[i]
  Stats <- ConfounderORChange(Input$Healthy, Input$Log2.3EH, Input[,var])
  Output[,(i)] <- Stats
  colnames(Output) <- paste0(RiskFactors,".AorA")
}
Confounder_Log2.3EH_Healthy <- as.data.frame(t(Output))
write.table(Confounder_Log2.3EH_Healthy,file="WHEALS_ConfounderORChange_Test.1213.ngmg_AorA.txt",sep="\t",quote=F,row.names=T)

#No Potential Confounders Found for Log2.3EH by 10% OR Change#
