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

#Add Health Status and Condition#
HumanDB$Healthy <- ifelse(HumanDB$Diagnosis.AorA==1,"AorA","Healthy")
HumanDB$Healthy<- factor(HumanDB$Healthy, c("Healthy","AorA"))


#Reduce to only Columns of Interest#
colnames(HumanDB)
HumanDB_META <- HumanDB[,c(1,24,3,22,23,14:21)]
HumanDB_DiHOME <- HumanDB_META[,c(1:5)]
head(HumanDB_DiHOME)
HumanDB_META$Study <- factor(HumanDB_META$Study, c("WHEALS","TIPS"))

#Test the significance of the relationship between 12,13 DiHOME and Health Status using logistical regression#
Input<-HumanDB_META
head(Input)

#Plot Data by Study#
plot <- ggplot(data=HumanDB_META, aes(x=Study, y=Test.1213.ngmg))+
  geom_boxplot()+
  geom_jitter(aes(shape=HumanDB_META$Healthy), size=5)+
  theme_classic() +scale_fill_brewer(palette = "Greys", direction=1)+
  scale_shape_manual(values = c(18,17,15,16))+
  guides(fill = guide_legend(reverse=FALSE))+
  labs(title="WHEALS_1213_Summary", y="pg.mg")
plot
ggsave(filename="GRAPH-1213DiHOME_byStudy.pdf", plot=plot, useDingbats=FALSE,width=8.5, height=11, units="in")

model = glm(HumanDB_META$Study ~ HumanDB_META$Test.1213.ngmg, 
            data = Input, 
            family = binomial(link="logit"))
summary<-summary(model)
summary
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

