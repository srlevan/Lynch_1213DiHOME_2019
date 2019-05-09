#Clear the Environment#
rm(list=ls())

#Load Required Packages#
library(lmerTest)
library(visreg)
library(ggplot2)
library(plyr)
library(gridExtra)
library(reshape2)
library(dplyr)

#Set Working Directory#
setwd("~/Box Sync/Lynch Lab/Code Respositories/Figure3_20190301/")

#Load Files#
dihome=read.csv("INPUT-MouseSerumD22.csv")
dihome

#Simplify the Data Frame#
colnames(dihome)
Input <- dihome[,c(1:2,21,23)]

#Calculate p-values compared to 3EH_CRA#
Input$Group <- factor(Input$Group, c("3EH_CRA","CRA","3EH_PBS","PBS"))
x <- matrix(data = NA, nrow =  4, ncol = 3, dimnames=list(1:4, c("Estimate", "Standard Error", "P")))
colnames(Input)

for (i in 3:4){
  model <- lmer(Input[[i]] ~ Group + (1|Assay), data = Input)
  y <- summary(model)
  x[,1] <- y$coefficients[,1]
  x[,2] <- y$coefficients[,2]
  x[,3] <- y$coefficients[,5]
  print(colnames(Input[i]))
  print(x)
}

#Calculate p-values compared to PBS#
Input$Group <- factor(Input$Group, c("PBS","3EH_PBS","CRA","3EH_CRA"))
x <- matrix(data = NA, nrow =  4, ncol = 3, dimnames=list(1:4, c("Estimate", "Standard Error", "P")))
colnames(Input)

for (i in 3:4){
  model <- lmer(Input[[i]] ~ Group + (1|Assay), data = Input)
  y <- summary(model)
  x[,1] <- y$coefficients[,1]
  x[,2] <- y$coefficients[,2]
  x[,3] <- y$coefficients[,5]
  print(colnames(Input[i]))
  print(x)
}


#Visualize Data with LME#
colnames(Input)
DC<-names(Input)
Input$Group <- factor(Input$Group, c("PBS","3EH_PBS","CRA","3EH_CRA"))


#Model 12,13 DiHOME in Serum#
var <- DC[[3]]
model<- lmer(Input[[3]]~Group+(1|Assay),data=Input)
y.model<-summary(model)
y.model
coef(model)
plot<-visreg(model,"Group",type="conditional")
resfit<-plot$res
resfit["Assay"] <- Input$Assay[c(11:38)]
resfit
levels(resfit$Assay) <- c("18","16","17")
write.table(resfit, file=paste0(var,"_lme_graphpoints.txt"),sep="\t",quote=F,row.names=F)

  
#Model 9,10 DiHOME in Serum#
var <- DC[[4]]
model<- lmer(Input[[4]]~Group+(1|Assay),data=Input)
y.model<-summary(model)
y.model
coef(model)
plot<-visreg(model,"Group",type="conditional")
resfit<-plot$res
resfit["Assay"] <- Input$Assay[c(11:38)]
resfit
levels(resfit$Assay) <- c("18","16","17")
write.table(resfit, file=paste0(var,"_lme_graphpoints.txt"),sep="\t",quote=F,row.names=F)

for (i in 3:4) {  
  var <- DC[[i]]
  lmepoint<-read.table(file=paste0(var,"_lme_graphpoints.txt"),header=T,sep="\t")
  lmepoint$Group <- factor(lmepoint$Group, c("PBS","3EH_PBS","CRA","3EH_CRA"))
  lmepoint
  #Generate plot space with lmepoint data#
  Plot <- ggplot(data=lmepoint, aes(x=Group, y=visregRes)) + 
    geom_boxplot()+
    geom_jitter(data=lmepoint, aes(x=Group, y=visregRes), position=position_jitter(0.2), size=5, shape=lmepoint$Assay)+
    theme_classic() + scale_fill_brewer(palette = "Greys", direction=1)+
    guides(fill = guide_legend(reverse=FALSE))
  Plot
  ggsave(filename=paste0(var,".pdf"), plot=Plot,  useDingbats=FALSE, width=8.5, height=11, units="in")
}

  
