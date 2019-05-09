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

#Load File#
dihome=read.csv("INPUT-Mouse_Tregs.csv")
dihome

#Simplify Input#
Input <- dihome[,c(1:2,5,7,9,11)]

#Calculate p-values for compared to CRA_3EH#
Input$Treatment <- factor(Input$Treatment, c("CRA_3EH","CRA_Glycerol","PBS_3EH","PBS_Glycerol"))
x <- matrix(data = NA, nrow =  4, ncol = 3, dimnames=list(1:4, c("Estimate", "Standard Error", "P")))
for (i in 3:6){
  model <- lmer(Input[[i]] ~ Treatment + (1|Assay), data = Input)
  y <- summary(model)
  x[,1] <- y$coefficients[,1]
  x[,2] <- y$coefficients[,2]
  x[,3] <- y$coefficients[,5]
  print(colnames(Input[i]))
  print(x)
}

#Calculate p-values for compared to PBS_GlycerolH#
Input$Treatment <- factor(Input$Treatment, c("PBS_Glycerol","PBS_3EH","CRA_Glycerol","CRA_3EH"))
x <- matrix(data = NA, nrow =  4, ncol = 3, dimnames=list(1:4, c("Estimate", "Standard Error", "P")))
for (i in 3:6){
  model <- lmer(Input[[i]] ~ Treatment + (1|Assay), data = Input)
  y <- summary(model)
  x[,1] <- y$coefficients[,1]
  x[,2] <- y$coefficients[,2]
  x[,3] <- y$coefficients[,5]
  print(colnames(Input[i]))
  print(x)
}



#Model Live#
colnames(Input)
DC<-names(Input)
var <- DC[[3]]
model<- lmer(Input[[3]]~Treatment+(1|Assay),data=Input)
y.model<-summary(model)
y.model
coef(model)
plot<-visreg(model,"Treatment",type="conditional")
resfit<-plot$res
resfit["Assay"] <- Input$Assay
resfit
levels(resfit$Assay) <- c("18","16","17")
write.table(resfit, file=paste0(var,"_lme_graphpoints.txt"),sep="\t",quote=F,row.names=F)

#Model Tregs#
colnames(Input)
DC<-names(Input)
var <- DC[[4]]
model<- lmer(Input[[4]]~Treatment+(1|Assay),data=Input)
y.model<-summary(model)
y.model
coef(model)
plot<-visreg(model,"Treatment",type="conditional")
resfit<-plot$res
resfit["Assay"] <- Input$Assay
resfit
levels(resfit$Assay) <- c("18","16","17")
write.table(resfit, file=paste0(var,"_lme_graphpoints.txt"),sep="\t",quote=F,row.names=F)


#Model Resident#
colnames(Input)
DC<-names(Input)
var <- DC[[5]]
model<- lmer(Input[[5]]~Treatment+(1|Assay),data=Input)
y.model<-summary(model)
y.model
coef(model)
plot<-visreg(model,"Treatment",type="conditional")
resfit<-plot$res
resfit["Assay"] <- Input$Assay
resfit
levels(resfit$Assay) <- c("18","16","17")
write.table(resfit, file=paste0(var,"_lme_graphpoints.txt"),sep="\t",quote=F,row.names=F)


#Model CD3 Resident#
colnames(Input)
DC<-names(Input)
var <- DC[[6]]
model<- lmer(Input[[6]]~Treatment+(1|Assay),data=Input)
y.model<-summary(model)
y.model
coef(model)
plot<-visreg(model,"Treatment",type="conditional")
resfit<-plot$res
resfit["Assay"] <- Input$Assay
resfit
levels(resfit$Assay) <- c("18","16","17")
write.table(resfit, file=paste0(var,"_lme_graphpoints.txt"),sep="\t",quote=F,row.names=F)

for (i in 3:6) {  
  var <- DC[[i]]
  lmepoint<-read.table(file=paste0(var,"_lme_graphpoints.txt"),header=T,sep="\t")
  lmepoint$Treatment= factor(lmepoint$Treatment, c("PBS_Glycerol","PBS_3EH","CRA_Glycerol","CRA_3EH"))
  lmepoint

#plot the table in a dot plot#
#Generate plot space with lmepoint data#
Plot <- ggplot(data=lmepoint, aes(x=Treatment, y=visregRes)) + 
  geom_boxplot()+
  geom_jitter(data=lmepoint, aes(x=Treatment, y=visregRes), position=position_jitter(0.2), size=5, shape=lmepoint$Assay)+
    theme_classic() + scale_fill_brewer(palette = "Greys", direction=1)+
    guides(fill = guide_legend(reverse=FALSE))
Plot
ggsave(filename=paste0(var,".pdf"), plot=Plot,  useDingbats=FALSE, width=8.5, height=11, units="in")
}

