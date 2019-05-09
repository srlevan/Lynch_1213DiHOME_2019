#Clear Envorinment#
rm(list=ls())

#load required packages#
library(lmerTest)
library(visreg)
library(ggplot2)
library(plyr)
library(gridExtra)
library(reshape2)

#Set the Working Directory
setwd("~/Box Sync/Lynch Lab/Code Respositories/Figure2_20190301/")

#Import csv with effects of 12,13 DiHOME-treated DCs on Tcell Subsets#
dihome=read.csv("INPUT_OtherLigands.csv")
dihome

#Input only single reagents#
Input<-droplevels(subset(dihome, Treatment!="9,10 DiHOME + 12,13 epOME" & Treatment!="9-HODE + 13-HODE"))
Input

Input$Treatment <- factor(Input$Treatment, c("No treat", "9,10 DiHOME- 50 uM", "12,13 epOME - 50 uM", 
                                             "9-HODE-100uM","13-HODE-100 uM"))

OtherTGraph <- ggplot() + 
  geom_boxplot(data=Input, aes(x=Treatment, y=Treg, fill=Treatment))+
  geom_jitter(data=Input, aes(x=Treatment, y=Treg, fill=Treatment), position=position_dodge(width=0.8), size=8)+
  theme_classic()+scale_fill_brewer(palette = "Greys", direction=1) + 
  theme(legend.position="none")+
  scale_shape_manual(values = c(18,16,15,17))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  ylim(0,3)+
  labs(title="OtherT_SubsetSamples", x="Treatment", y="OtherT Frequency")
OtherTGraph

ggsave(filename="GRAPH-OtherLigands.pdf", plot=OtherTGraph,  useDingbats=FALSE,width=11, height=6, units="in")

