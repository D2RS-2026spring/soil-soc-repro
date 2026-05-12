rm(list=ls())
## comparison of molecular C groups among samples
library(ggplot2)
library(multcomp)
library(gridExtra)
library(cowplot)
library(here)

# load NMR data 
data<-read.csv(here("data/raw/Hall_Ye_etal_NEON_NMR_biogeo.csv"),header=TRUE)
head(data)

###First, look at the relative abundance of C groups among all the samples
# make a vector with all the C groups determined by NMR data and the Baldock molecular mixing model, expressed as fractional abundance
test<-c(data$Carbohydrate_fraction,data$Protein_fraction,data$Lignin_fraction,data$Lipid_fraction,data$Char_fraction,data$Carbonyl_fraction)

length(data$Carbohydrate_fraction)

# make a vector to label each C group in the previous vector
trts<-c(rep("    Carbohydrate",42),rep("Protein",42),rep("Lignin",42),rep("Lipid",42),rep("Char",42),rep("Carbonyl",42))
trts <- as.factor(trts)
#make a data frame with the C group and label vectors
SOM<-data.frame(test,trts)
##test whether C groups  significantly differ among samples
lm1<-lm(test~trts, SOM)
anova(lm1)
#yes, there are differences
 h1<-glht(lm1, linfct = mcp(trts = "Tukey"))
# cld(h1)
#only carbonyl C is lower than the other C groups; the others are all similar

 #arrange factors
 SOM$trts<-factor(SOM$trts,levels=c("    Carbohydrate","Lignin","Lipid","Protein","Char","Carbonyl"))
 SOM$trts
#quartz(,4.5,4)
#windows(width=4.5, height=4)
plot_SOM<-ggplot(SOM, aes(trts,test))+
  geom_boxplot(outlier.shape = NA)+ # omit outliers so they don't show up twice with jitter
  geom_jitter(position=position_jitter(0.1), color="grey50")+
  theme_classic()+
  scale_x_discrete("C molecule")+
  scale_y_continuous("Fractional abundance")+
  theme(axis.text.x=element_text(angle = +45, hjust = 1))

# next, make a vector with all the underlying C moieties determined by NMR, expressed as fractional abundance
test1<-c(data$Alkyl_fraction, data$N_Alkyl_Methoxyl_fraction,	data$O_Alkyl_fraction,	data$Di_O_Alkyl_fraction, data$Aromatic_fraction,	data$Phenolic_fraction, data$Amide_Carboxyl_fraction)
# make a vector to label each C group in the previous vector
trts1<-c(rep("Alkyl",42),rep("N-Alkyl/Methoxyl",42),rep("O-Alkyl",42),rep("Di-O-Alkyl",42),rep("Aromatic",42),rep("Phenolic",42),rep("Amide/Carboxyl",42))
#make a data frame with the C group and label vectors
SOM1<-data.frame(test1,trts1)
##test whether C groups  significantly differ among samples
lm11<-lm(test~trts, SOM1)
anova(lm11)
#yes, there are differences
# h11<-glht(lm11, linfct = mcp(trts = "Tukey"))
# cld(h11)
head(SOM1)
tapply(SOM1$test1, SOM1$trts1, mean)
levels(SOM1$trts1)

SOM1$trts1<-factor(SOM1$trts,levels=c("O-Alkyl","Alkyl","Aromatic","Amide/Carboxyl","N-Alkyl/Methoxyl","Phenolic","Di-O-Alkyl"))

SOM2<-data.frame(SOM,rep(data$forest_vs_non),rep(data$Prescribed_fire))

plot_SOM_comp<-ggplot(SOM1, aes(trts1,test1))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitter(0.1), color="grey50")+
  theme_classic()+
  scale_x_discrete("C functional group")+
  scale_y_continuous("Fractional abundance")+
  theme(axis.text.x=element_text(angle = +45, hjust = 1))

###Make Figure 1
#quartz(,6.25,3.5)
windows(width=6.25, height=3.5)
# grid.arrange(plot_SOM_comp,plot_SOM,nrow=1)
p1<-plot_grid(plot_SOM_comp,plot_SOM,nrow=1,labels=c("a","b"))
print(p1)


names(SOM2)[names(SOM2)=="rep.data.forest_vs_non."]<-"Vegetation"
SOM2$Vegetation[SOM2$Vegetation==0]<-"Other"
SOM2$Vegetation[SOM2$Vegetation==1]<-"Forest"
SOM2$Vegetation<-factor(SOM2$Vegetation,levels=c("Forest","Other"))

names(SOM2)[names(SOM2)=="rep.data.Prescribed_fire."]<-"Fire"
SOM2$Fire[SOM2$Fire==0]<-"No"
SOM2$Fire[SOM2$Fire==1]<-"Yes"
SOM2$Fire<-factor(SOM2$Fire,levels=c("No","Yes"))

####
##Make Extended Data Fig 4 (revision)
#first panel
#quartz(,5,8)
#windows(width=5, height=8)
S4_a<-ggplot(SOM2, aes(trts,test,fill=Vegetation))+
  geom_boxplot(outlier.shape = NA)+
  scale_fill_brewer(palette="Pastel1")+
  geom_point(position=position_jitterdodge(),alpha=0.5)+ 
  #geom_jitter(position=position_jitter(0.1), color="grey50")+
  theme_classic()+
  scale_x_discrete("C functional group")+
  scale_y_continuous("Fractional abundance")+
  theme(axis.text.x=element_text(angle = +45, hjust = 1))

####
#second panel

S4_b<-ggplot(SOM2, aes(trts,test,fill=Fire))+
  geom_boxplot(outlier.shape = NA)+
  scale_fill_brewer(palette="Set2")+
  geom_point(position=position_jitterdodge(),alpha=0.5)+ 
  #geom_jitter(position=position_jitter(0.1), color="grey50")+
  theme_classic()+
  scale_x_discrete("C functional group")+
  scale_y_continuous("Fractional abundance")+
  theme(axis.text.x=element_text(angle = +45, hjust = 1))

##
#quartz(,4,8)
windows(width=4, height=8)

pED4<-plot_grid(S4_a,S4_b,nrow=2,labels=c("a","b"))
print(pED4)
