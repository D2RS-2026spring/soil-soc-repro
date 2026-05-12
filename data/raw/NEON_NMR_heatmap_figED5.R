#make Extended Data Fig 8
source(here("data/raw/NEON_NMR_overview_1.R"))

library(ggplot2)
library(multcomp)
library(gridExtra)
library(cowplot)
library(reshape)
library(viridis)
library(Hmisc)
library(GGally)
library(here)
## first run NEON_NMR_overview_1.R

# make a vector with all the C groups determined by NMR data and the Baldock molecular mixing model, expressed as fractional abundance
test<-c(data$Carbohydrate_fraction,data$Protein_fraction,data$Lignin_fraction,data$Lipid_fraction,data$Char_fraction,data$Carbonyl_fraction)

length(data$Carbohydrate_fraction)

# make a vector to label each C group in the previous vector
trts<-c(rep("    Carbohydrate",42),rep("Protein",42),rep("Lignin",42),rep("Lipid",42),rep("Char",42),rep("Carbonyl",42))
#make a data frame with the C group and label vectors
SOM<-data.frame(test,trts)

# next, make a vector with all the underlying C moieties determined by NMR, expressed as fractional abundance
test1<-c(data$Alkyl_fraction, data$N_Alkyl_Methoxyl_fraction,	data$O_Alkyl_fraction,	data$Di_O_Alkyl_fraction, data$Aromatic_fraction,	data$Phenolic_fraction, data$Amide_Carboxyl_fraction)
# make a vector to label each C group in the previous vector
trts1<-c(rep("Alkyl",42),rep("N-Alkyl/Methoxyl",42),rep("O-Alkyl",42),rep("Di-O-Alkyl",42),rep("Aromatic",42),rep("Phenolic",42),rep("Amide",42))
#make a data frame with the C group and label vectors
SOM1<-data.frame(test1,trts1)


### now plot heatmap of correlations among components

#C components
comps<-data1[,c("Carbohydrate_fraction","Protein_fraction","Lignin_fraction","Lipid_fraction","Char_fraction","Carbonyl_fraction")]
head(comps)
names(comps)<-c("Carbohydrate","Protein","Lignin","Lipid","Char","Carbonyl")

NMR_cor_vars<-as.matrix(cor(comps,use="pairwise.complete.obs"))
NMR_cor_mat<-melt(NMR_cor_vars)
names(NMR_cor_mat)[3]<-"r"

get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

NMR_cor_vars_upper<-get_upper_tri(NMR_cor_vars)
NMR_cor_mat<-melt(NMR_cor_vars_upper, na.rm=TRUE)
names(NMR_cor_mat)[3]<-"r"
str(NMR_cor_mat)
NMR_cor_mat$X2<-factor(NMR_cor_mat$X2, levels=c("Carbonyl","Char","Lipid","Lignin","Protein","Carbohydrate"))
NMR_cor_mat$X1<-factor(NMR_cor_mat$X1, levels=c("Carbonyl","Char","Lipid","Lignin","Protein","Carbohydrate"))

# NMR_cor_mat<-NMR_cor_mat[order(NMR_cor_mat$X1, NMR_cor_mat$X2),]
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

cormat<-reorder_cormat(NMR_cor_vars)
cormat_upper<-get_upper_tri(cormat)
NMR_cormat_upper<-melt(cormat_upper, na.rm=TRUE)
names(NMR_cormat_upper)[3]<-"r"

NMR_cormat_upper$X2<-factor(NMR_cormat_upper$X2, levels=c("Carbonyl","Lignin","Char","Lipid","Protein","Carbohydrate"))
NMR_cormat_upper$X1<-factor(NMR_cormat_upper$X1, levels=c("Carbonyl","Lignin","Char","Lipid","Protein","Carbohydrate"))


#get sig values
Pmat<-as.data.frame(rcorr(as.matrix(comps))$P)
r<-as.data.frame(rcorr(as.matrix(comps))$r)

str(Pmat)
#get the same format as for the matrices above
Pmat_red<-Pmat[1:6,1:6]

#replace non sig values with NA; threshold is p = 0.05/15
Pmat_red[Pmat_red>0.003 & Pmat_red<1]<-NA

#use Bonferroni correction (15 tests); replace sig values with 1 or 2
Pmat_red<-Pmat_red*15
# Pmat_red[Pmat_red< 0.0013 & Pmat_red> 0.001]<-3
Pmat_red[Pmat_red< 0.05 & Pmat_red> 0.01]<-4
Pmat_red[Pmat_red< 0.01 & Pmat_red> 0.001]<-3
Pmat_red[Pmat_red< 0.001 & Pmat_red> 0.0001]<-2
Pmat_red[Pmat_red< 0.0001 ]<-1


#now replace 1 and 2 with stars
Pmat_red[Pmat_red==1]<-"****"
Pmat_red[Pmat_red==2]<-"***"
Pmat_red[Pmat_red==3]<-"**"
Pmat_red[Pmat_red==4]<-"*"
#now eliminate the redundant stars from the list
Pmat_red[5,1]<-NA
Pmat_red[3,2]<-NA
Pmat_red[4,2]<-NA
Pmat_red[6,3]<-NA

#melt the matrix to get the same format as above
Pmat_red_melt<-melt(as.matrix(Pmat_red))
r_melt<-melt(as.matrix(r))
names(r_melt)[3]<-"r"
#paste the correlation vector to the p value data frame
#for color labels
Pmat_red_melt$Corr<-r_melt$r
Pmat_red_melt$Corr[Pmat_red_melt$Corr<0]<- -1
Pmat_red_melt$Corr[Pmat_red_melt$Corr>0]<- 1
Pmat_red_melt$Corr[Pmat_red_melt$Corr==1]<-"black"
Pmat_red_melt$Corr[Pmat_red_melt$Corr==-1]<-"white"

#round correlation values

NMR_cormat_upper$r<-round(NMR_cormat_upper$r,2)
NMR_cormat_upper$Corr<-NMR_cormat_upper$r
NMR_cormat_upper$Corr[NMR_cormat_upper$Corr<0]<- -1
NMR_cormat_upper$Corr[NMR_cormat_upper$Corr>0]<- 1
NMR_cormat_upper$Corr[NMR_cormat_upper$Corr==1]<-"black"
NMR_cormat_upper$Corr[NMR_cormat_upper$Corr==-1]<-"white"

#plot the corr matrix 

plot_corr_map<-ggplot(data=NMR_cormat_upper)+
  geom_tile(aes(x=X2,y=X1,fill=r))+
  #scale_fill_viridis(option="plasma")+
  scale_fill_viridis()+
  scale_y_discrete("")+
  scale_x_discrete("")+
  theme(axis.text.x=element_text(angle = +45, hjust = 1, size=14),axis.title.x=element_blank(),
        axis.text.y=element_text(size=14))+ 
  geom_text(data=Pmat_red_melt,                                           mapping=aes(x=X2,y=X1,label=value),
          color=Pmat_red_melt$Corr,
          size=10,
          nudge_y=0.2)+
  geom_text(data=NMR_cormat_upper,                                                 mapping=aes(x=X2,y=X1,label=r),
            color=NMR_cormat_upper$Corr)

#quartz(,6,5)
windows(width=6, height=5)
print(plot_corr_map)









