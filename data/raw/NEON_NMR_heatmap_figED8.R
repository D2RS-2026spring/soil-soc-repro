# heatmap of  NMR C functional group correlations with predictors
# Extended Data Fig 8
library(viridis)
library(Hmisc)
library(GGally)
library(reshape)
# first run NEON_NMR_1.R, NEON_NMR_stats_RC1.R, NEON_NMR_boxplot_figS2.R
head(datared)
#untransformed data
datared1<-datared[,-1]
# 
# #log-transformed data (selected vars)
# datared_log<-data1[,c("Oxalate_Fe_mg_g_log10",
#                       "Oxalate_Al_mg_g_ICP_log10",
#                       "Fe_crystalline_log10",
#                       "phH2o",
#                       "Ca_Na2SO4_mg_g_log10",
#                       "live_root_lessthan_4mm_approx_mg_cm3",
#                       "wgt_live_root_lessthan_4mm_CN_ratio_0_30cm",
#                       "MAP_minus_PET_mm",
#                       "MAT_C",
#                       "fungi_meanCopyNumber_ITS_log10",
#                       "bact_meanCopyNumber_16S_log10",
#                       "fung_bact_ratio_log10",
#                       "clayTotal_percent",
#                       "estimatedOC_g_kg",
#                       "OC_N_ratio",
#                       "nitrogenTot_g_kg")]
# names(datared_log)<-c("Feo","Alo","Fed-o","pH","Ca","Fine roots","Fine root C:N","MAP-PET","MAT","ITS","16S","ITS:16S","Clay","SOC","SOC:N","N")
datared_log<-data1[,c("Oxalate_Fe_mg_g_log10",
                      "Oxalate_Al_mg_g_ICP_log10",
                      # "Oxalate_Al_Fe_umol_g_log10",
                      "Fe_crystalline_log10",
                      "phH2o",
                      #"Ca_Na2SO4_mg_g_log10",
                      "Ca_plus_Mg_umol_log10",
                      "live_root_lessthan_4mm_approx_mg_cm3",
                      "wgt_live_root_lessthan_4mm_CN_ratio_0_30cm",
                      "MAP_minus_PET_mm",
                      "MAT_C",
                      "fungi_meanCopyNumber_ITS_log10",
                      "bact_meanCopyNumber_16S_log10",
                      "fung_bact_ratio_log10",
                      "clayTotal_percent",
                      "estimatedOC_g_kg",
                      "OC_N_ratio",
                      "nitrogenTot_g_kg",
                      "foliage_ligninPercent",
                      "foliage_CtoN",
                      "litter_ligninPercent",
                      "litter_CtoN",
                      "horizonTopDepth",
                      "horizon_depth",
                      "forest_vs_non",
                      "Prescribed_fire"
)]

names(datared_log)<-c("Feo",
                      "Alo",
                      # "Alo+Feo",
                      "Fed-o",
                      "pH",
                      #"Ca",
                      "Ca+Mg","Fine roots","Fine root C:N","MAP-PET","MAT","ITS","16S","ITS:16S","Clay","SOC","SOC:N","Soil N","Foliar lignin","Foliar C:N","Litter lignin","Litter C:N","O horizon (cm)","Horizon depth (cm)",
                      "Forest",
                      "Fire")




comps<-data1[,c("Alkyl_fraction","N_Alkyl_Methoxyl_fraction","O_Alkyl_fraction","Di_O_Alkyl_fraction","Aromatic_fraction","Phenolic_fraction","Amide_Carboxyl_fraction")]
head(comps)
names(comps)<-c("Alkyl","N-alkyl/methoxyl","O-alkyl","Di-O-alkyl","Aromatic","Phenolic","Amide/carboxyl")

# NMR_cor_vars<-as.matrix(cor(comps,datared,use="pairwise.complete.obs"))
# NMR_cor_mat<-melt(NMR_cor_vars)
# names(NMR_cor_mat)[3]<-"r"

NMR_cor_vars_log<-as.matrix(cor(comps,datared_log,use="pairwise.complete.obs"))
NMR_cor_mat_log<-melt(NMR_cor_vars_log)
names(NMR_cor_mat_log)[3]<-"r"

#make matrix of significance of correlations, using rcorr function in Hmisc
#extract the P-value matrix
Pmat<-as.data.frame(rcorr(as.matrix(comps),as.matrix(datared_log))$P)
str(Pmat)
#get the same format as for the matrices above
Pmat_red<-Pmat[1:7,8:length(Pmat)]


#replace non sig values with NA; threshold is p = 0.1/(7*16)
Pmat_red[Pmat_red>0.0009 & Pmat_red<1]<-NA
Pmat_red<-Pmat_red*56
#use Bonferroni correction; replace sig values with 1 or 2 or 3
Pmat_red[Pmat_red< 0.05 & Pmat_red> 0.01]<-4
Pmat_red[Pmat_red< 0.01 & Pmat_red> 0.001]<-3
Pmat_red[Pmat_red< 0.001 & Pmat_red> 0.0001]<-2
Pmat_red[Pmat_red< 0.0001 ]<-1

#now replace with stars
Pmat_red[Pmat_red==1]<-"****"
Pmat_red[Pmat_red==2]<-"***"
Pmat_red[Pmat_red==3]<-"**"
Pmat_red[Pmat_red==4]<-"*"



#melt the matrix to get the same format as above
Pmat_red_melt<-melt(as.matrix(Pmat_red))

#order factor levels for plotting
levels(NMR_cor_mat_log$X2)

# NMR_cor_mat_log$X2<-factor(NMR_cor_mat_log$X2, levels=c(
#   "MAT",
#   "MAP-PET",
#   "pH",
#   "Clay",
#   "SOC",
#   "N",
#   "SOC:N",
#   "Ca",
#   "Alo",
#   "Feo",
#   "Fed-o",
#   "Fine roots",
#   "Fine root C:N",
#   "16S",
#   "ITS",
#   "ITS:16S"
# ))

NMR_cor_mat_log$X2<-factor(NMR_cor_mat_log$X2, levels=c(
  "MAT",
  "MAP-PET",
  "Forest",
  "Fire",
  "pH",
  "Clay",
  "SOC",
  "Fine roots",
  "Soil N",
  "SOC:N",
  "Fine root C:N",
  "Foliar C:N",
  "Litter C:N",
  "Foliar lignin",
  "Litter lignin",
  "O horizon (cm)",
  "Horizon depth (cm)",
  # "Ca",
  "Ca+Mg",
  "Alo",
  "Feo",
  "Alo+Feo",
  "Fed-o",
  "16S",
  "ITS",
  "ITS:16S"
))
#paste the correlation vector to the p value data frame
#for color labels
Pmat_red_melt$Corr<-NMR_cor_mat_log$r
Pmat_red_melt$Corr[Pmat_red_melt$Corr<0]<- -1
Pmat_red_melt$Corr[Pmat_red_melt$Corr>0]<- 1
Pmat_red_melt$Corr[Pmat_red_melt$Corr==1]<-"black"
Pmat_red_melt$Corr[Pmat_red_melt$Corr==-1]<-"white"

#plot the corr matrix (version with several log-transformed vars)

#Supplementary Fig 10

quartz(,11,3.5)
ggplot(data=NMR_cor_mat_log)+
  geom_tile(aes(x=X2,y=X1,fill=r))+
  scale_fill_viridis(option="viridis",limits=c(-0.63,0.72))+
  scale_x_discrete("", labels=c(
    "MAT",
    "MAP-PET",
    "Forest",
    "Fire",
    "pH",
    "Clay",
    "SOC",
    "Fine roots",
    "Soil N",
    "SOC:N",
    "Fine root C:N",
    "Foliar C:N",
    "Litter C:N",
    "Foliar lignin",
    "Litter lignin",
    "O horizon thickness",
    "Horizon thickness",
    #expression(Ca[s]),
    expression(Ca[s]+Mg[s]),
    expression(Al[o]),
    expression(Fe[o]),
    # expression(Al[o]+Fe[o]),
    expression(Fe[d-o]),
    "16S",
    "ITS",
    "16S:ITS"
  ))+
  theme(axis.text.x=element_text(angle = +45, hjust = 1),axis.title.x=element_blank(),
        axis.title.y=element_blank())+geom_text(data=Pmat_red_melt,
            mapping=aes(x=X2,y=X1,label=value),
            color=Pmat_red_melt$Corr,
            size=7)

