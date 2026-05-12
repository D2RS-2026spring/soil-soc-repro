#Extended Data Fig 2
source(here("data/raw/NEON_NMR_overview_1.R"))
#run NEON_NMR_overview_1.R first

library(ggplot2)
library(gridExtra)
library(cowplot)
library(reshape)
library(here)

#get data frame with subset for graphing
data1$Fe_crystalline<-data1$DithioniteFe_g_kg-data1$OxalateFe_g_kg
data1$Fe_crystalline_log10<-log10(data1$Fe_crystalline)

datared<-data1[,c("siteID",
                  "Oxalate_Fe_mg_g_ICP",
                  "Oxalate_Al_mg_g_ICP",
                 "Fe_crystalline",
                 "phH2o",
                # "Ca_Na2SO4_mg_g",
                "Ca_plus_Mg_umol_g",
                 "live_root_lessthan_4mm_approx_mg_cm3",
                 "wgt_live_root_lessthan_4mm_CN_ratio_0_30cm",
                 "MAP_minus_PET_mm",
                 "MAT_C",
                 "fungi_meanCopyNumber_ITS",
                 "bacteria_and_archaea_meanCopyNumber_16S",
                 "fung_bact_ratio",
                 "clayTotal_percent",
                 "estimatedOC_g_kg",
                 "OC_N_ratio",
                 "nitrogenTot_g_kg",
                 "foliage_ligninPercent",
                 "foliage_CtoN",
                 "litter_ligninPercent",
                 "litter_CtoN"
)]

#rename for plotting names
names(datared)<-c("Site","Alo (mg g-1)","Feo (mg g-1)","Fed-o (mg g-1)","pH","Ca (mg g-1)","Fine roots (mg cm-3)","Fine root C:N", "MAP-PET (mm)","MAT (C)","ITS","16S","ITS:16S","Clay (%)","SOC (mg g-1)","SOC:N","Soil N (mg g-1)","Foliar lignin (%)","Foliar C:N","Litter lignin (%)","Litter C:N")



datared_melt<-melt(datared, id.vars="Site")
head(datared_melt)

datared_melt$ordering<-as.character(datared_melt$variable)
unique(datared_melt$ordering)

datared_melt$ordering[datared_melt$ordering=="MAT (C)"]<-1
datared_melt$ordering[datared_melt$ordering=="MAP-PET (mm)"]<-2
datared_melt$ordering[datared_melt$ordering=="pH"]<-3
datared_melt$ordering[datared_melt$ordering=="SOC (mg g-1)"]<-4
datared_melt$ordering[datared_melt$ordering=="Soil N (mg g-1)"]<-5
datared_melt$ordering[datared_melt$ordering=="SOC:N"]<-6
datared_melt$ordering[datared_melt$ordering=="Fine roots (mg cm-3)"]<-7
datared_melt$ordering[datared_melt$ordering=="Fine root C:N"]<-8
datared_melt$ordering[datared_melt$ordering=="Foliar C:N"]<-9
datared_melt$ordering[datared_melt$ordering=="Litter C:N"]<-10
datared_melt$ordering[datared_melt$ordering=="Foliar lignin (%)"]<-11
datared_melt$ordering[datared_melt$ordering=="Litter lignin (%)"]<-12

datared_melt$ordering[datared_melt$ordering=="Ca (mg g-1)"]<-13
datared_melt$ordering[datared_melt$ordering=="Alo (mg g-1)"]<-14
datared_melt$ordering[datared_melt$ordering=="Feo (mg g-1)"]<-15
datared_melt$ordering[datared_melt$ordering=="Fed-o (mg g-1)"]<-16
datared_melt$ordering[datared_melt$ordering=="Clay (%)"]<-17
datared_melt$ordering[datared_melt$ordering=="16S"]<-18
datared_melt$ordering[datared_melt$ordering=="ITS"]<-19
datared_melt$ordering[datared_melt$ordering=="ITS:16S"]<-20

datared_melt$ordering<-as.numeric(datared_melt$ordering)

datared_melt<-datared_melt[order(datared_melt$ordering),]

datared_melt$variable<-factor(datared_melt$variable, levels=c("MAT (C)",
                                                              "MAP-PET (mm)",
                                                              "pH",
                                                              "SOC (mg g-1)",
                                                              "Soil N (mg g-1)",
                                                              "SOC:N",
                                                              "Fine roots (mg cm-3)",
                                                              "Fine root C:N",
                                                              "Foliar C:N",
                                                              "Litter C:N",
                                                              "Foliar lignin (%)",
                                                              "Litter lignin (%)",
                                                              "Ca (mg g-1)",
                                                              "Alo (mg g-1)",
                                                              "Feo (mg g-1)",
                                                              "Fed-o (mg g-1)",
                                                              "Clay (%)",                                                         
                                                              "16S",
                                                              "ITS",
                                                              "ITS:16S"))

#plot Supplementary Fig 2
#quartz(,8,8)
windows(width=8, height=8)
ggplot(datared_melt, aes(1,value))+
  geom_boxplot(outlier.shape = NA)+ # omit outliers so they don't show up twice with jitter
  geom_jitter(position=position_jitter(0.1), color="grey50")+
  facet_wrap(~variable, scales="free_y")+
  theme_classic()+
  theme(axis.text.x=element_text(angle = +45, hjust = 1))+
  scale_x_discrete("")+
  scale_y_continuous("")


#make special annotations for panel labels
panel_names<-c(
  "MAT (C)" = "MAT*~(degree*C)" ,             
  "MAP-PET (mm)" = "MAP-PET*~(mm)",        
  "pH"="pH"           ,        
  "SOC (mg g-1)" = "SOC*~(mg*~g^-1)"  ,      
  "Soil N (mg g-1)" = "SON*~(mg*~g^-1)",      
  "SOC:N" = "SOC:N" ,    
  "Fine roots (mg cm-3)" = "Fine*~roots*~(mg*~cm^-3)", 
  "Fine root C:N"="Fine*~root*~C:N",
  "Foliar C:N"="Foliar*~C:N" , 
  "Litter C:N" ="Litter*~C:N",
  "Foliar lignin (%)" ="Foliar*~lignin*~('%')"   ,
  "Litter lignin (%)"="Litter*~lignin*~('%')"    ,
  "Ca (mg g-1)" ="Ca*~+Mg*~(mu*mol*~g^-1)", 
  "Alo (mg g-1)" ="Al[o]*~(mg*~g^-1)" ,
  "Feo (mg g-1)"  ="Fe[o]*~(mg*~g^-1)" ,
  "Fed-o (mg g-1)"="Fe[d-o]*~(mg*~g^-1)" ,
  "Clay (%)"="Clay*~('%')",
  "16S"="16*S",
  "ITS"="ITS",
  "ITS:16S"="ITS*':'*16*S"        )


print(ggplot(datared_melt, aes(1,value))+
  geom_boxplot(outlier.shape = NA)+ # omit outliers so they don't show up twice with jitter
  geom_jitter(position=position_jitter(0.1), color="grey50")+
  facet_wrap(~variable, scales="free_y",
             labeller = labeller(variable=as_labeller(panel_names,label_parsed)))+
  theme_classic()+
  theme(axis.text.x=element_text(angle = +45, hjust = 1))+
  scale_x_discrete("")+
  scale_y_continuous(""))



