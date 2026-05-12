#run the "NEON_NMR_overview_1" script first to load data

library(olsrr)
library(texreg)

#make some new predictors
#weathering index of Chittleborough 2007 Australian Journal of Earth Sciences
#WR = [(CaO + MgO + Na2O)/ZrO2] 
#Mjelm = mg / kg
#convert from metal to oxide
data1$WR_Zr<-(data1$caMjelm_mg_kg*((40.1+16)/40.1)+data1$mgMjelm_mg_kg*((24.3+16)/24.3)+data1$naMjelm_mg_kg*(23*2+16)/(23*2))/(data1$zrMjelm_mg_kg*(91.2+16*2)/91.2)

#compare with another weathering ratio
data1$Al_Si_total_ratio<-data1$alMjelm_mg_kg/data1$siMjelm_mg_kg
data1$Al_Si_total_ratio_log10<-log10(data1$alMjelm_mg_kg/data1$siMjelm_mg_kg)

#calculate horizon depth
data1$horizon_depth<-data1$horizonBottomDepth-data1$horizonTopDepth

#log transform some predictors
data1$WR_Zr_log10<-log10(data1$WR_Zr)
data1$siMjelm_log10<-log10(data1$siMjelm_mg_kg)
data1$Oxalate_Al_mg_g_ICP_log10<-log10(data1$Oxalate_Al_mg_g_ICP)
data1$Oxalate_Fe_mg_g_log10<-log10(data1$OxalateFe_g_kg)
data1$Ca_Na2SO4_mg_g_log10<-log10(data1$Ca_Na2SO4_mg_g)
data1$Fe_crystalline<-data1$DithioniteFe_g_kg-data1$OxalateFe_g_kg
data1$Fe_crystalline_log10<-log10(data1$Fe_crystalline)
data1$Dithionite_Fe_mg_g_log10<-log10(data1$DithioniteFe_g_kg)
data1$Fe_ox_di_ratio<-data1$Oxalate_Fe_mg_g_ICP/data1$DithioniteFe_g_kg
data1$alMjelm_mg_kg_log10<-log10(data1$alMjelm_mg_kg)
data1$feMjelm_mg_kg_log10<-log10(data1$feMjelm_mg_kg)
data1$caMjelm_mg_kg_log10<-log10(data1$caMjelm_mg_kg)
data1$live_total_roots_mg_cm3_log10<-log10(data1$live_total_roots_mg_cm3)
data1$total_roots_mg_cm3<-data1$dead_total_roots_mg_cm3+data1$live_total_roots_mg_cm3
data1$total_roots_mg_cm3_log10<-log10(data1$total_roots_mg_cm3)
data1$estimatedOC_g_kg_log10<-log10(data1$estimatedOC_g_kg)
data1$wgt_live_root_lessthan_4mm_CN_ratio_0_30cm_log10<-log10(data1$wgt_live_root_lessthan_4mm_CN_ratio_0_30cm)
data1$ti_si_ratio<-data1$tiMjelm_mg_kg/data1$siMjelm_mg_kg
data1$ti_si_ratio_log10<-log10(data1$tiMjelm_mg_kg/data1$siMjelm_mg_kg)
data1$zr_si_ratio_log10<-log10(data1$zrMjelm_mg_kg/data1$siMjelm_mg_kg)
data1$Al_oxalate_Si_total_ratio<-data1$Oxalate_Al_mg_g_ICP/data1$siMjelm
data1$Al_oxalate_Si_total_ratio_log10<-log10(data1$Al_oxalate_Si_total_ratio)
data1$Al_log10_Si_log10<-log10(data1$Oxalate_Al_mg_g_ICP)/log10(data1$siMjelm)

data1$fungi_meanCopyNumber_ITS_log10<-log10(data1$fungi_meanCopyNumber_ITS)
data1$bact_meanCopyNumber_16S_log10<-log10(data1$bacteria_and_archaea_meanCopyNumber_16S)
data1$fung_bact_ratio_log10<-log10(data1$fung_bact_ratio)

data1$Ca_plus_Mg_umol_log10<-log10(data1$Ca_plus_Mg_umol_g)
data1$Oxalate_Al_Fe_umol_g_log10<-log10(data1$Oxalate_Al_Fe_umol_g) 


#export "data1" for subsequent SEM, etc.
# write.csv(data1,"data1.csv", row.names=FALSE)


#make df with the predictors of interest for the global model
data1_global<-as.data.frame(scale(as.matrix(data1[,c("WR_Zr_log10","Fe_crystalline_log10","Oxalate_Fe_mg_g_log10","Oxalate_Al_mg_g_ICP_log10","Ca_Na2SO4_mg_g_log10","Ca_plus_Mg_umol_log10","phH2o","live_root_lessthan_4mm_approx_mg_cm3","wgt_live_root_lessthan_4mm_CN_ratio_0_30cm","MAP_minus_PET_mm","MAT_C","siMjelm_mg_kg","fungi_meanCopyNumber_ITS_log10","bact_meanCopyNumber_16S_log10","fung_bact_ratio_log10","litter_CtoN","litter_ligninPercent","foliage_CtoN","foliage_ligninPercent","horizonTopDepth","horizon_depth","forest_vs_non","Prescribed_fire")])))
data1_global$siteID<-data1$siteID
#bind the RC vectors onto this DF
data1_global$RC1<-data1$RC1
data1_global$RC2<-data1$RC2
data1_global$RC3<-data1$RC3

#only include samples from NEON
data1_global_NEON<-data1_global[!is.na(data1_global$wgt_live_root_lessthan_4mm_CN_ratio_0_30cm),]

#only include samples with qpcr
data1_global_NEON_qpcr<-data1_global_NEON[!is.na(data1_global_NEON$bact_meanCopyNumber_16S_log10),]

#only include forest

forest<-c("ABBY","BONA","CLBJ","GRSM","GUAN","HARV","JERC",
        "LENO","MLBS","ORNL","OSBS","RMNP","SCBI","SOAP","TALL",
        "TEAK","UKFS","UNDE","WREF","calz","carr","bayo","vanc",
        "reac", "samt", "marc", "icac","elve","redb")
data1_forest<-data1_global[data1_global$siteID %in% forest,]

#only include grass/shrublands
grass_shrub<-setdiff(data1_global$siteID,forest)
data1_grass_shrub<-data1_global[data1_global$siteID %in% grass_shrub,]

#we will analyze different datasets
length(data1_global_NEON$Oxalate_Fe_mg_g_log10)
length(data1_forest$Oxalate_Fe_mg_g_log10)
length(data1_global_NEON_qpcr$Oxalate_Fe_mg_g_log10)
length(data1_global$Oxalate_Fe_mg_g_log10)

### some linear models
#global model
lmX<-lm(RC1~Fe_crystalline_log10+Oxalate_Fe_mg_g_log10+Oxalate_Al_mg_g_ICP_log10+Ca_plus_Mg_umol_log10+phH2o+MAT_C+MAP_minus_PET_mm+horizon_depth+horizonTopDepth+forest_vs_non+Prescribed_fire, data1_global, na.action='na.fail')

ols_coll_diag(lmX) # remove Feo

lmX_noFeo<-lm(RC1~Fe_crystalline_log10+Oxalate_Al_mg_g_ICP_log10+Ca_plus_Mg_umol_log10+phH2o+MAT_C+MAP_minus_PET_mm+horizon_depth+horizonTopDepth+forest_vs_non+Prescribed_fire, data1_global, na.action='na.fail')

ols_coll_diag(lmX_noFeo) # remove pH

lmX_noFeo_nopH<-lm(RC1~Fe_crystalline_log10+Oxalate_Al_mg_g_ICP_log10+Ca_plus_Mg_umol_log10+MAT_C+MAP_minus_PET_mm+horizon_depth+horizonTopDepth+forest_vs_non+Prescribed_fire, data1_global, na.action='na.fail')

ols_coll_diag(lmX_noFeo_nopH) # OK

RC1_global <- ols_step_backward_p(lmX_noFeo_nopH, prem = 0.05)
RC1_global$model

RC1_global<-lm(RC1~Oxalate_Al_mg_g_ICP_log10+Ca_plus_Mg_umol_log10+forest_vs_non+Prescribed_fire,data=data1_global)
summary(RC1_global)

#try more conservative model selection approach
RC1_global_cons <- ols_step_backward_p(lmX_noFeo_nopH, prem = 0.01)
RC1_global_cons$model #just Alo

RC1_global_cons<-lm(RC1~Oxalate_Al_mg_g_ICP_log10,data=data1_global)
summary(RC1_global_cons)


#now, do backward selection, adding variables only present in the NEON samples
lmX_NEON<-lm(RC1~Oxalate_Al_mg_g_ICP_log10+Ca_plus_Mg_umol_log10+forest_vs_non+Prescribed_fire+fungi_meanCopyNumber_ITS_log10+bact_meanCopyNumber_16S_log10+wgt_live_root_lessthan_4mm_CN_ratio_0_30cm+live_root_lessthan_4mm_approx_mg_cm3,data=data1_global)

RC1_NEON <- ols_step_backward_p(lmX_NEON, prem = 0.05)

RC1_NEON<-lm(RC1~Oxalate_Al_mg_g_ICP_log10+Ca_plus_Mg_umol_log10+Prescribed_fire+fungi_meanCopyNumber_ITS_log10, data1_global)
#ITS now significant
summary(RC1_NEON)

#more conservative model

# RC1_NEON_cons <- ols_step_backward_p(lmX_NEON, prem = 0.01)

