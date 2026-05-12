
#first run the NEON_NMR_overview_1.R and NEON_NMR_stats_RC1.R

##global model

lmX<-lm(RC3~Fe_crystalline_log10+Oxalate_Fe_mg_g_log10+Oxalate_Al_mg_g_ICP_log10+Ca_plus_Mg_umol_log10+phH2o+MAT_C+MAP_minus_PET_mm+horizonTopDepth+horizon_depth+forest_vs_non+Prescribed_fire, data1_global, na.action='na.fail')

ols_coll_diag(lmX) #remove Feo

lmX_noFe<-lm(RC3~Fe_crystalline_log10+Oxalate_Al_mg_g_ICP_log10+Ca_plus_Mg_umol_log10+phH2o+MAT_C+MAP_minus_PET_mm+horizonTopDepth+horizon_depth+forest_vs_non+Prescribed_fire, data1_global, na.action='na.fail')

ols_coll_diag(lmX_noFe) #remove MAP

lmX_noFe_noMAP<-lm(RC3~Fe_crystalline_log10+Oxalate_Al_mg_g_ICP_log10+Ca_plus_Mg_umol_log10+MAT_C+phH2o+horizonTopDepth+horizon_depth+forest_vs_non+Prescribed_fire, data1_global, na.action='na.fail')

ols_coll_diag(lmX_noFe_noMAP) #remove pH

lmX_noFe_noMAP_nopH<-lm(RC3~Fe_crystalline_log10+Oxalate_Al_mg_g_ICP_log10+Ca_plus_Mg_umol_log10+MAT_C+horizonTopDepth+horizon_depth+forest_vs_non+Prescribed_fire, data1_global, na.action='na.fail')
ols_coll_diag(lmX_noFe_noMAP_nopH) #OK


RC3_global <- ols_step_backward_p(lmX_noFe_noMAP_nopH, prem = 0.05)
RC3_global$model

RC3_global<-lm(RC3~Ca_plus_Mg_umol_log10+horizonTopDepth,data1_global)
summary(RC3_global)# note that Ca + Mg and pH are strongly related; fit model with pH instead

RC3_global_alt<-lm(RC3~phH2o+horizonTopDepth,data1_global)
summary(RC3_global_alt)

##model selection on global model

lmX<-lm(RC3~Fe_crystalline_log10+Oxalate_Fe_mg_g_log10+Oxalate_Al_mg_g_ICP_log10+Ca_plus_Mg_umol_log10+phH2o+MAT_C+MAP_minus_PET_mm+horizonTopDepth+horizon_depth+forest_vs_non+Prescribed_fire, data1_global, na.action='na.fail')

RC3_global_alt1 <- ols_step_backward_p(lmX, prem = 0.05)
RC3_global_alt1$model

RC3_global_alt1<-lm(RC3~phH2o+horizonTopDepth+MAP_minus_PET_mm,data1_global)
summary(RC3_global_alt1)
ols_coll_diag(RC3_global_alt1) #OK

##conservative model

RC3_global_alt1_cons <- ols_step_backward_p(lmX, prem = 0.01)
RC3_global_alt1_cons$model

RC3_global_alt1_cons<-lm(RC3~phH2o+horizonTopDepth,data1_global)
summary(RC3_global_alt1_cons)

###NEON model
lmX_NEON<-lm(RC3~phH2o+horizonTopDepth+MAP_minus_PET_mm+fungi_meanCopyNumber_ITS_log10+bact_meanCopyNumber_16S_log10+wgt_live_root_lessthan_4mm_CN_ratio_0_30cm+live_root_lessthan_4mm_approx_mg_cm3,data1_global)

RC3_NEON<-ols_step_backward_p(lmX_NEON, prem = 0.05)

RC3_NEON<-lm(RC3~phH2o+horizonTopDepth+wgt_live_root_lessthan_4mm_CN_ratio_0_30cm,data1_global)

summary(RC3_NEON)



