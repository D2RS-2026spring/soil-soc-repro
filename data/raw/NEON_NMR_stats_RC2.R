#run the NEON NMR overview and PCA script first

#run RC1 analysis script first

##global dataset
lmX<-lm(RC2~Fe_crystalline_log10+Oxalate_Fe_mg_g_log10+Oxalate_Al_mg_g_ICP_log10+Ca_plus_Mg_umol_log10+phH2o+MAT_C+MAP_minus_PET_mm+horizonTopDepth+horizon_depth+forest_vs_non+Prescribed_fire, data1_global, na.action='na.fail')

#evaluate collinearity
ols_coll_diag(lmX) #eliminate Feo

lmX_noFeo<-lm(RC2~Fe_crystalline_log10+Oxalate_Al_mg_g_ICP_log10+Ca_plus_Mg_umol_log10+phH2o+MAT_C+MAP_minus_PET_mm+horizonTopDepth+horizon_depth+forest_vs_non+Prescribed_fire, data1_global, na.action='na.fail')

#evaluate collinearity
ols_coll_diag(lmX_noFeo) #eliminate pH

lmX_noFeo_nopH<-lm(RC2~Fe_crystalline_log10+Oxalate_Al_mg_g_ICP_log10+Ca_plus_Mg_umol_log10+MAP_minus_PET_mm+MAT_C+horizonTopDepth+horizon_depth+forest_vs_non+Prescribed_fire, data1_global, na.action='na.fail')

#evaluate collinearity
ols_coll_diag(lmX_noFeo_nopH) #OK

#model selection
RC2_global<-ols_step_backward_p(lmX_noFeo_nopH, prem = 0.05)
RC2_global$model

RC2_global<-lm(RC2~Fe_crystalline_log10+horizon_depth,data=data1_global)
summary(RC2_global)

ols_coll_diag(RC2_global) #OK
ols_plot_diagnostics(RC2_global) #OK

#model selection
RC2_global_cons<-ols_step_backward_p(lmX_noFeo_nopH, prem = 0.01)
RC2_global_cons$model

RC2_global_cons<-lm(RC2~Fe_crystalline_log10,data=data1_global)
summary(RC2_global_cons)

#add predictors from NEON 
#now, do backward selection, adding variables only present in the NEON samples
lmX_NEON<-lm(RC2~Fe_crystalline_log10+horizon_depth+fungi_meanCopyNumber_ITS_log10+bact_meanCopyNumber_16S_log10+wgt_live_root_lessthan_4mm_CN_ratio_0_30cm+live_root_lessthan_4mm_approx_mg_cm3,data=data1_global)

RC2_NEON<-ols_step_backward_p(lmX_NEON, prem = 0.05)

RC2_NEON$model

RC2_NEON<-lm(RC2~Fe_crystalline_log10+horizon_depth+wgt_live_root_lessthan_4mm_CN_ratio_0_30cm+live_root_lessthan_4mm_approx_mg_cm3,data1_global)

summary(RC2_NEON)

# RC2_NEON_test<-lm(RC2~Fe_crystalline_log10+Ca_plus_Mg_umol_log10,data=data1_global_NEON)
# summary(RC2_NEON_test)
