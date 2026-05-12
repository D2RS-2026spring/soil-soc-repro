rm(list = ls())

library(piecewiseSEM) # Implements piecewise structural equation modeling from a single list of structural equations

data<-read.csv("Hall_Ye_2020_NEON_NMR_20200716.csv",header=TRUE)
str(data)
head(data)

#calculate depth of the sampled mineral horizon
data$horizon_depth<-data$horizonBottomDepth-data$horizonTopDepth

#note that "horizonTopDepth" is also actually the depth of the organic horizon, since we always sampled the mineral horizon directly below the organic horizon.


# evaluate controls on different C principal components
# select potential models a-priori based 1) on the multiple regression results (conservative model), and 2) based on the critiques from reviewer 1 (i.e., do climate, O horizon, and/or sampling depth influence our results)

#here, horizonTopDepth indicates the depth of the organic (O) horizon. horizon_depth indicates the depth increment of the mineral horizon that was sampled.

#log-transform Oxalate Al
###RC1 
data$Oxalate_Al_mg_g_ICP_log10<-log10(data$Oxalate_Al_mg_g_ICP)

#make sure that there are no NAs in our predictors
sum(is.na(data$Oxalate_Al_mg_g_ICP_log10))
sum(is.na(data$horizonTopDepth))
sum(is.na(data$MAP_minus_PET_mm))
sum(is.na(data$MAT_C))
sum(is.na(data$horizon_depth))
sum(is.na(data$phH2o))


#make a list of the individual models that underlie our hypothesized SEM.

SEM_RC1_models<-list(
  lm(RC1~Oxalate_Al_mg_g_ICP_log10+horizonTopDepth+horizon_depth, data),
  lm(Oxalate_Al_mg_g_ICP_log10~MAT_C+MAP_minus_PET_mm+horizon_depth+horizonTopDepth,data),
  lm(horizon_depth~MAT_C+MAP_minus_PET_mm+horizonTopDepth,data),
  lm(horizonTopDepth~MAT_C+MAP_minus_PET_mm,data))

SEM_RC1_1<-as.psem(SEM_RC1_models)

summary(SEM_RC1_1, .progressBar = F)


###
###RC2 
sum(is.na(data$Fe_crystalline_log10))
sum(is.na(data$Ca_Na2SO4_mg_g))

SEM_RC2_models<-list(
  lm(RC2~Fe_crystalline_log10 + horizonTopDepth + horizon_depth, data),
  lm(Fe_crystalline_log10~MAT_C+horizonTopDepth+MAP_minus_PET_mm+ horizon_depth,data),
  lm(horizon_depth~MAT_C+MAP_minus_PET_mm+horizonTopDepth,data),
  lm(horizonTopDepth~MAT_C+MAP_minus_PET_mm,data))

SEM_RC2_1<-as.psem(SEM_RC2_models)

summary(SEM_RC2_1, .progressBar = F)


###RC3 
sum(is.na(data$phH2o))

SEM_RC3_models<-list(
  lm(RC3~phH2o + horizonTopDepth + horizon_depth+MAP_minus_PET_mm, data),
  lm(phH2o~MAT_C + MAP_minus_PET_mm + horizon_depth + horizonTopDepth,data),
  lm(horizon_depth~MAT_C+MAP_minus_PET_mm+horizonTopDepth,data),
  lm(horizonTopDepth~MAT_C+MAP_minus_PET_mm,data))

SEM_RC3_1<-as.psem(SEM_RC3_models)

summary(SEM_RC3_1, .progressBar = F)


###now, prune the models

###RC1 
#make a list of the individual models that underlie our hypothesized SEM.



  RC1_0<-psem(lm(RC1~Oxalate_Al_mg_g_ICP_log10+horizonTopDepth+horizon_depth, data),
  lm(Oxalate_Al_mg_g_ICP_log10~MAT_C+MAP_minus_PET_mm+horizon_depth+horizonTopDepth,data),
  lm(horizon_depth~MAT_C+MAP_minus_PET_mm+horizonTopDepth,data),
  lm(horizonTopDepth~MAT_C+MAP_minus_PET_mm,data))
  
  summary(RC1_0, .progressBar = F)
  
  RC1_1<-update(RC1_0, horizon_depth~MAT_C+horizonTopDepth)

  summary(RC1_1, .progressBar = F)
  AIC(RC1_0,RC1_1)
  
  RC1_2<-update(RC1_1, RC1~Oxalate_Al_mg_g_ICP_log10+horizon_depth)
  
  summary(RC1_2, .progressBar = F)
  AIC(RC1_1,RC1_2)
  
  RC1_3<-update(RC1_2, horizon_depth~horizonTopDepth)
  
  summary(RC1_3, .progressBar = F)
  AIC(RC1_2,RC1_3) #removing MAT from horizon_depth increased AIC;

  RC1_4<-update(RC1_2, Oxalate_Al_mg_g_ICP_log10~MAT_C+MAP_minus_PET_mm+horizon_depth)
  
  summary(RC1_4, .progressBar = F)
  AIC(RC1_2,RC1_4) 
  
  RC1_5<-update(RC1_2, RC1~Oxalate_Al_mg_g_ICP_log10)

  #keep RC1_2 as the optimal model
  
  summary(RC1_2, .progressBar = F)
  
  
  
###
###RC2 


RC2_0<-psem(
  lm(RC2~Fe_crystalline_log10 + horizonTopDepth + horizon_depth, data),
  lm(Fe_crystalline_log10~MAT_C+horizonTopDepth+MAP_minus_PET_mm+ horizon_depth,data),
  lm(horizon_depth~MAT_C+MAP_minus_PET_mm+horizonTopDepth,data),
  lm(horizonTopDepth~MAT_C+MAP_minus_PET_mm,data))

summary(RC2_0, .progressBar = F)

RC2_1<-update(RC2_0,Fe_crystalline_log10~MAT_C+horizonTopDepth+MAP_minus_PET_mm)

summary(RC2_1, .progressBar = F)

AIC(RC2_0,RC2_1)

RC2_2<-update(RC2_1,horizon_depth~MAT_C+horizonTopDepth)

summary(RC2_2, .progressBar = F)

AIC(RC2_1,RC2_2)

RC2_3<-update(RC2_2,horizon_depth~horizonTopDepth)

summary(RC2_3, .progressBar = F)

AIC(RC2_2,RC2_3)

RC2_4<-update(RC2_2,RC2~Fe_crystalline_log10 + horizon_depth)

summary(RC2_4, .progressBar = F)

AIC(RC2_2,RC2_5)

RC2_5<-update(RC2_2,Fe_crystalline_log10~horizonTopDepth+MAP_minus_PET_mm+ horizon_depth)

summary(RC2_5, .progressBar = F)

AIC(RC2_2,RC2_5)

#keep RC2_2


###RC3 

RC3_0<-psem(
  lm(RC3~phH2o + horizonTopDepth + horizon_depth+MAP_minus_PET_mm, data),
  lm(phH2o~MAT_C + MAP_minus_PET_mm + horizon_depth + horizonTopDepth,data),
  lm(horizon_depth~MAT_C+MAP_minus_PET_mm+horizonTopDepth,data),
  lm(horizonTopDepth~MAT_C+MAP_minus_PET_mm,data))

summary(RC3_0, .progressBar = F)

RC3_1<-update(RC3_0,phH2o~MAT_C + MAP_minus_PET_mm + horizonTopDepth)

summary(RC3_1, .progressBar = F)

AIC(RC3_0,RC3_1)

RC3_2<-update(RC3_1,horizon_depth~MAT_C+horizonTopDepth)

AIC(RC3_1,RC3_2)

summary(RC3_2, .progressBar = F)

RC3_3<-update(RC3_2,phH2o~MAP_minus_PET_mm + horizonTopDepth)

summary(RC3_3, .progressBar = F)

AIC(RC3_2,RC3_3)

RC3_4<-update(RC3_3,RC3~phH2o + horizonTopDepth +MAP_minus_PET_mm)

summary(RC3_4, .progressBar = F)

AIC(RC3_3,RC3_4)

#RC3_3 is best


