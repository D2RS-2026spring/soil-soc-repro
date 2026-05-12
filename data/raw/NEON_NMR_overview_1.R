#updated 20200522

#this script loads the "master" spreadsheet with NEON soil NMR spectra and other soil variables. Then, it does a rotated PCA, where the scores are used for downstream analyses.


#remove old data
rm(list=ls())

#some packages to load
library(psych)
library(MASS)
library(Hmisc)
#library(pca3d)
library(factoextra)

# load NMR data
data<-read.csv("Hall_Ye_etal_NEON_NMR_biogeo.csv",header=TRUE)
head(data)
str(data)

# make a data frame with the C molecules expressed as fractional abundance
mat1<-data.frame(data[,c("Lignin_fraction","Carbohydrate_fraction","Lipid_fraction","Protein_fraction","Char_fraction", "Carbonyl_fraction")])
#label the sample names in the data frame
row.names(mat1)<-data$siteID

#look at correlations among different functional groups
 cor(mat1)
 pairs(mat1)
 rcorr(as.matrix(mat1))

#do the same analyses with the C functional groups
 mat1_1<-data.frame(data[,c("Alkyl_fraction","N_Alkyl_Methoxyl_fraction","O_Alkyl_fraction","Di_O_Alkyl_fraction","Aromatic_fraction", "Phenolic_fraction","Amide_Carboxyl_fraction")])
  
 #label the sample names in the data frame
 row.names(mat1_1)<-data$Site
 rcorr(as.matrix(mat1_1))

# do a regular PCA based on the covariance matrix of the modeled C component data
mat2.1<-prcomp(mat1, scale.=FALSE)
# three principal components explain 97% of the data!
summary(mat2.1) # but  components don't map nicely onto the axes... reason to do rotation, shown later

# do a regular PCA based on the correlation matrix
mat2.2<-prcomp(mat1, scale.=TRUE)
# three principal components explain 90% of the data!
summary(mat2.2)
# the main difference here between the unscaled and scaled matrices is the carbonyl C

# #confirm that three axes are best here 
 # fviz_eig(mat2.1)

# #look at relationships of C components with sites; we see that a rotation might help better match vectors with axes
 # fviz_pca_biplot(mat2.1, show.labels=TRUE)
# #another version of the 2-d plot:
#   pca2d(mat2.1, biplot=TRUE, show.labels=FALSE)
# # #can look in 3d space now; confusing but interesting
#  pca3d(mat2.1, biplot=TRUE, show.labels=FALSE)
#  pca3d(mat2.1, biplot=TRUE, show.labels=TRUE)
# # 

#next, try a rotated PCA using varimax rotation, function from the psych package. Use the correlation matrix
 mat2.3<-principal(mat1, nfactors=3, rotate='varimax', covar=FALSE)

  #three vactors explain 90% of the variance, just as above in the regular PCA (which we expect). The main difference here are the correlations with the axes.
   
 #try rotated PCA with the covariance matrix
 mat2.4<-principal(mat1, nfactors=3, rotate='varimax', covar=TRUE)
 mat2.4
 # standardized loadings are very similar to the version above with correlation matrix. Therefore, stick with 2.3
 
 #extract component scores for correlating with other soil variables 
 data1<-cbind(data,mat2.3$scores) #RC1,RC2,RC3 are the scores for the rotated components 1,2,3
 
 
 
