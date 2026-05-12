#make Figure 2, plot of ordinations

#run NEON_NMR_overview_1.R first
source(here("data/raw/NEON_NMR_overview_1.R"))
library(psych)
library(here)
#rotated principal components with correlation matrix
#rename columns to remove "fraction"
colnames(mat1)<-c("Lignin","Carbohyd","Lipid","Protein","Char","Carbonyl")
mat2.3<-principal(mat1, nfactors=3, rotate='varimax', covar=FALSE)

# mat2.3
# mat2.3$loadings[1:18]
# weights<-data.frame(mat2.3$weights)

scores<-as.data.frame(mat2.3$scores)


#quartz(,12,4)
windows(width=12, height=4)
par(mfrow=c(1,3),mar=c(3,3,1.5,1.5), tck = -0.02, mgp=c(1.5,0.25,0),cex=1.2)
biplot(mat2.3, choose=c(1,2),
       xlim.s=c(-3,3),ylim.s=c(-3,3),
       xlim.f=c(-1.15,1.15),ylim.f=c(-1.15,1.15),
       pch=20,
       col=c(rgb(0,0,0,max=255,alpha=75),"darkgreen"))
mtext("a",font=2, at=c(-1.5),cex=1.3,line=0.5)
biplot(mat2.3, choose=c(2,3),
       xlim.s=c(-3,3),ylim.s=c(-3,3),
       xlim.f=c(-1.15,1.15),ylim.f=c(-1.15,1.15),
       pch=20,
       col=c(rgb(0,0,0,max=255,alpha=75),"darkgreen"))
mtext("b",font=2, at=c(-1.5),cex=1.3,line=0.5)
biplot(mat2.3, choose=c(1,3),
       xlim.s=c(-3,3),ylim.s=c(-3,3),
       xlim.f=c(-1.15,1.15),ylim.f=c(-1.15,1.15),
       pch=20,
       col=c(rgb(0,0,0,max=255,alpha=75),"darkgreen"))
mtext("c",font=2, at=c(-1.5),cex=1.3,line=0.5)




# 
# quartz(,12,4)
# par(mfrow=c(1,3),mar=c(3,3,1.5,1.5), tck = -0.02, mgp=c(1.5,0.25,0),cex=1.2)
# plot(RC2~RC1,
#      scores,
#      type='n',
#      xlim=c(-7,7),
#      ylim=c(-7,7))
# text(scores$RC1,scores$RC2, row.names(scores))
# plot(RC3~RC2,
#      scores,
#      type='n',
#      xlim=c(-7,7),
#      ylim=c(-7,7))
# text(scores$RC2,scores$RC3, row.names(scores))
# plot(RC3~RC1,
#      scores,
#      type='n',
#      xlim=c(-7,7),
#      ylim=c(-7,7))
# text(scores$RC1,scores$RC3, row.names(scores))
# 



