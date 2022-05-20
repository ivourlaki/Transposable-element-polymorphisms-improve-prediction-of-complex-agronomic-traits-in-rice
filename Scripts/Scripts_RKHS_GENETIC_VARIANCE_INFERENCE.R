
######################## Genetic Variance Inference ##################################################
######################## SNPs are used here ###########################################################


library(AGHmatrix)
library(data.table)
library(genio)
library(BEDMatrix)
library(BGLR)
library(gridExtra)      

##--------------------- 1. SNPS ------------------------------------------------------
nIter=100000

load("Additive_Matrix.RData")
##-- Additive_Matrix.Rdata contans the below matrices:
##-- where G_VanRadenPine1, the additive matrix for SNPs
##-- where G_VanRadenPine2, the additive matrix for MITE/DTX
##-- where G_VanRadenPine3, the additive matrix for RLX/RIX


all.phenotypes.snps<- NULL
all.phenotypes.snps.origin  <- NULL
all.phenotypes.snps.pca  <- NULL
snps.origin.b.sd  <- NULL
snps.pca.b.sd  <- NULL




## input phenotypes
all.phenotypes <- read.csv("Accessions_Traits.csv", header=TRUE)




data<- cbind(all.phenotypes[,5:15])
names<-names(data)
list.phen<-list(all.phenotypes[,5:15])

fm= list()
for (i in 1:11) {
  
  #---- SCALE PHENOTYPES
  y<-(unlist(list.phen[[1]][i]))
  y<-scale(y,center=TRUE,scale=TRUE)
  
  
  
  
  fm1 = BGLR(y=y,ETA=list(ETA1=list(K=G_VanRadenPine1,model='RKHS')), nIter=nIter,saveAt=sprintf("./fmgSNPS_y%.0f_",i),verbose=T)
  
  
  varG1=scan(sprintf("./fmgSNPS_y%.0f_ETA_ETA1_varU.dat",i))
  varE=scan(sprintf("./fmgSNPS_y%.0f_varE.dat",i))
  h2G1 = mean(varG1/(varG1+varE))
  print (h2G1) ; (var(varE)) 
  all.values.h2 <-cbind(h2G1,var(varE))
  
  
  
  
  snps.origin.b.sd <-c(fm1$ETA$ETA0$b,fm1$ETA$ETA0$SD.b )
  
  
  
  
  snps.pca.b.sd <- c(fm1$ETA$ETA01$b,fm1$ETA$ETA01$SD.b )
  
  
  
  # gather all the phenotypes together 
  all.phenotypes.snps<-rbind(all.phenotypes.snps,all.values.h2)
  
  
  all.phenotypes.snps.origin<- rbind(all.phenotypes.snps.origin,snps.origin.b.sd)
  
  all.phenotypes.snps.pca<- rbind(all.phenotypes.snps.pca,snps.pca.b.sd)
  
  fm [i] =list(fm1)
  
}




row.names(all.phenotypes.snps.origin)<-names
colnames(all.phenotypes.snps.origin)<- c("ADM b","ARO b","AUS b","IND b","JAP b","ADM sd b","ARO sd b","AUS sd b","IND sd b","JAP sd b")
png(file=sprintf("./all.phenotypes.snps.origin.png"), height = 100*nrow( all.phenotypes.snps.origin), width = 150*ncol(all.phenotypes.snps.origin))
grid.table(all.phenotypes.snps.origin)

dev.off()

colnames(all.phenotypes.snps)<- c("h2","var(varE)")
row.names(all.phenotypes.snps)<-names
png(file=sprintf("./all.phenotypes.snps.png"), height = 50*nrow( all.phenotypes.snps), width = 200*ncol(all.phenotypes.snps))
grid.table(all.phenotypes.snps)
dev.off()


row.names(all.phenotypes.snps.pca)<-names
colnames(all.phenotypes.snps.pca)<- c("PCA1 b","PCA2 b","PCA3 b","PCA4 b","PCA1 sb","PCA2 sb","PCA3 sb","PCA4 sb")
png(file=sprintf("./all.phenotypes.snps.pca.png"), height = 100*nrow(all.phenotypes.snps.pca), width = 150*ncol(all.phenotypes.snps.pca))
grid.table(all.phenotypes.snps.pca)


dev.off()

