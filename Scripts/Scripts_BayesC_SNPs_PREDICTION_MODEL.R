

#        R script for running prediction across populations (1st part) and within population (2st part)

#        BayesC method is used here provided by BGLR package



###                          PREDICTION USING SNPs MARKERS                              ###
########################### 1st part: PREDICTION ACORSS POPULATION ########################
library(data.table)
library(genio)
library(BEDMatrix)
library(BGLR)
library(gridExtra) 


nIter=100000
pi = 0.01
p0 = 5
all.phenotypes.COR<- NULL


#--- input phenotypes
all.phenotypes <- read.csv("Accessions_Traits.csv", header=TRUE)



#---- LIST WITH ALL PHENOTYPES FOR LOOP
data<- cbind(all.phenotypes[,5:15])
names<-names(data)
list.phen<-list(all.phenotypes[,5:15])

somePDFPath = "BAYESCADMARO_snps.pdf"
pdf(file=somePDFPath, onefile=TRUE)  
par(mfrow=c(2,2))


fm= list()
for (i in 1:11) {
  
  #---- SCALE PHENOTYPES
  y<-(unlist(list.phen[[1]][i]))
  y<-scale(y,center=TRUE,scale=TRUE)
  
  yNA=y
  tst=which(all.phenotypes$Group %in% c("ADM", "ARO") & y != "NA")
  tstADM=which(all.phenotypes$Group %in% c("ADM") & y != "NA")
  tstARO=which(all.phenotypes$Group %in% c("ARO") & y != "NA")
  tstNA=which(is.na(y)) #NO USED
  yNA[tst]=NA
  
  
  fm1 = BGLR(y=yNA,ETA=list(ETA1=list(X=snps, model='BayesC', probIn=pi, counts=p0)), nIter=nIter,saveAt=sprintf("./fmgBAYESCaroadmsnps_y%.0f_",i),verbose=F)
  
  corel1=cor(fm1$yHat[tst], y[tst])
  
  
  yHat<-fm1$yHat 
  y1=as.numeric(y)
  tmp<-range(c(yHat,y),na.rm = T) 
  plot(yHat[tstADM]~ y1[tstADM], col="red", xlim=tmp,ylim=c(-2,2),xlab="",ylab="")
  par(new=TRUE)
  plot(yHat[tstARO]~y1[tstARO], col="BLUE", xlim=tmp,ylim=c(-2,2),xlab="Observed",ylab="Predicted",main=c(names[i],sprintf("ARO,ADM snps BayesC, cor%.02f", corel1)))
  par(new=TRUE)
  tstAROADM<-c(tstARO,tstADM)
  abline(lm(yHat[tstAROADM]~y1[tstAROADM]))
  
  all.phenotypes.COR<- rbind( all.phenotypes.COR,corel1)
  
  fm [i] =list(fm1)
}

dev.off() 



row.names(all.phenotypes.COR)<-names
colnames(all.phenotypes.COR)<- c("correlation")
png(file=sprintf("./corr_BAYESC_snps_AROADM.png"))
grid.table(all.phenotypes.COR)

dev.off()



save.image("BAYESCsnpsAROADM.RData")
rm(list = ls())




########################### 2st part: PREDICTION WITHIN POPULATION ########################
library(AGHmatrix)
library(data.table)
library(genio)
library(BEDMatrix)
library(BGLR)
library(gridExtra)   


nIter=100000
pi = 0.01
p0 = 5
all.phenotypes.COR<- NULL



## input new varieties with phenotypes
all.phenotypes <- read.csv("Accessions_Traits.csv", header=TRUE)



#---- FOR PHENOTYPES
data<- cbind(all.phenotypes[,5:15])
names<-names(data)
list.phen<-list(all.phenotypes[,5:15])


somePDFPath = "bayesCIND_snps.pdf"
pdf(file=somePDFPath, onefile=TRUE)  
par(mfrow=c(2,2))

fm= list()
for (i in 1:11) {
  
  #---- SCALE PHENOTYPES
  y<-(unlist(list.phen[[1]][i]))
  y<-scale(y,center=TRUE,scale=TRUE)
  yNA=y
  tst=which(all.phenotypes$Group == "IND" & all.phenotypes$Status == "I" &  y != "NA") 
  tstNA=which(is.na(y))
  yNA[tst]=NA
  
  fm2 = BGLR(y=yNA,ETA=list(ETA01=list(X=V1,model="FIXED"),ETA1=list(X=snps, model='BayesC', probIn=pi, counts=p0)), nIter=nIter,saveAt=sprintf("./fmgBAYESCINDsnps_y%.0f_",i),verbose=F)
  
  corel2=cor(fm2$yHat[tst], y[tst])
  
  yHat<-fm2$yHat 
  y1=as.numeric(y)
  tmp<-range(c(yHat,y),na.rm = T) 
  plot(yHat[tst]~ y1[tst], col="red", xlim=tmp,ylim=c(-2,2),xlab="Observed",ylab="Predicted",main=c(names[i],sprintf("IND snps BayesC, cor%.02f", corel2)))
  par(new=TRUE)
  abline(lm(yHat[tst]~y1[tst]))
  
  all.phenotypes.COR<- rbind( all.phenotypes.COR,corel2)
  
  fm [i] =list(fm2)
  
}

dev.off() 

row.names(all.phenotypes.COR)<-names
colnames(all.phenotypes.COR)<- c("correlation")
png(file=sprintf("./corr_BAYESc_snps_IND.png"))
grid.table(all.phenotypes.COR)

dev.off()


save.image("BayesCsnpsIND.RData")
rm(list = ls())
