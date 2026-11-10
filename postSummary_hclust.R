rm(list=ls())
library(waveslim)
library(stringr)
library(cluster)
setwd("~/fPDPM/")

dist_cluVec<-function(cluVec1, cluVec2, wtVec=NULL){
  if(is.null(wtVec)){
    wtVec<-rep(1, length(cluVec1))
  }
  if(length(cluVec1) != length(cluVec2)| length(cluVec2)!=length(wtVec)){
    stop("Wrong Length")
  }
  res<-1-sum(as.numeric(cluVec1==cluVec2)*wtVec)/sum(wtVec)
  return(res)
}

args <- commandArgs(TRUE)
# args<- c("diagErr_highSNR_nLoc32_nSmp300_b3","3")
print(args)
dfName<-trimws(args[1])
K_max<-as.integer(args[2])
nLoc<-32
nSmp<-300

nJ<-log(nLoc,base = 2)
thinIn<-3

f.path<-"./MCMC_Res/"
f.name<-paste0(dfName,"_Kmax",K_max,".RDS")
fpdpm<-readRDS(paste0(f.path,f.name))
fpdpm.clu<-fpdpm$res$clu
rm(fpdpm)
n.pstSmp<-dim(fpdpm.clu)[3]
selIdx<-seq(0,n.pstSmp,thinIn)[-1]
pairDist.arr<-array(NA,dim=c(nSmp,nSmp,length(selIdx)))
wtVec<-1/(1:nJ[1])
wtIdx<-rep(2^(2* (1: nJ[1])-1) < nSmp)
wtVec[!wtIdx]<-wtVec[!wtIdx]/2
wtVec<-rep(wtVec, each=3)
for(sIdx in 1:length(selIdx)){
  ite<-selIdx[sIdx]
  # print(ite)
  diag(pairDist.arr[,,sIdx])<-0
  for(i1 in 1:(nSmp-1)){
    for(i2 in (i1+1):nSmp){
      pairDist.arr[i1,i2,sIdx]<-pairDist.arr[i2,i1,sIdx]<-
        dist_cluVec(cluVec1=fpdpm.clu[i1,-dim(fpdpm.clu)[2],ite], 
                    cluVec2 = fpdpm.clu[i2,-dim(fpdpm.clu)[2],ite], 
                    wtVec=wtVec)
    }
  }
}

saveRDS(pairDist.arr,file=paste0(f.path,"pairDist_",f.name))


nClu_maxSil<-30
distMat<-apply(pairDist.arr,c(1,2),mean)
silMat<-matrix(NA,nr=nClu_maxSil-1,nc=3)
tr_pdpm<-hclust(d = as.dist(distMat))
silVec_pdpm<-rep(NA,29)
for(k in 2:nClu_maxSil){
  tmp<-silhouette(cutree(tr_pdpm,k),dist = as.dist(distMat))
  tmp_sum<-summary(tmp)
  silVec_pdpm[k-1]<-tmp_sum$avg.width
}
nClu_sil<-which.max(silVec_pdpm)+1
cutree(tree = tr_pdpm,k=nClu_sil)
plot(tr_pdpm)
