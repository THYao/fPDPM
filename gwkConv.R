rm(list=ls())
library(waveslim)
library(ggplot2)
library(stringr)
library(coda)
setwd("~/fPDPM/")

args <- commandArgs(TRUE)
# args<- c("diagErr_highSNR_nLoc32_nSmp300_b3","3","1")
print(args)
dfName<-trimws(args[1])
K_max<-as.integer(args[2])
chIdx<-as.integer(args[3])

pdpm_muMat<-function(cluMat_, coefLs, wf_="haar"){
  n_<-nrow(cluMat_)
  nL_<-length(coefLs)
  D_type.vec<-c("LH","HL","HH")
  D_Num<-length(D_type.vec)
  nJ_ <- nL_ /D_Num
  res<-matrix(NA,nr=n_, ncol=2^(nJ_*2))
  for(i in 1:n_){
    clu_i<-cluMat_[i,]
    tmp.dwt<-dwt.2d(matrix(0,nr=2^nJ_,nc=2^nJ_), wf=wf_, J=nJ_)
    for(j in (1:nJ_)-1){
      for(d in 1:D_Num){
        l<-j *D_Num +d
        D_type<-D_type.vec[d]
        if(D_type!="HH"){
          corSgn<-(-1)
        }else{
          corSgn<-(1)
        }
        tmp.dwt[[paste0(D_type.vec[d],nJ_-j)]]<-matrix(coefLs[[l]][clu_i[l],]*corSgn,nrow=2^j)
      }
    }
    res[i,]<-idwt.2d(tmp.dwt)
  }
  return(res)
}

nLoc<-32
nSmp<-300
thinIn<-3
f.path<-"./MCMC_Res/"
paste0(dfName,"_Kmax",K_max,"_ch",chIdx,".RDS")
fpdpm.fit<-readRDS(file=paste0(f.path,f.name))$res
fpdpm.clu<-fpdpm.fit$clu
n_pstNum<-dim(fpdpm.clu)[3]
selIdx<-seq(0,n_pstNum,thinIn)[-1]
muArr<-array(NA,dim=c(nSmp,nLoc^2,length(selIdx)))
for(sIdx in 1:length(selIdx)){
  ite<-selIdx[sIdx]
  muArr[,,sIdx]<-pdpm_muMat(cluMat_ = fpdpm.clu[,,ite], coefLs = fpdpm.fit$coef[[ite]])
}
gwkRes.arr<-array(NA,dim=c(nSmp,nLoc^2))
for(i in 1:nSmp){
  print(i)
  gwkRes<-geweke.diag(t(muArr[i,,]))
  gwkRes.arr[i,]<-gwkRes$z
}

fdr_p_Prop<-rep(NA,nSmp)
for(i in 1:nSmp){
  fdr_p_Prop[i]<-sum(p.adjust(2*(1-pnorm(abs(gwkRes.arr[i,]))),method="BH")>0.05)/(nLoc^2)
}
summary(fdr_p_Prop)




