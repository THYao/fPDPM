rm(list=ls())
library(stringr)
setwd("~/fPDPM/")
source("./Code/PDPM_dwt2d_LRD_func.R")

args <- commandArgs(TRUE)
# args<- c("diagErr_highSNR_nLoc32_nSmp300_b3","3")
print(args)
dfName<-trimws(args[1])
K_max<-as.integer(args[2])

obsDf<-readRDS(file=paste0("./Data/",dfName,".RDS"))$df

iteNum<-3000
burnIn<-iteNum*0.9
print(iteNum)
seed_= iteNum +K_max*12
ch1.t<-system.time({ch1.res<-fPDPM_mcmc(obsDf = obsDf, iteNum = iteNum, wf_="haar",
                                        burnIn=burnIn, seed_=seed_, K.max = K_max,
                                        KcutOff = 0.7, K.tol=1e-1, saveCov = F)})
plot(x=1:iteNum, y=ch1.res$llh)
print(ch1.t)
# file.exists(paste0("./MCMC_Res/",dfName,"_Kmax",K_max,".RDS"))
saveRDS(list(t=ch1.t,res=ch1.res),
        file=paste0("./MCMC_Res/",dfName,"_Kmax",K_max,".RDS"))

