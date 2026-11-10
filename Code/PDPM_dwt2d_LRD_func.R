# Load required libraries
library(MASS)
library(extraDistr)
library(MCMCpack)
library(waveslim)
Rcpp::sourceCpp("./Code/dmv_LRD_woodbury.cpp")

# GenPsi_Haar.2d<-function(nLoc_, rslnIdx_, digit=6){
#   nJ_<-log2(nLoc_)
#   D_type.vec<-c("LH","HL","HH")
#   D_Num<-length(D_type.vec)
#   tmp.dwt<-dwt.2d(matrix(0,nr=nLoc_[1],nc=nLoc_[2]),wf='haar', J=nJ_[1])
#   n_<-prod(nLoc_)
#   tmp.coef<-matrix(rnorm(n_^2),nr=n_)
#   df<-matrix(NA,nr=n_,nc=n_)
#   for(i in 1:n_){
#     coef_i<-tmp.coef[i,]
#     tmp.dwt[[paste0("LL",nJ_[1])]]<-matrix(coef_i[rslnIdx_==-1])
#     for(j in (1:nJ_[1])-1){
#       for(d in 1:D_Num){
#         l<-j *D_Num +d
#         D_type<-D_type.vec[d]
#         if(D_type!="HH"){
#           corSgn<-(-1)
#         }else{
#           corSgn<-(1)
#         }
#         tmp.dwt[[paste0(D_type.vec[d],nJ_[1]-j)]]<-
#           matrix(coef_i[rslnIdx_==l]*corSgn,nrow=2^j)
#         df[i,]<-as.vector(idwt.2d(tmp.dwt))
#       }
#     }
#   }
#   Psi<-t(solve(tmp.coef) %*% df)
#   Psi<-round(Psi,digit)
#   return(Psi)
# }

# C = W W^t + s I_p
# C^{-1} = (1/s) I_p + (1/s) W M^{-1} W^t; M = W^t W + s I_k
# det(C) = det(s I_p +W W^t) = s^p det(I_p + 1/s W W^t) = s^p det(I_k + 1/s W^t W) = s^(p-k) det(s I_k + W^t W)
# log_dmvn_wdbInv<-function(xVec_, mu_, W_, s2_){
#   L_<-nrow(W_)
#   d_<-ncol(W_)
#   detCov<-log(det(s2_ * diag(d_) + t(W_) %*% W_)) + (L_-d_) * log(s2_)
#   M = t(W_) %*% W_ + s2_ * diag(d_)
#   M_inv<-solve(M)
#   WM_inv <- W_ %*% M_inv
#   # covInv<- (1/s2) *diag(p) - (1/s2) * WM_inv %*% t(W)
#   # max(abs(covInv - solve(S)))
#   aVec <- matrix(xVec_ - mu_, nc=1)
#   expNum<- sum(aVec^2) - t(aVec) %*% WM_inv %*% t(W_) %*% aVec
#   expNum <- expNum/s2_
#   res<- detCov + expNum + L_*log(2*pi)
#   res <- res * (-0.5)
#   return(as.numeric(res))
# }

update_u <- function(nSmp, weightLs, cluMat) {
  nL<-length(weightLs)
  u.res<-matrix(0,nr=nSmp,nc=nL)
  for(l in 1:nL){
    weightVec<-weightLs[[l]]
    clu_l<-cluMat[,l]
    for (i in 1:nSmp) {
      u.res[i,l] <- runif(1,0,weightVec[clu_l[i]])
    }
  }
  return(u.res)
}

update_stick<-function(cnctrnAlpha, uMat, stickLs, cluMat, weightLs){
  nSmp<-length(uMat)
  nL<-length(weightLs)
  weightRes<-stickRes<-list()
  for(l in 1:nL){
    stick_l<-stickLs[[l]]
    weight_l<-weightLs[[l]]
    nClu_l<-length(weight_l)
    clu_l<-cluMat[,l]
    u_l<-uMat[,l]
    for(clu in 1:nClu_l){
      cluIdx<-clu_l==clu
      bigCluIdx<-clu_l>clu
      if(any(cluIdx)){
        if(clu==1){
          alpha_clu<-max(u_l[cluIdx])
        }else{
          alpha_clu<-max(u_l[cluIdx]/prod(1-stick_l[1:(clu-1)]))
        }
      }else{
        alpha_clu<-0
      }
      if(any(bigCluIdx)){
        denom<-unlist(lapply(clu_l[bigCluIdx], function(x) stick_l[x]*prod(1-stick_l[setdiff(1:(x-1),clu)])))
        beta_clu<-1-max(u_l[bigCluIdx]/denom)
      }else{
        beta_clu<-1
      }
      stick_l[clu]<-truncdist::rtrunc(n=1, spec = "beta",a = alpha_clu, b = beta_clu, shape1=1, shape2=cnctrnAlpha)
      if(clu==1){
        weight_l[clu]<-stick_l[clu]
      }else{
        weight_l[clu]<-stick_l[clu] * prod(1-stick_l[1:(clu-1)])
      }
      # if(!all(weightRes[clu_l]>uVec)){
      #   stop(print(clu))
      # }
    }
    stickRes[[l]]<-stick_l
    weightRes[[l]]<-weight_l
  }
  return(list(stick=stickRes,weight=weightRes))
}

chk_nClu<-function(uVec, weightVec, stickVec, cnctrnAlpha){
  minU<-min(uVec)
  accuSumWeight<-cumsum(weightVec)
  idxCand<-accuSumWeight>1-min(uVec)
  while(!any(idxCand)){
    newStick<-rbeta(1,1,cnctrnAlpha)
    newWeight<-newStick*prod(1-stickVec)
    stickVec<-c(stickVec,newStick)
    weightVec<-c(weightVec,newWeight)
    accuSumWeight<-cumsum(weightVec)
    idxCand<-accuSumWeight>1-min(uVec)
  }
  return(list(stick=stickVec,weight=weightVec))
}

spawn_newClu_coef<-function(nClu, nCoef){
  res<-list()
  coefMat<-matrix(rnorm(nClu*nCoef),nrow=nClu,nc=nCoef)
  tauDiag<-matrix(1,nrow=nClu,nc=nCoef)
  res$coef<-coefMat
  res$tau<-tauDiag
  return(res)
}

spawn_newClu_diagVar<-function(nClu){
  res<-rep(1,nClu)
  return(res)
}

spawn_newClu_LRD<-function(nClu, K_, n_, L_, LRD_hyper){
  kappaLs<-deltaLs<-phiLs<-Lambda<-list()
  xiVec<-rep(NA, nClu)
  for(cIdx in 1:nClu){
    Lambda_r<-matrix(0, nrow=L_,ncol=K_)
    phi_r<-matrix(rgamma(L_*K_,shape = LRD_hyper['nu1']/2, rate = LRD_hyper['nu2']/2),nrow=L_)
    deltaVec_r<-c(rgamma(1, LRD_hyper['a1'], rate=1), rgamma(K_-1,LRD_hyper['a2'], rate=1))
    kappaVec_r<-cumprod(deltaVec_r)
    xiVec[cIdx] <- rgamma(1, LRD_hyper['ax'], rate=LRD_hyper['bx'])
    
    Lambda[[cIdx]]<-Lambda_r
    phiLs[[cIdx]]<-phi_r
    deltaLs[[cIdx]]<-deltaVec_r
    kappaLs[[cIdx]]<-kappaVec_r
  }
  return(list(Lambda = Lambda, phi=phiLs, 
              delta = deltaLs, kappa=kappaLs, xi=xiVec))
}

update_tauVar<-function(coefBeta, lambda=1){
  coefLen<-length(coefBeta)
  # invVar<-LaplacesDemon::rinvgaussian(coefLen,
  #                                     mu=sqrt(1/(lambda * coefBeta^2)),
  #                                     lambda = lambda)
  invVar<-rwald(coefLen,
                mu=abs(lambda/coefBeta),
                lambda = lambda^2)
  newVar<-1/invVar
  return(newVar)
}

update_coef<-function(obsDf, cluMat, paramList, wf_){
  nSmp<-dim(obsDf)[3]
  nLoc<-dim(obsDf)[1:2]
  nJ<-log2(nLoc)
  nL<-ncol(cluMat)
  D_type.vec<-c("LH","HL","HH")
  D_Num<-length(D_type.vec)
  coefParam<- paramList$coef
  tauDiag<-paramList$tau
  varParam<-paramList$var
  facMean<-gen_facMean(Lambda_ = paramList$LRD$Lambda, eta_ = paramList$LRD$eta, cluVec = cluMat[,nL])
  facMean<-array(facMean,dim=c(nLoc,nSmp))
  # tmp.wd<-wd(rep(0,nT), filter.number = 1, family="DaubExPhase")
  # Psi<-GenW(n=nT, filter.number = 1, family="DaubExPhase")
  # waveMat<-t(Psi)
  for(j in (1:nJ[1])-1){
    for(d in 1:D_Num){
      l<-j *D_Num +d
      D_type<-D_type.vec[d]
      if(D_type!="HH"){
        corSgn<-(-1)
      }else{
        corSgn<-(1)
      }
      cluVec_l<-cluMat[,l]
      uq.clu<-sort(unique(cluVec_l))
      for(uqIdx in 1:length(uq.clu)){
        clu<-uq.clu[uqIdx]
        lc_Idx<-cluVec_l==clu
        n_lc<-sum(lc_Idx)
        y_lc<-array(obsDf[,,lc_Idx],dim=c(nLoc,n_lc))
        facMean_lc<-array(facMean[,,lc_Idx],dim=c(nLoc,n_lc))
        # yBarjc<-apply(df_lc,c(1:2),mean)
        var_lc<-varParam[cluMat[lc_Idx,nL]]
        varArr_lc<-array(rep(var_lc,each=prod(nLoc)),dim=c(nLoc,n_lc))
        aArr_lc<-(y_lc - facMean_lc)/varArr_lc
        aMat_lc<-apply(aArr_lc, c(1:2),sum)
        
        tau_diag_lc<-1/(sum(1/var_lc)+1/tauDiag[[l]][clu,])
        # newHyperMeanjc<- as.numeric(matrix((n_jc/(sdSigma^2)) * tau_diag_jc * t(t(Psi_[,rslnIdx ==j]) %*% yBarjc),nrow=1))
        # tmp.dwt<-dwt.2d(yBarjc, wf='haar', J=nJ[1])
        tmp.dwt<-dwt.2d(aMat_lc, wf=wf_, J=nJ[1])
        tmp.dwt_l<-as.numeric(tmp.dwt[[paste0(D_type,nJ[1]-j)]] * corSgn)
        coef_lc<- as.numeric(matrix(tau_diag_lc * tmp.dwt_l ,nrow=1))
        
        newCoef_lc<-rnorm(2^(2*j),mean=coef_lc, sd= sqrt(tau_diag_lc))
        coefParam[[l]][clu,]<-newCoef_lc
        
        newVarTau<-update_tauVar(coefBeta = newCoef_lc)
        tauDiag[[l]][clu,]<-newVarTau
      }
    }
  }
  return(list(coef=coefParam, tau=tauDiag, facMean = facMean))
}

update_varDiag<-function(obsDf, cluMat, paramList, wf_, facMean_ = NULL,a_s=3,b_s=3){
  nSmp<-dim(obsDf)[3]
  nLoc<-dim(obsDf)[1:2]
  nJ<-log2(nLoc)
  nL<-ncol(cluMat)
  D_type.vec<-c("LH","HL","HH")
  D_Num<-length(D_type.vec)
  coefParam<- paramList$coef
  resVar<-paramList$var
  rsdu.vec<-rep(NA,nSmp)
  tmp.dwt<-dwt.2d(matrix(0,nr=nLoc[1],nc=nLoc[1]), wf=wf_, J=nJ[1])
  # facMean<-paramList$LRD$eta %*% t(paramList$LRD$Lambda)
  if(is.null(facMean_)){
    facMean_ <- gen_facMean(Lambda_ = paramList$LRD$Lambda, eta_ = paramList$LRD$eta, cluVec = cluMat[,nL])
    facMean_ <- array(facMean_,dim=c(nLoc,nSmp))
  }
  muMat_ <- gen_muMat(cluMat_ = cluMat, coefLs = coefParam, wf_=wf_)
  # muMat_<-matrix(NA,nr=nSmp, nc = prod(nLoc))
  for(i in 1:nSmp){
    clu_i<-cluMat[i,]
    # for(j in (1:nJ[1])-1){
    #   for(d in 1:D_Num){
    #     l<-j *D_Num +d
    #     D_type<-D_type.vec[d]
    #     if(D_type!="HH"){
    #       corSgn<-(-1)
    #     }else{
    #       corSgn<-(1)
    #     }
    #     tmp.dwt[[paste0(D_type.vec[d],nJ[1]-j)]]<-matrix(coefParam[[l]][clu_i[l],]*corSgn,nrow=2^j)
    #   }
    # }
    # muMat_[i,]<-idwt.2d(tmp.dwt)
    theta_i<-muMat_[i,] + as.numeric(facMean_[,,i])
    rsdu_i<-sum(( as.numeric(obsDf[,,i]) - theta_i)^2)
    rsdu.vec[i]<-rsdu_i
  }
  cluVec_s<-cluMat[,nL]
  uq.clu<-sort(unique(cluVec_s))
  for(uqIdx in 1:length(uq.clu)){
    clu<-uq.clu[uqIdx]
    sc_Idx<-cluVec_s==clu
    n_sc<-sum(sc_Idx)
    new_a<- a_s + (n_sc*prod(nLoc))/2
    new_b<- b_s + sum(rsdu.vec[sc_Idx])/2
    newVar_sc<-1/rgamma(1,shape = new_a, rate=new_b)
    resVar[clu]<-newVar_sc
  }
  return(list(var=resVar, muMat = muMat_))
}

update_clu<-function(obsDf, weightLs, stickLs, uMat, LRD_hyper, wf_,
                     cluMat, paramList, cnctrnAlpha){
  nSmp<-dim(obsDf)[3]
  nLoc<-dim(obsDf)[1:2]
  nJ<-log2(nLoc)
  nL<-length(weightLs)
  D_type.vec<-c("LH","HL","HH")
  D_Num<-length(D_type.vec)
  
  newStick<-stickLs
  newWeight<-weightLs
  newCoefLs<-paramList$coef
  newTau<-paramList$tau
  newVar<-paramList$var
  newLRD<-paramList$LRD
  # facMean<-paramList$LRD$eta %*% t(paramList$LRD$Lambda)
  # set.seed(1)
  for(l in 1:nL){
    tmp.stk<-chk_nClu(uVec=uMat[,l], weightVec = weightLs[[l]],
                      stickVec = stickLs[[l]], cnctrnAlpha = cnctrnAlpha)
    newStick[[l]]<-tmp.stk$stick
    newWeight[[l]]<-tmp.stk$weight
    new_nClu<-length(newStick[[l]])
    if(l!=nL){
      j<-floor((l-1)/D_Num)
      if(new_nClu != length(weightLs[[l]])){
        add.CoefLs<-spawn_newClu_coef(nClu = new_nClu - length(weightLs[[l]]), nCoef = 2^(2*j))
        newCoefLs[[l]]<-rbind(paramList$coef[[l]],add.CoefLs$coef)
        newTau[[l]]<-rbind(paramList$tau[[l]],add.CoefLs$tau)
      }else{
        newCoefLs[[l]]<-paramList$coef[[l]]
        newTau[[l]]<-paramList$tau[[l]]
      }
    }else{
      if(new_nClu != length(weightLs[[l]])){
        add.VarLs<-spawn_newClu_diagVar(nClu = new_nClu - length(weightLs[[l]]))
        newVar<-c(newVar,add.VarLs)
        add.LRD<-spawn_newClu_LRD(nClu = new_nClu - length(weightLs[[l]]), 
                                  K_ = LRD_hyper['Kmax'], n_ = nSmp, L_ = prod(nLoc),
                                  LRD_hyper = LRD_hyper)
        newLRD$Lambda <- c(paramList$LRD$Lambda, add.LRD$Lambda)
        newLRD$phi<- c(paramList$LRD$phi, add.LRD$phi)
        newLRD$delta <- c(paramList$LRD$delta, add.LRD$delta)
        newLRD$kappa <- c(paramList$LRD$kappa, add.LRD$kappa)
        newLRD$xi <- c(paramList$LRD$xi, add.LRD$xi)
      }
    }
  }
  # subMeanMat.tmp<-cppFun_gen_subMeanMat(incpt = paramList$int.mean, cluMat_ = cluMat, meanMatLs_ = newCoefLs)
  # cluRes.tmp<-cppFun_update_clu(obsDf_ = obsDf,sigma_ = sigma_,
  #                               cluMat_ = cluMat, subMeanMat_ = subMeanMat.tmp,
  #                               meanMatLs_ = newCoefLs, uMat_ = uMat,
  #                               weightLs_ = newWeight, rslnIdx_ = rsolnIdx)
  
  # coefMat_n<-matrix(NA,nrow=nrow(obsDf),ncol=nL)
  cluRes<-cluMat
  newLambda<-newLRD$Lambda
  etaLs <- newLRD$eta
  tmp.dwt<-dwt.2d(matrix(0, nr=nLoc[1],nc=nLoc[2]),wf=wf_, J=nJ[1])
  # facMean <- gen_facMean(Lambda_ = newLRD$Lambda, eta_ = newLRD$eta, cluVec = cluMat[,nL])
  for(i in 1:nSmp){
    # print(i)
    clu_i<-cluMat[i,]
    eta_i<-etaLs[[i]]
    tmp.dwt[[paste0("LL",nJ[1])]]<-matrix(paramList$incpt[i])
    for(j in (1:nJ[1])-1){
      for(d in 1:D_Num){
        l<-j*D_Num+d
        D_type<-D_type.vec[d]
        if(D_type!="HH"){
          tmp.dwt[[paste0(D_type,nJ[1]-j)]]<-matrix(newCoefLs[[l]][clu_i[l],]*(-1),nrow=2^j)
        }else{
          tmp.dwt[[paste0(D_type,nJ[1]-j)]]<-matrix(newCoefLs[[l]][clu_i[l],],nrow=2^j)
        }
      }
    }
    # for(tmpL in 1:n){
    #     coefVec_i<-c(coefVec_i, newCoefLs[[tmpJ]][clu_i[tmpJ],])
    # }
    for(l in 1:nL){
      j<-floor((l-1)/D_Num)
      d<-l-j*D_Num
      D_type<-D_type.vec[d]
      if(D_type!="HH"){
        corSgn<-(-1)
      }else{
        corSgn<-(1)
      }
      idxCand<-newWeight[[l]]>uMat[i,l]
      # new_nClu_l<-length(idxCand)
      # prob_i<-rep(-Inf,new_nClu_l)
      candClu_l<-which(idxCand)
      tmpEta_i<-list()
      if(length(candClu_l)==1){
        clu_il<-clu_i[l]<-candClu_l
        if(l==nL){
          Lambda_ik<- newLambda[[clu_il]]
          K_i <- ncol(Lambda_ik)
          if(K_i > nrow(eta_i)){
            eta_ik<-matrix(c(eta_i, rep(0, K_i - nrow(eta_i))),nc=1)
          }else if(K_i < nrow(eta_i)){
            eta_ik<-matrix(eta_i[1:K_i,],nc=1)
          }else{
            eta_ik <- eta_i
          }
          tmpEta_i[[1]]<-eta_ik
        }
      }else{
        prob_i<-rep(0,length(candClu_l))
        for(kIdx in 1:length(candClu_l)){
          k<-candClu_l[kIdx]
          if(l!=nL){
            tmp.dwt[[paste0(D_type,nJ[1]-j)]]<-matrix(newCoefLs[[l]][k,]*corSgn,nrow=2^j)
            # mean_i <- as.numeric(idwt.2d(tmp.dwt)) + facMean[i,]
            # prob_i[kIdx]<-sum(dnorm(as.numeric(obsDf[,,i]),mean=mean_i, sd=sqrt(newVar[clu_i[nL]]), log=T))
            mean_ik <- as.numeric(idwt.2d(tmp.dwt))
            prob_i[kIdx]<-log_dmvn_wdb(xVec_ = as.numeric(obsDf[,,i]), mu_ = mean_ik, 
                                          W_ = newLambda[[clu_i[nL]]], s2_ = newVar[clu_i[nL]])
          }else{
            Lambda_ik<- newLambda[[k]]
            K_i <- ncol(Lambda_ik)
            if(K_i > nrow(eta_i)){
              eta_ik<-matrix(c(eta_i, rep(0, K_i - nrow(eta_i))),nc=1)
            }else if(K_i < nrow(eta_i)){
              eta_ik<-matrix(eta_i[1:K_i,],nc=1)
            }else{
              eta_ik <- eta_i
            }
            tmpEta_i[[kIdx]]<-eta_ik
            # facMean_ik<- Lambda_ik %*% eta_ik
            # mean_ik <- as.numeric(idwt.2d(tmp.dwt)) + facMean_ik
            # prob_i[kIdx]<-sum(dnorm(as.numeric(obsDf[,,i]),mean=mean_ik, sd=sqrt(newVar[k]), log=T))
            # facMean_ik<- Lambda_ik %*% eta_ik
            mean_ik <- as.numeric(idwt.2d(tmp.dwt))
            # var_ik<- Lambda_ik %*% t(Lambda_ik) + newVar[k] * diag(prod(nLoc))
            # prob_i[kIdx]<- sum(mvnfast::dmvn(as.numeric(obsDf[,,i]), mu = mean_ik, sigma = var_ik, log = T))
            prob_i[kIdx]<-log_dmvn_wdb(xVec_ = as.numeric(obsDf[,,i]), mu_ = mean_ik, W_ = Lambda_ik, s2_ = newVar[k])
          }
        }
        prob_i<-prob_i-matrixStats::logSumExp(prob_i)
        # clu_i[j+1]<-clu_ij<-which.max(exp(prob_i))
        clu_il<-clu_i[l]<-sample(candClu_l,1,prob=exp(prob_i))
      }
      if(l!=nL){
        tmp.dwt[[paste0(D_type,nJ[1]-j)]]<-matrix(newCoefLs[[l]][clu_il,]*corSgn,nrow=2^j)
      }else{
        etaLs[[i]]<-tmpEta_i[[which(candClu_l==clu_i[l])]]
      }
    }
    cluRes[i,]<-clu_i
  }
  newLRD$eta <- etaLs
  return(list(clu=cluRes, var=newVar, LRD=newLRD,
              coef=newCoefLs, tau=newTau,
              stick=newStick, weight=newWeight))
}

gen_muMat<-function(cluMat_, coefLs, wf_){
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

gen_facMean<-function(Lambda_, eta_, cluVec){
  n_<-length(cluVec)
  L_<-dim(Lambda_[[1]])[1]
  res<-matrix(NA,nr=n_,nc=L_)
  for(i in 1:n_){
    clu_i<-cluVec[i]
    Lambda_i<-Lambda_[[clu_i]]
    # K_i<-ncol(Lambda_i)
    res[i,]<- Lambda_i %*% matrix(eta_[[i]],nc=1)
  }
  return(res)
}

update_lambda<-function(dfMat_, muMat_, diagVar_, cluVec, LRD_){
  etaLs<- LRD_$eta
  phiLs<- LRD_$phi 
  kapLs<- LRD_$kappa
  xiVec <- LRD_$xi
  # K_<- ncol(etaMat)
  L_<- ncol(dfMat_)
  n_<- nrow(dfMat_)
  res<-LRD_$Lambda
  # sigma2Vec<-diagVar_[cluVec]
  # eta2_sigma<-matrix(0, nr=K_, nc=K_)
  # S_inv<- phiMat_ * matrix(rep(kappaVec_, each=L_),nrow=L_)
  # res<-matrix(NA, nr=L_, nc=K_)
  # for(i in 1:n_){
  #   eta2_sigma <- eta2_sigma + matrix(etaMat_[i,],nc=1) %*% matrix(etaMat_[i,],nr=1)/sigma2Vec[i]
  # }
  # for(l in 1:L_){
  #   if(length(S_inv[l,])==1){
  #     precMat <- S_inv[l,] + eta2_sigma
  #   }else{
  #     precMat <- diag(S_inv[l,]) + eta2_sigma
  #   }
  #   varMat <- solve(precMat)
  #   aVec<-rep(0, K_)
  #   for(i in 1:n_){
  #     aVec <- aVec + etaMat_[i,] * (dfMat_[i,l] - muMat_[i,l])/sigma2Vec[i]
  #   }
  #   meanVec <- varMat %*% aVec
  #   res[l,]<- MASS::mvrnorm(1, mu=meanVec, Sigma = varMat)
  # }
  
  uq.clu<-sort(unique(cluVec))
  for(uqIdx in 1:length(uq.clu)){
    clu<-uq.clu[uqIdx]
    iIdx_s<-cluVec==clu
    n_s<-sum(iIdx_s)
    sigma2_s<-diagVar_[clu]
    K_s <- ncol(LRD_$Lambda[[clu]])
    etaMat_s <- t(matrix(unlist(etaLs[iIdx_s]),nc=n_s))
    eta_sigmaInv_s<-etaMat_s/sigma2_s
    eta2_sigmaInv_s = t(eta_sigmaInv_s) %*% etaMat_s
    Ss_inv<- phiLs[[clu]] * matrix(rep(kapLs[[clu]], each=L_),nrow=L_) *xiVec[clu]
    aMat_s<- t(matrix(dfMat_[iIdx_s,] -muMat_[iIdx_s,],nr=n_s)) %*% eta_sigmaInv_s
    res_s<-matrix(NA,nr=L_,nc=K_s)
    for(l in 1:L_){
      Ssl_inv<-Ss_inv[l,]
      if(length(Ssl_inv)==1){
        precMat= Ssl_inv + eta2_sigmaInv_s
      }else{
        precMat=diag(Ssl_inv) + eta2_sigmaInv_s
      }
      aVec_s<- matrix(aMat_s[l,],nc=1)
      prec_chol<-chol(precMat)
      sqrtVar_aVec<-solve(t(prec_chol), aVec_s)
      meanVec <- solve(prec_chol, sqrtVar_aVec)
      z<-matrix(rnorm(K_s), ncol = 1)
      sqrtVar_z <- solve(prec_chol, z)
      res_s[l,]<-meanVec + sqrtVar_z
    }
    res[[clu]]<-res_s
  }
  # sigma2Mat<-matrix(rep(sigma2Vec, each=K_),byrow=T,ncol=K_)
  # eta_sigmaInv<-etaMat_/sigma2Mat
  # eta2_sigma = t(eta_sigmaInv) %*% etaMat_
  # S_inv<- phiMat_ * matrix(rep(kappaVec_, each=L_),nrow=L_)
  # aMat<- t(dfMat_ -muMat_) %*% eta_sigmaInv
  # res<-matrix(NA,nr=L_,nc=K_)
  # for(l in 1:L_){
  #   Sl_inv<-S_inv[l,]
  #   if(length(Sl_inv)==1){
  #     precMat= Sl_inv + eta2_sigma
  #   }else{
  #     precMat=diag(Sl_inv) + eta2_sigma
  #   }
  #   aVec<- matrix(aMat[l,],nc=1)
  #   prec_chol<-chol(precMat)
  #   sqrtVar_aVec<-solve(t(prec_chol), aVec)
  #   meanVec <- solve(prec_chol, sqrtVar_aVec)
  #   z<-matrix(rnorm(K_), ncol = 1)
  #   sqrtVar_z <- solve(prec_chol, z)
  #   res[l,]<-meanVec + sqrtVar_z
  # }
  return(res)
}

update_eta<-function(LRD_, diagVar_, dfMat_, muMat_, cluVec){
  n_<-nrow(dfMat_)
  res<-LRD_$eta
  # K_<-ncol(res)
  LamLs<-LRD_$Lambda
  sigma2Vec<-diagVar_[cluVec]
  uq.clu<-sort(unique(cluVec))
  varMatLs<-varMat_qrLs<-list()
  for(uqIdx in 1:length(uq.clu)){
    clu<-uq.clu[uqIdx]
    sigma2_s<-diagVar_[clu]
    Lam_s<- LamLs[[clu]]
    K_s <- ncol(Lam_s)
    Lam2_s <- t(Lam_s) %*% Lam_s
    precMat= diag(K_s) + Lam2_s/sigma2_s
    prec_chol<-chol(precMat)
    prec_qr<-qr(prec_chol)
    varMat_qr <- solve(prec_qr$qr)
    varMat <- varMat_qr %*% t(varMat_qr)
    varMatLs[[uqIdx]]<-varMat
    varMat_qrLs[[uqIdx]]<-varMat_qr
  }
  # Lam2<-t(Lambda_) %*% Lambda_
  # for(i in 1:n_){
  #   sigma2_i<-sigma2Vec[i]
  #   precMat= diag(K_) + Lam2/sigma2_i
  #   varMat <- solve(precMat)
  #   meanVec<-varMat %*% t(Lambda_) %*% (dfMat_[i,] - muMat_[i,]) /sigma2_i
  #   res[i,]<- MASS::mvrnorm(1, mu=meanVec, Sigma = varMat)
  # }
  for(i in 1:n_){
    sigma2_i<-sigma2Vec[i]
    clu_i<-cluVec[i]
    # Lam_i<-LamLs[[cluVec[i]]]
    # Lam2<-t(Lam_i) %*% Lam_i
    # precMat= diag(K_) + Lam2/sigma2_i
    # prec_chol<-chol(precMat)
    # prec_qr<-qr(prec_chol)
    # varMat_qr <- solve(prec_qr$qr)
    sIdx<-which(uq.clu==clu_i)
    varMat <- varMatLs[[sIdx]]
    K_i <- nrow(varMat)
    meanVec<-varMat %*% t(LamLs[[clu_i]]) %*% (dfMat_[i,] - muMat_[i,]) /sigma2_i
    res[[i]]<-meanVec + varMat_qrLs[[sIdx]] %*% matrix(rnorm(K_i),nc=1)
  }
  return(res)
}

update_phi<-function(nu_1, nu_2, LRD_, cluVec){
  LamLs <- LRD_$Lambda
  kapLs <- LRD_$kappa
  xiVec <- LRD_$xi
  res<-LRD_$phi
  L_<-nrow(LamLs[[1]])
  uq.clu<-sort(unique(cluVec))
  for(uqIdx in 1:length(uq.clu)){
    clu<-uq.clu[uqIdx]
    kappaMat<-matrix(rep(kapLs[[clu]], each=L_),nrow=L_)
    Lam_s<-LamLs[[clu]]
    res_s<-rgamma(prod(dim(Lam_s)), nu_1+0.5, rate=nu_2 + (xiVec[clu]* as.vector(kappaMat * Lam_s^2))/2)
    res_s<-matrix(res_s, nr=L_)
    res[[clu]]<-res_s
  }
  return(res)
  # kappaMat<-matrix(rep(kappaVec_, each=L_),nrow=L_)
  # res<-rgamma(prod(dim(Lambda_)), (nu_+1)/2, rate=(nu_ + as.vector(kappaMat * Lambda_^2))/2)
  # return(matrix(res,nrow=L_))
}

update_delta_kappa<-function(a_1, a_2, LRD_, cluVec){
  phiLs <- LRD_$phi
  LamLs <- LRD_$Lambda
  xiVec <- LRD_$xi
  res.delta<-LRD_$delta
  res.kappa<-LRD_$kappa
  
  L_<-nrow(LamLs[[1]])
  uq.clu<-sort(unique(cluVec))
  for(uqIdx in 1:length(uq.clu)){
    clu<-uq.clu[uqIdx]
    tmpMat<-  phiLs[[clu]] * LamLs[[clu]]^2
    resDelta_s<-res.delta[[clu]]
    resKappa_s <- cumprod(resDelta_s)
    K_s<- ncol(LamLs[[clu]])
    for(k in 1:K_s){
      if(k==1){
        new_a<-a_1+L_*K_s/2
        new_b<-1+0.5* xiVec[clu] * sum(resKappa_s * colSums(tmpMat)/resDelta_s[1])
      }else{
        new_a<-a_2+L_*(K_s-k+1)/2
        new_b<-1+0.5* xiVec[clu]* sum(resKappa_s[k:K_s] * colSums(matrix(tmpMat[, k:K_s],nrow=L_))/resDelta_s[k])
      }
      resDelta_s[k] <- rgamma(1, new_a, rate=new_b)
      resKappa_s <- cumprod(resDelta_s)
    }
    res.delta[[clu]]<-resDelta_s
    res.kappa[[clu]]<-resKappa_s
  }
  # tmpMat<- phiMat_ * Lambda_^2
  # resKappa<-cumprod(deltaVec_)
  # resDelta<-deltaVec_
  # for(k in 1:K_){
  #   if(k==1){
  #     new_a<-a_1+L_*K_/2
  #     new_b<-1+0.5* sum(resKappa * colSums(tmpMat)/resDelta[1])
  #   }else{
  #     new_a<-a_2+L_*(K_-k+1)/2
  #     new_b<-1+0.5* sum(resKappa[k:K_] * colSums(matrix(tmpMat[, k:K_],nrow=L_))/resDelta[k])
  #   }
  #   resDelta[k] <- rgamma(1, new_a, rate=new_b)
  #   resKappa <- cumprod(resDelta)
  # }
  return(list(kappa=res.kappa,delta=res.delta))
}

update_xi<-function(LRD_, a_x, b_x, cluVec){
  phiLs <- LRD_$phi
  LamLs <- LRD_$Lambda
  kapLs <- LRD_$kappa
  res.xiVec<-LRD_$xi
  L_<-nrow(LamLs[[1]])
  
  uq.clu<-sort(unique(cluVec))
  for(uqIdx in 1:length(uq.clu)){
    clu<-uq.clu[uqIdx]
    tmpMat<- phiLs[[clu]] * LamLs[[clu]]^2
    K_s<- ncol(LamLs[[clu]])
    new_a<-a_x+L_*K_s/2
    new_b<-b_x+0.5* sum(kapLs[[clu]] * colSums(tmpMat))
    res.xiVec[clu]<-rgamma(1, new_a, rate=new_b)
  }
  return(res.xiVec)
}

K_adaption<-function(b0 = 1, b1 = 0.0005, tol, propCutOff,
                     LRD_, ite_, cluVec, Kmax){
  LamLs <- LRD_$Lambda
  L_ <- nrow(LamLs[[1]])
  # K_ <- ncol(Lambda_)
  prob <- 1 / exp(b0 + b1 * ite_)  # probability of adapting
  uq.clu<-sort(unique(cluVec))
  K_res<-list()
  change_Flg<-rep(F,length(uq.clu))
  for(uqIdx in 1:length(uq.clu)){
    clu<-uq.clu[uqIdx]
    u_s <- runif(1)
    K_s <- ncol(LamLs[[clu]])
    lind_s <- colSums(abs(LamLs[[clu]]) < tol) / L_  # proportion of elements in each column less than eps in magnitude
    rdntCol_s <- lind_s >= propCutOff
    rdntNum_s <- sum(rdntCol_s)  # number of redundant columns
    if (u_s < prob) {
      if (ite_ > 20 && rdntNum_s == 0 && all(lind_s < 0.995)) {
      # if (ite_ > 20 && rdntNum_s == 0) {
        if(K_s < Kmax){
          K_s <- K_s + 1
          change_Flg[uqIdx]<-T
        }
      } else if (rdntNum_s > 0) {
        # nonRdnt <- setdiff(1:K_, which(rdntCol))
        # K_s <- max(K_s - rdntNum_s, 1)
        K_s <- K_s - rdntNum_s
        change_Flg[uqIdx]<-T
        # if(K_s==0){
        #   rdntCol_s<-c(F,rep(T, ncol(LamLs[[clu]])-1))
        # }
      }
    }
    K_res[[uqIdx]]<-list(k=K_s,rdntIdx=rdntCol_s)
  }
  return(list(uqClu=uq.clu, K_info=K_res, changeFlg=change_Flg))
  # u <- runif(1)
  # lind <- colSums(abs(Lambda_) < tol) / L_  # proportion of elements in each column less than eps in magnitude
  # rdntCol <- lind >= propCutOff
  # rdntNum <- sum(rdntCol)  # number of redundant columns
  # if (u < prob) {
  #   if (ite_ > 20 && rdntNum == 0 && all(lind < 0.995)) {
  #     K_ <- K_ + 1
  #   } else if (rdntNum > 0) {
  #     # nonRdnt <- setdiff(1:K_, which(rdntCol))
  #     K_ <- max(K_ - rdntNum, 1)
  #   }
  # }
  # return(list(k=K_,rdntIdx=rdntCol))
}

update_LRD<-function(obsDf_, cluVec, muMat_, LRD, diagVar_, LRD_hyper , ite_){
  n_<-length(cluVec)
  nLoc<-dim(obsDf_)[1:2]
  obsDf_<-matrix(obsDf_, ncol=prod(dim(obsDf_)[1:2]))
  res.LRD<-LRD
  
  etaLs <- update_eta(LRD_ = res.LRD, diagVar_ = diagVar_, cluVec= cluVec,
                      dfMat_ = obsDf_, muMat_ = muMat_)
  res.LRD$eta <- etaLs
  rm(etaLs)
  
  Lambda <- update_lambda(dfMat_=obsDf_, muMat_=muMat_, cluVec = cluVec,
                          diagVar_=diagVar_, LRD_ = res.LRD)
  res.LRD$Lambda <- Lambda
  rm(Lambda)
  
  phi <- update_phi(LRD_ = res.LRD, cluVec=cluVec,
                    nu_1=LRD_hyper['nu1'], nu_2 = LRD_hyper['nu2'])
  res.LRD$phi<-phi
  rm(phi)
  
  tmp <- update_delta_kappa(LRD_ = res.LRD, cluVec = cluVec,
                            a_1 = LRD_hyper['a1'], a_2 = LRD_hyper['a2'])
  res.LRD$delta <- tmp$delta
  res.LRD$kappa<- tmp$kappa
  rm(tmp)
  
  xi <- update_xi(LRD_ = res.LRD, cluVec= cluVec,
                  a_x = LRD_hyper['ax'], b_x = LRD_hyper['bx'])
  res.LRD$xi <- xi
  rm(xi)
  
  k_res<-K_adaption(LRD_ = res.LRD, ite_ = ite_, cluVec = cluVec, 
                    Kmax = LRD_hyper['Kmax'],tol = LRD_hyper['Ktol'], 
                    propCutOff = LRD_hyper['KcutOff'])
  uq.clu<-k_res$uqClu
  for(uqIdx in 1:length(uq.clu)){
    if(!k_res$changeFlg[uqIdx]){
      next
    }
    clu<-uq.clu[uqIdx]
    iIdx_s<-cluVec==clu
    n_s<-sum(iIdx_s)
    LamMat_s<-res.LRD$Lambda[[clu]]
    phiMat_s<-res.LRD$phi[[clu]]
    delVec_s<-res.LRD$delta[[clu]]
    kapVec_s<-res.LRD$kappa[[clu]]
    K_s <- k_res$K_info[[uqIdx]]$k
    etaMat_s <- t(matrix(unlist(res.LRD$eta[iIdx_s]),nc=n_s))
    rdntCol_s <- k_res$K_info[[uqIdx]]$rdntIdx
    if(K_s > ncol(LamMat_s)){
      LamMat_s <- cbind(LamMat_s, matrix(0, nrow = prod(nLoc), ncol = 1))
      # etaMat <- cbind(etaMat, matrix(rnorm(nSmp), ncol = 1))
      phiMat_s <- cbind(phiMat_s, rgamma(prod(nLoc), shape = LRD_hyper['nu1']/2, rate =LRD_hyper['nu2']/2))
      delVec_s <- c(delVec_s, rgamma(1, LRD_hyper['a2'], rate=1))
      kapVec_s <- cumprod(delVec_s)
      etaMat_s <- cbind(etaMat_s, rnorm(n_s))
      res.LRD$eta[iIdx_s]<-lapply(seq_len(nrow(etaMat_s)), function(i) matrix(etaMat_s[i,],nc=1))
    }else if(K_s<ncol(LamMat_s)){
      if(K_s==0){
        LamMat_s <- matrix(0,nr=nrow(LamMat_s))
        # etaMat <- matrix(etaMat[, nonRdnt],nr=nrow(etaMat))
        nonRdnt_s <- 1
      }else{
        nonRdnt_s <- setdiff(1:ncol(LamMat_s), which(rdntCol_s))
        LamMat_s <- matrix(LamMat_s[, nonRdnt_s],nr=nrow(LamMat_s))
        # etaMat <- matrix(etaMat[, nonRdnt],nr=nrow(etaMat))
      }
      phiMat_s <- matrix(phiMat_s[, nonRdnt_s],nr=nrow(phiMat_s))
      delVec_s <- delVec_s[nonRdnt_s]
      kapVec_s <- cumprod(delVec_s)
      etaMat_s <- matrix(etaMat_s[,nonRdnt_s],nr=n_s)
      res.LRD$eta[iIdx_s]<-lapply(seq_len(nrow(etaMat_s)), function(i) matrix(etaMat_s[i,],nc=1))
    }
    res.LRD$Lambda[[clu]] <- LamMat_s
    res.LRD$phi[[clu]] <- phiMat_s
    res.LRD$delta[[clu]] <- delVec_s
    res.LRD$kappa[[clu]] <- kapVec_s
  }
  return(res.LRD)
}

llh_cal<-function(obsDf, cluMat, weightLs, paramList, wf_){
  nLoc<-dim(obsDf)[1:2]
  nSmp<-dim(obsDf)[3]
  nL<-length(weightLs)
  varVec<-paramList$var
  LambdaLs <- paramList$LRD$Lambda
  # facMean<-paramList$LRD$eta %*% t(paramList$LRD$Lambda)
  # facMean<- gen_facMean(Lambda_ = paramList$LRD$Lambda, eta_ = paramList$LRD$eta, cluVec = cluMat[,nL])
  muMat <- gen_muMat(cluMat_ = cluMat, coefLs = paramList$coef, wf_=wf_)
  # res<-0
  llh_res<-matrix(NA,nr=nSmp, nc=nL+1)
  for(i in 1:nSmp){
    clu_i<-cluMat[i,]
    # w_i<-0
    mu_i<-muMat[i,]
    # facMean_i <-facMean[i,]
    for(l in 1:nL){
      # w_i<-w_i+log(weightLs[[l]][clu_i[l]])
      llh_res[i,l]<-log(weightLs[[l]][clu_i[l]])
    }
    # mean_i<-mu_i + facMean_i
    # llh_i<-sum(dnorm(as.numeric(obsDf[,,i]), mean= mean_i, sd=sqrt(varVec[clu_i[nL]]), log=T))
    # llh_res[i,nL+1]<-sum(dnorm(as.numeric(obsDf[,,i]), mean= mean_i, sd=sqrt(varVec[clu_i[nL]]), log=T))
    llh_res[i,nL+1]<-log_dmvn_wdb(xVec_ = as.numeric(obsDf[,,i]), mu_ = mu_i, 
                                  W_ = LambdaLs[[clu_i[nL]]], s2_ = varVec[clu_i[nL]])
    # res<-res+ sum(llh_i,w_i)
  }
  # return(res)
  return(sum(llh_res))
}

fPDPM_mcmc<-function(obsDf, iteNum, burnIn, seed_=1, wf_,
                     K.max = NULL, KcutOff = NULL, K.tol=1e-1,
                     init.nClu=3, hyperAlpha=3, saveCov=T){
  set.seed(seed_)
  nSmp<-dim(obsDf)[3]
  nLoc<-dim(obsDf)[1:2]
  nJ<-log2(nLoc)
  D_type.vec<-c("LH","HL","HH")
  nL<-nJ[1]*length(D_type.vec)+1
  
  # rslnIdx<--1
  # for(j in (1:nJ[1])-1){
  #   for(d in 1:length(D_type.vec)){
  #     l <- j*length(D_type.vec)+d
  #     rslnIdx<-c(rslnIdx, rep(l,2^(2*j)))
  #   }
  # }
  # Psi<-GenPsi_Haar.2d(nLoc_ = nLoc, rslnIdx_ = rslnIdx, digit=5)
  
  # hyperParam init
  asg<-2.5
  bsg<-3
  ax <- nu1<-3
  bx <- nu2<-2
  a1 <- 2.1
  a2 <- 3.1
  if(is.null(K.max)){
    K.max<-floor(log(prod(nLoc)) * 3)
  }
  if(is.null(KcutOff)){
    KcutOff <- 1
  }
  LRD_hyper=c(asg, bsg, nu1, nu2, a1, a2, ax, bx, K.max, KcutOff, K.tol)
  names(LRD_hyper)=c("asg", "bsg", "nu1", "nu2", "a1", "a2", 
                     "ax", "bx", "Kmax", "KcutOff", "Ktol")
  
  # Initiation of the weights and clu of PDPM
  weight.ls<-stick.ls<-list()
  initClu<-clu.mat<-matrix(NA,nrow=nSmp, ncol=nL)
  for(l in 1:nL){
    j<-floor((l-1)/length(D_type.vec))
    d<-l-j*length(D_type.vec)
    D_type<-D_type.vec[d]
    stick.vec<-rbeta(init.nClu,1,hyperAlpha)
    weight.vec<-rep(NA,length(stick.vec))
    for(cIdx in 1:init.nClu){
      if(cIdx==1){
        weight.vec[cIdx]<-stick.vec[cIdx]
      }else{
        weight.vec[cIdx]<-stick.vec[cIdx] * prod(1-stick.vec[1:(cIdx-1)])
      }
    }
    stick.ls[[l]]<-stick.vec
    weight.ls[[l]]<-weight.vec
    initClu[,l]<-clu.mat[,l]<-sample(1:init.nClu,nSmp,replace = T)
    if(l!=nL){
      names(weight.ls)[l]<-names(stick.ls)[l]<-paste0(D_type,nJ[1]-j)
    }else{
      names(weight.ls)[nL]<-names(stick.ls)[nL]<-"Var"
    }
  }
  colnames(clu.mat)<-names(weight.ls)
  
  # Initiation of PDPM coef
  for(l in 1:nL){
    j<-floor((l-1)/length(D_type.vec))
    d<-l-j*length(D_type.vec)
    D_type<-D_type.vec[d]
    if(l==1){
      incpt<-rep(0,nSmp)
      coef.ls<-list()
      tau.ls<-list()
    }
    if(l!=nL){
      newCoef.ls<-spawn_newClu_coef(nClu = init.nClu, 2^(2*j))
      coef.ls[[l]]<- newCoef.ls$coef
      tau.ls[[l]]<-newCoef.ls$tau
      names(tau.ls)[l]<-names(coef.ls)[l]<-paste0(D_type,nJ[1]-j)
    }else{
      # Initiation of LRD
      newVar.vec<-spawn_newClu_diagVar(nClu=init.nClu)
      LRD<-spawn_newClu_LRD(nClu = init.nClu, K_=K.max, n_=nSmp, L_ = prod(nLoc), LRD_hyper = LRD_hyper)
      etaMat<-matrix(rnorm(nSmp*K.max),ncol=K.max)
      etaLs<- lapply(seq_len(nrow(etaMat)), function(i) matrix(etaMat[i,],nc=1))
      LRD$eta<-etaLs
    }
  }
  
  param.list<-list()
  param.list$coef<-coef.ls
  param.list$tau<-tau.ls
  param.list$incpt<-incpt
  param.list$var<-newVar.vec
  param.list$LRD<-LRD
  
  # storage
  LRDResList<-diagVarList<-coefResList<-list()
  cluRes.arr<-array(NA,dim=c(nSmp,nL,iteNum-burnIn))
  llh.vec<-rep(NA,iteNum)
  rm(list=c("coef.ls", "tau.ls", "incpt","nJ","nL", "K.max","KcutOff",
            "newCoef.ls", "D_type.vec", "D_type","newVar.vec",
            "LRD", "etaMat","etaLs"))
  gc()
  for(ite in 1:iteNum){
    if(ite %% 1==0){print (ite)}
    u.mat<-update_u(nSmp = nSmp, weightLs = weight.ls, cluMat = clu.mat)
    
    stk.tmp<-update_stick(cnctrnAlpha = hyperAlpha, uMat = u.mat,
                          stickLs = stick.ls, cluMat = clu.mat, weightLs = weight.ls)
    stick.ls<-stk.tmp$stick
    weight.ls<-stk.tmp$weight
    rm(stk.tmp)
    
    coef.tmp<-update_coef(obsDf=obsDf, cluMat = clu.mat, paramList = param.list, wf_=wf_)
    param.list$coef<-coef.tmp$coef
    param.list$tau<-coef.tmp$tau
    facMean <- coef.tmp$facMean
    rm(coef.tmp)
    
    var.tmp<-update_varDiag(obsDf=obsDf, cluMat = clu.mat, paramList = param.list,wf_=wf_,
                            facMean_ = facMean, a_s=asg, b_s=bsg)
    param.list$var<-var.tmp$var
    muMat<- var.tmp$muMat
    rm(var.tmp)
    
    LRD.tmp<-update_LRD(obsDf_ = obsDf, cluVec=clu.mat[,ncol(clu.mat)], muMat_ = muMat, 
                        LRD =param.list$LRD, diagVar_ = param.list$var, 
                        LRD_hyper = LRD_hyper, ite_ = ite)
    param.list$LRD<-LRD.tmp
    rm(LRD.tmp)
    
    clu.tmp<-update_clu(obsDf = obsDf, uMat = u.mat, wf_=wf_,
                        weightLs = weight.ls, stickLs = stick.ls,
                        cluMat = clu.mat, paramList = param.list,
                        cnctrnAlpha = hyperAlpha, LRD_hyper = LRD_hyper)
    clu.mat<-clu.tmp$clu
    param.list$coef<-clu.tmp$coef
    param.list$tau<-clu.tmp$tau
    param.list$var<-clu.tmp$var
    stick.ls<-clu.tmp$stick
    weight.ls<-clu.tmp$weight
    param.list$LRD <- clu.tmp$LRD
    rm(clu.tmp)

    llh.vec[ite]<-llh_cal(obsDf = obsDf, paramList = param.list, wf_=wf_,
                          cluMat = clu.mat, weightLs = weight.ls)
    if(ite>burnIn){
      # LRDResList[[ite - burnIn]]<-param.list$LRD
      cluRes.arr[,,(ite-burnIn)]<-clu.mat
      coefResList[[(ite-burnIn)]]<-param.list$coef
      diagVarList[[ite - burnIn]]<-param.list$var
      if(saveCov){
        LRDResList[[ite - burnIn]]<-param.list$LRD
      }
    }
  }
  return(list(initClu=initClu, llh=llh.vec, clu=cluRes.arr,
              coef=coefResList, LRD=LRDResList, var = diagVarList))
}

