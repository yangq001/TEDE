

TEDE_Sc=function(effect_GX,se_GX=c(),effect_GY,se_GY=c(),n1,n2,LDcov,correlated_snps=TRUE,GX_joint=FALSE,GY_joint=FALSE){
  w=effect_GX
  r=effect_GY
  sxk=se_GX
  syk=se_GY

  if(length(se_GX)==0){
    sxk=w
  }


  #LDA MR-Egger
  if(correlated_snps==FALSE){
    LDcov=diag(diag(LDcov))
  }
  #LDA MR-Egger


  ### build joint
  #G->X
  if(GX_joint==FALSE){
    ss1=JointSum(w,sxk,N=n1,XX=LDcov,YY0=diag(1),adj_Y=0)
    w2=ss1$beta[,1]
    sxk2=sqrt(diag(ss1$cov))
    rcx=ss1$cov
  }else{
    w2=w
    sxk2=sxk
    rcx=diag(sxk2^2)
  }
  #G->Y
  if(GY_joint==FALSE){
    ss2=JointSum(r,syk,N=n2,XX=LDcov,YY0=diag(1),adj_Y=0)
    r2=ss2$beta[,1]
    syk2=sqrt(diag(ss2$cov))
    rcy=ss2$cov
  }else{
    r2=r
    syk2=syk
    rcy=diag(syk2^2)
  }
  ###





  p=length(w2)

  #aSPU summary
  B1=B2=r
  S1=S2=syk
  N=n2
  lala=YYX5b(B1,S1,B2,S2,N,XX=LDcov,YY0=diag(1))
  GG=lala$GG
  GY=lala$GY
  YY=lala$YY

  XY=w2%*%GY
  XX=t(w2)%*%GG%*%w2

  beta_TWAS2=c(XY/XX)  #checked
  sigY2=c(sqrt((YY+beta_TWAS2^2*XX-2*beta_TWAS2*XY)/(N-1)))

  #sigY2=sigY2+beta_TWAS2^2*(sum(sxk2^2))   #wrong conservative

  UU2=GY-c(beta_TWAS2)*t(t(w2)%*%GG)
  UUnew2=UU2/sigY2^2   #correct scores
  #GX=t(G2)%*%Xstar
  GX=as.matrix(w*c(XX))
  cov_UUnew2=matrix(nrow=length(UUnew2),ncol=length(UUnew2))
  J=p
  for(j in 1:J){
    for(k in 1:J){
      cov_UUnew2[j,k]=GG[j,k]
    }
  }
  cov_UUnew2=cov_UUnew2/sigY2^2

  wal=t(UUnew2)%*%ginv(cov_UUnew2)%*%UUnew2 #summary
  #1-pchisq(wal,df=J-1)

  1-pchisq(wal,df=qr(cov_UUnew2)$rank)

  control=1
  if(control==1){
    cov_UU=GG/sigY2^2   #old
    qusi=GG%*%rcx%*%GG
    cov_UU2=cov_UU+beta_TWAS2^2*qusi/sigY2^4   #consider var(w)
    UUnew2=UU2/sigY2^2
    cov_UUnew2=cov_UU2
  }

  wal2=t(UUnew2)%*%ginv(cov_UUnew2)%*%UUnew2

  mow=c(1-pchisq(wal,df=qr(cov_UUnew2)$rank),1-pchisq(wal2,df=qr(cov_UUnew2)$rank))

  names(mow)=c("TEDE-Sc","TEDE-Sc2")

  if(correlated_snps==FALSE){
    pQ=CRQ(w,sxk,r,syk)
    mow=c(mow,pQ)
    names(mow)=c("TEDE-Sc","TEDE-Sc2","Q")
  }

  if(length(se_GX)==0){
    return(mow[1])
  }else{
    return(mow)
  }
}


TEDE_aSPU=function(effect_GX,se_GX,effect_GY,se_GY,n1,n2,LDcov,correlated_snps=TRUE,n.perm=1000,distribution_based=FALSE,GX_joint=FALSE,GY_joint=FALSE){
  w=effect_GX
  r=effect_GY
  sxk=se_GX
  syk=se_GY

  if(length(se_GX)==0){
    sxk=w
  }

  aSPU_trans=0

  #LDA MR-Egger
  if(correlated_snps==FALSE){
    LDcov=diag(diag(LDcov))
  }
  #LDA MR-Egger


  ### build joint
  #G->X
  if(GX_joint==FALSE){
    ss1=JointSum(w,sxk,N=n1,XX=LDcov,YY0=diag(1),adj_Y=0)
    w2=ss1$beta[,1]
    sxk2=sqrt(diag(ss1$cov))
    rcx=ss1$cov
  }else{
    w2=w
    sxk2=sxk
    rcx=diag(sxk2^2)
  }
  #G->Y
  if(GY_joint==FALSE){
    ss2=JointSum(r,syk,N=n2,XX=LDcov,YY0=diag(1),adj_Y=0)
    r2=ss2$beta[,1]
    syk2=sqrt(diag(ss2$cov))
    rcy=ss2$cov
  }else{
    r2=r
    syk2=syk
    rcy=diag(syk2^2)
  }
  ###

  p=length(w2)

  #aSPU summary
  B1=B2=r
  S1=S2=syk
  N=n2
  lala=YYX5b(B1,S1,B2,S2,N,XX=LDcov,YY0=diag(1))
  GG=lala$GG
  GY=lala$GY
  YY=lala$YY

  XY=w2%*%GY
  XX=t(w2)%*%GG%*%w2

  beta_TWAS2=c(XY/XX)  #checked
  sigY2=c(sqrt((YY+beta_TWAS2^2*XX-2*beta_TWAS2*XY)/(N-1)))

  #sigY2=sigY2+beta_TWAS2^2*(sum(sxk2^2))   #wrong conservative

  UU2=GY-c(beta_TWAS2)*t(t(w2)%*%GG)
  UUnew2=UU2/sigY2^2   #correct scores
  #GX=t(G2)%*%Xstar
  GX=as.matrix(w*c(XX))
  cov_UUnew2=matrix(nrow=length(UUnew2),ncol=length(UUnew2))
  J=p
  for(j in 1:J){
    for(k in 1:J){
      cov_UUnew2[j,k]=GG[j,k]
    }
  }
  cov_UUnew2=cov_UUnew2/sigY2^2

  wal=t(UUnew2)%*%ginv(cov_UUnew2)%*%UUnew2 #summary
  #1-pchisq(wal,df=J-1)

  1-pchisq(wal,df=qr(cov_UUnew2)$rank)

  aa=UUnew2/sqrt(diag(cov_UUnew2)) #standardize -> Z
  bb=cov2cor(cov_UUnew2)

  if(aSPU_trans==1){
    rua=chol(ginv(bb))
    aa=rua%*%aa
    bb=rua%*%bb%*%t(rua)
  }
  if(aSPU_trans==2){
    rua=ginv(bb)
    aa=rua%*%aa
    bb=rua%*%bb%*%t(rua)
  }

  if(distribution_based==FALSE){
    resultA0b=aSPUs(aa,bb,pow=c(1,2,4,8,Inf),n.perm=n.perm,prune=FALSE)
    ka=resultA0b$pvs[length(resultA0b$pvs)]
  }else{
    resultA0b=aSPUsD(aa,bb)
    ka=resultA0b[length(resultA0b)]
  }

  control=1
  if(control==1){
    cov_UU=GG/sigY2^2   #old
    qusi=GG%*%rcx%*%GG
    cov_UU2=cov_UU+beta_TWAS2^2*qusi/sigY2^4   #consider var(w)
    UUnew2=UU2/sigY2^2
    cov_UUnew2=cov_UU2
  }

  wal2=t(UUnew2)%*%ginv(cov_UUnew2)%*%UUnew2

  aa=UUnew2/sqrt(diag(cov_UUnew2)) #standardize -> Z
  bb=cov2cor(cov_UUnew2)

  if(aSPU_trans==1){
    rua=chol(ginv(bb))
    aa=rua%*%aa
    bb=rua%*%bb%*%t(rua)
  }
  if(aSPU_trans==2){
    rua=ginv(bb)
    aa=rua%*%aa
    bb=rua%*%bb%*%t(rua)
  }

  if(distribution_based==FALSE){
    resultA0c=aSPUs(aa,bb,pow=c(1,2,4,8,Inf),n.perm=n.perm,prune=FALSE)
    ka2=resultA0c$pvs[length(resultA0c$pvs)]
  }else{
    resultA0c=aSPUsD(aa,bb)
    ka2=resultA0c[length(resultA0c)]
  }

  #t(UUnew2/sd(UUnew2))%*%ginv(cov_UUnew2/sd(UUnew2)/sd(UUnew2))%*%(UUnew2/sd(UUnew2))  #no use


  #mow=c(resultA0b$pvs[length(resultA0b$pvs)],resultA0c$pvs[length(resultA0c$pvs)],1-pchisq(wal,df=J-1))
  mow=c(ka,ka2)

  names(mow)=c("TEDE-aSPU","TEDE-aSPU2")

  if(correlated_snps==FALSE){
    pQ=CRQ(w,sxk,r,syk)
    mow=c(mow,pQ)
    names(mow)=c("TEDE-aSPU","TEDE-aSPU2","Q")
  }

  if(length(se_GX)==0){
    return(mow[1])
  }else{
    return(mow)
  }
}




YYX5b=function(B1,S1,B2,S2,N,XX=diag(1,nrow=1),YY0){ #x first
  B1=as.matrix(B1)
  B2=as.matrix(B2)
  S1=as.matrix(S1)
  S2=as.matrix(S2)
  n=max(N)
  ny=ncol(B2)
  xx=XX[1,1]
  s1=S1[1]
  b1=B1[1]
  S21=S2[1,]
  B21=B2[1,]

  nrx=1
  #nrx=nrow(B1)
  N=t(N)


  for(rx in 1:nrx){
    ny=ncol(B2)
    xx=XX[rx,rx]
    s1=S1[rx]
    b1=B1[rx]
    S21=S2[rx,]
    B21=B2[rx,]

    if(rx>nrow(N)){
      rx2=1
    }else{
      rx2=rx
    }

    if(rx==1){
      yy1=N[1,rx2]/n*(n-1)*xx*s1^2+xx*b1^2
      YY2=N[-1,rx2]/n*(n-1)*xx*S21^2+xx*B21^2
    }
    else{
      yy1=yy1+N[1,rx2]/n*(n-1)*xx*s1^2+xx*b1^2
      YY2=YY2+N[-1,rx2]/n*(n-1)*xx*S21^2+xx*B21^2
    }
  }
  yy1=yy1/nrx
  YY2=YY2/nrx


  YYself=c(yy1,YY2)
  xxd=diag(XX)
  XY1=xxd*B1       #vector
  XY2=t(t(B2)*xxd)   #matrix
  nx=ncol(XX)
  A=matrix(nrow=nx,ncol=nx)
  A[1:nx,1:nx]=XX

  B=c(XY1)

  if(det(A)<1e-20){
    sA=MASS::ginv(A)
  }else{
    sA=MASS::ginv(A)
  }

  YY=YYself
  GY=XY1
  GG=XX

  beta=sA%*%B
  sigma2=(yy1-t(beta)%*%B)/(n-(nx))
  se=sigma2[1,1]*sA   #cov

  pvalue=c()
  for(i in 1:length(beta)){
    pvalue[i]=2*pnorm(abs(beta[i]/sqrt(se[i,i])),lower.tail=FALSE)
  }
  return(list(beta=beta,cov=se,pvalue=pvalue,sigma2=sigma2,YY=YY,GY=XY1,GG=XX))
}



prune_pv=function(LD,pv,chrs,cutoff=0.5){
  killed=c()
  for(i in 1:(nrow(LD)-1)){
    for(j in (i+1):nrow(LD)){
      if(i %in% killed){
        break
      }else{
        if(j %in% killed){
          next
        }else{
          haha=LD[i,j]*(chrs[i,1]==chrs[j,1])
          if(abs(haha) > cutoff){
            if(pv[i]>pv[j]){
              killed=c(killed,i)
            }else{
              killed=c(killed,j)
            }
          }
        }
      }
    }
  }
  return(killed)
}



exe1=function(a){
  return(a[1])
}
exe2=function(a){
  return(a[2])
}



CRQ=function(w,sxk,r,syk,rcx=diag(sxk)^2,rcy=diag(syk)^2){
  ra=r/w   #ratio/beta
  sra=syk/w      #SE of ratio

  S=diag(abs(sra))
  cov_ra=S%*%cov2cor(rcy)%*%S  #cov(ratios)

  vinv=ginv(cov_ra)

  l1=rep(1,length(w))
  wra=c(t(l1)%*%vinv/c((t(l1)%*%vinv%*%l1)))  #weights: equal to (w/syk)^2/sum((w/syk)^2) for independent SNPs
  beta0=c(t(wra)%*%ra)     #equal to IVW estimate for independent SNPs

  S2=t(t(diag(length(ra)))-wra)
  vinv2=ginv(S2%*%cov_ra%*%S2)
  Q2=t(ra-beta0)%*%vinv2%*%(ra-beta0)
  pvalue=1-pchisq(abs(Q2),df=qr(vinv)$rank)
  return(pvalue)
}




LDA_Egger=function(w,sxk,r,syk,n1,n2,LDcov,correlated_snps=TRUE,flip=TRUE,weights=FALSE){

  if(correlated_snps==FALSE){
    LDcov=diag(diag(LDcov)) #If we know for sure the SNPs are not correlated, we can set the LDcov matrix to be diagonal
  }

  ### build joint models
  #G->X
  if(weights==FALSE){
    ss1=JointSum(w,sxk,N=n1,XX=LDcov,YY0=diag(1),adj_Y=0)
    w2=ss1$beta[,1]
    sxk2=sqrt(diag(ss1$cov))
    rcx=ss1$cov
  }else{    #weights = TRUE means the "w" we input already represents joint effects, so there is no need to build a joint model
    w2=w
    sxk2=sxk
    rcx=diag(sxk2^2)
  }
  #G->Y
  ss2=JointSum(r,syk,N=n2,XX=LDcov,YY0=diag(1),adj_Y=0)
  r2=ss2$beta[,1]
  syk2=sqrt(diag(ss2$cov))
  rcy=ss2$cov
  ###

  if(flip==TRUE){   #flip the alleles to make G->X effects positive, as done in MR-Egger
    ho=sign(w2)
    ha=diag(sign(w2))
    r2=r2*ho
    w2=w2*ho
    rcx=ha%*%rcx%*%ha
    rcy=ha%*%rcy%*%ha
  }

  #estimate coefficients and covariance
  l=rep(1,length(w2))
  S_inv=ginv(rcy)   #Sigma^(-1)
  a1=t(l)%*%S_inv%*%l
  a2=t(l)%*%S_inv%*%w2
  a4=t(w2)%*%S_inv%*%w2
  SS=matrix(c(a1,a2,a2,a4),nrow=2)
  a5=c(t(l)%*%S_inv%*%r2,t(w2)%*%S_inv%*%r2)

  ab_cov=ginv(SS)
  ab_est=ab_cov%*%a5

  ma2=cbind(ab_est,sqrt(diag(ab_cov)),ab_est/sqrt(diag(ab_cov)),2*pnorm(-abs(ab_est/sqrt(diag(ab_cov)))))

  rownames(ma2)=c("Intercept","Slope")
  colnames(ma2)=c("Estimate","Std. Error","t","Pr(>|t|)")

  return(ma2)
}

