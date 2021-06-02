#conditional model (call)
library(MASS)

JointSum=function(B1,S1,B2=0,S2=0,N,XX=diag(1,nrow=1),YY0,adj_Y=1,lam=0){
  if(adj_Y==1){
    return(YYX4(B1,S1,B2,S2,N,XX,YY0,lam))
  }
  else{
    if(length(B2)<2){
    B2=B1
    S2=S1
    N2=cbind(N,N)
    YY0=diag(2)
    return(YYX5(B1,S1,B2,S2,N2,XX,YY0))
    }
    else{
      return(YYX5(B1,S1,B2,S2,N,XX,YY0))
    }
  }
}

#conditional model (ind)

YYX5=function(B1,S1,B2,S2,N,XX=diag(1,nrow=1),YY0){ #x first
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


  beta=sA%*%B
  sigma2=(yy1-t(beta)%*%B)/(n-(nx))
  se=sigma2[1,1]*sA   #cov

  pvalue=c()
  for(i in 1:length(beta)){
    pvalue[i]=2*pnorm(abs(beta[i]/sqrt(se[i,i])),lower.tail=FALSE)
  }
  return(list(beta=beta,cov=se,pvalue=pvalue,sigma2=sigma2))
}

#conditional model (sum)


YYX4=function(B1,S1,B2,S2,N,XX=diag(1,nrow=1),YY0,lam=0){ #x first
  B1=as.matrix(B1)
  B2=as.matrix(B2)
  S1=as.matrix(S1)
  S2=as.matrix(S2)
  n=max(N)
  
  N=t(N)

  nrx=1
  for(rx in 1:nrx){
    ny=ncol(B2)
    xx=XX[rx,rx]
    s1=S1[rx]
    b1=B1[rx]
    S21=S2[rx,]
    B21=B2[rx,]
    if(rx==1){
      yy1=N[1,rx]/n*(n-1)*xx*s1^2+xx*b1^2
      YY2=N[-1,rx]/n*(n-1)*xx*S21^2+xx*B21^2
    }
    else{
      yy1=yy1+N[1,rx]/n*(n-1)*xx*s1^2+xx*b1^2
      YY2=YY2+N[-1,rx]/n*(n-1)*xx*S21^2+xx*B21^2
    }
  }
  yy1=yy1/nrx
  YY2=YY2/nrx

  YYself=c(yy1,YY2)
  xxd=diag(XX)
  XY1=xxd*B1       #vector
  XY2=t(t(B2)*xxd)   #matrix

  YYcross=matrix(nrow=ny+1,ncol=ny+1)
  for(i1 in 1:(ny+1)){
    for(j1 in 1:(ny+1)){
      if(i1==j1){
        YYcross[i1,j1]=YYself[i1]
      }
      else{
        YYcross[i1,j1]=YY0[i1,j1]*sqrt(YYself[i1]*YYself[j1])
      }
    }
  }

  nx=ncol(XX)
  A=matrix(nrow=ny+nx,ncol=ny+nx)
  A[1:nx,1:nx]=XX
  A[(nx+1):(ny+nx),1:nx]=t(XY2)
  A[1:nx,(nx+1):(ny+nx)]=XY2

  A[(nx+1):(ny+nx),(nx+1):(ny+nx)]=YYcross[2:(ny+1),2:(ny+1)]

  B=c(XY1,YYcross[1,2:(ny+1)])

  O=matrix(0,ncol=(ny+nx),nrow=(ny+nx))
  for(i in 1:(ny+nx)){
    O[i,i]=1
  }
  A=(A+lam*O)/(1+lam)
  
  
  beta=MASS::ginv(A)%*%B
  sigma2=(yy1-t(beta)%*%B)/(n-(nx+ny))
  se=sigma2[1,1]*MASS::ginv(A)   #cov

  if(sum(diag(se)<=0)==0){
    pvalue=c()
    for(i in 1:length(beta)){
      pvalue[i]=2*pnorm(abs(beta[i]/sqrt(se[i,i])),lower.tail=FALSE)
    }
  }else{
    diag(se)=abs(diag(se))
    pvalue=beta*0-1
  }


  return(list(beta=beta,cov=se,pvalue=pvalue,sigma2=sigma2))
}
