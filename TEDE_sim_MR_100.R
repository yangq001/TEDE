
#Simulations in the MR context (independent SNPs, p=100)

library("mr.raps") 
library("jointsum")
library("MASS")
library("MendelianRandomization")
library("bindata")
library("aSPU")
library("data.table")
library("PMR")
library("TEDE")

####### SIMULATIONS ########



  
#################################
  
p=100 #number of SNPs (30 | 100)
rho=0 #correlation = 0 for MR
maf=0.3
MAF=rep(maf,p) #MAFs of SNPs
  
#Define the correlation and covariance matrices
corsnps=matrix(nrow=p,ncol=p)
for(i in 1:p){
  for(j in 1:p){
    corsnps[i,j]=rho^(abs(i-j))
  }
}
varsnps=diag(MAF*(1-MAF))
covsnps=sqrt(varsnps)%*%corsnps%*%sqrt(varsnps)

n1=10000 #sample size for the first dataset (G, X)
n1=50000 #sample size for the first dataset (G, X)
n2=10000 #sample size for the second dataset (G, Y)
n=n1+n2 #total sample size

iter=1000 #number of replications

size_XY=0 #X -> Y effect (0 | 0.2)

scenario = 1 #(1 | 2 | 3)

results=matrix(nrow=7,ncol=4)   #each column corresponds to one value of %invalid
colnames(results)=c("0","10%","30%","50%")
row.names(results)=c("Q","egger","score","aSPU","score2","aSPU2","PMR")

for(pinvalid in c(0, 0.1, 0.3, 0.5)){           #proportion of IVs that are invalid
  if(pinvalid == 0 & scenario > 2){
    next                   #when pinvalid = 0, it is the same for all scenarios
  }
  
  p_Q=p_egger=p_PMR=c()
  p_aspu=p_aspu2=p_score=p_score2=c()

  it=1
  for(it in 1:iter){
    set.seed(it)
    
    #Generate G->X effect sizes
    trunk=0.08                   
    size_GX=rnorm(p*90,sd=0.15) 
    size_GX=size_GX[size_GX>trunk | size_GX< -trunk]
    size_GX=size_GX[1:p]    
    
    #Generate G->Y direct effect sizes and G->U effect sizes
    if(pinvalid == 0){
      size_GY=rep(0,p)
      size_GU=rep(0,p)
    }
    if(pinvalid != 0 & scenario==1){
      size_GY=rnorm(p)*sqrt(0.075)
      size_GU=rep(0,p)
    }
    if(pinvalid != 0 & scenario==2){
      size_GY=0.1+rnorm(p)*sqrt(0.025)
      size_GU=rep(0,p)
    }
    if(pinvalid != 0 & scenario==3){
      size_GY=0.1+rnorm(p)*sqrt(0.025)
      size_GU=0.1*runif(p)
    }
    size_GY=size_GY*sign(size_GX)
    
    #set valid IV's size_GU and size_GY to 0
    nvalid=round(p*(1-pinvalid))
    vav=sample(1:p,nvalid)
    size_GY[vav]=0
    size_GU[vav]=0
    
    #Generate G and U
    G=matrix(rbinom(n*p,2,prob=maf),ncol=p)
    U=G%*%size_GU+rnorm(nrow(G))
    LDcov=cov(G)
    LDcor=cov2cor(LDcov)

    #Generate X (effect sizes are modified to control proportion of X explained by G)
    explained=0.2  #proportion of X explained by G
    if(explained>0){
    sdga=sd(G%*%size_GX)
    sii=sqrt(2*explained/(1-explained))/sdga
    size_GX=size_GX*sii
    }
    X=G%*%size_GX + U + rnorm(nrow(G))
    var(G%*%size_GX)/var(X)    
    
    #Generate Y (effect sizes are modified to control proportion of Y explained by X and G)
    bexp=0.02   #proportion of Y explained by X
    exa=0.003   #proportion of Y explained by G's direct effects
    if(sum(abs(size_GY))>0){
    sdga=sd(G%*%size_GY)
    sii2=sqrt(2*exa/(1-exa))/sdga
    size_GY=size_GY*sii2
    }
    if(size_XY!=0){
    tua=sd(size_XY*X)
    sii3=sqrt(2*bexp/(1-bexp))/tua
    size_XY2=size_XY*sii3   #resize size_XY
    }else{
      size_XY2=0
    }
    Y=G%*%size_GY + size_XY2*X  + U + rnorm(nrow(G)) 
    
    #Split the data into two datasets and center variables at 0
    part1=1:n1
    part2=(n1+1):nrow(G)
    G1=G[part1,]
    G2=G[part2,]
    G1=t(t(G1)-colMeans(G1))
    G2=t(t(G2)-colMeans(G2))
    X1=X[part1]-mean(X[part1])
    Y2=Y[part2]-mean(Y[part2])
    
    #Get MARGINAL summary statistics (effect size, SE) for Y~SNP (r, syk) and X~SNP (w, sxk)
    w=r=syk=sxk=c()

    gg=colSums(G2^2)
    r=c(t(G2)%*%Y2)/c(gg)
    Y2_pred=t(t(G2)*r)
    Y2_res=-(Y2_pred-Y2)
    sr=sqrt(colSums(Y2_res^2)/(length(part2)-2))
    syk=sr/sqrt(gg)
    
    gg=colSums(G1^2)
    w=c(t(G1)%*%X1)/c(gg)
    X1_pred=t(t(G1)*w)
    X1_res=-(X1_pred-X1)
    sr=sqrt(colSums(X1_res^2)/(length(part1)-2))
    sxk=sr/sqrt(gg)
    
    #Analysis using summary statistics

    #MR-Egger
    mr_object<-mr_input(bx = w, bxse = sxk, by = r, byse = syk)
    su1=mr_egger(mr_object)
    p_egger[it]=su1$Pvalue.Int
    
    #TEDE (building joint models from marginal is included in the function)
    su3=TEDE(w,sxk,r,syk,n1=n1,n2=n2,LDcov=LDcov,correlated_snps=FALSE,method="aSPU")
    su4=TEDE(w,sxk,r,syk,n1=n1,n2=n2,LDcov=LDcov,correlated_snps=FALSE,method="score")
    p_aspu[it]=su3[1]
    p_aspu2[it]=su3[2]
    p_score[it]=su4[1]
    p_score2[it]=su4[2]
    
    #Q
    p_Q[it]=su3[3]
    
    #PMR
    su2=PMR_summary_Egger(w/sxk,r/syk,LDcor,LDcor,n1,n2)
    p_PMR[it]=su2$pleiotropy_pvalue
    
    if(it%%10==0){
      print(it)  
    }
  }
  
  #show results for one setting
  data.frame(p,rho,pinvalid,size_XY,scenario)
  ev=cbind(p_Q,p_egger,p_score,p_aspu,p_score2,p_aspu2,p_PMR)
  print(colMeans(ev<0.05))
  
  bexp;n1;n2
  
  coll=which(c(0, 0.1, 0.3, 0.5)==pinvalid)
  results[,coll]=colMeans(ev<0.05)
}
 
results


#save results
setwd("/Users/dengy/Dropbox/Data Copy/TEDE")  
nam=paste("TEDE_MR_sc",scenario,"_p",p,"_rho",rho,"_beta",size_XY,"_nx",n1,"_ny",n2,"_hxy",bexp,".csv",sep="")
write.csv(results,file=nam)
  
  


  