
TEDE=function(effect_GX,se_GX=c(),effect_GY,se_GY,GX_joint=FALSE,GY_joint=FALSE,n1=1,n2=1,LDcov,correlated_snps=TRUE,method="aSPU",distribution_based=FALSE,n.perm=1000){
  if(length(se_GY)==0){
    stop("se_GY is required")
  }

  if(method == "aSPU"){
    result=TEDE_aSPU(effect_GX=effect_GX,se_GX=se_GX,effect_GY=effect_GY,se_GY=se_GY,n1=n1,n2=n2,LDcov=LDcov,correlated_snps=correlated_snps,n.perm=n.perm,distribution_based=distribution_based,GX_joint=GX_joint,GY_joint=GY_joint)
    return(result)
  }else{
    result=TEDE_Sc(effect_GX=effect_GX,se_GX=se_GX,effect_GY=effect_GY,se_GY=se_GY,n1=n1,n2=n2,LDcov=LDcov,correlated_snps=correlated_snps,GX_joint=GX_joint,GY_joint=GY_joint)
    return(result)
  }
}





