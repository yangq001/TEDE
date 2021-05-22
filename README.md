# TEDE
TEsting Direct Effects for MR (Mendelian Randomization) or TWAS (Transcriptome-Wide Association Studies)

Please intall the package (TEDE_0.19.tar.gz) and check the help document
?TEDE

Required packages: aSPU, jointsum

Usage:
TEDE(effect_GX,se_GX=c(),effect_GY,se_GY,GX_joint=FALSE,GY_joint=FALSE,n1,n2,LDcov,correlated_snps=TRUE,method="aSPU",distribution_based=FALSE,n.perm=1000)

effect_GX: a vector containing the effect sizes of p SNPs on X.<br />
se_GX: a vector containing the standard errors of effect_GX. Do not specify this if this information is unavailable, in which case TEDE-Sc2 and TEDE-aSPU2 will not be performed.
effect_GY: a vector containing the effect sizes of p SNPs on Y.
se_GY: a vector containing the standard errors of effect_GY.
GX_joint: FALSE: the G->X effects are based on marginal models; TRUE: the G->X effects are based on joint models.
GY_joint: FALSE: the G->Y effects are based on marginal models; TRUE: the G->Y effects are based on joint models.
n1: the sample size used to get effect_GX. No need to specify this if GX_joint is TRUE.
n2: the sample size used to get effect_GY. No need to specify this if GY_joint is TRUE.
LDcov: a covariance matrix of the p SNPs estimated from a reference panel.
correlated_snps: whether the SNPs are correlated. Usually MR uses uncorrelated SNPs, and TWAS uses correlated SNPs.
method: "aSPU": TEDE-aSPU (and TEDE-aSPU2); "score": TEDE-Sc (and TEDE-Sc2).
distribution_based: FALSE: apply the standard aSPU test with summary statistics; TRUE: apply the distribution-based aSPU test. No need to specify this if method is "score".
n.perm: the number of iterations for the aSPU test. No need to specify this if distribution_based is TRUE.

result: a table containing the p-values.
