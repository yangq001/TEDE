# TEDE
**TEsting Direct Effects** for MR (Mendelian Randomization) or TWAS (Transcriptome-Wide Association Studies)<br />
Authors: Yangqing Deng and Wei Pan

Please intall the package (TEDE_0.31.tar.gz) and check the help document
?TEDE

Required packages: aSPU, jointsum


<br /><br />

Usage:
TEDE(effect_GX,se_GX=c(),effect_GY,se_GY,GX_joint=FALSE,GY_joint=FALSE,n1,n2,LDcov,correlated_snps=TRUE,method="aSPU",distribution_based=FALSE,n.perm=1000)

Input:<br />
**effect_GX**: a vector containing the effect sizes of p SNPs on X.<br />
**se_GX**: a vector containing the standard errors of effect_GX. Do not specify this if this information is unavailable, in which case TEDE-Sc2 and TEDE-aSPU2 will not be performed.<br />
**effect_GY**: a vector containing the effect sizes of p SNPs on Y.<br />
**se_GY**: a vector containing the standard errors of effect_GY.<br />
**GX_joint**: FALSE: the G->X effects are based on marginal models; TRUE: the G->X effects are based on joint models.<br />
**GY_joint**: FALSE: the G->Y effects are based on marginal models; TRUE: the G->Y effects are based on joint models.<br />
**n1**: the sample size used to get effect_GX. No need to specify this if GX_joint is TRUE.<br />
**n2**: the sample size used to get effect_GY. No need to specify this if GY_joint is TRUE.<br />
**LDcov**: a covariance matrix of the p SNPs estimated from a reference panel.<br />
**correlated_snps**: whether the SNPs are correlated. Usually MR uses uncorrelated SNPs, and TWAS uses correlated SNPs.<br />
**method**: "aSPU": TEDE-aSPU (and TEDE-aSPU2); "score": TEDE-Sc (and TEDE-Sc2).<br />
**distribution_based**: FALSE: apply the standard aSPU test with summary statistics; TRUE: apply the distribution-based aSPU test. No need to specify this if method is "score".<br />
**n.perm**: the number of iterations for the aSPU test. No need to specify this if distribution_based is TRUE.

Output:<br />
**result**: a table containing the p-values.
