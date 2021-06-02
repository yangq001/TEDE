\name{TEDE}
\alias{TEDE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
TEDE
}
\description{
TEsting Direct Effects (for MR or TWAS)
}
\usage{
TEDE(effect_GX,se_GX=c(),effect_GY,se_GY,GX_joint=FALSE,GY_joint=FALSE,n1,n2,LDcov,correlated_snps=TRUE,method="aSPU",distribution_based=FALSE,n.perm=1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{effect_GX}{
    a vector containing the effect sizes of p SNPs on X.
  }
  \item{se_GX}{
    a vector containing the standard errors of effect_GX. Do not specify this if this information is unavailable, in which case TEDE-Sc2 and TEDE-aSPU2 will not be performed.
  }
  \item{effect_GY}{
    a vector containing the effect sizes of p SNPs on Y.
  }
  \item{se_GY}{
    a vector containing the standard errors of effect_GY.
  }
  \item{GX_joint}{
    FALSE: the G->X effects are based on marginal models; TRUE: the G->X effects are based on joint models.
  }
  \item{GY_joint}{
    FALSE: the G->Y effects are based on marginal models; TRUE: the G->Y effects are based on joint models.
  }
  \item{n1}{
    the sample size used to get effect_GX. No need to specify this if GX_joint is TRUE.
  }
  \item{n2}{
    the sample size used to get effect_GY. No need to specify this if GY_joint is TRUE.
  }
  \item{LDcov}{
    a covariance matrix of the p SNPs estimated from a reference panel.
  }
  \item{correlated_snps}{
    whether the SNPs are correlated. Usually MR uses uncorrelated SNPs, and TWAS uses correlated SNPs.
  }
  \item{method}{
    "aSPU": TEDE-aSPU (and TEDE-aSPU2); "score": TEDE-Sc (and TEDE-Sc2).
  }
  \item{distribution_based}{
    FALSE: apply the standard aSPU test with summary statistics; TRUE: apply the distribution-based aSPU test. No need to specify this if method is "score".
  }
  \item{n.perm}{
    the number of iterations for the aSPU test. No need to specify this if distribution_based is TRUE.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{result}{
    a table containing the p-values.
  }
}
\references{
%Deng, Y., Pan, W. (2017). Conditional analysis of multiple quantitative traits based on marginal GWAS summary statistics. Genet Epidemiol. doi: 10.1002/gepi.22046.
}
\author{
Yangqing Deng and Wei Pan.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library(aSPU)
library(TEDE)

##### Testing horizontal pleiotropy for MR
set.seed(1)
p=20
n1=n2=100
effect_GX=rnorm(p)
se_GX=rep(1,p)
effect_GY=rnorm(p)
se_GY=rep(1,p)
LDcov=diag(p)

#TEDE-Sc and TEDE-Sc2
TEDE(effect_GX=effect_GX,se_GX=se_GX,effect_GY=effect_GY,se_GY=se_GY,n1=n1,n2=n2,LDcov=LDcov,correlated_snps=FALSE,method="score")

#TEDE-aSPU and TEDE-aSPU2
TEDE(effect_GX=effect_GX,se_GX=se_GX,effect_GY=effect_GY,se_GY=se_GY,n1=n1,n2=n2,LDcov=LDcov,correlated_snps=FALSE,method="aSPU",distribution_based=TRUE)


##### Testing horizontal pleiotropy for TWAS
set.seed(1)
p=20
n1=n2=100
A=diag(p)
A[1,2]=A[2,1]=0.3
effect_GX=c(rmvnorm(n=1,sigma=A))
se_GX=rep(1,p)
effect_GY=c(rmvnorm(n=1,sigma=A))
se_GY=rep(1,p)
LDcov=A

#TEDE-Sc and TEDE-Sc2
TEDE(effect_GX=effect_GX,se_GX=se_GX,effect_GY=effect_GY,se_GY=se_GY,n1=n1,n2=n2,LDcov=LDcov,correlated_snps=TRUE,method="score")

#TEDE-aSPU and TEDE-aSPU2
TEDE(effect_GX=effect_GX,se_GX=se_GX,effect_GY=effect_GY,se_GY=se_GY,n1=n1,n2=n2,LDcov=LDcov,correlated_snps=TRUE,method="aSPU",distribution_based=TRUE)


##### Testing horizontal pleiotropy for TWAS using weights
set.seed(1)
p=20
n1=n2=100
A=diag(p)
A[1,2]=A[2,1]=0.3
effect_GY=c(rmvnorm(n=1,sigma=A))
se_GY=rep(1,p)
LDcov=A

# We assume effect_GX are weights from previous studies based on joint models with variable selection (e.g. eNet) and se(effect_GX) is not available
effect_GX=c(rmvnorm(n=1,sigma=ginv(A)))

#TEDE-Sc only (TEDE-Sc2 cannot be applied since se(effect_GX) is not available)
TEDE(effect_GX=effect_GX,effect_GY=effect_GY,se_GY=se_GY,GX_joint=TRUE,n1=n1,n2=n2,LDcov=LDcov,correlated_snps=TRUE,method="score")

#TEDE-aSPU only (TEDE-aSPU2 cannot be applied since se(effect_GX) is not available)
TEDE(effect_GX=effect_GX,effect_GY=effect_GY,se_GY=se_GY,GX_joint=TRUE,n1=n1,n2=n2,LDcov=LDcov,correlated_snps=TRUE,method="aSPU",distribution_based=TRUE)


}
