\name{JointSum}
\alias{JointSum}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
JointSum
}
\description{
Use summary statistics to build a joint linear model for one trait vs. L SNPs or one trait vs. L SNPs + (K-1) traits.
}
\usage{
JointSum(B1,S1,B2=0,S2=0,N,XX=diag(1,nrow=1),YY0,adj_Y=1,lam=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{B1}{
    a L*1 matrix containing marginal effects on the trait treated as response.
  }
  \item{S1}{
    a L*1 matrix containing standard errors for B1.
  }
  \item{B2}{
    a L*(K-1) matrix containing marginal effects on the (K-1) traits to adjust for. If K=1, do not specify this.
  }
  \item{S2}{
    a L*(K-1) matrix containing standard errors for B2. If K=1, do not specify this.
  }
  \item{N}{
    a L*K matrix containing sample sizes for each coefficients in B1, B2.
  }
  \item{XX}{
    a L*L estimated covariance matrix for the L SNPs.
  }
  \item{YY0}{
    a K*K estimated correlation matrix for the K traits.
  }
  \item{adj_Y}{
    whether traits should be adjusted for. If it is 0, adjust for SNPs only. Otherwise adjust for both SNPs and traits.
  }
  \item{lam}{
    a modifying parameter in [0,1). It is used only if adj_Y=1.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{beta}{
    coefficient estimates (SNPs first).
  }
  \item{cov}{
    the covariance matrix for coefficients.
  }
  \item{pvalue}{
    p-values for coefficients.
  }
  \item{sigma2}{
    estimated mean squared error.
  }
}
\references{
Deng, Y., Pan, W. (2017). Conditional analysis of multiple quantitative traits based on marginal GWAS summary statistics. Genet Epidemiol. doi: 10.1002/gepi.22046.
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
#2 correlated SNPs, 2 traits
set.seed(13)
x1=rbinom(1000,1,0.3)
x2=c(x1[1:300],rbinom(1000-300,1,0.3))
y2=rnorm(1000)+x1
y1=rnorm(1000)+y2/2

#standardization
x1=x1-mean(x1)
x2=x2-mean(x2)
y1=y1-mean(y1)
y2=y2-mean(y2)

#summary statistics
a=summary(lm(y1~x1-1))$coefficients
b=summary(lm(y1~x2-1))$coefficients
c=summary(lm(y2~x1-1))$coefficients
d=summary(lm(y2~x2-1))$coefficients

B1=as.matrix(c(a[1],b[1]))
S1=as.matrix(c(a[2],b[2]))
B2=as.matrix(c(c[1],d[1]))
S2=as.matrix(c(c[2],d[2]))
N=matrix(1000,nrow=2,ncol=2)
XX=cov(cbind(x1,x2))
YY0=cor(cbind(y1,y2))

#model Y1 ~ X1 + X2
JointSum(B1,S1,B2,S2,N,XX,YY0,adj_Y=0)
#or
JointSum(B1,S1,N=N[,1],XX=XX,YY0=diag(1),adj_Y=0)

#model Y1 ~ X1 + X2 + Y2
JointSum(B1,S1,B2,S2,N,XX,YY0,adj_Y=1)

#may compare with joint models using individual level data
summary(lm(y1~x1+x2-1))
summary(lm(y1~x1+x2+y2-1))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
