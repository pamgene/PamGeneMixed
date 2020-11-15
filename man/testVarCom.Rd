\name{testVarCom}
\alias{testVarCom}
\title{Wald Test to check variances are equal or not }


\description{
Testing whether the group specific variances are the same or not. \cr


}

\usage{
testVarCom(ObjPamFit)
}


\arguments{
  \item{ObjPamFit}{Output object from the \code{\link{PamGeneMix}}.} 
  }


\value{
  \item{A list}{A list with p-value(s) of the test. And estimated variances and their standard errors for each group. If there are no groups then NA will be returned.}

}


\references{
Engle, Robert F. (1983). "Wald, Likelihood Ratio, and Lagrange Multiplier Tests in Econometrics". In Intriligator, M. D.; and Griliches, Z.. Handbook of Econometrics. II. Elsevier. pp. 796-801. ISBN 978-0-444-86185-6.
}

\author{Pushpike Thilakarathne, Ziv Shkedy and Dan Lin}


\seealso{\code{\link{PamGeneMix}}}


\examples{
data(TestPepModelData)

\donttest{
PTx<-TestPepModelData


# log2 transform the response 
PTx$y<-log2(PTx[,ncol(PTx)]) 

# we now create interaction between ResState and Treatment as follows.

n.groups<-length(unique(PTx[,c("ResState")]))*length(unique(PTx[,c("TreatName")]))
PTx$ResTrt<-0
xhlp<-lm(y~-1+ResState:TreatName,data=PTx,x=TRUE)$x
for (i in 1:n.groups) PTx$ResTrt[xhlp[,i]==1]<-i
PTx$ResTrt<-as.factor(PTx$ResTrt)
levels(PTx$ResTrt)<-colnames(xhlp)

#unique cell lines
cellLines<-levels(PTx$CellName)
PTx$ResState<-as.factor(PTx$ResState)

# create interaction between Treatment and cell lines
PTx$CellLineResTrt<-0
xhlp<-lm(y~-1+CellName:TreatName,data=PTx,x=TRUE)$x
ncols<-ncol(xhlp)
for (i in 1:ncols) PTx$CellLineResTrt[xhlp[,i]==1]<-i
PTx$CellLineResTrt<-as.factor(PTx$CellLineResTrt)
levels(PTx$CellLineResTrt)<-colnames(xhlp)

M3gamm<-PamGeneMix(formula=y~-1+ResTrt+s(time,by=ResTrt,bs="tp",m=3),Weights=varIdent(form=~1|ResTrt) ,
PTx=PTx,Random.structure=list(ArrayNum=~-1,CellLineResTrt=~1+time+time2,ID=~1+time+time2),
Control.list=list(maxIter=200, msMaxIter=250 ,msMaxEval=1000,apVar=TRUE))

#we can now use testVarCom function to test the hypothesis that 
#H0: Var(group_R_i) - Var(group_R_j) = 0 and  Var(group_NR_i) - Var(group_NR_j) = 0. 
#That is testing for treatment effects conditional on the responsive statues 

testVarCom(M3gamm)  
}
}

\keyword{Wald,likelihood-ratio test}
