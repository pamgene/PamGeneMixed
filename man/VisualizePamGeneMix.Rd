\name{VisualizePamGeneMix}
\alias{VisualizePamGeneMix}
\title{Visualize PamGeneMix Output Objects}


\description{
This function visualizes following plots using the fitted PamGeneMix object.




}

\usage{
VisualizePamGeneMix(PamObject,
plot.type=c("velocity","velocityCI","smooth.fit","cellline","replicate","SubSpeVelocity"),
name.cell=NULL)
}


\arguments{
  \item{PamObject}{Object out from the \code{\link{PamGeneMix}}.}
  \item{plot.type}{Required plot types: 
  
   \enumerate{
        \item velocity: Produce first order derivatives for each group in single plot.
        \item velocityCI: Produce first order derivative along with the 95\% confidence intervals for each group separately.
        \item smooth.fit: Produce smoothing fit for each group severalty with observed profiles as background.
        \item cellline: Cell line specific fitted smoothing curves for each group separately.
        \item replicate: Replicate specific smoothing fit for a given cell line. Cell line name should be given.
        \item SubSpeVelocity: Replicate specific first order derivatives of the fitted smoothing curves for each group separately.
   }
    
  } 
  \item{name.cell}{Name of the cell line that needs to be visualized and this is only applicable when replicate specific smoothing fit needs to be visualized.} 
}



\value{
No value returns
}


\author{Pushpike Thilakarathne, Ziv Shkedy and Dan Lin}


\seealso{\code{\link{ObservedRepProfiles}}, \code{\link{ObservedGroupProfiles}}, \code{\link{ObservedCellProfiles}}}


\examples{
#--------

\donttest{
data(TestPepModelData)
head(TestPepModelData)

PTx<-TestPepModelData

# log2 transform the response 
PTx$y<-log2(PTx[,ncol(PTx)]) 

PTx<-PTx[PTx[,c("ResState")]=="R",]

n.groups<-length(unique(PTx[,c("ResState")]))*length(unique(PTx[,c("TreatName")]))
#n.groups

xhlp<-lm(y~-1+TreatName,data=PTx,x=TRUE)$x
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

#--------------- the model -----------------
M1gamm<-PamGeneMix(formula=y~-1+ResTrt+s(time,by=ResTrt,bs="tp",m=3),
        PTx=PTx,Random.structure=list(ArrayNum=~-1,CellLineResTrt=~1+time+time2,ID=~1+time+time2),
        Control.list=list(maxIter=200, msMaxIter=250 ,msMaxEval=1000,apVar=TRUE))
        
#----------------------------------------------------------------------        
        
VisualizePamGeneMix(M1gamm,plot.type="smooth.fit")

VisualizePamGeneMix(M1gamm,plot.type="cellline")

VisualizePamGeneMix(M1gamm,plot.type="velocityCI")

VisualizePamGeneMix(M1gamm,plot.type="replicate",name.cell= "A431" )

VisualizePamGeneMix(M1gamm,plot.type="replicate",name.cell="DU145" )

VisualizePamGeneMix(M1gamm, plot.type="SubSpeVelocity")

 }
}

\keyword{plot,smoothing, non-linear, velocity, confidence intervals}
