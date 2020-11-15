\name{ObservedRepProfiles}
\alias{ObservedRepProfiles}
\title{Observed replicate specific profiles }


\description{
Observed replicate specific profiles for a given cell line and peptide. 
}

\usage{
ObservedRepProfiles(pepname,DataMat,cells,Trt.Group=c("Treatment","Control"),
Res="R",log.true=FALSE,cell.mean=T)
}


\arguments{
  \item{pepname}{Name of the peptide to be used in visualization.}
  \item{DataMat}{A data frame of PamChip.} 
  \item{cells}{Name of the cell that needs to be visualized.} 
  \item{Trt.Group}{Treatment group for which visualization is required. Default is for both "Treatment" and "Control" groups. } 
  \item{Res}{Indicator for Responsive (R) or non-responsive (NR) group.} 
   \item{log.true}{Default \code{FALSE} and no log2 transformation is performed for intensities.  }
   \item{cell.mean}{A boolean parameter and default is \code{TRUE}. If \code{TRUE} cell line specific mean curve is drawn.}
}


\value{
none.
}


\references{
#add reference..
}

\author{Pushpike Thilakarathne, Ziv Shkedy and Dan Lin}


\seealso{\code{\link{plot}},  \code{\link{ObservedCellProfiles}}, \code{\link{ObservedGroupProfiles}}}


\examples{
  data(PamChipData)

\donttest{
    ObservedRepProfiles(pepname="Pep1",DataMat=PamChipData,
    cells="SKBR3",Trt.Group=c("Treatment","Control"),Res="R",log.true=FALSE)
    
    ObservedRepProfiles(pepname="Pep3",DataMat=PamChipData,
    cells="SKBR3",Trt.Group=c("Treatment","Control"),Res="R",log.true=FALSE)
    
    ObservedRepProfiles(pepname="Pep130",DataMat=PamChipData,
    cells="SUM159",Trt.Group="Treatment",Res="R",log.true=FALSE)
    
    ObservedRepProfiles(pepname="Pep130",DataMat=PamChipData,
    cells="SUM159",Trt.Group="Control",Res="R",log.true=FALSE)
    
}

}

\keyword{Plot}
