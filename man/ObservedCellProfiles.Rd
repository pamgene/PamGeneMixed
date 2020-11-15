\name{ObservedCellProfiles}
\alias{ObservedCellProfiles}
\title{Observed Cell Line Specific Profiles}


\description{
This function visualizes the observed Cell line specific profiles for a given peptide of the PamChip microarray data.
}

\usage{
ObservedCellProfiles(pepname,DataMat,
Trt.Group=c("Treatment","Control"),Res,log.true=FALSE)
}


\arguments{
  \item{pepname}{Name of the peptide to be used in visualization.}
  \item{DataMat}{A data frame of PamChip.} 
  \item{Trt.Group}{Treatment group for which visualization is required. Default is for both "Treatment" and "Control" groups. } 
  \item{Res}{Indicator for Responsive (R) or non-responsive (NR) group.} 
  \item{log.true}{Default \code{FALSE} and no log2 transformation is performed for intensities.  }
}



\value{
none
}


\author{Pushpike Thilakarathne, Ziv Shkedy and Dan Lin}

\seealso{\code{\link{ObservedRepProfiles}}, \code{\link{ObservedGroupProfiles}}}


\examples{
#--------

\donttest{

for (i in 1:148) {
    par(ask=TRUE)
    ObservedCellProfiles(pepname=paste("Pep",i,sep=""),DataMat=PamChipData,
    Trt.Group=c("Treatment","Control"),Res=c("R","NR"),log.true=FALSE)
    cat(i)
}

ObservedCellProfiles(pepname="Pep58",DataMat=PamChipData,
Trt.Group=c("Treatment","Control"),Res=c("R","NR"),log.true=FALSE)

ObservedCellProfiles(pepname="Pep9",DataMat=PamChipData,
Trt.Group=c("Treatment","Control"),Res=c("R","NR"),log.true=FALSE)

#------------log2 transformed -----------------------------
ObservedCellProfiles(pepname="Pep9",DataMat=PamChipData,
Trt.Group=c("Treatment","Control"),Res=c("R","NR"),log.true=TRUE)

ObservedCellProfiles(pepname="Pep26",DataMat=PamChipData,
Trt.Group=c("Treatment","Control"),Res=c("R","NR"),log.true=TRUE)



 }
}

\keyword{plot}
