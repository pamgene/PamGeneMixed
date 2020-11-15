\name{ObservedGroupProfiles}
\alias{ObservedGroupProfiles}
\title{Observed Group Specific Profiles}


\description{
This function estimate and visualize observed average profiles for each group for a given peptide.   
}

\usage{
ObservedGroupProfiles(pepname="Pep5",DataMat,
Trt.Group=c("Treatment","Control"),Res="R",log.true=FALSE)
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

\seealso{\code{\link{ObservedRepProfiles}}, \code{\link{ObservedCellProfiles}}}



\examples{
  data(PamChipData)


  \donttest{  
ObservedGroupProfiles(pepname="Pep142",DataMat=PamChipData,
Trt.Group=c("Treatment","Control"),Res=c("R","NR"),log.true=FALSE)

ObservedGroupProfiles(pepname="Pep142",DataMat=PamChipData,
Trt.Group=c("Treatment","Control"),Res=c("R","NR"),log.true=TRUE)

ObservedGroupProfiles(pepname="Pep108",DataMat=PamChipData,
Trt.Group=c("Treatment","Control"),Res=c("R","NR"),log.true=FALSE)

ObservedGroupProfiles(pepname="Pep108",DataMat=PamChipData,
Trt.Group=c("Treatment","Control"),Res=c("R","NR"),log.true=TRUE)

  }

}

\keyword{plot}
