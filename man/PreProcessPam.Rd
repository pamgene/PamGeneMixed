\name{PreProcessPam}
\alias{PreProcessPam}
\title{Pre-processing Cell Line specific data for a given Peptide}


\description{
This function removes Cell Line specific extreme observations, profiles and negative intensity measurements of a given peptide.  Negative intensity measurements within a profile are replaced by half of
the positive minimum intensity measurement of that profile. Extreme profiles have been removed via residuals based on one-way ANOVA (Response for the ANOVA is AUC: Average(residuals)
plus or minus 2 times SE(residuals)). Extreme observations within a profile were deleted if they differ more than plus or minus 2 times standard deviations of the profile-specific mean. 
}

\usage{
PreProcessPam(p=108, cel =4, d=2, PamS, plotting=TRUE )
}


\arguments{
  \item{p}{Peptide to be considered. And should be supplied as numeric value that corresponds to the peptide index.}
  \item{cel}{Cell Line index.}
  \item{d}{Flexibility for the extreme observations to be removed. Eg.  2 or 3 times standard deviations.  Default is 2.}
  \item{PamS}{A PamChipData data frame.}
  \item{plotting}{A Boolean parameter. Default is \code{TRUE} and observed and filtered data are visualized for a given cell line of a peptide.}
}

\value{
A list of objects of class \code{\link{CleanedPeptide}}.
}


\author{Pushpike Thilakarathne, Ziv Shkedy and Dan Lin}



\seealso{\code{\link{PreProcessAllPeptides}}}







\examples{
  data(PamChipData)

\donttest{
#-----------
PreProcessPam(p   =20, cel =4,   d   =2, PamS=PamChipData, plotting=T )

#-----------
PreProcessPam(p   =6, cel =4,   d   =2, PamS=PamChipData, plotting=T )

#-----------
PreProcessPam(p   =78, cel =14,   d   =2, PamS=PamChipData, plotting=T)

#-----------
PreProcessPam(p   =143, cel =2,   d   =2, 
PamS=PamChipData[PamChipData[,c("ResState")]=="R",],  plotting=T )

#----------- Sequentially visualize the preprocessing steps -----------
for (i in 1:148){
    for (cc in 1:20){
    par(ask=T)
    PreProcessPam(p   =i, cel =cc,   d   =2, PamS=PamChipData, plotting=T )
    }
}

    }
}

\keyword{anova, datasets}
