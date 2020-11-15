\name{PreProcessAllPeptides}
\alias{PreProcessAllPeptides}
\title{Pre-processing PamChip Array data}


\description{ This is a wrapper function for \code{\link{PreProcessPam}}. And auotmates the preprocessing of list of peptides and Cell lines. 
This function removes peptide specific extreme observations, profiles and negative intensity measurements.  Negative intensity measurements within a profile are replaced by half of
the positive minimum intensity measurement of that profile. Extreme profiles have been removed via residuals based on one-way ANOVA (Response for the ANOVA is AUC: Average(residuals)
plus or minus 2 times SE(residuals)). Extreme observations within a profile were deleted if they differ more than plus or minus 2 times standard deviations of the profile-specific mean. 
}

\usage{
PreProcessAllPeptides(pep.names,PamSig,PathOutPut)
}

\arguments{
  \item{pep.names}{List of peptide names to be considered.}
  \item{PamSig }{A PamChipData data frame.}
  \item{PathOutPut}{Preprocessed RData files are created for each peptide and  by default save in the current working directory. Otherwise path to which files can be saved should be given. }
}

\value{A list of objects of class \code{\link{CleanedPamData}}

}


\author{Pushpike Thilakarathne, Ziv Shkedy and Dan Lin}


\seealso{\code{\link{PreProcessPam}}}



\examples{
data(PamChipData)

#-- Peptides Names --
PeptidesNames<-setdiff(colnames(PamChipData),
c("ID", "ResState", "ArrayNum", "CellName", "TreatName","Time"))


\donttest{
#--preprocess five peptides ----
temp<-PreProcessAllPeptides(pep.names=PeptidesNames[1:5],PamSig=PamChipData)

show(temp)
#-------- PreProcessed List of Peptides --------
#Class         : CleanedPamData 
#No. Peptides  : 5 
#No. Cell lines: 20 
#Columns names for a certain peptide:
# ID ResState ArrayNum CellName TreatName Time Pep1 
#Unique Time points: 0 5 10 15 20 25 30 35 40 45 50 55 60 
#PreProcessed Data frames are saved as RData Objects in : C:/Users/

#--preprocess three peptides only for responsive cell lines ----
temp2<-PreProcessAllPeptides(pep.names=PeptidesNames[1:3],
PamSig=PamChipData[PamChipData[,c("ResState")]=="R",])


#--preprocess three peptides only for NON-responsive cell lines ----
temp3<-PreProcessAllPeptides(pep.names=PeptidesNames[1:3],
PamSig=PamChipData[PamChipData[,c("ResState")]=="NR",])


    }
}

\keyword{anova, datasets}
