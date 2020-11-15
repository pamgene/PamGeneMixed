\name{VolcanoPam}
\alias{VolcanoPam}
\title{Visualize  Output Object from AutoPamGeneMix}


\description{
This function visualizes  Volcano plots of p-values returns from the   \code{\link{AutoPamGeneMix}}.  
It also returns Top K Peptides according to the test at t=0, t=end, t=any time point, and t=entire profile. Corresponding p-values are also returned.\cr


}

\usage{
VolcanoPam(ObjectAutoPam,Topp=5,plotting=TRUE)
}


\arguments{
  \item{ObjectAutoPam}{Object out from the \code{\link{AutoPamGeneMix}}.}
  \item{Topp}{Number of top K peptides to be considered. Default is five.}
  \item{plotting}{A boolean parameter indicating whether the plots to be shown.}
}



\value{
  \item{PP}{A data frame with p-values.}
    \item{TopList}{Top K peptides at time points where the tests are performed.}

}


\author{Pushpike Thilakarathne, Ziv Shkedy and Dan Lin}


\seealso{\code{\link{AutoPamGeneMix}}}


\examples{
#--------

\donttest{
    data(PamChipData)
    
    PeptidesNames<-setdiff(colnames(PamChipData),
        c("ID", "ResState", "ArrayNum", "CellName", "TreatName","Time"))


 #------- Four groups in the data ---
    # Apply preprocessing steps for each peptide
    temp2<-PreProcessAllPeptides(pep.names=PeptidesNames,PamSig=PamChipData,PathOutPut="C:/Temp")
    
    #-----------------------------------
    #  Fit the model for each peptide : This will take a while !!
    ResultsAll<-AutoPamGeneMix(formula=y~-1+ResTrt+s(time,by=ResTrt,bs="tp",m=3),
    Weights=varIdent(form=~1|ResTrt),
    Random.structure=list(ArrayNum=~1,CellLineResTrt=~1+time+time2,ID=~1+time+time2),
    temp2,TestAt=30,PathOutPut="C:/FittedGamm")
    
    show(ResultsAll)
    
    
    #---------------------------------
    VolcanoPam(ResultsAll,Topp=10,plotting=FALSE)
    
    #---------------------------------
    VolcanoPam(ResultsAll,Topp=5,plotting=TRUE)
 }
}

\keyword{plot}
