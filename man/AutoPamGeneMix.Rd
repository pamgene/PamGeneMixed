\name{AutoPamGeneMix}
\alias{AutoPamGeneMix}
\title{Simultaneously Fitting Semi-Parametric Mixed Models to PamChip Data}

\description{
This is a wrapper function for \code{\link{PamGeneMix}} and can be used to fit peptide-specific semi-parametric mixed effects 
models to cleaned PamChip microarray data from function \code{\link{PreProcessAllPeptides}}. The structure of the function is more or less the same as \code{\link{PamGeneMix}}.
However, it automates the analysis for list of peptides.\cr

Function can handel complex mean structures and complex random effects structures depending on the experimental design of the PamChip arrays.
It estimates mean smoothing function and first order derivative, the velocity, of the fitted curve along with point-wise 95\% confidence intervals. 
The estimated group specific velocities are compared and tested at first, end time points as well as for entire velocity profile.  
However, velocities can be compared at any time point within the available time range. All these tests are performed for each peptide separately.\cr

This function can handel two or more grouping structures in the data.  Function returns list of data frames in which p-values are given. 
And corresponsing fitted gamm output ojects for each peptide is saved as Rdata file. \cr

For S4 class information, s. \code{\link{PamChipAllPepResults}}.
}

\usage{
AutoPamGeneMix(formula=y~-1+ResTrt+s(time,by=ResTrt,bs="tp",m=3),
Weights=varIdent(form=~1|ResTrt),
Random.structure=list(ArrayNum=~1,CellLineResTrt=~1+time+time2,ID=~1+time+time2),
CleanedPamData,
TestAt=NA,
PathOutPut)
}

\arguments{

\item{formula}{A GAM formula (see also \code{\link{formula.gam}} and \code{\link{PamGeneMix}}). This is like the formula for a glm except that smooth terms (\code{\link{s}} and \code{\link{te}}) can be added to the right hand side of the formula. Note that ids for smooths and fixed smoothing parameters are not supported.}
\item{Weights}{In the generalized case, weights with the same meaning as \code{\link{glm}} weights. An \code{\link{lme}} type weights argument may only be used in the identity link gaussian case, with no offset (see documentation for \code{\link{lme}} for details of how to use such an argument)}
\item{Random.structure}{ The (optional) random effects structure as specified in a call to \code{\link{lme}}: only the list form is allowed, to facilitate manipulation of the random effects structure within \code{\link{gamm}} in order to deal with smooth terms. See example below.}
\item{CleanedPamData}{A data frame with cleaned PamChip array data from \code{\link{PreProcessAllPeptides}}. }
\item{TestAt}{Time point at which test should be performed, possible value should be  0<t<end time. By default the estimated group specific velocities will be compared and tested at first, end time points as well as for entire velocity profile}
\item{PathOutPut}{Direct path where the output object from PamGeneMix can be saved. If missing objects are saved to current working directory.}
}

\value{ A list of objects of class \code{\link{PamChipAllPepResults}}. 
   
}



\references{
Bowman,A.W. and Azzalini,A. (1997) Applied Smoothing Techniques for Data Analysis.
\emph{Oxford Science Publications}, New York.

Fitzmaurice,G. et al. (eds) (2008) \emph{Longitudinal Data Analysis}. Chapman & Hall/CRC,
New York

Pushpike J. Thilakarathne, Lieven Clement, Dan Lin,Ziv Shkedy, Adetayo Kasim, Willem Talloen, Matthias Versele, and Geert Verbeke. 
The use of semiparametric mixed models to analyze PamChip peptide array data: an application to an oncology experiment. 
\emph{Bioinformatics} (2011) 27(20): 2859-2865 first published online August 16, 2011 doi:10.1093/bioinformatics/btr475. 

Versele,M. et al. (2009) Response prediction to a multitargeted kinase inhibitor in cancer cell lines and xenograft tumors using high-content 
tyrosine peptide arrays with a kinetic readout. \emph{Mol. Cancer Therap.}, 8, 1846-1855.

Wood,S.N. (2006) \emph{Generalized Additive Models}: An Introduction with R . Chapman &
Hall/CRC, New York.



Wood, S.N. (2003) Thin plate regression splines. \emph{J.R.Statist.Soc.B} 65(1):95-114 
}


\author{Pushpike Thilakarathne, Ziv Shkedy and Dan Lin}


\seealso{\code{\link{PamGeneMix}}, \code{\link{PreProcessAllPeptides}}, \code{\link{PreProcessPam}}, \code{\link{gamm}}, \code{\link{lme}}, \code{\link{gam}}, \code{\link{glm}} }






\examples{

\donttest{
    data(PamChipData)
    
    PeptidesNames<-setdiff(colnames(PamChipData),
    c("ID", "ResState", "ArrayNum", "CellName", "TreatName","Time"))
    
    #------- only two groups in the data -------------------
    temp1<-PreProcessAllPeptides(pep.names=PeptidesNames[1:5],
    PamSig=PamChipData[PamChipData[,c("ResState")]=="R",],PathOutPut="C:/Temp")
       
    
    # Case I -----------------------------------
    # Testing the group specific velocities at t=0,t=end and for entire profile.
    
    Results1<-AutoPamGeneMix(formula=y~-1+ResTrt+s(time,by=ResTrt,bs="tp",m=3),
    Weights=varIdent(form=~1|ResTrt),
    Random.structure=list(ArrayNum=~1,CellLineResTrt=~1+time+time2,ID=~1+time+time2),
    temp1,PathOutPut="C:/FittedGamm")
    
    #---- use show 
    show(Results1)
    #- output -
    #-------- Results for set of Peptides --------
    #Class        : PamChipAllPepResults 
    #No. Peptides : 5 
    #Total Run time (in seconds) : 133.89 
    #Number of Groups : 2 
    #Fitted Gamm Objects are saved as RData Objects in :
    # C:/FittedGamm 
    #-----------------------------------------    

    #use --
    VolcanoPam(Results1,Topp=5,plotting=TRUE)
    
    # Case II-----------------------------------  
    # Testing the group specific velocities at t=0,t=end, t=30, and for entire profile.                      
    Results2<-AutoPamGeneMix(formula=y~-1+ResTrt+s(time,by=ResTrt,bs="tp",m=3),
    Weights=varIdent(form=~1|ResTrt),
    Random.structure=list(ArrayNum=~1,CellLineResTrt=~1+time+time2,ID=~1+time+time2),
    temp1,TestAt=30,PathOutPut="C:/FittedGamm")
    
    #---- use show 
    #show(Results2)
    #VolcanoPam(Results2,Topp=5,plotting=TRUE)
         
    #------- Four groups in the data ---
    # Apply preprocessing steps for each peptide
    temp2<-PreProcessAllPeptides(pep.names=PeptidesNames,
    PamSig=PamChipData,PathOutPut="C:/Temp")
    
    #-----------------------------------
    #  Fit the model for each peptide
    #  This would take approximately 6 hrs to run !!! and depends on your system speed.
    ResultsAll<-AutoPamGeneMix(formula=y~-1+ResTrt+s(time,by=ResTrt,bs="tp",m=3),
    Weights=varIdent(form=~1|ResTrt),
    Random.structure=list(ArrayNum=~1,CellLineResTrt=~1+time+time2,ID=~1+time+time2),
    temp2,TestAt=30,PathOutPut="C:/FittedGamm")
    
}

}

\keyword{Mixed models, Smoothing, Generalized Additive Models, Thin-plate regression splines}
