\name{PamGeneMix}
\alias{PamGeneMix}
\title{Fitting Semi-Parametric Mixed Models to PamChip Data}


\description{
This is a wrapper function for \code{\link{gamm}} and can be used to fit semi-parametric mixed effects model to PamChip microarray data. Function can handel complex mean structures and complex random effects structures depending on the experimental design of the PamChip arrays.
It estimates mean smoothing function and first order derivative, the velocity, of the fitted curve along with point-wise 95\% confidence intervals. 
The estimated group specific velocities will be compared and tested at first, end time points as well as for entire velocity profile.  However, velocities can be compared at any time point within the available time range.  }

\usage{
PamGeneMix(formula,Correlation=NULL,Weights=varIdent(form=~1),
PTx=list(),Random.structure=NULL,
Control.list=list(maxIter=200, msMaxIter=250 ,msMaxEval=1000,apVar=TRUE),
test.at=NA )
}

\arguments{

\item{formula}{A GAM formula (see also \code{\link{formula.gam}} and \code{\link{gam.models}}). This is like the formula for a glm except that smooth terms (\code{\link{s}} and \code{\link{te}}) can be added to the right hand side of the formula. Note that ids for smooths and fixed smoothing parameters are not supported.}
\item{Random.structure}{ The (optional) random effects structure as specified in a call to \code{\link{lme}}: only the list form is allowed, to facilitate manipulation of the random effects structure within \code{\link{gamm}} in order to deal with smooth terms. See example below.}
\item{Correlation}{An optional corStruct object (see \code{\link{corClasses}}) as used to define correlation structures in \code{\link{lme}}. Any grouping factors in the formula for this object are assumed to be nested within any random effect grouping factors, without the need to make this explicit in the formula (this is slightly different to the behaviour of \code{\link{lme}}).}
\item{Weights}{In the generalized case, weights with the same meaning as \code{\link{glm}} weights. An \code{\link{lme}} type weights argument may only be used in the identity link gaussian case, with no offset (see documentation for \code{\link{lme}} for details of how to use such an argument)}
\item{ Control.list}{A list of fit control parameters for \code{\link{lme}} to replace the defaults returned by \code{\link{lmeControl}}. Note the setting for the number of EM iterations used by \code{\link{lme}}: smooths are set up using custom pdMat classes, which are currently not supported by the EM iteration code. If you supply a list of control values, it is advisable to include niterEM=0, as well, and only increase from 0 if you want to perturb the starting values used in model fitting (usually to worse values!). The optimMethod option is only used if your version of R does not have the nlminb optimizer function.}

\item{PTx}{A data frame with PamChip array data. }
\item{test.at}{Time point at which test should be performed, possible value should be  0<t<end time. By default the estimated group specific velocities will be compared and tested at first, end time points as well as for entire velocity profile}
}


\value{A list of objects of class \code{\link{PamChipMixed}}.  The results can be visualized using the method \code{\link{VisualizePamGeneMix}}.
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


\seealso{ \code{\link{VisualizePamGeneMix}}, \code{\link{AutoPamGeneMix}}, \code{\link{testVarCom}},\code{\link{gamm}}, \code{\link{lme}}, \code{\link{gam}}}

\examples{
#--------------

\donttest{
# Examples of using PamGeneMix
data(TestPepModelData)
head(TestPepModelData)

#---------------- section I -------------------
#In this section, we use only Responisive cell lines data such that we 
#have only two groups: treatment and control for the analysis. 
#
#In this particular model, the cell line-specific random intercept is
#considered to capture correlation of the intensity measurement over time
#within the cell line. We assumed cell line-specific random slopes for
#linear as well as for quadratic time effects to capture different evolution
#of kinase activity over time. Moreover, we allow these cell lines-specific
#random structures to be different for each group. We assumed group-specific
#random structure for cell lines.

PTx<-TestPepModelData

# log2 transform the response 
PTx$y<-log2(PTx[,ncol(PTx)]) 

PTx<-PTx[PTx[,c("ResState")]=="R",]

n.groups<-length(unique(PTx[,c("ResState")]))*length(unique(PTx[,c("TreatName")]))
n.groups

xhlp<-lm(y~-1+TreatName,data=PTx,x=TRUE)$x
for (i in 1:n.groups) PTx$ResTrt[xhlp[,i]==1]<-i
PTx$ResTrt<-as.factor(PTx$ResTrt)
levels(PTx$ResTrt)<-colnames(xhlp)

#--unique cell lines
cellLines<-levels(PTx$CellName)
PTx$ResState<-as.factor(PTx$ResState)

#--create interaction between Treatment and cell lines
PTx$CellLineResTrt<-0
xhlp<-lm(y~-1+CellName:TreatName,data=PTx,x=TRUE)$x
ncols<-ncol(xhlp)
for (i in 1:ncols) PTx$CellLineResTrt[xhlp[,i]==1]<-i
PTx$CellLineResTrt<-as.factor(PTx$CellLineResTrt)
levels(PTx$CellLineResTrt)<-colnames(xhlp)

#------------------ Model I --------------------
M1gamm<-PamGeneMix(formula=y~-1+ResTrt+s(time,by=ResTrt,bs="tp",m=3),
        PTx=PTx,Random.structure=list(ArrayNum=~1,CellLineResTrt=~1+time+time2,ID=~1+time+time2),
        Control.list=list(maxIter=200, msMaxIter=250 ,msMaxEval=1000,apVar=TRUE))

show(M1gamm)
#output
#-------- Fitted PamGeneMix Model --------
#Class: PamChipMixed 
#Note:  Model is successfully fitted 
#Number of Groups:  2 
#P-values for comparing group specific velocities 
#Test at t=0 :  3.34599e-12 
#Test at t=max(t) :  0.1044876 
#Test for entire profile :  0 
#-----------------------------------------

M1gammt30<-PamGeneMix(formula=y~-1+ResTrt+s(time,by=ResTrt,bs="tp",m=3),
        PTx=PTx,Random.structure=list(ArrayNum=~1,CellLineResTrt=~1+time+time2,ID=~1+time+time2),
        Control.list=list(maxIter=200, msMaxIter=250 ,msMaxEval=1000,apVar=TRUE),test.at=30)

show(M1gammt30)
#output
#-------- Fitted PamGeneMix Model --------
#Class: PamChipMixed 
#Note:  Model is successfully fitted 
#Number of Groups:  2 
#P-values for comparing group specific velocities 
#Test at t=0 :  3.34599e-12 
#Test at t= 30  :  0.12664 
#Test at t=max(t) :  0.1044876 
#Test for entire profile :  0 
#-----------------------------------------


#------- summary ----------
summary(M1gammt30)


#------------------ Model II --------------------

#In this particular model is more or less the same as M1.  
#However, correlation structure for the random effects are now common for 
#treatment and control groups. See CellName in Random.structure.


M2gamm<-PamGeneMix(formula=y~-1+ResTrt+s(time,by=ResTrt,bs="tp",m=3),
        PTx=PTx,Random.structure=list(ArrayNum=~1,CellName=~1+time+time2,ID=~1+time+time2),
        Control.list=list(maxIter=200, msMaxIter=250 ,msMaxEval=1000,apVar=TRUE))

show(M2gamm)

#------- summary ----------
summary(M2gamm)

#---------------- section II -------------------
#In this section, we use all Responisive and Non-responsive cell lines data such that we 
#have only two groups: treatment and control for the analysis. 

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


#------------Model III ------------------

#The more complex model for which we assumed that group-specific smoothing parameter and 
#group-specific variance covariance structure. we used thin plate regression splines with third 
#order derivative. In this model nested random effects structure is assumed for replicates within
#the cell lines. And for the cell lines - specific random effects are assumed to be realized 
#from each group separately.

#This model assuming the heteroscedastic variances for responsive and non-responsive cell lines.

M3gamm<-PamGeneMix(formula=y~-1+ResTrt+s(time,by=ResTrt,bs="tp",m=3),Weights=varIdent(form=~1|ResTrt) ,
PTx=PTx,Random.structure=list(ArrayNum=~-1,CellLineResTrt=~1+time+time2,ID=~1+time+time2),
Control.list=list(maxIter=200, msMaxIter=250 ,msMaxEval=1000,apVar=TRUE))

show(M3gamm)

#we can now use testVarCom function to test the hypothesis that 
#H0: Var(group_R_i) - Var(group_R_j) = 0 and  Var(group_NR_i) - Var(group_NR_j) = 0. 
#That is testing for treatment effects conditional on the responsive statues 

testVarCom(M3gamm)  


# we refit the model with common variance for responsive and non-responsive cell lines

M4gamm<-PamGeneMix(formula=y~-1+ResState*TreatName+s(time,by=ResTrt,bs="tp",m=3),Weights=varIdent(form=~1|ResState) ,
PTx=PTx,Random.structure=list(ArrayNum=~-1,CellLineResTrt=~1+time+time2,ID=~1+time+time2),
Control.list=list(maxIter=200, msMaxIter=250 ,msMaxEval=1000,apVar=TRUE))

show(M4gamm)

#we can now use testVarCom function to test the hypothesis that 
#H0: Var(group_R_i) - Var(group_R_j) = 0 and  Var(group_NR_i) - Var(group_NR_j) = 0. 
#That is testing for treatment effects conditional on the responsive statues 

testVarCom(M4gamm) 


}

}

\keyword{Mixed models, Smoothing, Generalized Additive Models, Thin-plate regression splines}
