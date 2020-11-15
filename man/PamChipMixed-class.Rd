\name{PamChipMixed-class}
\docType{class}
\alias{PamChipMixed-class}
\alias{PamChipMixed}
\alias{PamChipMixed-method}
\alias{show,PamChipMixed-method}

\title{"PamChipMixed"}
\description{Object returned by the method \code{\link{PamGeneMix}}.}
\section{Slots}{
     \describe{
  \item{res.t1}{A vector of length two in which \emph{p-value} for comparing group-specific velocities at t=0 along with test statistic are returned.} 
  \item{res.tend}{A vector of length two in which \emph{p-value} for comparing group-specific velocities at t=max(t) along with test statistic are returned.}
  \item{res.tany}{If test.at is given then a vector of length two in which \emph{p-value} for comparing group-specific velocities at t=tgiven along with test statistic are returned. Otherwise NA is returned.}
  \item{test.at}{Time point at which additional test needs to be performed.} 
  \item{PValProfile}{A vector of size two where first gives p-value for comparing group-specific velocities for entire profile and second value gives the Ch-Squared test statistic.}
  \item{t.mesh}{Fine grid of time points that has been used to produce smooth curves.} 
  \item{End.time}{Maximum time point in the Dataset.} 
  \item{n.g}{Number of groups in the data} 
  \item{gnames}{Group names.}
  \item{v}{Estimated velocity. And it is a vector of size equals to length of the t.mesh times number of groups.}
  \item{v.sd}{Estimated variances at each time point.  And it is a vector of size equals to length of the t.mesh times number of groups.}
  \item{Vall}{Full variance covariance matrix of the estimated velocities.}
  \item{yfitfixef}{Fitted mean function. Returns as a vector of size equals to \emph{length(t.mesh)} times \emph{n.g}}
  \item{yfit}{Estimated \emph{BLUP}.}
  \item{gammran3}{Fitted \code{\link{gamm}} object.}
  \item{PTx}{A data frame used to fit the model.}
  }
}
\section{Methods}{
  \describe{
    \item{show}{Use \code{show(PamChipMixed-object)} for brief information.}
    \item{summary}{Use \code{summary(PamChipMixed-object)} to apply the
                   classic \code{summary()} function to the  \code{PamChipMixed-object)},
                   s. \code{\link{summary,PamChipMixed-method}}}
     }
}

\author{Pushpike Thilakarathne, Ziv Shkedy and Dan Lin}

\seealso{\code{\link{PamGeneMix}}}
\keyword{Mixed models, Smoothing, Generalized Additive Models, Thin-plate regression splines}
