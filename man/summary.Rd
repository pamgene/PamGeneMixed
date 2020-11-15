\name{summary}
\alias{summary,PamChipMixed-method}
\title{Summarize Fitted Gamm Object for a given Peptide}
\description{This method principally does nothing more than
applying the pre-implemented \code{summary()} function to the slot \code{gammran3}
of an object of class \code{\link{PamChipMixed}}. One then obtains estimated group specific velocities at t=0, t=max(t), t=any, and t=end.  And also fitted \code{\link{lme}} and \code{\link{gam}} objects are shown.}

\arguments{
  \item{object}{An object of class \code{\link{PamChipMixed}}.}
}

\value{No return.}

\author{Pushpike Thilakarathne, Ziv Shkedy and Dan Lin}

\seealso{\code{\link{VisualizePamGeneMix}}, \code{\link{gamm}}.}
\keyword{Semi-Parametric Mixed Models}
