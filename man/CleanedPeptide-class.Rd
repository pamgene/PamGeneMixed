\name{CleanedPeptide-class}
\docType{class}
\alias{CleanedPeptide-class}
\alias{CleanedPeptide}
\alias{show,CleanedPeptide-method}

\title{"CleanedPeptide"}
\description{Object returned by \code{\link{PreProcessPam}}.}
\section{Slots}{
	 \describe{
    \item{\code{Profiles.Set.ID}:}{A data frame of cleaned cell line data of a peptide.}
    \item{\code{p}:}{Index of the Peptide.}
    \item{\code{cel}:}{Index of the cell line used.}
    
  }
}
\section{Methods}{
  \describe{
    \item{show}{Use \code{show(CleanedPeptide-object)} for brief information.}
	 }
}


\author{Pushpike Thilakarathne, Ziv Shkedy and Dan Lin}

\seealso{\code{\link{PreProcessAllPeptides}}}

\keyword{anova, AUC}
