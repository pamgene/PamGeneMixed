\name{CleanedPamData-class}
\docType{class}
\alias{CleanedPamData-class}
\alias{CleanedPamData}
\alias{show,CleanedPamData-method}

\title{"CleanedPamData"}
\description{Object returned by \code{\link{PreProcessAllPeptides}}.}
\section{Slots}{ 
	 \describe{
    \item{\code{PeptideData}:}{A list of data frames as many as peptides analyzed.}
    \item{\code{np}:}{Number of peptides analyzed.}
    \item{\code{nc}:}{Number of  cell lines used.}
    \item{\code{path}:}{Path where the Cleaned Peptide specific RData have been saved.}
  }
}
\section{Methods}{
  \describe{
    \item{show}{Use \code{show(CleanedPamData-object)} for brief information.}
	 }
}


\author{Pushpike Thilakarathne, Ziv Shkedy and Dan Lin}

\seealso{\code{\link{PreProcessPam}}}

\keyword{anova, AUC}
