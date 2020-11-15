\name{PamChipAllPepResults-class}
\docType{class}
\alias{PamChipAllPepResults-class}
\alias{PamChipAllPepResults}
\alias{show,PamChipAllPepResults-method}

\title{"PamChipAllPepResults"}
\description{Object returned by \code{\link{AutoPamGeneMix}}.}
\section{Slots}{ 
	 \describe{  
    \item{\code{AllPepRes}:}{A list of data frames as many as peptides analyzed.}
    \item{\code{RunTime}:}{Total run time (in seconds) to analyze the all peptides.}
    \item{\code{NoOfPeptides}:}{Number of peptides analyzed.}
    \item{\code{TestAt}:}{Additional time point at which group-specific velocities are compared.}
    \item{\code{n.groups}:}{Number of  groups in the data.}
    \item{\code{path}:}{Path where the Cleaned Peptide specific RData have been saved.}
  }
}
\section{Methods}{
  \describe{
    \item{show}{Use \code{show(PamChipAllPepResults-object)} for brief information.}
	 }
}


\author{Pushpike Thilakarathne, Ziv Shkedy and Dan Lin}

\seealso{\code{\link{PamGeneMix}}, \code{\link{AutoPamGeneMix}}}

\keyword{Mixed Model, Smoothing, Thin Plate regression}
