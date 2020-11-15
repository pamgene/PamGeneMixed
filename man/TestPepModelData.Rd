\name{TestPepModelData}

\alias{TestPepModelData}

\docType{data}


\title{PamChip Peptide Specific Cleaned Data }


\description{
This is a cleaned PamChip Peptide specific dataset  where there are columns for ArrayID, Cell line names, Responsive Statuses, Treatment indicator, Time points, and Replicate ID. 
Intensity measurements for the peptide is given in the last column.
}


\usage{data(TestPepModelData)}


\format{
The format is a data frame which has following structure.

\tabular{ccccccc}{
ID	\tab	ResState	\tab	ArrayNum	\tab	CellName	\tab	TreatName	\tab	Time	\tab	Peptide	\cr
6	\tab	R	\tab	0559C80	\tab	H3255	\tab	Control	\tab	11	\tab	2706	\cr
6	\tab	R	\tab	0559C80	\tab	H3255	\tab	Control	\tab	16	\tab	4442	\cr
6	\tab	R	\tab	0559C80	\tab	H3255	\tab	Control	\tab	21	\tab	5045	\cr
6	\tab	R	\tab	0559C80	\tab	H3255	\tab	Control	\tab	26	\tab	5491	\cr
6	\tab	R	\tab	0559C80	\tab	H3255	\tab	Control	\tab	31	\tab	5950	\cr
6	\tab	R	\tab	0559C80	\tab	H3255	\tab	Control	\tab	36	\tab	6124	\cr
}

}


\seealso{\code{\link{PamChipData}}}

\author{Pushpike Thilakarathne, Ziv Shkedy and Dan Lin}

\examples{
#-- load Data --
data(TestPepModelData)

#-- Check Data --
TestPepModelData[1:5,]
  
#-- Unique Time Points --
UTime<-unique(TestPepModelData[,c("Time")])
UTime

#-- Responsive and Non-Responsive Cell Lines --
RnRcellines<-unique(TestPepModelData[,c("CellName","ResState")])
rownames(RnRcellines)<-1:20
RnRcellines

#-- Number of Plates (96 well plates)---
Plates<-unique(TestPepModelData[,c("ArrayNum")])
Plates


}


\keyword{datasets}
