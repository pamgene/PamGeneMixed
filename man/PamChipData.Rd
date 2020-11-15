\name{PamChipData}

\alias{PamChipData}

\docType{data}


\title{PamChip Microarray Data }


\description{
This is PamChip Data frame where there are columns for ArrayID, Cell line names, Responsive Statuses, Treatment indicator, Time points, and Replicate ID. And there are columns as many as peptides.
}


\usage{data(PamChipData)}


\format{
The format is a data frame which has following structure.

\tabular{ccccccccc}{
 ID 	\tab	 ResState 	\tab	 ArrayNum 	\tab	 CellName 	\tab	 TreatName 	\tab	 Time 	\tab	 Pep1 	\tab	 Pep2 	\tab	 Pep3 	\cr
1	\tab	 NR 	\tab	 0559C80 	\tab	 MDA231PAR 	\tab	 Control 	\tab	11	\tab	2237	\tab	2104	\tab	287	\cr
1	\tab	 NR 	\tab	 0559C80 	\tab	 MDA231PAR 	\tab	 Control 	\tab	16	\tab	3746	\tab	3696	\tab	391	\cr
1	\tab	 NR 	\tab	 0559C80 	\tab	 MDA231PAR 	\tab	 Control 	\tab	21	\tab	3979	\tab	4051	\tab	460	\cr
1	\tab	 NR 	\tab	 0559C80 	\tab	 MDA231PAR 	\tab	 Control 	\tab	26	\tab	4361	\tab	4392	\tab	394	\cr
1	\tab	 NR 	\tab	 0559C80 	\tab	 MDA231PAR 	\tab	 Control 	\tab	31	\tab	4218	\tab	4452	\tab	393	\cr
7	\tab	 R 	\tab	 0559C80 	\tab	 N87 	\tab	 Control 	\tab	36	\tab	4679	\tab	4855	\tab	832	\cr
7	\tab	 R 	\tab	 0559C80 	\tab	 N87 	\tab	 Control 	\tab	41	\tab	4726	\tab	4777	\tab	920	\cr
7	\tab	 R 	\tab	 0559C80 	\tab	 N87 	\tab	 Control 	\tab	46	\tab	5072	\tab	5072	\tab	952	\cr
7	\tab	 R 	\tab	 0559C80 	\tab	 N87 	\tab	 Control 	\tab	51	\tab	5147	\tab	5344	\tab	1039	\cr
7	\tab	 R 	\tab	 0559C80 	\tab	 N87 	\tab	 Control 	\tab	56	\tab	5091	\tab	5204	\tab	999	\cr
7	\tab	 R 	\tab	 0559C80 	\tab	 N87 	\tab	 Control 	\tab	61	\tab	5060	\tab	5172	\tab	1131	

}

}

\seealso{\code{\link{TestPepModelData}}}

\author{Pushpike Thilakarathne, Ziv Shkedy and Dan Lin}

\references{
Thilakarathne,P.J. et al. (2011) The use of semiparametric mixed models to analyze PamChip peptide array data: an application to an oncology experiment Bioinformatics 27(20): 2859-2865 doi:10.1093/bioinformatics/btr475. 

Versele,M. et al. (2009) Response prediction to amultitargeted kinase inhibitor in cancer cell lines and xenograft tumors using high-content tyrosine peptide arrays with a kinetic readout. Mol. Cancer Therap., 8, 1846-1855.
}

\examples{
#-- load Data --
data(PamChipData)

#-- Check Data --
PamChipData[1:5,1:10]
  
#-- Unique Time Points --
UTime<-unique(PamChipData[,c("Time")])
UTime

#-- Responsive and Non-Responsive Cell Lines --
RnRcellines<-unique(PamChipData[,c("CellName","ResState")])
rownames(RnRcellines)<-1:20
RnRcellines

#-- Number of Plates (96 well plates)---
Plates<-unique(PamChipData[,c("ArrayNum")])
Plates


}


\keyword{datasets}
