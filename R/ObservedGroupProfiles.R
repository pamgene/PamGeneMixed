
# Function 6:  Group-specific observed profiles (average over cell lines and replicates)
##################################################################################################################3

ObservedGroupProfiles<-function(pepname="Pep5",DataMat,Trt.Group=c("Treatment","Control"),Res="R",log.true=FALSE){

options(warn=-1)
if (missing(pepname)) stop("Peptide name should be given: eg. 'Pep1' ")
if (missing(DataMat)) stop("A data frame should be given")

	if (log.true) { pep<-log2(DataMat[ ,is.element(colnames(DataMat),pepname)])
		pep[is.na(pep)]<-1
		Ylabname<-"log2[Signal Intensity]"
	} else {
		pep<-DataMat[ ,is.element(colnames(DataMat),pepname)]
		Ylabname<-"Signal Intensity"
	}
	
cells<-unique(DataMat[,c("CellName")])
t.points<-unique(DataMat[,c("Time")])

# Replicate Specific ------
par(mfrow=c(length(Trt.Group),length(Res)))
  for (Ri in Res){
     for (Ti in Trt.Group){
           plot(c(min(t.points),max(t.points)),c(min(pep),max(pep)),
           type="n",xlab="Time",ylab="Signal Intensity",main=paste(Ri,Ti,sep=":"))
           for (Ci in cells) {
           index.id<-DataMat[,c("ResState")]==Ri&DataMat[,c("TreatName")]==Ti&DataMat[,c("CellName")]==Ci
           id<-unique(DataMat[index.id,c("ID")])
           temp<-cbind(DataMat[index.id,c("ID")],DataMat[index.id,c("Time")],pep[index.id])
            for (i in id)  lines(temp[temp[,1]==i,2],temp[temp[,1]==i,3],col="gray")
            }
            index.group<-DataMat[,c("ResState")]==Ri&DataMat[,c("TreatName")]==Ti
            temp1<-cbind(DataMat[index.group,c("ID")],DataMat[index.group,c("Time")],pep[index.group])
            lines(t.points,sapply(1:length(t.points),function(i)mean(temp1[temp1[,2]==t.points[i],3])),lwd=3,lty=1)
      }
  }

}

# END OF Function 6: Cell line-specific observed profiles
##################################################################################################################3
