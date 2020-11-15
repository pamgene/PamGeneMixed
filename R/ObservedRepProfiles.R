# Function 4:  Repicate-specific observed profiles for a given cell line
##################################################################################################################3

ObservedRepProfiles<-function(pepname,DataMat,cells,Trt.Group=c("Treatment","Control"),Res="R",log.true=FALSE,cell.mean=T){

options(warn=-1)
if (missing(pepname)) stop("Peptide name should be given: eg. 'Pep1' ")
if (missing(DataMat)) stop("A data frame should be given")
if (missing(cells)) stop("Cell Name Should be Given")

if (log.true) { pep<-log2(DataMat[ ,is.element(colnames(DataMat),pepname)])
	            pep[is.na(pep)]<-1
                Ylabname<-"log2[Signal Intensity]"
} else {
                pep<-DataMat[ ,is.element(colnames(DataMat),pepname)]
                Ylabname<-"Signal Intensity"
}

t.points<-unique(DataMat[,c("Time")])

# Replicate Specific ------
par(mfrow=c(length(Res),length(Trt.Group)))
  for (Ri in Res){
     for (Ti in Trt.Group){
         for (Ci in cells) {
           index.id<-DataMat[,c("ResState")]==Ri&DataMat[,c("TreatName")]==Ti&DataMat[,c("CellName")]==Ci
           id<-unique(DataMat[index.id,c("ID")])
           temp<-cbind(DataMat[index.id,c("ID")],DataMat[index.id,c("Time")],pep[index.id])
           plot(c(min(t.points),max(t.points)),c(min(temp[,3],na.rm=T),max(temp[,3],na.rm=T)),
           type="n",xlab="Time",ylab=Ylabname,main=paste(Ri,":",Ti,"\n Cell Name:",cells,sep=" "))
           for (i in id)  lines(temp[temp[,1]==i,2],temp[temp[,1]==i,3])
            }
         if (cell.mean)   lines(t.points,sapply(1:length(t.points),function(i) mean(temp[temp[,2]==t.points[i],3])),col="red",lwd=3,lty=2)
      }
  }

}

# END OF Function 4: Repicate-specific observed profiles for a given cell line
##################################################################################################################3
