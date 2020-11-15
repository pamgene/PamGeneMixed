####**********************************************************************
####**********************************************************************
####
####  PreProcessAllPeptides
####
####  Copyright 2012, I-BioStat, Uhasselt
####
####  This program is free software; you can redistribute it and/or
####  modify it under the terms of the GNU General Public License
####  as published by the Free Software Foundation; either version 2
####  of the License, or (at your option) any later version.
####
####  This program is distributed in the hope that it will be useful,
####  but WITHOUT ANY WARRANTY; without even the implied warranty of
####  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
####  GNU General Public License for more details.
####
####  ----------------------------------------------------------------
####  Written by:
####    --------------------------------------------------------------
####    Pushpike J Thilakarathne, PhD.
####    Interuniversity Institute for Biostatistics and statistical Bioinformatics (I-BioStat) 
####    University of Hasselt and Catholic university of Leuven
####    3000 Leuven
####    Belgium
####
####    email:  pushpike@gmail.com
####    URL:    http://www.uhasselt.be/fiche?voornaam=Pushpike&naam=Thilakarathne
####
####    
####**********************************************************************
####**********************************************************************

########################################################################
# Primary R function for Pre-processing PamChip Array data. Peptide specific extreme observations, profiles and negative intensitiy measurments will be removed. 
########################################################################


##################################################################################################################3

PreProcessAllPeptides<-function(
pep.names,
PamSig,
PathOutPut
){

  if (missing(pep.names)) stop("Peptide names should be given...")
  if (missing(PamSig))    stop("PamGene Data Object should be given ...")
  if (missing(PathOutPut)) {  path<-getwd() 
    } else {
    
    if (file.exists(PathOutPut)){
                        #setwd(file.path(PathOutPut))
                        path<-PathOutPut 
                    } else {
                        dir.create(file.path(PathOutPut))
                        #setwd(file.path(PathOutPut))
                        path<-PathOutPut 
                    
                    }
                              
     }
   
   
   
    Removed.info<-as.data.frame(matrix(NA,ncol=2,nrow=1))
    colnames(Removed.info)<-c("Cell","Peptide")
    t.points<-unique(PamSig[,c("Time")])

    PeptideData<-list()
    
    cellines.unique<-unique(PamSig[,c("CellName")])
    
    WrapPreProcessPam<-function(p){
    #for (p in 1:length(pep.names)) {

        Good.P<-as.data.frame(matrix(rep(1,7),nrow=1,ncol=7) )
        colnames(Good.P)<-c("ID","ResState","ArrayNum","CellName","TreatName","Time",pep.names[p])
                 for (cc in 1:length(cellines.unique)) {
                      CellProfiles<-PreProcessPam(p=p,cel= cc, d =2, PamS=PamSig,plotting=FALSE)
                      Good.Profiles<-CellProfiles@Profiles.Set.ID 
                      colnames(Good.Profiles)<-c("ID","ResState","ArrayNum","CellName","TreatName","Time",pep.names[p])
                                if (nrow(Good.Profiles)>=length(t.points)) {
                                        Good.P<-rbind(Good.P,Good.Profiles)
                                      } else {
                                        Removed.info<-rbind(Removed.info,Good.Profiles)
                                      }
                        }
                    Filttered<-Good.P[-1,];Filttered<-Filttered[!(Filttered$ResState==0),];Filttered<-Filttered[!(Filttered$CellName==0),]
                    Filttered<-Filttered[!(Filttered$TreatName==0),];Filttered<-Filttered[!(Filttered$ArrayNum==0),]
                    Filttered<-Filttered[!(Filttered$Time==0),]
                    Filttered<-Filttered[!is.na(Filttered[,7]),]
                    #write.csv(Filttered,file=paste(P2out,"\\",SPL[p,2],".csv",sep=""),row.names = F)

    
    n.profiles=nrow(Filttered)/length(t.points)
    TempData<-list(Filttered=Filttered,n.profiles=n.profiles)
    
    

    save(TempData,file=paste(path,"/PepProcessed",pep.names[p],".RData",sep=""))
    return(TempData)
    cat("Peptide Index ",p,"\n")
    }
np<-length(pep.names)  
nc<-length(cellines.unique)  
ppp<-1:np
PeptideData<-lapply(ppp,WrapPreProcessPam)

#class(PeptideData)<-"PamGenMixed"
#cat("Results Data frames are saved as RData Objects in ", path)
return(new("CleanedPamData",PeptideData=PeptideData,np=np,nc=nc,path=path))
}
# END OF PreProcessAllPeptides
#################################################
