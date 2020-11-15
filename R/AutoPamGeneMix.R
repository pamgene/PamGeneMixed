####**********************************************************************
####**********************************************************************
####
####  AutoPamGeneMix
####
####  Copyright 2012, CenStat, Uhasselt
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
# Primary R function for fitting semi-parametric mixed effects model for PamChip array data. This is a wrapper for gamm
########################################################################


AutoPamGeneMix<-function(
formula=y~-1+ResTrt+s(time,by=ResTrt,bs="tp",m=3),
Weights=varIdent(form=~1|ResTrt),
Random.structure=list(ArrayNum=~1,CellLineResTrt=~1+time+time2,ID=~1+time+time2),
CleanedPamData,
TestAt=NA,
PathOutPut
){

  if (missing(CleanedPamData)) stop("Cleaned PamData  should be given...")
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
   
#if (is.na(TestAt))  {TestAt<--1}

NoOfPeptides<-length(CleanedPamData@PeptideData)


PTData<-CleanedPamData@PeptideData[[1]]$Filttered
              
n.groups<-length(unique(PTData[,c("ResState")]))*length(unique(PTData[,c("TreatName")]))



ptm <- proc.time()
AllPepRes<-lapply(1:NoOfPeptides,function(pIndex) {

#                            if (missing(PathInPut)) { load(file=paste("PepProcessed",pep.names[i],".RData",sep="")) 
#                                 } else {
#                                                      load(file=paste(PathInPut,"PepProcessed",pep.names[i],".RData",sep="")) 
#                             }
#    PTx<-TempData

#i=1                 
TempGamm<-NULL
                        PTx<-CleanedPamData@PeptideData[[pIndex]]$Filttered
                        #colnames(PTx)
                        
                        PTx$y<-log2(PTx[,ncol(PTx)]) 
                        
                        nRNR<-length(unique(PTx[,c("ResState")]))
                        
                        n.groups<-length(unique(PTx[,c("ResState")]))*length(unique(PTx[,c("TreatName")]))
                        #n.groups
                        
                       
                        PTx$TreatName<-as.factor(PTx$TreatName)
                        PTx$CellName<-as.factor(PTx$CellName)
                        PTx$ResState<-as.factor(PTx$ResState)               
                        
                        
                        if (nRNR>=2) {
                                    PTx$ResTrt<-0
                                    xhlp<-lm(y~-1+ResState:TreatName,data=PTx,x=TRUE)$x
                                    for (k in 1:n.groups) PTx$ResTrt[xhlp[,k]==1]<-k
                                    PTx$ResTrt<-as.factor(PTx$ResTrt)
                                    levels(PTx$ResTrt)<-colnames(xhlp)
                                    
                                    #unique cell lines
                                    cellLines<-levels(PTx$CellName)
                                    
                                    
                                    # create interaction between Treatment and cell lines
                                    PTx$CellLineResTrt<-0
                                    xhlp<-lm(y~-1+CellName:TreatName,data=PTx,x=TRUE)$x
                                    ncols<-ncol(xhlp)
                                    for (j in 1:ncols) PTx$CellLineResTrt[xhlp[,j]==1]<-j
                                    PTx$CellLineResTrt<-as.factor(PTx$CellLineResTrt)
                                    levels(PTx$CellLineResTrt)<-colnames(xhlp)
                                    
                          TempGamm<-PamGeneMix(formula=formula,Weights=Weights ,PTx=PTx,Random.structure=Random.structure,
                                    Control.list=list(maxIter=200, msMaxIter=250 ,msMaxEval=1000,apVar=TRUE),test.at=TestAt)
                                    
                          if ((is.na(TempGamm@gammran3))[1])  {   
                                                  Random.structure<-Random.structure[-1]   
                                                  TempGamm<-PamGeneMix(formula=formula,Weights=Weights ,PTx=PTx,Random.structure=Random.structure,
                                                            Control.list=list(maxIter=200, msMaxIter=250 ,msMaxEval=1000,apVar=TRUE),test.at=TestAt)  
                                                  } 
                                    
                        }
                        
                     if (nRNR==1) {
                            n.groups<-length(unique(PTx[,c("ResState")]))*length(unique(PTx[,c("TreatName")]))
                            #n.groups
                            
                            xhlp<-lm(y~-1+TreatName,data=PTx,x=TRUE)$x
                            for (i in 1:n.groups) PTx$ResTrt[xhlp[,i]==1]<-i
                            PTx$ResTrt<-as.factor(PTx$ResTrt)
                            levels(PTx$ResTrt)<-colnames(xhlp)
                            
                            
                            #unique cell lines
                            cellLines<-levels(PTx$CellName)
                            PTx$ResState<-as.factor(PTx$ResState)
                            
                            # create interaction between Treatment and cell lines
                            PTx$CellLineResTrt<-0
                            xhlp<-lm(y~-1+CellName:TreatName,data=PTx,x=TRUE)$x
                            ncols<-ncol(xhlp)
                            for (l in 1:ncols) PTx$CellLineResTrt[xhlp[,l]==1]<-i
                            PTx$CellLineResTrt<-as.factor(PTx$CellLineResTrt)
                            levels(PTx$CellLineResTrt)<-colnames(xhlp)
                            
                            
                            TempGamm<-PamGeneMix(formula=formula,Weights=Weights ,PTx=PTx,Random.structure=Random.structure,
                                    Control.list=list(maxIter=200, msMaxIter=250 ,msMaxEval=1000,apVar=TRUE),test.at=TestAt)
                          
                          if ((is.na(TempGamm@gammran3))[1])  {   
                                                  Random.structure<-Random.structure[-1]   
                                                  TempGamm<-PamGeneMix(formula=formula,Weights=Weights ,PTx=PTx,Random.structure=Random.structure,
                                                            Control.list=list(maxIter=200, msMaxIter=250 ,msMaxEval=1000,apVar=TRUE),test.at=TestAt)  
                                                  } 
                     }
                        
         
          save(TempGamm,file=paste(path,"/Gamm3FitPep",pIndex,".RData",sep="") ) 
          
          if (is.na(TestAt)) {  PMat<-data.frame(TempGamm@res.t1,TempGamm@res.tend,TempGamm@PValProfile)
                               } else { 
                                  PMat<-data.frame(TempGamm@res.t1,TempGamm@res.tend,TempGamm@res.tany,TempGamm@PValProfile)
          }
          #print(PMat)
          return(as.data.frame(PMat))
        
         }
      )
EndT<-as.vector(proc.time() - ptm )
RunTime<-EndT[1]
#cat("Total run time (in seconds): ", EndT[1] ,"\n","Number Peptides: ",NoOfPeptides,"\n Fitted Gamm Objects are saved in ", path,"\n")

#class(AllPepRes)<-"PamChipAllPepResults"
#if (is.null(TestAt)) rownames(AllPepRes)<-c("t0","tEnd","Tpro")
#if (!is.null(TestAt)) rownames(AllPepRes)<-c("t0","tEnd","tAny","Tpro")
return(new("PamChipAllPepResults",AllPepRes=AllPepRes,RunTime=RunTime,NoOfPeptides=NoOfPeptides,TestAt=TestAt,n.groups=n.groups,path=path))
}
