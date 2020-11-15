####**********************************************************************
####**********************************************************************
####
####  VolcanoPam
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

#ObjectAutoPam<-Results2

#setGeneric("VolcanoPam", function(ObjectAutoPam,Topp=5,plotting=TRUE,...) standardGeneric("VolcanoPam"))
#
#setMethod("VolcanoPam", signature(ObjectAutoPam = "PamChipAllPepResults"),
      
VolcanoPam<-function(ObjectAutoPam,Topp=5,plotting=TRUE){

        if (class(ObjectAutoPam)!="PamChipAllPepResults") stop("Argument shoud be of class 'PamChipAllPepResults'...")
        
        Npep<-length(ObjectAutoPam@AllPepRes)
        
        PepName<-paste("Pep",1:Npep,sep="")
        
        tempP<-sapply(1:Npep,function(i)  as.numeric(ObjectAutoPam@AllPepRes[[i]][1,]))
        tenpStat<-sapply(1:Npep,function(i)  as.numeric(ObjectAutoPam@AllPepRes[[i]][2,]))
        
        TestAt<-ObjectAutoPam@TestAt
        
        
#        dotsCall <- substitute(list(...))
#        ll <- eval(dotsCall)
#        if(!hasArg(xlab)) ll$xlab <- "Test Statistics"
#        if(!hasArg(ylab)) ll$ylab <- "-log10[p-value]"
#        #if(!hasArg(ylim)) ll$ylim <- c(0,1)
#        if(!hasArg(pch)) ll$pch<-20
#        if(!hasArg(col)) ll$col<-"lightblue"

                        if ((is.na(TestAt))&(length(TestAt)==1)) {  J=3 
                         M<-c("T=0","T=End","Entire Profile")
                          if (plotting) par(mfrow=c(1,3))
                        } else  {
                        J=4
                        M<-c("T=0","T=End",paste("T=",TestAt,sep=""),"Entire Profile \n Chi-Squared")
                       
                        if (plotting) par(mfrow=c(2,2))
                        }  
                        TopPep<-list()
                        Pvalues<-list()
                        
                        for (j in 1:J){
                                TempP<-tempP[j,]
                                TempS<-tenpStat[j,]
                                names(TempP)<-PepName
                                Pvalues[[j]]<-TempP
                                
                                TempP<-TempP[!is.na(TempP)]
                                TempS<-TempS[!is.na(TempS)]
                                PepName<-PepName[!is.na(TempP)]
                                xx<-TempS
                                xx[xx>100]<-100
                                yy<--log(TempP)
                                
                                inde<-order(yy,decreasing = TRUE)
                               
                               if (plotting) { 
                                #plot(xx,yy,main=M[j],xlab=ll$xlab,ylab=ll$ylab,pch=ll$pch,col=ll$col)
                                plot(xx,yy,main=M[j],xlab="Test Statistics",ylab="-log10[p-value]",pch=20,col="lightblue")
                                points(xx,yy,pch=1,col="blue")
                                text(xx[inde][1:Topp],yy[inde][1:Topp],PepName[inde][1:Topp],col="red",cex=0.6)
                                }
                                
                                TopPep[[j]]<-PepName[inde][1:Topp]
                        }
                        
                        if ((is.na(TestAt))&(length(TestAt)==1)) {  PP<-data.frame(PvalT0=Pvalues[[1]],PvalTEnd=Pvalues[[2]],PvalTProfile=Pvalues[[3]])
                                                TopList<-data.frame(PvalT0=TopPep[[1]],PvalTEnd=TopPep[[2]],PvalTProfile=TopPep[[3]])
                                             } else {
                                               PP<-data.frame(PvalT0=Pvalues[[1]],PvalTEnd=Pvalues[[2]],PvalTAny=Pvalues[[3]],PvalTProfile=Pvalues[[4]])
                                               TopList<-data.frame(PvalT0=TopPep[[1]],PvalTEnd=TopPep[[2]],PvalTAny=TopPep[[3]],PvalTProfile=TopPep[[4]])
                         }


return(list(PP,TopList))
}
#)
