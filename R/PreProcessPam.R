####**********************************************************************
####**********************************************************************
####
####  PreProcessPam
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




#                                   removing outline and negative profiles 
#--------------------------------------------------------------------------------------------------------------
# Function 7:   removing outline and negative profiles 
##################################################################################################################3

# Trapezoidal rule
# t values need be equally spaced No missing values
# d=t2-t1
Trapezoid <- function(dix,y) (y[1]+2*sum(y[seq(2,(length(y)-1))])+y[length(y)])*dix/2
#  Function 1.1: Trapezoidal rule
# t values need not be equally spaced and missing values are allows
# if missing, linear interpolation is used
Trapezoid.Miss <- function(d,y) {
    index<-which(!is.na(y))
    n<-length(index)
    AUC<-d*sum((y[index[1:n-1]]+y[index[2:n]])*(index[2:n]-index[1:n-1]))/2
   return(AUC)
}

# END OF Function 7:   removing outline and negative profiles 
##################################################################################################################3



# Function 8: remove observation if it is extrem by looking at the individual profiles
##################################################################################################################3

  Extrem.pset<-function(pset,IDnum){
  for (ID in IDnum){
        Temp.P<-pset[pset[,1]==ID,7]
        Temp.P[Temp.P<=0]<-0.01
        TempSet<-log2(Temp.P)
        pset.lower<-mean(TempSet)-3*sd(TempSet)
        pset.upper<-mean(TempSet)+3*sd(TempSet)
        Temp.P[(Temp.P<(2^pset.lower))|((2^pset.upper)<Temp.P)]<-NA 
        Temp.P[1]<-(pset[pset[,1]==ID,7])[1]
       # iinde<-c(T,T,T,T,F,F,F,F,F,T,T,T,T)
       # Temp.P[iinde]<-(pset[pset[,1]==ID,7])[iinde] # only spikes will be removed
       pset[(pset[,1]==ID),7]<-Temp.P
       }
   return(pset)
  }
# END OF Function 8: remove observation if it is extrem by looking at the individual profiles
##################################################################################################################3



# Function 9:  estmiate upper and lower limits to remove outliers Based on anova model with AUC
##################################################################################################################3

AUC.trt.anova<-function(propset,d=2){
           
           Check.Res<-function(anv,sdauc,i){
                upper.anova<-mean(anv[anv[,2]==i,1])+3*sdauc
                lower.anova<-mean(anv[anv[,2]==i,1])-3*sdauc
                if (any((upper.anova<anv[anv[,2]==i,1])|(lower.anova>anv[anv[,2]==i,1]))) {RM.ID<-anv[((upper.anova<anv[anv[,2]==i,1])|(lower.anova>anv[anv[,2]==i,1])),3]
                } else {RM.ID<-NA}
                return(RM.ID)
            }
            #out.ID=NA; out.ID.obs=NA
            Pprof.ID<-unique(propset[,1])
            Trt<-sapply(Pprof.ID,function(ind) unique(propset[(propset[,1]==ind),5]))
            AUC<-sapply(Pprof.ID,function(ind) Trapezoid.Miss(5,propset[(propset[,1]==ind),7]))
            ID.AUC<-data.frame(Pprof.ID,AUC,Trt); colnames(ID.AUC)<-c("PosID","AUC","Trt")
            
           if (length(unique(Trt))>1){
                                        # fit ANOVA with AUC values
                                        lm.auc<-lm(AUC~as.factor(Trt),data=as.data.frame(ID.AUC))
                                        anova.auc<-aov(AUC~as.factor(Trt),data=as.data.frame(ID.AUC))
                                        sdauc<-sqrt(sum(resid(anova.auc)^2)/anova.auc$df)
                                        anv<-cbind(lm.auc$residuals,Trt,ID.AUC)
                                        TC<-unique(Trt)
                                        if (!is.na(sdauc)){# check EMS
                                            
                                            for (j in TC) { #removing outline profiles 
                                                            t1.ID<-Check.Res(anv,sdauc,j)
                                                            if (any(!is.na(t1.ID))) for (i in t1.ID) propset<-propset[propset[,1]!=i,]
                                                           }
                                        } 
            }
            # use AUC to remove outliers if at least two profiles per Trt present
            if ((length(unique(propset[propset[,5]==TC[2],1]))>2)&(length(unique(propset[propset[,5]==TC[1],1]))>2)){
                    upper.AUC<-mean(AUC)+d*sqrt(var(AUC))
                    lower.AUC<-mean(AUC)-d*sqrt(var(AUC))
                    # if there are outliers then remove them
                    if ((any((lower.AUC>ID.AUC[,2])|(ID.AUC[,2]>upper.AUC)))&(dim(ID.AUC)[1]>1)) { out.ID<-ID.AUC[(lower.AUC>ID.AUC[,2])|(ID.AUC[,2]>upper.AUC),1]
                                                                                                     for (i in out.ID) propset<-propset[propset[,1]!=i,]     }
             }
             
return(propset)
}
# END OF Function 9:  estmiate upper and lower limits to remove outliers Based on anova model with AUC
##################################################################################################################3



# Function 10:  ploting filltering steps
##################################################################################################################3

Plot.Good.Profiles<-function(Profile,RedProfiles,Peptide,cc,p,OutMe="Without Outliers",Res,Trt){
       
        # rescaled time points
        tt<-unique(Profile[,c("Time")])
        if (tt[1]!=0) Time<-tt-tt[1]
         
       if (dim(Profile)[1]==dim(RedProfiles)[1]) { Prof<-Profile; if (any(is.na(RedProfiles[,7]))) {Prof<-RedProfiles;
                                                                  } else { if (all(RedProfiles[,7]>0)) Prof<-RedProfiles;}
       } else {Prof<-RedProfiles}  
              
        # Number of Profiles per peptide in a given cell line :some peptides have more than six profiles        
        jEnd      <-(dim(Prof)[1])/length(Time) 
        
        if (unique(Profile[,c("ResState")])==Res[1]) {R<-"Responsive Cell line: "
                                         } else {R<-"NON-Responsive Cell line: "}
        
        MName=paste(OutMe,"\n Peptide ",p, " and Cell line ",cc,sep="")
        
        plot(c(0,12),c(min(Profile[,7])-20,max(Profile[,7])+20),type="n",xlab="Time ",ylab="Normalized Signal",main=MName,ask=T,axes=F,cex.main=1,cex.lab=1.2,cex.axis=1.1)
         axis(1,at=0:(length(Time)-1),Time[1:length(Time)]); axis(2); box(); 
        for (i in 1:jEnd)
            {coll="blue";ltyx=1; CHP=17;
             if (all(Prof[((i-1)*length(Time)+1):(i*length(Time)),c("TreatName")]==Trt[1])) {coll="red"; ltyx=1;CHP=19;
                                                     }  
             lines(Time/5,Prof[((i-1)*length(Time)+1):(i*length(Time)),7],col=coll,lty=ltyx,lwd=1) 
             points(Time/5,Prof[((i-1)*length(Time)+1):(i*length(Time)),7],pch=CHP,col=coll) 
            }   
            legend(0,(max(Profile[,7])+10),c("Control","Compound"),col=c("red","blue"),pch=c(17,19),bty="n",lwd=1, pt.cex=1.1,cex=1.1);
        
}
# END OF Function 10:  ploting filltering steps
##################################################################################################################3



# Function 11:  Filter profiles AUC and NEG 
# this function will call Function 1.1,2,3, and 4
# first 6 columns :"ID","ResState","ArrayNum","CellName","TreatName","Time"
##################################################################################################################3

PreProcessPam<-function(
        p   =108, # Peptide Index
        cel =4, #Cell line index
        d   =2, #  d*sd(AUC)
        PamS,
      plotting=TRUE
){

        LP<-setdiff(colnames(PamS),c("ID", "ResState", "ArrayNum", "CellName", "TreatName","Time")) 
        Cell<-as.vector(unique(PamS[,c("CellName")]))   
        n.t<-length(unique(PamS[,c("Time")]))
        Profiles.Set.ID<-PamS[(PamS[,c("CellName")]==Cell[cel]),c(1:6,which(LP[p]==colnames(PamS)))]
        Observed.Profiles<-Profiles.Set.ID
       
        Res<-PamS[(PamS[,c("CellName")]==Cell[cel]),c("ResState")]
        Trt<-unique(PamS[(PamS[,c("CellName")]==Cell[cel]),c("TreatName")])
        
        # plot two set of profiles: reduced and all time points    
        if (plotting) {par(mfrow=c(1,2)) 
        Plot.Good.Profiles(Profile=Observed.Profiles,RedProfiles=Profiles.Set.ID,Peptide=LP,cc=cel,p=p,OutMe="Observed Data",Res=Res,Trt=Trt)
        }
        #replacing negative values or zeros by min of the positive values 
        if (any(Profiles.Set.ID[,7]<=0)) { 
                        Neprof.index<-(Profiles.Set.ID[,7]<=0)
                        negprof.ID<-unique(Profiles.Set.ID[Neprof.index,1])
                            for (i in negprof.ID)
                            {  val.of.id<-Profiles.Set.ID[Profiles.Set.ID[,1]==i,7]
                                if (any(val.of.id>0)) {
                                    min.of.pos<-min(val.of.id[val.of.id>0])/2; if (min.of.pos==0) {min.of.pos<-0.01}
                                } else { min.of.pos<-0.01}
                                
                                val.of.id[val.of.id<=0]<-min.of.pos
                                Profiles.Set.ID[Profiles.Set.ID[,1]==i,7]<-val.of.id
                            }
        }
        
        
        # removing outline observations within a profile
        # recall Function Extrem.pset

          if (length(Trt)>1){ #apply function treatment vice
            h=0
            ddd=0
                    for (k in Trt) { h=h+1
                     eval(parse(text= paste("ddd<-Extrem.pset(pset=Profiles.Set.ID[Profiles.Set.ID[,5]=='",k,"',],IDnum=unique(Profiles.Set.ID[Profiles.Set.ID[,5]=='",k,"',1]))",sep="")))
                    if (h==1) dddall<-ddd
                    if (h>1) dddall<-rbind(dddall,ddd)
                    } 
          #ddd2<-Extrem.pset(pset=Profiles.Set.ID[Profiles.Set.ID[,5]=="DMSO",],IDnum=unique(Profiles.Set.ID[Profiles.Set.ID[,5]=="DMSO",1]))
          Profiles.Set.ID<-dddall
          #Profiles.Set.ID<-rbind(ddd1,ddd2)
          } else {  
          Pprof.ID       <-unique(Profiles.Set.ID[,1])
          Profiles.Set.ID<-Extrem.pset(pset=Profiles.Set.ID,IDnum=Pprof.ID)
          }
      
        # if at least one profiles are present
        if (dim(Profiles.Set.ID)[1]>=n.t) { 
        # recall Function AUC.trt.anova
             TempSet<-AUC.trt.anova(propset=Profiles.Set.ID,d=d)
             Profiles.Set.ID<-TempSet
             
        # recall Function Plot.Good.Profiles
         if (plotting)   Plot.Good.Profiles(Profile=Observed.Profiles,RedProfiles=Profiles.Set.ID,Peptide=LP,cc=cel,p=p,OutMe="Filttered Profiles",Res=Res,Trt=Trt)
            }
        # record completly removed peptides and cell lines combination
        removed.C.P<-matrix(NA,nrow=1,ncol=2)
        if (!(dim(Profiles.Set.ID)[1]>=n.t)){ removed.C.P<-as.data.frame(cbind(cel,p))
                    colnames(removed.C.P)<-c("Cell","Peptide")
        return(removed.C.P)           
        } else { # record Filtered profiles
        return(new("CleanedPeptide",Profiles.Set.ID=Profiles.Set.ID,p=p,cel=cel))
        }
} # end of PreProcessPam

# END OF Function 11:  Filter profiles AUC and NEG 
##################################################################################################################3
