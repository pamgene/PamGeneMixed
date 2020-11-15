####**********************************************************************
####**********************************************************************
####
####  testVarCom
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

#setGeneric("testVarCom", function(ObjPamFit) standardGeneric("testVarCom"))
#
#setMethod("testVarCom", signature(ObjPamFit = "PamChipMixed"), function(ObjPamFit){
# 
testVarCom<-function(ObjPamFit){

if (class(ObjPamFit)!="PamChipMixed") stop("Argument is not class of 'PamChipMixed'")

#ngroups<-ObjPamFit@n.g
gammOB <-ObjPamFit@gammran3

if ((is.na(gammOB))[1]) stop("Model is not Fitted ! or has failed to converge.")    
 
          #if (missing(ngroups)) stop("Number groups should be provided !!")
           
           if (!is.character(gammOB$lme$apVar)&!missing(gammOB)) {
                NameGroups<-attr(gammOB$lme$modelStruct$varStruct, "groupNames")
               
                NameGroups<-NameGroups[!is.na(NameGroups)]
                
                sigma2<-summary(gammOB$lme)$sigma^2 
                Sigma2<-(exp(attr(gammOB$lme$apVar,'Pars')))^2
                Var.groups<-sigma2*c(1,Sigma2[-c(1:(length(Sigma2)-length(NameGroups)),length(Sigma2))])
                names(Var.groups)<-NameGroups
                ngroups<-length(Var.groups)
                
                selog = sqrt(diag(gammOB$lme$apVar))
                senolog = exp(attr(gammOB$lme$apVar,'Pars')) * selog
                sevarnolog = 2 * exp(attr(gammOB$lme$apVar,'Pars')) * senolog
                
                sevg<-c(sevarnolog[-(1:(length(sevarnolog)-ngroups))][ngroups],sevarnolog[-(1:(length(sevarnolog)-ngroups))][ngroups]*sevarnolog[-(1:(length(sevarnolog)-ngroups))][1:(ngroups-1)])
                Test.Stat<-Var.groups/sevg
                
                names(sevg)<-NameGroups
                
                if (ngroups==4) { 
                            L1=c(1,-1,0,0)
                            p1<-(1-pnorm(abs(Test.Stat%*%L1/sqrt(2))))   
                            L2=c(0,0,1,-1)
                            p2<-(1-pnorm(abs(Test.Stat%*%L2/sqrt(2)))) 
                           p1=c(p1,p2)
                           return(list(p1=p1,Var.groups=Var.groups,sevg=sevg))
                           }
                if (ngroups==2) {L1=c(1,-1)
                             p1<-(1-pnorm(abs(Test.Stat%*%L1/sqrt(2)))) 
                             return(list(p1=p1,Var.groups=Var.groups,sevg=sevg))
                             }
                
                 } else      {
                stop("There is no grouping structure in the data")}
              
    }
 #)
