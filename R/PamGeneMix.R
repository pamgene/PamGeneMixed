####**********************************************************************
####**********************************************************************
####
####  PamGeneMix
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



PamGeneMix<-function(
    formula,
    Correlation=NULL,
    Weights=varIdent(form=~1),
    PTx=list(),
    Random.structure=NULL,
    Control.list=list(maxIter=200, msMaxIter=250 ,msMaxEval=1000,apVar=TRUE) ,
    test.at=NA ){


require(nlme,quietly=FALSE)
require(mgcv,quietly=FALSE)

unique.Res<-unique(PTx[,c("ResState")])
unique.Trt<-unique(PTx[,c("TreatName")])


if (!is.data.frame(PTx)) stop(" should be a data frame")

if (is.null(unique.Trt)) stop("treatments labels should be provided e.g. c('treatment','control')") 

if (is.null(unique.Res)) { n.groups<-length(unique.Trt) } else { n.groups<-length(unique.Trt)*length(unique.Res) }

if (is.null(Control.list))     Control.list <- list(niterEM=0,optimMethod="L-BFGS-B") 



if (min(PTx$Time)>0) { PTx$time<-PTx$Time-min(PTx$Time) 
    } else { PTx$time<-PTx$Time }
    PTx$time2<-(PTx$time)**2
    
#n.knots<-unique(PTx$time)

PTx$TreatName<-as.factor(PTx$TreatName)
PTx$CellName<-as.factor(PTx$CellName)
PTx$ResState<-as.factor(PTx$ResState)

n.NRR<-length(unique.Res)
n.trt<-length(unique.Trt)
              
End.time<-max(PTx$time,na.rm=T)
t.mesh<-0:End.time

n.g<-n.NRR*n.trt
          
res.t1<-NA
res.tend<-NA
res.tany<-NA

gammran3<-NA

options(warn=-1)

try(gammran3<-gamm(formula,data=PTx,random=Random.structure,correlation=Correlation,control=Control.list,weights=eval(Weights)), silent = TRUE )

if ((!is.na(gammran3))[1]) {   
              
             yfit<-predict(gammran3$lme) 
             
             
              if (n.g==4) {
                hlp<-as.data.frame( cbind(rep(seq(0,End.time,1),n.g),rep(sort(rep(c(1:n.NRR),(End.time+1))),n.NRR),sort(rep(rep(c(1:n.trt),(End.time+1)),n.trt)),sort(rep(1:n.g,(End.time+1)))) )
                names(hlp)<-c("time","ResState","TreatName","ResTrt")
                hlp$ResTrt<-as.factor(hlp$ResTrt)
                hlp$ResState<-as.factor(hlp$ResState)
                hlp$TreatName<-as.factor(hlp$TreatName)
                levels(hlp$ResTrt)<-levels(PTx$ResTrt)
                levels(hlp$ResState)<-levels(PTx$ResState)
                levels(hlp$TreatName)<-levels(PTx$TreatName)
              }
                    
              if (n.g==2) {
                hlp<-as.data.frame( cbind(rep(seq(0,End.time,1),n.trt),sort(rep(c(1:n.trt),(End.time+1))) ) )
                names(hlp)<-c("time","ResTrt")
                 hlp$ResTrt<-as.factor(hlp$ResTrt)
                 levels(hlp$ResTrt)<-levels(PTx$ResTrt)
                #hlp$TreatName<-as.factor(hlp$TreatName)
                #levels(hlp$TreatName)<-levels(PTx$TreatName)
                }      
                                                    
                yfitfixef<-predict(gammran3$gam,hlp)
                dim(yfitfixef)<-c((End.time+1),n.g)
                yfitfixef<-cbind(yfitfixef,0:End.time)
                                    
              if (is.null(weights)&(n.NRR>1)) {
                                    # common variance 
                                    yfit<-predict(gammran3$lme) #BLUP
                                    hlp<-as.data.frame(cbind(rep(seq(0,End.time,1),n.g),sort(rep(1:n.g,(End.time+1)))))
                                    names(hlp)<-c("time","ResTrt")
                                    hlp$ResTrt<-as.factor(hlp$ResTrt)
                                    levels(hlp$ResTrt)<-levels(PTx$ResTrt)
                                    yfitfixef<-predict(gammran3$gam,hlp)
                                    dim(yfitfixef)<-c((End.time+1),n.g)
                                    yfitfixef<-cbind(yfitfixef,0:End.time)
                }
        
                nCoef<-length(gammran3$gam$coef)
                nbasis<-(nCoef-n.g)/n.g
                
                if (n.g==4)  {L<-c(1,-1,-1,1)  
                H<-matrix(0,nrow=nbasis,ncol=nCoef)
                for (i in 1:(nbasis))  {  H[i, c(0,nbasis,2*(nbasis),3*(nbasis))+i+n.g]<-L }
                }
                
                if (n.g==2) {L<-c(1,-1)
                H<-matrix(0,nrow=nbasis,ncol=nCoef)
                for (i in 1:(nbasis))  {  H[i, c(0,nbasis)+i+n.g]<-L }
                } 
                
                #comparing entire profiles
                LH<-H%*%gammran3$gam$coef
                covLH<-H%*%gammran3$gam$Vp%*%t(H)
                stat<-t(LH)%*%solve(covLH)%*%LH
                
                P.profile<-1-pchisq(stat,nbasis)
                
                PValProfile<-c(P.profile,stat)
                names(PValProfile)<-c("P-Value","TestStat")
                
                delta<-1e-10
                #t.mesh<-0:End.time
                hlp$time<-rep(seq(0,End.time,1),n.g)-delta  
                X0<-predict(gammran3$gam,hlp,type="lpmatrix")
                hlp$time<-rep(seq(0,End.time,1),n.g)+delta 
                X1<-predict(gammran3$gam,hlp,type="lpmatrix")
                #Xp<-(X1-X0)/2/delta
                Xp<-(X1-X0)/delta
                v<-Xp%*%gammran3$gam$coef
                v.sd<-rowSums(Xp%*%gammran3$gam$Vp*Xp)^.5

                 Vall<-Xp%*%gammran3$gam$Vp%*%t(Xp) 
                 
                 
                 tp=1
                 if (n.g==4) {
                 #-------------------- t=0
                 ind.extract<-c(tp ,((tp + 1) + End.time),((tp + 2) + 2*End.time),((tp + 3) + 3*End.time))
                 VarVini<-Vall[ind.extract,ind.extract] # extract varaince covaraince matrix of the initial velocities t=0;
                 VintE<-c(v[ind.extract[1]],v[ind.extract[2]],v[ind.extract[3]],v[ind.extract[4]])
                 
                 Test.t<-(L%*%VintE)/sqrt(t(L)%*%VarVini%*%L)
                 P.initial<-(1-pnorm(abs(Test.t)))
                 res.t1  =c(P.initial,Test.t)
                 names(res.t1)<-c("P-Value","TestStat")
                 
                 #--------------------- t=max(t)
                 tp=End.time+1
                 ind.extract<-c(tp ,((tp + 1) + End.time),((tp + 2) + 2*End.time),((tp + 3) + 3*End.time))
                 VarVendT<-Vall[ind.extract,ind.extract] # extract varaince covaraince matrix of the initial velocities t=0;
                 Vend<-c(v[ind.extract[1]],v[ind.extract[2]],v[ind.extract[3]],v[ind.extract[4]])
                 
                 Test.t<-(L%*%Vend)/sqrt(t(L)%*%VarVendT%*%L)
                 P.last<-(1-pnorm(abs(Test.t)))
                 res.tend=c(P.last,Test.t)
                 names(res.tend)<-c("P-Value","TestStat")
                 }
                 
                 if (n.g==2) {
                 #-------------------- t=0
                 ind.extract<-c(tp ,((tp + 1) + End.time))
                 VarVini<-Vall[ind.extract,ind.extract] # extract varaince covaraince matrix of the initial velocities t=0;
                 VintE<-c(v[ind.extract[1]],v[ind.extract[2]])
                 
                 Test.t<-(L%*%VintE)/sqrt(t(L)%*%VarVini%*%L)
                 P.initial<-(1-pnorm(abs(Test.t)))
                 res.t1  =c(P.initial,Test.t)
                 names(res.t1)<-c("P-Value","TestStat")
                 
                 #--------------------- t=max(t)
                 tp=End.time+1
                 ind.extract<-c(tp ,((tp + 1) + End.time))
                 VarVendT<-Vall[ind.extract,ind.extract] # extract varaince covaraince matrix of the initial velocities t=0;
                 Vend<-c(v[ind.extract[1]],v[ind.extract[2]])
                 
                 Test.t<-(L%*%Vend)/sqrt(t(L)%*%VarVendT%*%L)
                 P.last<-(1-pnorm(abs(Test.t)))
                 res.tend=c(P.last,Test.t)
                 names(res.tend)<-c("P-Value","TestStat")
                 }

                 
                
                if (!is.na(test.at)) {
                tp=test.at+1
                        if (n.g==4) {
                         ind.extract<-c(tp ,((tp + 1) + End.time),((tp + 2) + 2*End.time),((tp + 3) + 3*End.time))
                         VarVanyT<-Vall[ind.extract,ind.extract] # extract varaince covaraince matrix of the initial velocities t=0;
                         VanyT<-c(v[ind.extract[1]],v[ind.extract[2]],v[ind.extract[3]],v[ind.extract[4]])
                         
                            Test.t<-(L%*%VanyT)/sqrt(t(L)%*%VarVanyT%*%L)
                            P.test.at<-(1-pnorm(abs(Test.t)))
                            res.tany<-c(P.test.at,Test.t)
                            names(res.tany)<-c("P-Value","TestStat")
                         }
                         
                         if (n.g==2) {
                         ind.extract<-c(tp ,((tp + 1) + End.time))
                         VarVanyT<-Vall[ind.extract,ind.extract] # extract varaince covaraince matrix of the initial velocities t=0;
                         VanyT<-c(v[ind.extract[1]],v[ind.extract[2]])
                         
                            Test.t<-(L%*%VanyT)/sqrt(t(L)%*%VarVanyT%*%L)
                            P.test.at<-(1-pnorm(abs(Test.t)))
                            res.tany<-c(P.test.at,Test.t)
                            names(res.tany)<-c("P-Value","TestStat")
                         }
    
                 } 
                 
# if model is fitted then return output -------              
return(new("PamChipMixed",
            res.t1  =res.t1,
            res.tend=res.tend,
            res.tany=res.tany,
            test.at=test.at,
            PValProfile=PValProfile,
            t.mesh=t.mesh,
            End.time=End.time,
            n.g=n.g,
            gnames=levels(PTx$ResTrt),
            v=v,
            v.sd=v.sd,
            Vall=Vall,
            yfitfixef=yfitfixef,
            yfit=yfit,
            gammran3=gammran3,
            PTx=PTx)
         )
# end of gamm fit
    }  
    
# if model is NOT fitted then return output ------- 
if((is.na(gammran3))[1])  return(new("PamChipMixed",
                            res.t1  =rep(NA,2),
                            res.tend=rep(NA,2),
                            res.tany=rep(NA,2),
                            test.at=NA,
                            PValProfile=rep(NA,2),
                            t.mesh=t.mesh,
                            End.time=End.time,
                            n.g=n.g,
                            gnames=levels(PTx$ResTrt),
                            v=as.vector(1),
                            v.sd=as.vector(1),
                            Vall=as.matrix(1),
                            yfitfixef=as.matrix(1),
                            yfit=1,
                            gammran3=NA,
                            PTx=PTx)
                       )
    
    
}  
 
