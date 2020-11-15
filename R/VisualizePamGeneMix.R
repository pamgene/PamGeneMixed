####**********************************************************************
####**********************************************************************
####
####  VisualizePamGeneMix
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
# Primary R function for visualizing the results from PamGeneMix object
########################################################################
#VisualizePamGeneMix


# if true results are graphically visualized. fitted average curves, cell line-specific and replicate-specific, and first order derivative together with its pointwise 95% will be produced separately.  

VisualizePamGeneMix<-function(PamObject,plot.type=c("velocity","velocityCI","smooth.fit","cellline","replicate","SubSpeVelocity"),name.cell=NULL){

if (missing(PamObject)) stop("PamObject object is missing!")

if(class(PamObject)!="PamChipMixed") stop("Invalid Class...! and should be a class of PamChipMixed")

cellLines<-levels(PamObject@PTx$CellName)

plot.type <- match.arg(plot.type)
if(!is.element(plot.type, eval(formals(VisualizePamGeneMix)$plot.type)))
stop("Invalid 'plot.type' specified \n")


if (PamObject@n.g==1)  {Nrow=1;Ncol=1}
if (PamObject@n.g==2)  {Nrow=1;Ncol=2}
if (PamObject@n.g==4)  {Nrow=2;Ncol=2}


if (plot.type=="velocity") {
        #windows("First order derivative over time")
        plot(c(0,(PamObject@End.time+1)),c(min(PamObject@v),max(PamObject@v)),type="n",main=expression(frac(dS(t), dt)),ylab="velocity",xlab="time")
        for (i in 1:PamObject@n.g) lines(PamObject@v[((i-1)*PamObject@End.time+i):(i*PamObject@End.time+i)],col=i,lty=i)
        legend(PamObject@End.time/3,max(PamObject@v),PamObject@gnames,col=1:PamObject@n.g,lty=1:PamObject@n.g)
}

if (plot.type=="velocityCI") {       
       # windows("First order derivative over time with CI")
        par(mfrow=c(Nrow,Ncol))
   
          for (i in 1:PamObject@n.g)
        { 
         plot(c(0:(PamObject@End.time)),ylim=range(c(PamObject@v+2*PamObject@v.sd+0.025,PamObject@v-2*PamObject@v.sd)),type="n",
         main=PamObject@gnames[i],ylab="velocity",xlab="time",cex.main=1.2,cex.lab=1,cex.axis=1)
         polygon(c(PamObject@t.mesh,rev(PamObject@t.mesh)), 
                 c(PamObject@v[((i-1)*PamObject@End.time+i):(i*PamObject@End.time+i)]+2*PamObject@v.sd[((i-1)*PamObject@End.time+i):(i*PamObject@End.time+i)],
                 rev(PamObject@v[((i-1)*PamObject@End.time+i):(i*PamObject@End.time+i)]-2*PamObject@v.sd[((i-1)*PamObject@End.time+i):(i*PamObject@End.time+i)])),
          col="gray",border = NA) 
         lines(PamObject@t.mesh,PamObject@v[((i-1)*PamObject@End.time+i):(i*PamObject@End.time+i)],col=2)
        }
}


if (plot.type=="SubSpeVelocity") {    

 NcolR<-dim(PamObject@gammran3$lme$coef$random$ID)[2]
 if (NcolR==3)  XDer<-cbind(0,1,2*PamObject@t.mesh) 
 if (NcolR==2)  XDer<-cbind(0,1) 
 
 
   DerRan<-PamObject@gammran3$lme$coef$random$ID%*%t(XDer)   

        #windows("Subject Specific First order derivatives")
        par(mfrow=c(Nrow,Ncol))
   
        for (i in 1:PamObject@n.g)
        { 
        
         plot(c(0:(PamObject@End.time)),ylim=range(c(PamObject@v+2*PamObject@v.sd+0.025,PamObject@v-2*PamObject@v.sd)),type="n",
         xlab="time",ylab="First Order Derivative",main=PamObject@gnames[i],cex.main=1.2,cex.lab=1,cex.axis=1)
         
         for (j in 1:nrow(DerRan)) lines(PamObject@t.mesh,PamObject@v[((i-1)*PamObject@End.time+i):(i*PamObject@End.time+i)]+DerRan[j,],col="red")
               
         lines(PamObject@t.mesh,PamObject@v[((i-1)*PamObject@End.time+i):(i*PamObject@End.time+i)],col=1,lwd=2)
        }
}






if (plot.type=="smooth.fit") {    

        
        par(mfrow=c(Nrow,Ncol))
        teller<-0
        for (k in PamObject@gnames)
        {
        teller<-teller+1
        plot(PamObject@PTx[,c("time","y")],ylim=range(PamObject@PTx$y),col=0,main=PamObject@gnames[teller],xlab="Time[min]",ylab="log2[ Intensity ]")
        for (j in cellLines)
            for (i in unique(PamObject@PTx$ID)) 
            {
            lines(PamObject@PTx[PamObject@PTx$ID==i&PamObject@PTx$CellName==j&PamObject@PTx$ResTrt==k,c("time","y")],ylim=range(PamObject@PTx$y),pch=20,col="darkgray",cex = 0.75,lty=1,)
            }
            lines(PamObject@yfitfixef[,c(ncol(PamObject@yfitfixef),teller)],lwd=4)
        }

}



if (plot.type=="cellline"){

        par(mfrow=c(Nrow,Ncol))
        teller<-0
         for (k in PamObject@gnames)
        {
        teller<-teller+1
        plot(PamObject@PTx[,c("time","y")],ylim=range(PamObject@PTx$y),col=0,main=PamObject@gnames[teller],xlab="Time[min]",ylab="log2[ Intensity ]")
        for (j in cellLines)
                    {
                    aaa<-cbind(PamObject@PTx[PamObject@PTx$CellName==j&PamObject@PTx$ResTrt==k,c("time")],PamObject@yfit[PamObject@PTx$CellName==j&PamObject@PTx$ResTrt==k])
                    lines(seq(0,PamObject@End.time,5),sapply(seq(0,PamObject@End.time,5), function(r) mean(aaa[aaa[,1]==r,2])),
                    ylim=range(PamObject@PTx$y),col="blue",lty=2)
                    points(seq(0,PamObject@End.time,5),sapply(seq(0,PamObject@End.time,5), function(r) mean(aaa[aaa[,1]==r,2])),
                    ylim=range(PamObject@PTx$y),col="blue",pch=which(cellLines==j))
                    }
                    lines(PamObject@yfitfixef[,c(ncol(PamObject@yfitfixef),teller)],lwd=3)
        }

}


if (plot.type=="replicate") {
     
if (is.null(name.cell)) stop("Cell name should be provided")
if (length(name.cell)>1) stop("A single Cell name should be given !!!")
if (!is.element(name.cell,unique(PamObject@PTx[,c("CellName")]))) stop("Unknown Cell line Name")

cell.group<-unique(PamObject@PTx[PamObject@PTx[,c("CellName")]==name.cell,c("ResTrt")])

par(mfrow=c(1,length(cell.group)))
        for (k in cell.group) {
        teller<-which(k==PamObject@gnames)
                plot(PamObject@PTx[PamObject@PTx$ID==1,c("time","y")],ylim=range(PamObject@PTx$y),col=0,
                main=paste("Cell Line:",name.cell,"  \n",PamObject@gnames[teller],sep=""),xlab="Time[min]",ylab="log2[ Intensity ]",cex.main=1.2,cex.lab=1,cex.axis=1)
                for (j in name.cell)
                            for (i in unique(PamObject@PTx$ID)) 
                            {
                            points(PamObject@PTx[PamObject@PTx$ID==i&PamObject@PTx$CellName==j&PamObject@PTx$ResTrt==k,c("time","y")],ylim=range(PamObject@PTx$y),col="blue",pch=20)
                            lines(PamObject@PTx[PamObject@PTx$ID==i&PamObject@PTx$CellName==j&PamObject@PTx$ResTrt==k,c("time")],
                            PamObject@yfit[PamObject@PTx$ID==i&PamObject@PTx$CellName==j&PamObject@PTx$ResTrt==k],
                            ylim=range(PamObject@PTx$y),col="black",lty=2)
                            }
                           
         }
        
}



}
