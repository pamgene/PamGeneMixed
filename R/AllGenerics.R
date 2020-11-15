#
#

# Author: Pushpike Thilakarathne
###############################################################################

#//////////////////////////////////////////////////////////////////////////////////
setMethod("show", signature="PamChipAllPepResults", function(object){
             
          cat("-------- Results for set of Peptides --------\n")           
          cat("Class        :",class(object),"\n")
          cat("No. Peptides :",object@NoOfPeptides,"\n")
          cat("Total Run time (in seconds) :",object@RunTime,"\n")
          cat("Number of Groups :",object@n.groups,"\n")
          cat("Fitted Gamm Objects are saved as RData Objects in :\n",object@path,"\n")
          cat("-----------------------------------------\n")
          }
)


#//////////////////////////////////////////////////////////////////////////////////
setMethod("show", signature="CleanedPeptide", function(object){
             
          cat("-------- PreProcessed Peptide Data --------\n")           
          cat("Class        :",class(object),"\n")
          cat("Peptide Index:",object@p,"\n")
          cat("Cell Index   :",object@cel,"\n")
          tempt<-sort(unique((object@Profiles.Set.ID)[,c("Time")]))
          if (min(tempt)!=0) temp1<-tempt-min(tempt)
          cat("Unique Time points:",temp1,"\n")
          cat("-----------------------------------------\n")
          }
)

#//////////////////////////////////////////////////////////////////////////////////
setMethod("show", signature="CleanedPamData", function(object){
             
          cat("-------- PreProcessed List of Peptides --------\n")           
          cat("Class         :",class(object),"\n")
          cat("No. Peptides  :",object@np,"\n")
          cat("No. Cell lines:",object@nc,"\n")
          cat("Columns names for a certain peptide:\n",colnames(object@PeptideData[[1]]$Filttered),"\n")
          tempt<-sort(unique((object@PeptideData[[1]]$Filttered)[,c("Time")]))
          if (min(tempt)!=0) temp1<-tempt-min(tempt)
          cat("Unique Time points:",temp1,"\n")
          cat("PreProcessed Data frames are saved as RData Objects in \n:",object@path,"\n")
          cat("-----------------------------------------\n")
          }
)







#//////////////////////////////////////////////////////////////////////////////////
setMethod("show", signature="PamChipMixed", function(object){
         
          if(is.na(object@gammran3)) {Problem <- "Gamm Model failed to converge"   
          } else {
          Problem<-"Model is successfully fitted"
          }         
          cat("-------- Fitted PamGeneMix Model --------\n")           
          cat("Class:",class(object),"\n")
          cat("Note: ",Problem,"\n")
          
          if(!is.na(object@gammran3)) {
                  cat("Number of Groups: ",object@n.g,"\n")
                  cat("P-values for comparing group specific velocities \n")
                  cat("Test at t=0 : ",object@res.t1[1],"\n")
                  if(!is.na(object@test.at)) cat("Test at t=",object@test.at," : ",object@res.tany[1],"\n")
                  cat("Test at t=max(t) : ",object@res.tend[1],"\n")
                  cat("Test for entire profile : ",object@PValProfile[1],"\n")
                  }
          cat("-----------------------------------------\n")
          }
)
          

setGeneric("summary", function(object, ...) standardGeneric("summary"))
          
setMethod("summary", signature="PamChipMixed", function(object){
          if(is.na(object@gammran3)) {
                    cat("Model is failed to converge!\n")
          } else {
                  cat("--------- Estimated Initial Velocities (s.e.) ---- \n")
                   tp=1
                   ind.extract<-c(tp ,((tp + 1) + object@End.time),((tp + 2) + 2*object@End.time),((tp + 3) + 3*object@End.time))
                  for (i in 1:object@n.g) cat("group ",i," :",object@v[ind.extract[i]],"(",object@v.sd[ind.extract[i]],")\n")
                  cat("--------------------------------------------------------\n")
                  cat("-------- Fitted lme Object ---------\n")
                  cat("--------------------------------------------------------\n")
                  print(object@gammran3$lme)
                  cat("--------------------------------------------------------\n")
                  cat("-------- Fitted gam Object ---------\n")
                  cat("--------------------------------------------------------\n")
                  print(summary(object@gammran3$gam))
                  cat("--------------------------------------------------------\n")
          }
          
          }
)
