#
#
# Author: Pushpike Thilakarathne
###############################################################################





#--------set class for ouput of AutoPamGenMix -----------
setClass("PamChipAllPepResults", representation = representation(
           AllPepRes = "list",
           RunTime="numeric",
           NoOfPeptides="numeric",
           TestAt="vector",
           n.groups="numeric",
           path="character")
           )
           
setValidity("PamChipAllPepResults",
    function(object)
    {
        if (!is.list(slot(object, "AllPepRes")))
        {
            return("slot >AllPepRes< must be a list !")
        } else if  (!is.numeric(slot(object, "RunTime")))   {
            return("slot >RunTime< must be numeric !")
        } else if  (!is.numeric(slot(object, "NoOfPeptides")))   {
            return("slot >NoOfPeptides< must be numeric !")
        }  else if  (!is.vector(slot(object, "TestAt")))   {
            return("slot >TestAt< must be a vector !")
        }  else if  (!is.numeric(slot(object, "n.groups")))   {
            return("slot >n.groups< must be numeric !")
        }  else if  (!is.character(slot(object, "path")))   {
            return("slot >path< must be a character !")
        } 
    }
   )


#--------set class for ouput of PreProcessPam -----------
setClass("CleanedPeptide", representation = representation(
           Profiles.Set.ID = "data.frame",
           p="numeric",
           cel="numeric")
           )
           
setValidity("CleanedPeptide",
    function(object)
    {
        if (!is.data.frame(slot(object, "Profiles.Set.ID")))
        {
            return("slot >Profiles.Set.ID< must be a data frame !")
        } else if  (!is.numeric(slot(object, "p")))   {
            return("slot >p< must be numeric !")
        } else if  (!is.numeric(slot(object, "cel")))   {
            return("slot >cel< must be numeric !")
        } 
    }
   )

#--------set class for ouput of PreProcessAllPeptides -----------
setClass("CleanedPamData", representation = representation(
           PeptideData = "list",
           np="numeric",
           nc="numeric",
           path="character")
           )
           
setValidity("CleanedPamData",
    function(object)
    {
        if (!is.list(slot(object, "PeptideData")))
        {
            return("slot >PeptideData< must be a list !")
        } else if  (!is.numeric(slot(object, "np")))   {
            return("slot >np< must be numeric !")
        } else if  (!is.numeric(slot(object, "nc")))   {
            return("slot >nc< must be numeric !")
        } else if  (!is.character(slot(object, "path")))   {
            return("slot >path< must be a character !")
        } 
    }
   )



#--------set class for ouput of PamGeneMix -----------
setClass("PamChipMixed",
         representation = representation(
            res.t1  = "vector",
           res.tend = "vector",
           res.tany = "vector",
           test.at="vector",
           PValProfile = "vector",
           t.mesh = "vector",
           End.time = "numeric",
           n.g = "numeric",
           gnames = "vector",
           v = "vector",
           v.sd = "vector",
           Vall = "matrix",
           yfitfixef = "matrix",
           yfit = "numeric",
           gammran3 = "ANY",
           PTx = "data.frame")
           )


setValidity("PamChipMixed",
    function(object)
    {
        if (!is.vector(slot(object, "res.t1")))
        {
            return("slot >res.t1< must be a vector!")
        }
        else if (!is.vector(slot(object, "res.tend")))
        {
            return("slot >res.tend< must be a vector!")
        }
        else if (!is.vector(slot(object, "res.tany")))
        {
            return("slot >res.tany< must be a vector!")
        } else if (!is.vector(slot(object, "test.at")))
        {
            return("slot >test.at< must be a vector!")
        }        
        else if (!is.vector(slot(object, "PValProfile")))
        {
            return("slot >PValProfile< must be a vector!")
        }
        else if (!is.vector(slot(object, "t.mesh")))
        {
            return("slot >t.mesh< must be a vector!")
        }
        else if (!is.numeric(slot(object, "n.g")))
        {
            return("slot >n.g< must be a numeric!")
        }
        else if (!is.vector(slot(object, "gnames")))
        {
            return("slot >gnames< must be a vector!")
        }
        else if (!is.vector(slot(object, "v")))
        {
            return("slot >v< must be a v.sd!")
        }
        else if (!is.vector(slot(object, "v.sd")))
        {
            return("slot >v.sd< must be a vector!")
        }
        else if (!is.matrix(slot(object, "Vall")))
        {
            return("slot >Vall< must be a matrix!")
        }
        else if (!is.matrix(slot(object, "yfitfixef")))
        {
            return("slot >yfitfixef< must be a matrix!")
        }
        else if (!is.numeric(slot(object, "yfit")))
        {
            return("slot >yfit< must be a numeric vector!")
        }
#        else if (!is.list(slot(object, "gammran3")))
#        {
#            return("slot >gammran3< must be a list!")
#        }
        else if (!is.data.frame(slot(object, "PTx")))
        {
            return("slot >PTx< must be a data frame!")
        }
       
    }

 )
