# Init file for package PamGeneMixed
.onLoad <- function(lib, pkg)
{

#library.dynam("PamGeneMixed", pkg, lib)



#    if ((.Platform$OS.type == "windows") && (.Platform$GUI ==
#        "Rgui") && interactive()) {
#        vigFile = system.file("Meta", "vignette.rds", package = "PamGeneMixed")
#        if (!file.exists(vigFile)) {
#            warning(sprintf("PamGeneMixed vignette is missing, nothing is added to the menu bar"))
#        }
#        else {
#            vigMtrx = readRDS(vigFile)
#            vigs = file.path(chartr("\\", "/", .find.package("PamGeneMixed")), "doc", vigMtrx[,
#                "PDF"])
#            names(vigs) = vigMtrx[, "Title"]
#            if (!"Vignettes" %in% winMenuNames())
#                winMenuAdd("Vignettes")
#            pkgMenu = paste("Vignettes", "PamGeneMixed", sep = "/")
#            winMenuAdd(pkgMenu)
#            for (i in seq(along = vigs)) winMenuAddItem(pkgMenu,
#                names(vigs)[i], paste("shell.exec(\"", vigs[i],
#                  "\")", sep = ""))
#        }
#    }


  packageStartupMessage("+----------------------------+                                          \n",
                        "**********  WELCOME **********                                          \n",
                        "**********    TO    **********                                          \n",
                        "********  PAMGENEMIXED *******                                          \n",
                        "+----------------------------+                                          \n")

    version <- packageDescription("PamGeneMixed",fields="Version")
    packageStartupMessage( "Citation: Thilakarathne P.J. et al.,","\n",
      "PamGeneMixed: Semi-parametric mixed models for PamChipData","\n",
      "Based on Bioinformatics 27(20):2859-2865, 2011","\n",
      "PamGeneMixed Package Version ", version, "\n")
}

#.onUnload <- function(libpath)
#{
#    library.dynam.unload("PamGeneMixed", libpath)
#}
