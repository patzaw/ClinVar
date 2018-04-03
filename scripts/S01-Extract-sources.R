setwd("~/Shared/Data-Science/Data-Source-Model-Repository/clinVar/scripts/")

library(RJSONIO)
source("downloadSourceFiles.R")

desc <- readJSONStream("../DESCRIPTION.json")

sourceFiles <- desc$"source files"
urls <- unlist(lapply(
   sourceFiles,
   function(sf){
      toRet <- sf$"URL template"
      names(toRet) <- sf$"name"
      return(toRet)
   }
))
urls["ClinVarFullRelease.xml.gz"] <- sprintf(
   urls["ClinVarFullRelease.xml.gz"],
   format(Sys.Date(), "%Y-%m")
)
srcDir <- "../sources"

downloadSourceFiles(urls, srcDir)
