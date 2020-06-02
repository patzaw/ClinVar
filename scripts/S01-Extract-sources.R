library(RJSONIO)
library(here)
source(here("scripts/downloadSourceFiles.R"))

desc <- readJSONStream(here("DESCRIPTION.json"))

sourceFiles <- desc$"source files"
urls <- unlist(lapply(
   sourceFiles,
   function(sf){
      toRet <- sf$"URL template"
      names(toRet) <- sf$"name"
      return(toRet)
   }
))
urls["ClinVarFullRelease.xml.gz"] <- sub(
   "YYYY-MM",
   format(Sys.Date(), "%Y-%m"),
   urls["ClinVarFullRelease.xml.gz"]
)
srcDir <- here("sources")

downloadSourceFiles(urls, srcDir)
