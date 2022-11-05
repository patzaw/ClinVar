library(here)
library(XML)
library(tidyverse)
source(here("scripts/clinVar-Functions.R"))

# https://stackoverflow.com/questions/22643580/combine-values-in-huge-xml-files

cvp <- function(progress=1000) {
   res <- new.env(parent=emptyenv())   # for results
   it <- 0L                            # iterator -- nodes visited
   list(
      ClinVarSet=function(elt) {
         it <<- it + 1L
         if (it %% progress == 0L){
            message(it)
         }
         cvsId <- xmlAttrs(elt)["ID"]
         res[[as.character(it)]] <- parseCvs(elt)
      },
      getres = function() {
         ## retrieve the 'res' environment when done
         res
      }
   )
}

branches <- cvp()

message(Sys.time())
xmlEventParse(here("sources/test.xml"), handlers=NULL, branches=branches)
message(Sys.time())
a <- branches$getres() %>% as.list()
message(Sys.time())

