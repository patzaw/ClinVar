library(XML)
library(future.apply)
library(tidyverse)

CVS_TABLES <- c(
   "ReferenceClinVarAssertion", "rcvaAttributes", "observedIn",
   "measures", "measureNames", "measureAttributes",
   "measureCytogeneticLocations", "measureSequenceLocations",
   "measureXRef", "measureRelationships",
   "traits", "traitXRef", "traitNames",
   "ClinVarAssertions", "submitters", "cvaObservedIn"
)


###############################################################################@
#' Remove trailing spaces
#' 
#' @param x a character vector to clean
#' 
#' @return A character vector
#' 
rm_tr_spaces <- function(x){
   stringr::str_remove(x, stringr::regex("^ *")) %>%
      stringr::str_remove(stringr::regex(" *$"))
}

###############################################################################@
#' Extract non self nested elements from an XML connection
#' 
#' @param con connection to an XML file
#' @param tag the tag to extract
#' @param n the number of lines to read and
#' from which elements will be extracted. Negative values indicate that one
#' should read up to the end of input on the connection.
#' @param formerLines a character vector of lines to append in front of
#' those taken from the connection before extracting elements.
#' 
#' @details If elements are self nested in the document (e.g. <p><p></p></p>),
#' the function will not behave correctly. This function is intended to be used
#' on high level elements which can be seen as a list. For example:
#' <item id="1">...</item><item id="2">...</item><item id="3">...</item>...
#' 
#' @return A character vector in which each value corresponds to an extracted
#' XML element. The attribute "remaining" contains lines at the end of the 
#' connection chunk from which no elements could be extracted. Those lines
#' should be provided in the `formerLines` parameters during the following
#' call to this function. The attribute "read" contains the number of lines
#' that were actually read from the connection.
#' 
extract_nsn_xml_elements <- function(
   con, tag, n=-1L, formerLines=character(0)
){
   lines <- readLines(con, n=n)
   n <- length(lines)
   if(n==0){
      toRet <- character(0)
      attr(toRet, "remaining") <- formerLines
      attr(toRet, "n") <- n
      return(toRet)
   }
   lines <- c(formerLines, lines)
   
   
   stag <- sprintf("<%s[ >]", tag)
   etag <- sprintf("<[/]%s>", tag)
   starts <- grep(stag, lines)
   if(length(starts)==0){
      toRet <- character(0)
      attr(toRet, "remaining") <- character(0)
      attr(toRet, "n") <- n
      return(toRet)
   }
   ends <- grep(etag, lines)
   if(length(ends)==0){
      toRet <- character(0)
      attr(toRet, "remaining") <- c(formerLines, lines)
      attr(toRet, "n") <- n
      return(toRet)
   }
   
   starts <- starts[1:length(ends)]
   if(max(ends) < length(lines)){
      remaining <- lines[(max(ends)+1):length(lines)]
   }else{
      remaining <- character(0)
   }
   toRet <- apply(
      cbind(starts, ends), 1,
      function(x){
         sl <- lines[x[1]:x[2]]
         sl <- sub("^[[:blank:]]*", "", sl)
         sl[1] <- sub("^[^<]*", "", sl[1])
         sl[length(sl)] <- sub("[^>]*$", "", sl[length(sl)])
         paste(sl, collapse="\n")
      }
   )
   attr(toRet, "remaining") <- remaining
   attr(toRet, "n") <- n
   return(toRet)
}

###############################################################################@
#' Combine tables extracted form processed ClinVarSet elements from an XML file
#' 
#' @param dir a path to a directory where to read the data
#' 
#' @details 
#' 
#' @return A list of CVS_TABLES tibbles
#'
combine_ClinVarSets <- function(dir){
   combine_cvs_table <- function(dir, tn){
      toread <- list.files(dir, pattern=sprintf("^batch-.*-%s[.]rds$", tn))
      toreadi <- toread %>%
         stringr::str_remove(sprintf("-%s[.]rds$", tn)) %>% 
         stringr::str_remove("^batch-") %>% 
         as.numeric()
      toread <- toread[order(toreadi)]
      toreadi <- toreadi[order(toreadi)]
      lapply(file.path(dir, toread), readRDS) %>% 
         do.call(dplyr::bind_rows, .) %>% 
         dplyr::distinct() %>% 
         dplyr::as_tibble()
   }
   toRet <- list()
   for(tn in CVS_TABLES){
      toRet[[tn]] <- combine_cvs_table(dir, tn)
   }
   return(toRet)
}

###############################################################################@
#' Extract and process ClinVarSet elements from an XML file
#' 
#' @param file XML file
#' @param n_max maximum number of lines to read from the file
#' @param by number of lines to read and process together
#' 
#' @details 
#' 
#' @return A path to a directory where data are stored
#' 
process_ClinVarSets <- function(
      file, location, n_max=Inf, by=10^5, cores=15
){
   odir <- file.path(
      location,
      sprintf("CVS-batches-%s", format(Sys.time(), "%Y-%m-%d_%H%M%S"))
   )
   stopifnot(!file.exists(odir))
   dir.create(odir)
   toRet <- odir
   tag <- "ClinVarSet"
   by <- min(by, n_max)
   con <- file(file, "r")
   encoding <- readLines(con, n=1) %>% 
      sub('^[<][?]xml .*encoding="', "", .) %>% 
      sub('".*$', "", .)
   n <- 1
   remaining <- c()
   i <- 0
   while(TRUE){
      i <- i+1
      message(n, " lines read")
      if(n >= n_max){
         break()
      }
      elt_strings <- extract_nsn_xml_elements(
         con, tag=tag, n=by, formerLines=remaining
      )
      n <- n + attr(elt_strings, "n")
      if(
         length(elt_strings)==0 &&
         identical(attr(elt_strings, "remaining"), remaining)
      ){
         break()
      }
      remaining <- attr(elt_strings, "remaining")
      
      ## Parse nodes ----
      plan(multisession, workers=cores)
      on.exit(plan(sequential))
      batch_data <- future_lapply(
         elt_strings,
         function(x, e=encoding){
            parse_ClinVarSet(
               XML::xmlRoot(XML::xmlParse(x, encoding=e, asText=TRUE))
            )
         },
         future.globals=c(
            "parse_AttributeSet", "parse_ClinVarSet", 
            "parse_MeasureRelationship", "parse_MeasureSet", "parse_rcva_ObservedIn", 
            "parse_ReferenceClinVarAssertion", "parse_SequenceLocation", 
            "parse_XRef", "parse_TraitSet",
            "parse_ClinVarAssertions", "parse_cva_ObservedIn"
         )
      )
      plan(sequential)
      for(tn in CVS_TABLES){
         saveRDS(
            lapply(batch_data, function(x)x[[tn]]) %>%
               do.call(dplyr::bind_rows, .),
            file=file.path(odir, sprintf("batch-%s-%s.rds", i, tn)),
            compress=FALSE
         )
      }
      rm(batch_data)
      gc()

   }
   close(con)
   return(toRet)
}


###############################################################################@
#' Parse a ClinVarSet XML node
#' 
#' @param node the XML node to parse
#' 
#' @return A list of `CVS_TABLES` data.frames
#' 
parse_ClinVarSet <- function(node){
   
   ## Preparing the list to return ----
   toRet <- list()
   
   ## Just take current records ----
   if(XML::xmlValue(node[["RecordStatus"]]) != "current"){
      invisible(NULL)
   }
   
   
   ## ReferenceClinVarAssertion ----
   
   ### ClinVarSet identifier and title ----
   cvsId <- XML::xmlAttrs(node)["ID"]
   cvsTitle <- XML::xmlValue(node[["Title"]])
   ClinVarSet <- data.frame(
      cvs=cvsId,
      title=cvsTitle
   )
   
   ### RCVA  ----
   rcvaOutput <- parse_ReferenceClinVarAssertion(
      node[["ReferenceClinVarAssertion"]]
   )
   toRet$ReferenceClinVarAssertion <- cbind(
      ClinVarSet, rcvaOutput$ReferenceClinVarAssertion
   )
   toRet$rcvaAttributes <- rcvaOutput$rcvaAttributes
   toRet$observedIn <- rcvaOutput$observedIn

   toRet <- c(toRet, rcvaOutput$measureSet)
   toRet <- c(toRet, rcvaOutput$traitSet)

   ### ClinVarAssertions ----
   cvaOutput <- parse_ClinVarAssertions(node)
   toRet <- c(toRet, cvaOutput)
   toRet$ClinVarAssertions$cvs <- cvsId
   
   ##
   return(toRet)
   
}

###############################################################################@
#' Parse a ReferenceClinVarAssertion XML node
#' 
#' @param node the XML node to parse
#' 
#' @return A data.frame
#' 
parse_ReferenceClinVarAssertion <- function(node){
   
   ## Preparing the list to return ----
   toRet <- list()
   
   ## ReferenceClinVarAssertion ----
   rcvaId <- XML::xmlAttrs(node)["ID"]
   rcvAcc <- XML::xmlAttrs(node[["ClinVarAccession"]])["Acc"]
   assertion <- XML::xmlAttrs(node[["Assertion"]])["Type"]
   if("ClinicalSignificance" %in% names(node)){
      csn <- node[["ClinicalSignificance"]]
      rcvCSrev <- XML::xmlValue(csn[["ReviewStatus"]])
      rcvCSdes <- XML::xmlValue(csn[["Description"]])
      rcvCSexp <- XML::xmlValue(csn[["Explanation"]])
   }else{
      rcvCSrev <- rcvCSdes <- rcvCSexp <- NA
   }
   toRet$ReferenceClinVarAssertion <- data.frame(
      id=rcvaId,
      accession=rcvAcc,
      assertion=assertion,
      reviewStatus=rcvCSrev,
      clinicalSignificance=rcvCSdes,
      explanation=rcvCSexp
   )
   
   ## ReferenceClinVarAssertion AttributeSets ----
   rcvaAttributes <- parse_AttributeSet(node, c("Type", "integerValue"))
   if(!is.null(rcvaAttributes)){
      rcvaAttributes$rcvaId <- rcvaId
   }
   toRet$rcvaAttributes <- rcvaAttributes
   
   ## Observed In ----
   observedIn <- parse_rcva_ObservedIn(node)
   if(!is.null(observedIn)){
      observedIn$rcvaId <- rcvaId
   }
   toRet$observedIn <- observedIn

   ## MeasureSet ----
   toRet$measureSet <- parse_MeasureSet(node)
   toRet$measureSet$measures$rcvaId <- rcvaId

   # ## TraitSet ----
   toRet$traitSet <- parse_TraitSet(node)
   toRet$traitSet$traits$rcvaId <- rcvaId
   
   ##
   return(toRet)
   
}

###############################################################################@
#' Parse a list of AttributeSet from an XML node
#' 
#' @param node the XML node to parse
#' 
#' @return A data.frame
#' 
parse_AttributeSet <- function(
   node,
   attToGet=c("Type", "integerValue")
){
   toTake <- which(names(node)=="AttributeSet")
   if(length(toTake)==0){
      return(NULL)
   }else{
      toRet <- unique(do.call(rbind, lapply(
         node[toTake],
         function(x){
            av <- XML::xmlAttrs(x[["Attribute"]])
            toRet <- as.data.frame(
               t(av[attToGet])
            )
            colnames(toRet) <- attToGet
            toRet$value <- XML::xmlValue(x[["Attribute"]])
            return(toRet)
         }
      )))
      return(toRet)
   }
}

###############################################################################@
#' Parse a list of ObservedIn from a ReferenceClinVarAssertion XML node
#' 
#' @param node the XML node to parse
#' 
#' @return A data.frame
#' 
parse_rcva_ObservedIn <- function(node){
   toTake <- which(names(node)=="ObservedIn")
   if(length(toTake)==0){
      return(NULL)
   }else{
      toRet <- unique(do.call(rbind, lapply(
         node[toTake],
         function(x){
            sample <- x[["Sample"]]
            return(data.frame(
               origin= XML::xmlValue(sample[["Origin"]]),
               taxonomyId=if(is.null(
                  XML::xmlAttrs(sample[["Species"]])["TaxonomyId"]
               )){
                  "XML error"
               }else{
                  XML::xmlAttrs(sample[["Species"]])["TaxonomyId"]
               },
               species= XML::xmlValue(sample[["Species"]]),
               affectedStatus= XML::xmlValue(
                  sample[["AffectedStatus"]]
               ),
               numberTested= XML::xmlValue(
                  sample[["NumberTested"]]
               )
            ))
         }
      )))
      return(toRet)
   }
}

###############################################################################@
#' Parse MeasureSet element from a ReferenceClinVarAssertion XML node
#' 
#' @param node the XML node to parse
#' 
#' @return A data.frame
#' 
parse_MeasureSet <- function(node){
   if("MeasureSet" %in% names(node)){
      ms <- node[["MeasureSet"]]
   }else{
      ms <- node[["GenotypeSet"]][["MeasureSet"]]
   }
   toTake <- which(names(ms)=="Measure")
   measureSet <- lapply(
      ms[toTake],
      function(x){
         toRet <- list()
         
         ## Measure identification ----
         measure <- as.data.frame(
            t(XML::xmlAttrs(x)[c("Type", "ID")])
         )
         # colnames(measure) <- c("type", "id")
         toRet$measure <- measure
         
         ## Names of measure ----
         toTake <- which(names(x)=="Name")
         if(length(toTake)==0){
            measureNames <- NULL
         }else{
            measureNames <- unique(do.call(rbind, lapply(
               x[toTake],
               function(y){
                  ev <- y[["ElementValue"]]
                  return(data.frame(
                     name=XML::xmlValue(ev),
                     type=XML::xmlAttrs(ev)["Type"]
                  ))
               }
            )))
            measureNames$measureId <- measure$ID
         }
         toRet$measureNames <- measureNames
         
         ## Attributes ----
         toRet$measureAttributes <- parse_AttributeSet(
            x,
            attToGet=c("Type", "integerValue", "Change")
         )
         if(!is.null(toRet$measureAttributes)){
            toRet$measureAttributes$measureId <- measure$ID
         }
         
         ## Cytogenetic locations ----
         toTake <- which(names(x)=="CytogeneticLocation")
         if(length(toTake)==0){
            cytogeneticLocations <- NULL
         }else{
            cytogeneticLocations <- unique(do.call(rbind, lapply(
               x[toTake],
               function(y){
                  return(data.frame(
                     cytogenicLocation=XML::xmlValue(y)
                  ))
               }
            )))
            cytogeneticLocations$measureId <- measure$ID
         }
         toRet$cytogeneticLocations <- cytogeneticLocations
         
         ## Sequence locations ----
         sequenceLocations <- parse_SequenceLocation(x)
         if(!is.null(sequenceLocations)){
            sequenceLocations$measureId <- measure$ID
         }
         toRet$sequenceLocations <- sequenceLocations
         
         ## XREF ----
         xref <- parse_XRef(x)
         if(!is.null(xref)){
            xref$measureId <- measure$ID
         }
         toRet$xref <- xref
         
         ## MeasureRelationships ----
         measureRelationships <- parse_MeasureRelationship(x)
         if(!is.null(measureRelationships)){
            measureRelationships$measureId <- measure$ID
         }
         toRet$measureRelationships <- measureRelationships
         
         
         ## Return information
         return(toRet)
      }
   )
   
   ## Post processing
   toRet <- list()
   toRet$measures <- do.call(rbind, lapply(
      measureSet,
      function(x) x$measure
   ))
   
   toRet$measureNames <- do.call(rbind, lapply(
      measureSet,
      function(x) x$measureNames
   ))
   
   toRet$measureAttributes <- do.call(rbind, lapply(
      measureSet,
      function(x) x$measureAttributes
   ))
   
   toRet$measureCytogeneticLocations <- do.call(rbind, lapply(
      measureSet,
      function(x) x$cytogeneticLocations
   ))
   
   toRet$measureSequenceLocations <- do.call(rbind, lapply(
      measureSet,
      function(x) x$sequenceLocations
   ))
   
   toRet$measureXRef <- do.call(rbind, lapply(
      measureSet,
      function(x) x$xref
   ))
   
   toRet$measureRelationships <- do.call(rbind, lapply(
      measureSet,
      function(x) x$measureRelationships
   ))
   
   ##
   return(toRet)
}

###############################################################################@
#' Parse SequenceLocation element from a MeasureSet XML node
#' 
#' @param node the XML node to parse
#' 
#' @return A data.frame
#' 
parse_SequenceLocation <- function(
      node,
      attToGet=c(
         "Accession", "alternateAllele", "Assembly",
         "AssemblyAccessionVersion", "AssemblyStatus", "Chr",
         "display_start", "display_stop", "innerStart",
         "innerStop", "outerStart", "outerStop",
         "referenceAllele", "start", "stop",
         "Strand", "variantLength"
      )
){
   toTake <- which(names(node)=="SequenceLocation")
   if(length(toTake)==0){
      return(NULL)
   }else{
      toRet <- unique(do.call(rbind, lapply(
         node[toTake],
         function(x){
            xat <- XML::xmlAttrs(x)
            toRet <- as.data.frame(
               t(xat[attToGet])
            )
            colnames(toRet) <- attToGet
            return(toRet)
         }
      )))
      return(toRet)
   }
}

###############################################################################@
#' Parse XRef element from a MeasureSet XML node
#' 
#' @param node the XML node to parse
#' 
#' @return A data.frame
#'
parse_XRef <- function(node, attToGet=c("DB", "ID", "Type")){
   toTake <- which(names(node)=="XRef")
   if(length(toTake)==0){
      return(NULL)
   }else{
      toRet <- unique(do.call(rbind, lapply(
         node[toTake],
         function(x){
            xat <- XML::xmlAttrs(x)
            toRet <- as.data.frame(
               t(xat[attToGet])
            )
            colnames(toRet) <- attToGet
            return(toRet)
         }
      )))
      return(toRet)
   }
}

###############################################################################@
#' Parse XRef element from a MeasureSet XML node
#' 
#' @param node the XML node to parse
#' 
#' @return A data.frame
#'
parse_MeasureRelationship <- function(node){
   toTake <- which(names(node)=="MeasureRelationship")
   if(length(toTake)==0){
      return(NULL)
   }else{
      toRet <- unique(do.call(rbind, lapply(
         node[toTake],
         function(x){
            toRet <- data.frame(
               id="1",
               type=XML::xmlAttrs(x)["Type"]
            )
            ##
            toTake <- which(names(x)=="Name")
            if(length(toTake)==0){
               toAdd <- data.frame(
                  id="1",
                  name=NA,
                  # name.type=NA
               )
               # toAdd <- NULL
            }else{
               toAdd <- unique(do.call(rbind, lapply(
                  x[toTake],
                  function(y){
                     ev <- y[["ElementValue"]]
                     return(data.frame(
                        id="1",
                        name=XML::xmlValue(ev)
                        # name.type=XML::xmlAttrs(ev)["Type"]
                     ))
                  }
               )))
            }
            # toRet <- tidyr::crossing(toRet, toAdd)
            toRet <- merge(toRet, toAdd, by="id", all=T)
            ##
            toTake <- which(names(x)=="Symbol")
            if(length(toTake)==0){
               toAdd <- data.frame(
                  id="1",
                  symbol=NA,
                  # symbol.type=NA
               )
               # toAdd <- NULL
            }else{
               toAdd <- unique(do.call(rbind, lapply(
                  x[toTake],
                  function(y){
                     ev <- y[["ElementValue"]]
                     return(data.frame(
                        id="1",
                        symbol=XML::xmlValue(ev)
                        # symbol.type=XML::xmlAttrs(ev)["Type"]
                     ))
                  }
               )))
            }
            # toRet <- tidyr::crossing(toRet, toAdd)
            toRet <- merge(toRet, toAdd, by="id", all=T)
            toAdd <- parse_XRef(x)
            toAdd$id <- "1"
               # dplyr::filter(DB=="Gene") %>% 
               # select(ID)
            # toRet <- tidyr::crossing(toRet, toAdd)
            toRet <- merge(toRet, toAdd, by="id", all=T)
            return(toRet)
         }
      )))
      toRet <- dplyr::select(toRet, -"id")
      toRet <- toRet[which(toRet$DB=="Gene"),]
      return(toRet)
   }
}

###############################################################################@
#' Parse TraitSet element from a ReferenceClinVarAssertion XML node
#' 
#' @param node the XML node to parse
#' 
#' @return A data.frame
#' 
parse_TraitSet <- function(node){
   ts <- node[["TraitSet"]]
   toTake <- which(names(ts)=="Trait")
   traitSet <- lapply(
      ts[toTake],
      function(x){
         toRet <- list()
         ## Trait identification
         trait <- as.data.frame(
            t(XML::xmlAttrs(x)[c("Type", "ID")])
         )
         colnames(trait) <- c("type", "id")
         toRet$trait <- trait
         ## Direct xref
         xref <- parse_XRef(x)
         if(!is.null(xref)){
            xref$trait.id <- trait$id
         }
         toRet$xref <- xref
         ## Names
         toTake <- which(names(x)=="Name")
         detailedNames <- lapply(
            x[toTake],
            function(y){
               toRet <- list()
               toRet$name <- data.frame(
                  trait.id=trait$id,
                  trait.name.type=XML::xmlAttrs(y[["ElementValue"]])["Type"],
                  trait.name=XML::xmlValue(y[["ElementValue"]])
               )
               toRet$xref <- parse_XRef(y)
               return(toRet)
            }
         )
         toRet$traitNames <- unique(do.call(rbind, lapply(
            detailedNames,
            function(y){
               y$name
            }
         )))
         xrefToAdd <- unique(do.call(rbind, lapply(
            detailedNames,
            function(y){
               y$xref
            }
         )))
         if(!is.null(xrefToAdd)){
            xrefToAdd$trait.id <- trait$id
         }
         toRet$xref <- rbind(
            toRet$xref,
            xrefToAdd
         )
         ## Symbols
         toTake <- which(names(x)=="Symbol")
         xrefToAdd <- unique(do.call(rbind, lapply(
            x[toTake],
            parse_XRef
         )))
         if(!is.null(xrefToAdd)){
            xrefToAdd$trait.id <- trait$id
         }
         toRet$xref <- rbind(
            toRet$xref,
            xrefToAdd
         )
         
         ## Return information
         return(toRet)
      }
   )
   
   ## Post processing
   toRet <- list()
   toRet$traits <- do.call(rbind, lapply(
      traitSet,
      function(x) x$trait
   ))
   
   toRet$traitXRef <- do.call(rbind, lapply(
      traitSet,
      function(x) x$xref
   ))
   
   toRet$traitNames <- do.call(rbind, lapply(
      traitSet,
      function(x) x$traitNames
   ))
   
   ##
   return(toRet)
   
}

###############################################################################@
#' Parse a ClinVarAssertion elements from an XML node
#' 
#' @param node the XML node to parse
#' 
#' @return A list of data.frames
#' 
parse_ClinVarAssertions <- function(node){
   toTake <- which(names(node)=="ClinVarAssertion")
   toRet <- lapply(
      node[toTake],
      function(x){
         ClinVarAssertions <- data.frame(
            id=XML::xmlAttrs(x)["ID"],
            accession=XML::xmlAttrs(x[["ClinVarAccession"]])["Acc"],
            clinicalSignificance=XML::xmlValue(x[["ClinicalSignificance"]][["Description"]])
         )
         submitters <- data.frame(
            submitter=XML::xmlAttrs(x[["ClinVarSubmissionID"]])["submitter"],
            cvaId=ClinVarAssertions$id,
            primary=TRUE
         )
         if("AdditionalSubmitters" %in% names(x)){
            toAdd <- unlist(lapply(
               x[["AdditionalSubmitters"]][1:length(names(x[["AdditionalSubmitters"]]))],
               function(y){
                  return(XML::xmlAttrs(y)["SubmitterName"])
               }
            ))
            toAdd <- data.frame(
               submitter=toAdd,
               cvaId=ClinVarAssertions$id,
               primary=FALSE
            )
            submitters=rbind(submitters, toAdd)
         }
         cvaObservedIn <- parse_cva_ObservedIn(x)
         if(!is.null(cvaObservedIn)){
            cvaObservedIn$cvaId <- ClinVarAssertions$id
         }
         toRet <- list(
            ClinVarAssertions=ClinVarAssertions,
            submitters=submitters,
            cvaObservedIn=cvaObservedIn
         )
         return(toRet)
      }
   )
   ClinVarAssertions <- do.call(rbind, lapply(
      toRet,
      function(x) x$ClinVarAssertions
   ))
   submitters <- do.call(rbind, lapply(
      toRet,
      function(x) x$submitters
   ))
   cvaObservedIn <- do.call(rbind, lapply(
      toRet,
      function(x) x$cvaObservedIn
   ))
   return(list(
      ClinVarAssertions=ClinVarAssertions,
      submitters=submitters,
      cvaObservedIn=cvaObservedIn
   ))
}

###############################################################################@
#' Parse a list of ObservedIn from a ClinVarAssertion XML node
#' 
#' @param node the XML node to parse
#' 
#' @return A data.frame
#' 
parse_cva_ObservedIn <- function(node){
   toTake <- which(names(node)=="ObservedIn")
   if(length(toTake)==0){
      return(NULL)
   }else{
      toRet <- unique(do.call(rbind, lapply(
         node[toTake],
         function(x){
            return(data.frame(
               origin=XML::xmlValue(x[["Sample"]][["Origin"]]),
               species=XML::xmlValue(x[["Sample"]][["Species"]]),
               affectedStatus=XML::xmlValue(x[["Sample"]][["AffectedStatus"]])
            ))
         }
      )))
      return(toRet)
   }
}
