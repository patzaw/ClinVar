library(xml2)
library(future.apply)
library(tidyverse)

children_names <- function(node){
   xml2::xml_name(xml2::xml_children(node))
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
#' Extract and process ClinVarSet elements from an XML file
#' 
#' @param file XML file
#' @param n_max maximum number of lines to read from the file
#' @param by number of lines to read and process together
#' 
#' @details 
#' 
#' @return A list g.
#' 
process_ClinVarSets <- function(
      file, n_max=Inf, by=10^5
){
   tag <- "ClinVarSet"
   by <- min(by, n_max)
   con <- file(file, "r")
   encoding <- readLines(con, n=1) %>% 
      sub('^[<][?]xml .*encoding="', "", .) %>% 
      sub('".*$', "", .)
   n <- 1
   remaining <- c()
   toRet <- list()
   while(TRUE){
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
      batch_data <- future_lapply(
         elt_strings, function(x, encoding=encoding, ...){
            parse_ClinVarSet(read_xml(x, encoding=encoding))
         }
      )
      
      ## Combine tables ----
      toCombine <- c(
         "ReferenceClinVarAssertion", "rcvaAttributes", "observedIn",
         "measures", "measureNames", "measureAttributes",
         "measureCytogeneticLocations", "measureSequenceLocations",
         "measureXRef", "measureRelationships"
      )
      for(tc in toCombine){
         toRet[[tc]] <- dplyr::distinct(dplyr::bind_rows(
            toRet[[tc]],
            do.call(
               dplyr::bind_rows,
               lapply(batch_data, function(x) x[[tc]])
            )
         ))
      }
   }
   close(con)
   return(toRet)
}


###############################################################################@
#' Parse a ClinVarSet XML node
#' 
#' @param node the XML node to parse
#' 
#' @return A list with the following tibbles:
#' - ReferenceClinVarAssertion
#' - rcvaAttributes
#' - observedIn
#' 
parse_ClinVarSet <- function(node){
   
   ## Preparing the list to return ----
   toRet <- list()
   
   ## Just take current records ----
   if(xml2::xml_text(xml2::xml_child(node, "RecordStatus")) != "current"){
      invisible(NULL)
   }
   
   
   ## ReferenceClinVarAssertion ----
   
   ### ClinVarSet identifier and title ----
   cvsId <- xml2::xml_attr(node, "ID")
   cvsTitle <- xml2::xml_text(xml2::xml_child(node, "Title"))
   ClinVarSet <- dplyr::tibble(
      cvs=cvsId,
      title=cvsTitle
   )
   
   ### RCVA  ----
   rcvaOutput <- parse_ReferenceClinVarAssertion(
      xml2::xml_child(node, "ReferenceClinVarAssertion")
   )
   toRet$ReferenceClinVarAssertion <- dplyr::bind_cols(
      ClinVarSet, rcvaOutput$ReferenceClinVarAssertion
   )
   toRet$rcvaAttributes <- rcvaOutput$rcvaAttributes
   toRet$observedIn <- rcvaOutput$observedIn

   toRet <- c(toRet, rcvaOutput$measureSet)
   
   # toRet <- c(toRet, rcvaOutput$traitSet)
   # 
   # ## ClinVarAssertions ----
   # cvaOutput <- parse_ClinVarAssertions(node, cvsId)
   # toRet$ClinVarAssertions <- cvaOutput$cva
   # toRet$submitters <- cvaOutput$submitters
   # toRet$cvaObservedIn <- cvaOutput$observedIn
   
   ##
   return(toRet)
   
}


############################@
parse_ReferenceClinVarAssertion <- function(node){
   
   ## Preparing the list to return ----
   toRet <- list()
   
   ## ReferenceClinVarAssertion ----
   rcvaId <- xml2::xml_attr(node, "ID")
   rcvAcc <- xml2::xml_attr(xml2::xml_child(node, "ClinVarAccession"), "Acc")
   assertion <- xml2::xml_attr(xml2::xml_child(node, "Assertion"), "Type")
   if("ClinicalSignificance" %in% children_names(node)){
      csn <- xml2::xml_child(node, "ClinicalSignificance")
      rcvCSrev <- xml2::xml_text(xml2::xml_child(csn, "ReviewStatus"))
      rcvCSdes <- xml2::xml_text(xml2::xml_child(csn, "Description"))
      rcvCSexp <- xml2::xml_text(xml2::xml_child(csn, "Explanation"))
   }else{
      rcvCSrev <- rcvCSdes <- rcvCSexp <- NA
   }
   toRet$ReferenceClinVarAssertion <- dplyr::tibble(
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
   # toRet$traitSet <- parseTraitSet(node, rcvaId)
   
   ##
   return(toRet)
   
}

###############################################################################@
#' Parse a list of AttributeSet from an XML node
#' 
#' @param node the XML node to parse
#' 
#' @return A tibble
#' 
parse_AttributeSet <- function(
   node,
   attToGet=c("Type", "integerValue")
){
   toTake <- which(children_names(node)=="AttributeSet")
   if(length(toTake)==0){
      return(NULL)
   }else{
      toRet <- unique(do.call(dplyr::bind_rows, lapply(
         xml2::xml_children(node)[toTake],
         function(x){
            av <- xml2::xml_attrs(xml2::xml_child(x, "Attribute"))
            toRet <- dplyr::as_tibble(
               t(av[intersect(attToGet, names(av))])
            ) %>%
               dplyr::mutate(value=xml2::xml_text(xml2::xml_child(x, "Attribute")))
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
#' @return A tibble
#' 
parse_rcva_ObservedIn <- function(node){
   toTake <- which(children_names(node)=="ObservedIn")
   if(length(toTake)==0){
      return(NULL)
   }else{
      toRet <- unique(do.call(rbind, lapply(
         xml2::xml_children(node)[toTake],
         function(x){
            sample <- xml_child(x, "Sample")
            return(dplyr::tibble(
               origin= xml2::xml_text(xml2::xml_child(sample, "Origin")),
               taxonomyId=xml2::xml_attr(
                  xml_child(sample, "Species"),
                  "TaxonomyId"
               ),
               species= xml2::xml_text(xml2::xml_child(sample, "Species")),
               affectedStatus= xml2::xml_text(
                  xml2::xml_child(sample, "AffectedStatus")
               ),
               numberTested= xml2::xml_text(
                  xml2::xml_child(sample, "NumberTested")
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
#' @return A tibble
#' 
parse_MeasureSet <- function(node){
   if("MeasureSet" %in% children_names(node)){
      ms <- node %>% 
         xml2::xml_child("MeasureSet")
   }else{
      ms <- node %>% 
         xml2::xml_child("GenotypeSet") %>% 
         xml2::xml_child("MeasureSet")
   }
   toTake <- which(children_names(ms)=="Measure")
   measureSet <- lapply(
      xml2::xml_children(ms)[toTake],
      function(x){
         toRet <- list()
         
         ## Measure identification ----
         measure <- dplyr::as_tibble(
            t(xml2::xml_attrs(x)[c("Type", "ID")])
         )
         # colnames(measure) <- c("type", "id")
         toRet$measure <- measure
         
         ## Names of measure ----
         toTake <- which(children_names(x)=="Name")
         if(length(toTake)==0){
            measureNames <- NULL
         }else{
            measureNames <- unique(do.call(rbind, lapply(
               xml2::xml_children(x)[toTake],
               function(y){
                  ev <- xml2::xml_child(y, "ElementValue")
                  return(dplyr::tibble(
                     name=xml2::xml_text(ev),
                     type=xml2::xml_attr(ev, "Type")
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
         toTake <- which(children_names(x)=="CytogeneticLocation")
         if(length(toTake)==0){
            cytogeneticLocations <- NULL
         }else{
            cytogeneticLocations <- unique(do.call(rbind, lapply(
               xml2::xml_children(x)[toTake],
               function(y){
                  return(dplyr::tibble(
                     cytogenicLocation=xml2::xml_text(y)
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
   toRet$measures <- do.call(dplyr::bind_rows, lapply(
      measureSet,
      function(x) x$measure
   ))
   
   toRet$measureNames <- do.call(dplyr::bind_rows, lapply(
      measureSet,
      function(x) x$measureNames
   ))
   
   toRet$measureAttributes <- do.call(dplyr::bind_rows, lapply(
      measureSet,
      function(x) x$measureAttributes
   ))
   
   toRet$measureCytogeneticLocations <- do.call(dplyr::bind_rows, lapply(
      measureSet,
      function(x) x$cytogeneticLocations
   ))
   
   toRet$measureSequenceLocations <- do.call(dplyr::bind_rows, lapply(
      measureSet,
      function(x) x$sequenceLocations
   ))
   
   toRet$measureXRef <- do.call(dplyr::bind_rows, lapply(
      measureSet,
      function(x) x$xref
   ))
   
   toRet$measureRelationships <- do.call(dplyr::bind_rows, lapply(
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
#' @return A tibble
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
   toTake <- which(children_names(node)=="SequenceLocation")
   if(length(toTake)==0){
      return(NULL)
   }else{
      toRet <- unique(do.call(dplyr::bind_rows, lapply(
         xml2::xml_children(node)[toTake],
         function(x){
            xat <- xml2::xml_attrs(x)
            toRet <- dplyr::as_tibble(
               t(xat[intersect(attToGet, names(xat))])
            )
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
#' @return A tibble
#'
parse_XRef <- function(node, attToGet=c("DB", "ID", "Type")){
   toTake <- which(children_names(node)=="XRef")
   if(length(toTake)==0){
      return(NULL)
   }else{
      toRet <- unique(do.call(rbind, lapply(
         xml2::xml_children(node)[toTake],
         function(x){
            xat <- xml2::xml_attrs(x)
            toRet <- as.data.frame(
               t(xat[attToGet]),
               stringsAsFactors=FALSE
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
#' @return A tibble
#'
parse_MeasureRelationship <- function(node){
   toTake <- which(children_names(node)=="MeasureRelationship")
   if(length(toTake)==0){
      return(NULL)
   }else{
      toRet <- unique(do.call(rbind, lapply(
         xml2::xml_children(node)[toTake],
         function(x){
            toRet <- data.frame(
               id="1",
               type=xml2::xml_attr(x, "Type")
            )
            ##
            toTake <- which(children_names(x)=="Name")
            if(length(toTake)==0){
               toAdd <- data.frame(
                  id=1,
                  name=NA
                  # name.type=NA
               )
               # toAdd <- NULL
            }else{
               toAdd <- unique(do.call(rbind, lapply(
                  xml2::xml_children(x)[toTake],
                  function(y){
                     ev <- xml2::xml_child(y, "ElementValue")
                     return(data.frame(
                        id="1",
                        name=xml2::xml_text(ev)
                        # name.type=xml2::xml_attr(ev, "Type"),
                     ))
                  }
               )))
            }
            # toRet <- tidyr::crossing(toRet, toAdd)
            # toRet <- dplyr::full_join(toRet, toAdd, by="id")
            toRet <- merge(toRet, toAdd, by="id", all=T)
            ##
            toTake <- which(children_names(x)=="Symbol")
            if(length(toTake)==0){
               toAdd <- data.frame(
                  id="1",
                  symbol=NA,
                  # symbol.type=NA
               )
               # toAdd <- NULL
            }else{
               toAdd <- unique(do.call(rbind, lapply(
                  xml2::xml_children(x)[toTake],
                  function(y){
                     ev <- xml2::xml_child(y, "ElementValue")
                     return(data.frame(
                        id="1",
                        symbol=xml2::xml_text(ev)
                        # symbol.type=xml2::xml_attr(ev, "Type")
                     ))
                  }
               )))
            }
            # toRet <- tidyr::crossing(toRet, toAdd)
            # toRet <- dplyr::full_join(toRet, toAdd, by="id")
            toRet <- merge(toRet, toAdd, by="id", all=T)
            toAdd <- parse_XRef(x) %>% 
               dplyr::mutate(id="1")
               # dplyr::filter(DB=="Gene") %>% 
               # select(ID)
            # colnames(toAdd) <- c("xref.db", "xref.id", "xref.type")
            # toRet <- tidyr::crossing(toRet, toAdd)
            # toRet <- dplyr::full_join(toRet, toAdd, by="id") %>% 
            toRet <- merge(toRet, toAdd, by="id", all=T)
            return(toRet)
         }
      )))
      return(toRet %>% dplyr::select(-"id") %>% as_tibble())
   }
}
