library(xml2)
library(future.apply)
library(tidyverse)

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
#' Extract and process elements from an XML file
#' 
#' @param file XML file
#' @param tag the tag to extract
#' @param n_max maximum number of lines to read
#' @param by number of lines to read and process together
#' @param fun function to process XML elements
#' 
#' @details 
#' 
#' @return A list g.
#' 
process_xml_elements <- function(
      file, tag, n_max=Inf, by=10^5,
      fun=xml2::as_list
){
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
      toRet <- c(
         toRet,
         future_lapply(
            elt_strings, function(x, encoding=encoding, ...){
               fun(xml2::read_xml(x, encoding=encoding), ...)
            }
         )
      )
   }
   close(con)
   return(toRet)
}
