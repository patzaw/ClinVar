library(here)

#' Write a file with last resource update information
#'
#' @param who a character (default: Sys.info()["user"])
#' @param when a date (default: Sys.Date)
#' @param file the file in which the information should be written.
#' (default: here(".lastUpdate.txt"))
#' @param append a logical indicating if the information should be added at
#' the end of the file when it exists (default: TRUE)
#' @param overwrite a logical indicating if the file should be overwritten
#' when it exists and append is FALSE (default: TRUE)
#'
writeLastUpdate <- function(
   who=Sys.info()["user"],
   when=Sys.Date(),
   file=here(".lastestUpdates.txt"),
   append=TRUE, overwrite=TRUE
){
   stopifnot(
      length(who)==1,
      is.character(who),
      length(when)==1,
      inherits(when, "Date"),
      length(file)==1,
      is.character(file),
      length(append)==1,
      is.logical(append),
      length(overwrite)==1,
      is.logical(overwrite)
   )
   towrite <- data.frame(
      who=who,
      when=format(when, "%Y-%m-%d"),
      stringsAsFactors=FALSE
   )
   append <- file.exists(file) & append
   if(!file.exists(file) || append || overwrite){
      write.table(
         towrite,
         file=file,
         row.names=FALSE,
         col.names=!append,
         sep="\t",
         qmethod="double",
         append=append
      )
   }else{
      warning(
         "The file already exists. Nothing has been added to it",
         " and it has not been overwritten."
      )
   }
}
