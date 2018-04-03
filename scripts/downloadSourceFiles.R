library(RCurl)
downloadSourceFiles <- function(urls, directory){
   if(any(names(urls)=="ARCHIVES")){
      stop('"ARCHIVES" is a reserved file name: choose another one.')
   }
   if(!file.exists(directory)){
      stop(sprintf("The %s directory does not exist.", directory))
   }
   archDir <- file.path(directory, "ARCHIVES")
   archFile <- file.path(archDir, "ARCHIVES.txt")
   if(!file.exists(archDir)){
      dir.create(archDir)
   }
   if(file.exists(archFile)){
      archTable <- read.table(
         archFile,
         sep="\t", header=TRUE,
         stringsAsFactors=FALSE,
         check.names=FALSE
      )
   }else{
      archTable <- data.frame(
         file=character(),
         url=character(),
         current=character(),
         inUse=logical(),
         stringsAsFactors=FALSE
      )
   }
   archTable$current <- as.POSIXct(archTable$current, tz="GMT")
   for(name in names(urls)){
      url <- urls[name]
      ufile <- basename(url)
      current <- archTable$current[which(
         archTable$file==name & archTable$url==url
      )]
      protocol <- ifelse(
         length(grep("^ftp://", url, ignore.case=TRUE)==1), "FTP",
         ifelse(
            length(grep("^http://", url, ignore.case=TRUE)==1), "HTTP",
            ifelse(
               length(grep("^https://", url, ignore.case=TRUE)==1), "HTTP",
               NA
            )
         )
      )
      if(is.na(protocol)){
         stop("Unknown URL protocol")
      }
      if(protocol=="HTTP"){
         hg <- basicHeaderGatherer()
         httpHEAD(url, headerfunction=hg$update)
         if(hg$value()[["status"]]!="200"){
            stop(paste(url, "not available"))
         }
         rcurrent <- as.POSIXct(
            hg$value()["Last-Modified"],
            format="%a, %d %b %Y %T GMT",
            tz="GMT"
         )
         if(is.na(current) || (rcurrent - current) > 0){
            toSave <- TRUE
            message("Downloading ", ufile, "...")
            message("   ", Sys.time())
            content <- getBinaryURL(url, headerfunction=hg$update)
            rcurrent <- as.POSIXct(
               hg$value()["Last-Modified"],
               format="%a, %d %b %Y %T GMT",
               tz="GMT"
            )
            message("   ", Sys.time())
            message("...Done")
         }else{
            toSave <- FALSE
         }
      }
      if(protocol=="FTP"){
         udir <- paste0(dirname(url), "/")
         dirls <- unlist(strsplit(getURL(udir), split="\n"))
         finfo <- grep(paste0("[^>] ", ufile, "$"), dirls, value=T)
         if(length(finfo)==0){
            stop(paste(url, "not available"))
         }
         fdate <- sub(paste0(ufile, "$"), "", finfo)
         torm <- regexpr("^.* [[:alpha:]]", fdate)
         fdate <- substr(fdate, attr(torm, "match.length"), nchar(fdate))
         fdate <- sub("^ *", "", sub(" *$", "", fdate))
         rcurrent <- as.POSIXct(
            fdate,
            format="%b %d %H:%M",
            tz="GMT"
         )
         if(is.na(rcurrent)){
            rcurrent <- as.POSIXct(
               fdate,
               format="%b %d %Y",
               tz="GMT"
            )
         }
         if(is.na(rcurrent)){
            stop(sprintf("Could not find a date for file %s", name))
         }
         if(length(current)==0 || (rcurrent - current) > 0){
            toSave <- TRUE
            message("Downloading ", ufile, "...")
            message("   ", Sys.time())
            content <- getBinaryURL(url)
            message("   ", Sys.time())
            message("...Done")
         }else{
            toSave <- FALSE
         }
      }
      if(toSave){
         archTable[which(archTable$file==name), "inUse"] <- rep(
            FALSE,
            sum(archTable$file==name)
         )
         if(length(current)==0){
            archTable <- rbind(
               archTable,
               data.frame(
                  file=name,
                  url=url,
                  current=rcurrent,
                  inUse=TRUE,
                  stringsAsFactors=FALSE
               )
            )
         }else{
            archTable[
               which(archTable$file==name & archTable$url==url),
               "current"
            ] <- rcurrent
            archTable[
               which(archTable$file==name & archTable$url==url),
               "inUse"
            ] <- TRUE
         }
         destFiles <- c(
            file.path(directory, name),
            file.path(archDir, Sys.Date(), ufile)
         )
         for(f in destFiles){
            message("Writting ", f, "...")
            message("   ", Sys.time())
            dir.create(dirname(f), recursive=TRUE, showWarnings=FALSE)
            writeBin(content, con=f)
            message("   ", Sys.time())
            message("...Done")
         }
         success <- try(write.table(
            archTable,
            file=archFile,
            sep="\t",
            row.names=FALSE,
            col.names=TRUE,
            quote=FALSE
         ), silent=TRUE)
         if(inherits(success, "try-error")){
            for(f in destFiles){
               file.remove(f)
            }
            stop(as.character(success))
         }
      }
   }
}
