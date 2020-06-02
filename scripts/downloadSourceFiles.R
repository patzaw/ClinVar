library(RCurl)
downloadSourceFiles <- function(
   urls, # a named vector: names are the names of the final files
   directory,
   httpForce=FALSE,
   httpFoundValues=c("200"),
   httpRedirectValues=c("302", "301"),
   maxRedirections=5
){
   stopifnot(is.logical(httpForce), length(httpForce)==1, !is.na(httpForce))
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
         fVal <- c(httpFoundValues, httpRedirectValues)
         hg <- basicHeaderGatherer()
         urlList <- url
         httpHEAD(url, headerfunction=hg$update)
         r <- 0
         while(
            hg$value()["status"] %in% httpRedirectValues &&
            r <= maxRedirections
         ){
            r <- r + 1
            hgv <- hg$value()
            if("location" %in% names(hgv)){
               url <- hgv["location"]
            }else{
               if("Location" %in% names(hgv)){
                  url <- hgv["Location"]
               }else{
                  stop("No location information for redirection")
               }
            }
            urlList <- c(urlList, url)
            httpHEAD(url, headerfunction=hg$update)
         }
         if(!httpForce){
            if(!hg$value()["status"] %in% httpFoundValues){
               stop(paste(paste(urlList, collapse=" -> "), " : not available"))
            }
            rcurrent <- as.POSIXct(
               hg$value()["Last-Modified"],
               format="%a, %d %b %Y %T GMT",
               tz="GMT"
            )
            if(is.na(rcurrent)){
               rcurrent <- Sys.time()
               rcurrent <- as.POSIXct(
                  format(rcurrent, "%a, %d %b %Y %T GMT", tz="GMT"),
                  format="%a, %d %b %Y %T GMT",
                  tz="GMT"
               )
            }
         }
         if(httpForce || length(current)==0 || is.na(current) || (rcurrent - current) > 0){
            toSave <- TRUE
            message("Downloading ", ufile, "...")
            message("   ", Sys.time())
            content <- getBinaryURL(url, headerfunction=hg$update)
            message("   ", Sys.time())
            message("...Done")
            if(httpForce){
               rcurrent <- as.POSIXct(
                  hg$value()["Last-Modified"],
                  format="%a, %d %b %Y %T GMT",
                  tz="GMT"
               )
               if(is.na(rcurrent)){
                  rcurrent <- Sys.time()
                  rcurrent <- as.POSIXct(
                     format(rcurrent, "%a, %d %b %Y %T GMT", tz="GMT"),
                     format="%a, %d %b %Y %T GMT",
                     tz="GMT"
                  )
               }
            }
         }else{
            toSave <- FALSE
         }
         url <- urlList[1]
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
            rcurrent <- Sys.time()
            rcurrent <- as.POSIXct(
               format(rcurrent, "%a, %d %b %Y %T GMT", tz="GMT"),
               format="%a, %d %b %Y %T GMT",
               tz="GMT"
            )
         }
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
         if(length(current)==0 || is.na(current) || (rcurrent - current) > 0){
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
         # if(length(current)==0){
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
         # }else{
         #    archTable[
         #       which(archTable$file==name & archTable$url==url),
         #       "current"
         #    ] <- rcurrent
         #    archTable[
         #       which(archTable$file==name & archTable$url==url),
         #       "inUse"
         #    ] <- TRUE
         # }
         destFiles <- c(
            file.path(directory, name),
            file.path(archDir, Sys.Date(), ufile)
         )
         for(f in destFiles){
            message("Writing ", f, "...")
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
