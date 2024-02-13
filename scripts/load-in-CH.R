source("~/opt/KMT.R")
library(here)

fdb <- read_fileMDB(here())
rn <- db_info(fdb)$name
chdb <- get_MDB(.tkcon, rn, check=FALSE)

stopifnot(db_info(fdb)$name==db_info(chdb)$name)
tbkmspec <- read_fileMDB(here(sprintf("TBKM/TBKM specification for %s", rn)))
collibra <- read_fileMDB(here(sprintf("TBKM/Collibra metadata for %s", rn)))
tst <-  max(c(
   file.info(data_files(fdb)$dataFiles)$mtime,
   file.info(data_files(tbkmspec)$dataFiles)$mtime,
   file.info(data_files(collibra)$dataFiles)$mtime
))
fdb <- c(fdb, collibra, tbkmspec)
db_info(fdb)$timestamp <- tst

# message("In memo")
# message(Sys.time())
# mdb <- as_memoMDB(fdb, check=FALSE)
# message("In CH")
# message(Sys.time())
# chdb <- as_chMDB(mdb, .tkcon, timestamp=db_info(mdb)$timestamp, by=10^7)
# message(Sys.time())

newdb <- list(
   ori=select(
      chdb, -all_of(intersect(names(chdb), c(names(tbkmspec), names(collibra))))
   ),
   part=select(fdb, c(names(tbkmspec), names(collibra)))
) %>% set_names(c(rn, "part"))
db_info(newdb$part)$name <- "part"
newdb <- metaMDB(
   MDBs=newdb, relationalTables=NULL,
   dataModel=data_model(fdb),
   dbInfo=db_info(fdb),
   check=FALSE
)
chdb <- as_chMDB(newdb, .tkcon, timestamp=db_info(newdb)$timestamp, by=10^7)
