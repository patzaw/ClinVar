setwd("~/Shared/Data-Science/Data-Source-Model-Repository/clinVar/scripts/")

library(XML)
library(parallel)

source("clinVar-Functions.R")

##
mc.cores <- 55
sdir <- "../sources"
ddir <- "../data"

###############################################################################@
## Source information ----
###############################################################################@

sfi <- read.table(
   file.path(sdir, "ARCHIVES/ARCHIVES.txt"),
   sep="\t",
   header=T,
   stringsAsFactors=FALSE
)
ClinVar_sourceFiles <- sfi[which(sfi$inUse), c("url", "current")]

###############################################################################@
## Data from ClinVarFullRelease.xml.gz ----
###############################################################################@

## * Loading and parsing XML
message("Loading XML...")
message(Sys.time())
xmlFile <- file.path(sdir, "ClinVarFullRelease.xml.gz")
cvList <- readClinVar(xmlFile) #, n=400000) # total: >37210660
encoding <- attr(cvList, "encoding")
message(Sys.time())
message("... Done\n")

##
# x <- cvList[[36816]] # 18, 3
# node <- xmlRoot(xmlParse(x, encoding=encoding))
# # node <- node[["ReferenceClinVarAssertion"]]

message("Parsing XML...")
message(Sys.time())
scope <- c(
   "lrbind",
   "readClinVar",
   "parseCvs",
   "parseRcva",
   "getAttrSet",
   "parseObservedIn",
   "parseCvaObservedIn",
   "parseMeasureSet",
   "getSeqLoc",
   "getXRef",
   "getMeasureRelationship",
   "parseTraitSet",
   "parseCva",
   "encoding"
)
cl <- makeCluster(mc.cores)
clusterExport(cl, scope)
clusterEvalQ(cl, {library(XML)})
cvDbList <- parLapply(
   cl,
   cvList,
   function(x){
      parseCvs(xmlRoot(xmlParse(x, encoding=encoding)))
   }
)
stopCluster(cl)
rm(cl)
# cvDbList <- mclapply(
#   cvList,
#   function(x){
#     parseCvs(xmlRoot(xmlParse(x, encoding=encoding)))
#   },
#   mc.cores=mc.cores
# )
message(Sys.time())
message("... Done\n")

## * Generating tables ----
message("Generating tables...")
message(Sys.time())
toMerge <- unique(unlist(lapply(cvDbList, names)))
for(tm in toMerge){
   message(tm)
   assign(
      x=paste("ClinVar_", tm, sep=""),
      value=lrbind(cvDbList, tm)
   )
   message(Sys.time())
}
message(Sys.time())
message("... Done\n")

## * Organizing tables ----
message("Organizing tables...")
message(Sys.time())
##
ClinVar_ReferenceClinVarAssertion <- merge(
   ClinVar_ReferenceClinVarAssertion,
   ClinVar_ClinVarSet,
   by.x="cvs",
   by.y="id",
   all=T
)
rm(ClinVar_ClinVarSet)
##
ClinVar_rcvaInhMode <- ClinVar_rcvaAttributes[
   which(ClinVar_rcvaAttributes$Type=="ModeOfInheritance"),
   c("rcva.id", "value")
   ]
ClinVar_rcvaInhMode$value <- toupper(ClinVar_rcvaInhMode$value)
ClinVar_rcvaInhMode <- unique(ClinVar_rcvaInhMode)
colnames(ClinVar_rcvaInhMode) <- c("rcvaId", "inhMode")
otherAttTypes <- unique(setdiff(ClinVar_rcvaAttributes$Type, "ModeOfInheritance"))
if(length(otherAttTypes)>0){
   warning("Other RCVA attributes: ", paste(otherAttTypes, collapse=", "))
}
rm(ClinVar_rcvaAttributes)
##
ClinVar_rcvaObservedIn <- ClinVar_observedIn
ClinVar_rcvaObservedIn <- ClinVar_rcvaObservedIn[,c(
   "rcvaId", setdiff(colnames(ClinVar_rcvaObservedIn), "rcvaId")
)]
rm(ClinVar_observedIn)
##
# ClinVar_ClinVarAssertions unchanged
##
##
ClinVar_cvaSubmitters <- ClinVar_submitters[,c("cvaId", "submitter", "primary")]
rm(ClinVar_submitters)
##
ClinVar_cvaObservedIn <- ClinVar_cvaObservedIn[,c(
   "cvaId", setdiff(colnames(ClinVar_cvaObservedIn), "cvaId")
)]
##
ClinVar_rcvaTraits <- ClinVar_traits[,c("rcvaId", "id", "type")]
colnames(ClinVar_rcvaTraits) <- c("rcvaId", "t.id", "traitType")
rm(ClinVar_traits)
##
ClinVar_traits <- ClinVar_traitNames[
   which(ClinVar_traitNames$trait.name.type=="Preferred"),
   c("trait.id", "trait.name")
   ]
toAdd <- unique(setdiff(ClinVar_traitNames$trait.id, ClinVar_traits$trait.id))
if(length(toAdd) > 0){
   toAdd <- ClinVar_traitNames[
      which(
         ClinVar_traitNames$trait.id %in% toAdd &
            !duplicated(ClinVar_traitNames$trait.id)
      ),
      c("trait.id", "trait.name")
      ]
   ClinVar_traits <- rbind(ClinVar_traits, toAdd)
}
colnames(ClinVar_traits) <- c("id", "name")
##
ClinVar_traitNames <- ClinVar_traitNames[, c(
   "trait.id", "trait.name", "trait.name.type"
)]
colnames(ClinVar_traitNames) <- c("t.id", "name", "type")
##
ClinVar_traitCref <- ClinVar_traitXRef[,c("trait.id", "ID", "DB", "Type")]
colnames(ClinVar_traitCref) <- c("t.id", "id", "db", "type")
rm(ClinVar_traitXRef)
## Cleaning trait cross references
dbCleanTable <- read.table(
   "DB-ID-Cleaning-Table.txt",
   sep="\t", header=TRUE, stringsAsFactors=FALSE
)
for(i in 1:nrow(dbCleanTable)){
   oriDb <- dbCleanTable$ClinVar.DB[i]
   finalDb <- dbCleanTable$DB[i]
   prefix <- dbCleanTable$Prefix.to.remove[i]
   ClinVar_traitCref[which(ClinVar_traitCref$db==oriDb), "id"] <- sub(
      paste0("^", prefix), "",
      ClinVar_traitCref[which(ClinVar_traitCref$db==oriDb), "id"]
   )
   ClinVar_traitCref[which(ClinVar_traitCref$db==oriDb), "db"] <- finalDb
}
##
ClinVar_variants <- unique(ClinVar_measures[,c("id", "type")])
ClinVar_rcvaVariant <- ClinVar_measures[,c("id", "rcvaId")]
colnames(ClinVar_rcvaVariant) <- c("varId", "rcvaId")
rm(ClinVar_measures)
##
ClinVar_varNames <- ClinVar_measureNames[,c("measureId", "name", "type")]
colnames(ClinVar_varNames) <- c("varId", "name", "type")
rm(ClinVar_measureNames)
ClinVar_variants <- merge(
   ClinVar_variants,
   ClinVar_varNames,
   by.x="id", by.y="varId",
   all=T
)
ClinVar_variants <- ClinVar_variants[order(ClinVar_variants$type.y, decreasing=T),]
ClinVar_variants <- ClinVar_variants[
   which(!duplicated(ClinVar_variants$id)),
   c("id", "type.x", "name")
   ]
colnames(ClinVar_variants) <- c("id", "type", "name")
##
ClinVar_varEntrez <- unique(ClinVar_measureRelationships[
   which(
      # ClinVar_measureRelationships$type=="variant in gene" &
      ClinVar_measureRelationships$xref.db=="Gene"
   ),
   c("measureId", "xref.id", "type")
   ])
colnames(ClinVar_varEntrez) <- c("varId", "entrez", "type")
##
ClinVar_entrezNames <- unique(ClinVar_measureRelationships[
   which(ClinVar_measureRelationships$xref.db=="Gene"),
   c("xref.id", "name", "symbol")
   ])
colnames(ClinVar_entrezNames) <- c("entrez", "name", "symbol")
# otherRSTypes <- unique(setdiff(ClinVar_measureRelationships$type, "variant in gene"))
# if(length(otherRSTypes)>0){
#   warning("Other measureRelationShips types: ", paste(otherRSTypes, collapse=", "))
# }
rm(ClinVar_measureRelationships)
##
ClinVar_varCytoLoc <- ClinVar_measureCytogeneticLocations[,c("measureId", "cytogenicLocation")]
colnames(ClinVar_varCytoLoc) <- c("varId", "location")
rm(ClinVar_measureCytogeneticLocations)
##
ClinVar_varSeqLoc <- ClinVar_measureSequenceLocations[,c(
   "measureId", setdiff(colnames(ClinVar_measureSequenceLocations), "measureId")
)]
colnames(ClinVar_varSeqLoc) <- c(
   "varId", "Accession", "alternateAllele", "Assembly",
   "AssemblyAccessionVersion", "AssemblyStatus", "Chr",
   "display_start", "display_stop", "innerStart", "innerStop",
   "outerStart", "outerStop", "referenceAllele", "start", "stop",
   "Strand", "variantLength"
)
rm(ClinVar_measureSequenceLocations)
##
ClinVar_varXRef <- ClinVar_measureXRef[,c("measureId", "ID", "DB", "Type")]
colnames(ClinVar_varXRef) <- c("varId", "id", "db", "type")
rm(ClinVar_measureXRef)
##
ClinVar_varAttributes <- ClinVar_measureAttributes[
   ,
   c("measureId", setdiff(colnames(ClinVar_measureAttributes), "measureId"))
   ]
colnames(ClinVar_varAttributes) <- c("varId", "Type", "integerValue", "Change", "value")
rm(ClinVar_measureAttributes)
##
message(Sys.time())
message("... Done\n")

###############################################################################@
## Data from disease_names ----
###############################################################################@

# message("Parsing diseases_names...")
# message(Sys.time())
# ##
# rd <- readLines(file.path(sdir, "disease_names"))
# header <- unlist(strsplit(rd[1], split="\t"))
# rdlist <- strsplit(rd[-1], "\t")
# rdlist <- lapply(
#    rdlist, function(x){
#       x <- sub("^ ", "", sub(" $", "", x))
#       toRet <- x
#       if(length(x)>7){
#          toRet <- c(
#             paste(x[1:(length(x)-6)], collapse=" // "),
#             x[(length(x)-5):length(x)]
#          )
#       }
#       return(toRet)
#    }
# )
# ClinVar_diseaseNames <- as.data.frame(do.call(
#    rbind,
#    rdlist
# ), stringsAsFactors=F)
# colnames(ClinVar_diseaseNames) <- c(
#    "name", "source", "concept", "sourceID", "MIM", "LastModif", "Category"
# )
# rm(rd, rdlist, header)
# ##
# message(Sys.time())
# message("... Done\n")

###############################################################################@
## Custom information ----
###############################################################################@
ClinVar_clinSigOrder <- data.frame(
   label=c(
      "protective",
      "Benign",
      "Likely benign",
      "drug response",
      "confers sensitivity",
      "conflicting data from submitters", 
      "not provided",
      "Uncertain significance",
      "other",
      "association",
      "risk factor",
      "Likely pathogenic",
      "Pathogenic"
   ),
   order=1:13,
   stringsAsFactors=FALSE
)
ClinVar_revStatOrder <- data.frame(
   label=c(
      "not classified by submitter",
      "classified by single submitter", 
      "classified by multiple submitters",
      "reviewed by professional society",
      "reviewed by expert panel"
   ),
   order=1:5,
   stringsAsFactors=FALSE
)

###############################################################################@
## Writing tables ----
###############################################################################@

message("Writing tables...")
message(Sys.time())
toSave <- grep("^ClinVar[_]", ls(), value=T)
for(f in toSave){
   message(paste("   Writing", f))
   ## Ensure unicity
   toWrite <- get(f)
   write.table(
      get(f),
      file=file.path(ddir, paste(f, ".txt", sep="")),
      sep="\t",
      row.names=FALSE, col.names=TRUE,
      quote=TRUE,
      qmethod="double"
   )
}
message(Sys.time())
message("... Done\n")
