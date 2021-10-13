library(here)
library(parallel)
library(XML)
library(tidyverse)

source(here("scripts/writeLastUpdate.R"))

source(here("scripts/clinVar-Functions.R"))

rm_tr_spaces <- function(x){
   stringr::str_remove(x, stringr::regex("^ *")) %>%
      stringr::str_remove(stringr::regex(" *$"))
}

##
mc.cores <- 12
sdir <- here("sources")
ddir <- here("data")

###############################################################################@
## Source information ----
###############################################################################@

sfi <- read.table(
   file.path(sdir, "ARCHIVES/ARCHIVES.txt"),
   sep="\t",
   header=T,
   stringsAsFactors=FALSE
)
ClinVar_sourceFiles <- sfi[which(sfi$inUse), c("url", "current")] %>%
   as_tibble() %>%
   mutate(current=as.character(as.Date(current)))

###############################################################################@
## Data from ClinVarFullRelease.xml.gz ----
###############################################################################@

###############################################################################@
## _+ Loading and parsing XML ----
message("Loading XML...")
message(Sys.time())
xmlFile <- file.path(sdir, "ClinVarFullRelease.xml.gz")
cvList <- readClinVar(xmlFile) #, n=100000) # total: >37210660
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

###############################################################################@
## _+ Generating tables ----
message("Generating tables...")
message(Sys.time())
toMerge <- unique(unlist(lapply(cvDbList, names)))
for(tm in toMerge){
   message(tm)
   assign(
      x=paste0("ClinVar_", tm),
      value=lrbind(cvDbList, tm)
   )
   message(Sys.time())
}
save(
   list=paste0("ClinVar_", toMerge),
   file=paste0("tmp-tables-", format(Sys.time(), "%Y-%m-%d-%H%M%S"), ".rda")
)
message(Sys.time())
message("... Done\n")

###############################################################################@
## Organizing tables ----
###############################################################################@

message("Organizing tables...")
message(Sys.time())

ClinVar_rcvaAttributes <- ClinVar_rcvaAttributes %>%
   as_tibble() %>%
   unique()

###############################################################################@
## _+ ReferenceClinVarAssertion ----

ClinVar_ReferenceClinVarAssertion <- merge(
   ClinVar_ReferenceClinVarAssertion,
   ClinVar_ClinVarSet,
   by.x="cvs",
   by.y="id",
   all=T
)
rm(ClinVar_ClinVarSet)
ClinVar_ReferenceClinVarAssertion <- ClinVar_ReferenceClinVarAssertion %>%
   as_tibble() %>%
   mutate(
      cvs=as.integer(cvs),
      id=as.integer(id)
   ) %>%
   unique()

###############################################################################@
## _+ rcvaInhMode ----

ClinVar_rcvaInhMode <- ClinVar_rcvaAttributes %>%
   filter(Type=="ModeOfInheritance") %>%
   select(rcva.id, value) %>%
   mutate(
      value=toupper(value),
      rcva.id=as.integer(rcva.id)
   ) %>%
   unique() %>%
   rename("rcvaId"="rcva.id", "inhMode"="value")


otherAttTypes <- unique(setdiff(ClinVar_rcvaAttributes$Type, "ModeOfInheritance"))
if(length(otherAttTypes)>0){
   warning("Other RCVA attributes: ", paste(otherAttTypes, collapse=", "))
}
rm(ClinVar_rcvaAttributes)

###############################################################################@
## _+ rcvaObservedInd ----

ClinVar_rcvaObservedIn <- ClinVar_observedIn %>%
   as_tibble() %>%
   mutate(
      rcvaId=as.integer(rcvaId),
      numberTested=as.integer(numberTested)
   ) %>%
   select(
      rcvaId, origin, taxonomyId, species, affectedStatus, numberTested
   ) %>%
   unique()
rm(ClinVar_observedIn)

###############################################################################@
## _+ ClinVarAssertions ----

ClinVar_ClinVarAssertions <- ClinVar_ClinVarAssertions %>%
   as_tibble() %>%
   mutate(
      cvs=as.integer(cvs),
      id=as.integer(id)
   ) %>%
   unique()

###############################################################################@
## _+ cvaSubmitters ----

ClinVar_cvaSubmitters <- ClinVar_submitters %>%
   as_tibble() %>%
   select(cvaId, submitter, primary) %>%
   mutate(cvaId=as.integer(cvaId)) %>%
   unique()
rm(ClinVar_submitters)

###############################################################################@
## _+ cvaObservedIn ----

ClinVar_cvaObservedIn <- ClinVar_cvaObservedIn %>%
   as_tibble() %>%
   select(cvaId, origin, species, affectedStatus) %>%
   mutate(cvaId=as.integer(cvaId)) %>%
   unique()

###############################################################################@
## _+ rcvaTraits ----

ClinVar_rcvaTraits <- ClinVar_traits %>%
   as_tibble() %>%
   select(rcvaId, id, type) %>%
   rename("t.id"="id", "traitType"="type") %>%
   mutate(
      rcvaId=as.integer(rcvaId),
      t.id=as.integer(t.id)
   ) %>%
   unique()
rm(ClinVar_traits)

###############################################################################@
## _+ traits ----

ClinVar_traits <- ClinVar_traitNames %>%
   as_tibble() %>%
   filter(trait.name.type=="Preferred") %>%
   select(trait.id, trait.name)
toAdd <- unique(setdiff(ClinVar_traitNames$trait.id, ClinVar_traits$trait.id))
if(length(toAdd) > 0){
   ClinVar_traits <- ClinVar_traits %>%
      bind_rows(
         ClinVar_traitNames %>%
         as_tibble() %>%
         filter(trait.id %in% toAdd & !duplicated(trait.id)) %>%
         select(trait.id, trait.name)
      )
}
ClinVar_traits <- ClinVar_traits%>%
   rename("id"="trait.id", "name"="trait.name") %>%
   mutate(id=as.integer(id)) %>%
   filter(!duplicated(id))

###############################################################################@
## _+ traitNames ----

ClinVar_traitNames <- ClinVar_traitNames %>%
   as_tibble() %>%
   select(trait.id, trait.name, trait.name.type) %>%
   rename("t.id"="trait.id", "name"="trait.name", "type"="trait.name.type") %>%
   mutate(t.id=as.integer(t.id))

###############################################################################@
## _+ traitCref ----

ClinVar_traitCref <- ClinVar_traitXRef %>%
   as_tibble() %>%
   select(trait.id, ID, DB, Type) %>%
   rename("t.id"="trait.id", "id"="ID", "db"="DB", "type"="Type") %>%
   mutate(t.id=as.integer(t.id)) %>%
   unique()
rm(ClinVar_traitXRef)
## Cleaning trait cross references
dbCleanTable <- read_tsv(
   here("scripts/DB-ID-Cleaning-Table.txt")
)
ClinVar_traitCref$id <- rm_tr_spaces(ClinVar_traitCref$id)
for(i in 1:nrow(dbCleanTable)){
   oriDb <- dbCleanTable$ClinVar.DB[i]
   finalDb <- dbCleanTable$DB[i]
   prefix <- dbCleanTable$Prefix.to.remove[i]
   blanks <- dbCleanTable$Blanks[i]
   if(!is.na(prefix)){
      ClinVar_traitCref[which(ClinVar_traitCref$db==oriDb), "id"] <- str_remove(
         ClinVar_traitCref$id[which(ClinVar_traitCref$db==oriDb)],
         sprintf("^(%s)*", prefix)
      )
   }
   if(!blanks){
      ClinVar_traitCref[which(ClinVar_traitCref$db==oriDb), "id"] <- str_remove(
         ClinVar_traitCref$id[which(ClinVar_traitCref$db==oriDb)],
         "[[:space:]]+.*$"
      )
   }
   ClinVar_traitCref[which(ClinVar_traitCref$db==oriDb), "db"] <- finalDb
}

###############################################################################@
## _+ variants, rcvaVariant, varNames ----

##
ClinVar_variants <- ClinVar_measures %>%
   as_tibble() %>%
   select(id, type) %>%
   mutate(id=as.integer(id)) %>%
   unique()
dupMistake <- ClinVar_variants %>%
   filter(duplicated(id))
if(nrow(dupMistake) > 0){
   message(
      "Warning message: The following duplicated variants have been removed:"
   )
   print(dupMistake)
   ClinVar_variants <- ClinVar_variants %>%
      filter(!duplicated(id))
}
##
ClinVar_rcvaVariant <- ClinVar_measures %>%
   as_tibble() %>%
   select(id, rcvaId) %>%
   rename("varId"="id") %>%
   mutate(
      varId=as.integer(varId),
      rcvaId=as.integer(rcvaId)
   ) %>%
   unique()
rm(ClinVar_measures)
##
ClinVar_varNames <- ClinVar_measureNames %>%
   as_tibble() %>%
   select(measureId, name, type) %>%
   rename("varId"="measureId") %>%
   mutate(varId=as.integer(varId)) %>%
   unique()
rm(ClinVar_measureNames)
##
ClinVar_variants <- ClinVar_variants %>%
   left_join(
      ClinVar_varNames %>%
         arrange(desc(type)) %>%
         filter(!duplicated(varId)) %>%
         rename("id"="varId") %>%
         select(id, name)
   )

###############################################################################@
## _+ varEntrez, entrezNames ----

ClinVar_varEntrez <- ClinVar_measureRelationships %>%
   as_tibble()%>%
   filter(xref.db=="Gene") %>%
   select(measureId, xref.id, type) %>%
   rename("varId"="measureId", "entrez"="xref.id") %>%
   mutate(
      varId=as.integer(varId),
      entrez=as.integer(entrez)
   ) %>%
   unique()
##
ClinVar_entrezNames <- ClinVar_measureRelationships %>%
   as_tibble()%>%
   filter(xref.db=="Gene") %>%
   select(xref.id, name, symbol) %>%
   rename("entrez"="xref.id") %>%
   mutate(
      entrez=as.integer(entrez)
   ) %>%
   unique() %>% 
   distinct(entrez, .keep_all=TRUE)
rm(ClinVar_measureRelationships)

###############################################################################@
## _+ varCytoLoc ----

ClinVar_varCytoLoc <- ClinVar_measureCytogeneticLocations %>%
   as_tibble() %>%
   select(measureId, cytogenicLocation) %>%
   rename("varId"="measureId", "location"="cytogenicLocation") %>%
   mutate(varId=as.integer(varId)) %>%
   unique()
colnames(ClinVar_varCytoLoc) <- c("varId", "location")
rm(ClinVar_measureCytogeneticLocations)

###############################################################################@
## _+ varSeqLoc ----

ClinVar_varSeqLoc <- ClinVar_measureSequenceLocations %>%
   as_tibble(ClinVar_measureSequenceLocations) %>%
   select(
      measureId, Accession, alternateAllele, Assembly,
      AssemblyAccessionVersion, AssemblyStatus, Chr,
      display_start, display_stop, innerStart, innerStop,
      outerStart, outerStop, referenceAllele, start, stop,
      Strand, variantLength
   ) %>%
   rename("varId"="measureId") %>%
   mutate(
      varId=as.integer(varId),
      display_start=as.integer(display_start),
      display_stop=as.integer(display_stop),
      innerStart=as.integer(innerStart),
      innerStop=as.integer(innerStop),
      outerStart=as.integer(outerStart),
      outerStop=as.integer(outerStop),
      start=as.integer(start),
      stop=as.integer(stop)
   ) %>%
   unique()
rm(ClinVar_measureSequenceLocations)

###############################################################################@
## _+ varXRef ----

ClinVar_varXRef <- ClinVar_measureXRef %>%
   as_tibble() %>%
   select(measureId, ID, DB, Type) %>%
   rename("varId"="measureId", "id"="ID", "db"="DB", "type"="Type") %>%
   mutate(varId=as.integer(varId)) %>%
   unique()
rm(ClinVar_measureXRef)

###############################################################################@
## _+ varAttributes ----

ClinVar_varAttributes <- ClinVar_measureAttributes %>%
   as_tibble() %>%
   select(
      measureId, Type, integerValue, Change, value
   ) %>%
   rename("varId"="measureId") %>%
   mutate(
      varId=as.integer(varId),
      integerValue=as.integer(integerValue),
      value=ifelse(value=="", NA,value)
   )
rm(ClinVar_measureAttributes)
##
message(Sys.time())
message("... Done\n")

###############################################################################@
## Custom information ----
###############################################################################@
ClinVar_clinSigOrder <- c(
	"protective",
	"Benign",
	"Benign, other",
	"Benign/Likely benign",
	"Benign/Likely benign, other",
	"Likely benign",
	"Likely benign, other",
	"drug response",
	"confers sensitivity",
	"conflicting data from submitters", 
	"association not found",
	"not provided",
	"Uncertain significance",
	"Uncertain significance, other",
	"Uncertain significance, risk factor",
	"Conflicting interpretations of pathogenicity, other",
	"Conflicting interpretations of pathogenicity",
	"Conflicting interpretations of pathogenicity, association, risk factor",
	"Conflicting interpretations of pathogenicity, risk factor",
	"Conflicting interpretations of pathogenicity, association",
	"other",
	"other, risk factor",
	"association",
	"Benign, risk factor",
	"Benign/Likely benign, risk factor",
	"Likely benign, risk factor",
	"risk factor",
	"Uncertain significance, Affects",
	"Affects, other",
	"Affects",
	"Likely pathogenic, risk factor",
	"Pathogenic/Likely pathogenic, association",
	"Pathogenic/Likely pathogenic, risk factor",
	"Likely pathogenic",
	"Likely pathogenic, other",
	"Likely pathogenic, Affects",
	"Pathogenic/Likely pathogenic",
	"Pathogenic/Likely pathogenic, other",
	"Pathogenic, drug response",
	"Pathogenic, risk factor",
	"Pathogenic, association",
	"Pathogenic, other",
	"Pathogenic, Affects",
	"Pathogenic"
)
ClinVar_clinSigOrder <- tibble(
	label=ClinVar_clinSigOrder,
	order=1:length(ClinVar_clinSigOrder)
)
ClinVar_revStatOrder <- c(
	"no assertion provided",
	"not classified by submitter",
	"no assertion criteria provided",
	"classified by single submitter", 
	"criteria provided, single submitter",
	"criteria provided, conflicting interpretations",
	"classified by multiple submitters",
	"criteria provided, multiple submitters, no conflicts",
	"reviewed by professional society",
	"practice guideline",
	"reviewed by expert panel"
)
ClinVar_revStatOrder <- tibble(
	label=ClinVar_revStatOrder,
	order=1:length(ClinVar_revStatOrder)
)

###############################################################################@
## Writing tables ----
###############################################################################@

message("Writing tables...")
message(Sys.time())
file.rename(ddir, paste0(ddir, "_BCK_", Sys.Date()))
dir.create(ddir)
toSave <- grep("^ClinVar[_]", ls(), value=T)
for(f in toSave){
   message(paste("   Writing", f))
   tv <- get(f)
   for(cn in colnames(tv)){
      if(class(pull(tv, !!cn))=="character"){
         tv[, cn] <- rm_tr_spaces(pull(tv, !!cn))
      }
   }
   tv <- distinct(tv)
   write_tsv(
      tv,
      file=file.path(ddir, paste(f, ".txt", sep="")),
      quote="all", na="<NA>"
   )
}
message(Sys.time())
message("... Done\n")

writeLastUpdate()
