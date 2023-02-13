###############################################################################@
## Config ----
###############################################################################@

library(here)
library(ReDaMoR)

source(here("scripts/supporting-functions-XML.R"))
sdir <- here("sources")
ddir <- here("data")

dm <- read_json_data_model(here("model/ClinVar.json"))

###############################################################################@
## Parsing XML ----
###############################################################################@

xmlf <- file.path(sdir, "ClinVarFullRelease.xml.gz")

# print(Sys.time())
# con <- file(xmlf)
# a <- extract_nsn_xml_elements(con=con, n=5*10^6, tag="ClinVarSet")
# close(con)
# print(Sys.time())
# ax2 <- xmlRoot(xmlParse(a[[100]], asText=TRUE))

# ax <- read_xml(a[[14360]])
# node <- xml_child(ax, "ReferenceClinVarAssertion") %>% xml_child("MeasureSet") %>% xml_child("Measure")
# system.time(ar <- parse_MeasureRelationship(node))
# library(XML)
# bx <- xmlRoot(xmlParse(a[[14360]]))
# bxnode <- bx[["ReferenceClinVarAssertion"]][["MeasureSet"]][["Measure"]]
# system.time(br <- getMeasureRelationship(bxnode, "74470"))
# system.time(br2 <- parse_MeasureRelationship(bxnode))


message("Parsing XML...")
message(Sys.time())
cvs_dir <- process_ClinVarSets(
   file=xmlf, location=here("scripts"),
   # n_max=1*10^6,
   # by=3*10^5,
   n_max=Inf,
   by=5*10^7,
   cores=15
)
message(Sys.time())

message("Combining records...")
message(Sys.time())
# cvs_dir <- here("scripts/CVS-batches-2023-02-12_195208/")
d <- combine_ClinVarSets(cvs_dir)
message(Sys.time())


###############################################################################@
## Organizing tables ----
###############################################################################@

message("Organizing tables...")
message(Sys.time())

###############################################################################@
### ReferenceClinVarAssertion ----

ClinVar_ReferenceClinVarAssertion <- d$ReferenceClinVarAssertion %>% 
   mutate(
      cvs=as.integer(cvs),
      id=as.integer(id)
   )
d$ReferenceClinVarAssertion <- NULL

#### clinSigOrder (CUSTOM) ----
ClinVar_clinSigOrder <- c(
   "protective",
   "Benign",
   "Benign; association",
   "Benign, other",
   "Benign; other",
   "Benign/Likely benign",
   "Benign/Likely benign, other",
   "Benign/Likely benign; other",
   "Likely benign",
   "Likely benign, other",
   "Likely benign; other",
   "drug response",
   "confers sensitivity",
   "conflicting data from submitters", 
   "association not found",
   "not provided",
   "Uncertain significance",
   "Uncertain significance, other",
   "Uncertain significance; other",
   "Uncertain significance; other; risk factor",
   "Uncertain significance, risk factor",
   "Uncertain significance; risk factor",
   "Uncertain significance; Pathogenic/Likely pathogenic",
   "Uncertain significance; Affects",
   "Conflicting interpretations of pathogenicity, other",
   "Conflicting interpretations of pathogenicity; other",
   "Conflicting interpretations of pathogenicity",
   "Conflicting interpretations of pathogenicity, association, risk factor",
   "Conflicting interpretations of pathogenicity; association; risk factor",
   "Conflicting interpretations of pathogenicity, risk factor",
   "Conflicting interpretations of pathogenicity; risk factor",
   "Conflicting interpretations of pathogenicity, association",
   "Conflicting interpretations of pathogenicity; association",
   "Uncertain risk allele",
   "other",
   "other, risk factor",
   "other; risk factor",
   "association",
   "Benign, risk factor",
   "Benign; risk factor",
   "Benign/Likely benign, risk factor",
   "Benign/Likely benign; risk factor",
   "Likely benign; association",
   "Likely benign, risk factor",
   "Likely benign; risk factor",
   "Likely risk allele; risk factor",
   "risk factor",
   "Established risk allele",
   "Uncertain significance; association",
   "Uncertain significance, Affects",
   "Uncertain significance/Uncertain risk allele",
   "Affects, other",
   "Affects",
   "Conflicting interpretations of pathogenicity; other; risk factor", 
   "Likely pathogenic; risk factor",
   "Likely pathogenic, risk factor",
   "Pathogenic/Likely pathogenic, association",
   "Pathogenic/Likely pathogenic; association",
   "Pathogenic/Likely pathogenic, risk factor",
   "Pathogenic/Likely pathogenic; risk factor",
   "Likely risk allele",
   "Likely pathogenic, low penetrance",
   "Likely pathogenic",
   "Likely pathogenic, other",
   "Likely pathogenic; other",
   "Likely pathogenic; other; risk factor",
   "Likely pathogenic/Likely risk allele",
   "Likely pathogenic, Affects",
   "Likely pathogenic; Affects",
   "Pathogenic/Pathogenic, low penetrance",
   "Pathogenic/Likely pathogenic/Likely risk allele",
   "Pathogenic/Likely risk allele",
   "Pathogenic/Likely pathogenic/Established risk allele",
   "Pathogenic/Likely pathogenic/Established risk allele; risk factor",
   "Pathogenic/Likely pathogenic/Pathogenic, low penetrance",
   "Pathogenic/Likely pathogenic",
   "Pathogenic/Likely pathogenic, other",
   "Pathogenic/Likely pathogenic; other",
   "Pathogenic, drug response",
   "Pathogenic; drug response",
   "Pathogenic, risk factor",
   "Pathogenic; risk factor",
   "Pathogenic/Established risk allele",
   "Pathogenic, association",
   "Pathogenic; association",
   "Pathogenic, low penetrance" ,
   "Pathogenic, other",
   "Pathogenic; other",
   "Pathogenic, Affects",
   "Pathogenic; Affects",
   "Pathogenic"
)
missing <- setdiff(
   ClinVar_ReferenceClinVarAssertion$clinicalSignificance,
   ClinVar_clinSigOrder
)
if(length(missing) > 0){
   print(missing)
   stop("Missing clinical significance values in clinSigOrder")
}
ClinVar_clinSigOrder <- tibble(
   label=ClinVar_clinSigOrder,
   order=1:length(ClinVar_clinSigOrder)
)

#### revStatOrder (CUSTOM) ----
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
missing <- setdiff(
   ClinVar_ReferenceClinVarAssertion$reviewStatus,
   ClinVar_revStatOrder
)
if(length(missing) > 0){
   print(missing)
   stop("Missing review status values in revStatOrder")
}
ClinVar_revStatOrder <- tibble(
   label=ClinVar_revStatOrder,
   order=1:length(ClinVar_revStatOrder)
)

###############################################################################@
### rcvaInhMode ----

ClinVar_rcvaInhMode <- d$rcvaAttributes %>%
   filter(Type=="ModeOfInheritance") %>%
   select(rcvaId, value) %>%
   mutate(
      value=toupper(value),
      rcvaId=as.integer(rcvaId)
   ) %>%
   distinct() %>%
   rename("inhMode"="value")

otherAttTypes <- unique(setdiff(d$rcvaAttributes$Type, "ModeOfInheritance"))
if(length(otherAttTypes)>0){
   stop("Other RCVA attributes: ", paste(otherAttTypes, collapse=", "))
}
d$rcvaAttributes <- NULL

###############################################################################@
### rcvaObservedInd ----

ClinVar_rcvaObservedIn <- d$observedIn %>%
   mutate(
      rcvaId=as.integer(rcvaId),
      numberTested=as.integer(numberTested)
   ) %>%
   select(
      rcvaId, origin, taxonomyId, species, affectedStatus, numberTested
   ) %>%
   distinct()
d$observedIn <- NULL

###############################################################################@
### ClinVarAssertions ----

ClinVar_ClinVarAssertions <- d$ClinVarAssertions %>%
   mutate(
      cvs=as.integer(cvs),
      id=as.integer(id)
   )
d$ClinVarAssertions <- NULL

###############################################################################@
### cvaSubmitters ----

ClinVar_cvaSubmitters <- d$submitters %>%
   select(cvaId, submitter, primary) %>%
   mutate(cvaId=as.integer(cvaId)) %>%
   distinct()
d$submitters <- NULL

###############################################################################@
### cvaObservedIn ----

ClinVar_cvaObservedIn <- d$cvaObservedIn %>%
   select(cvaId, origin, species, affectedStatus) %>%
   mutate(cvaId=as.integer(cvaId)) %>%
   distinct()
d$cvaObservedIn <- NULL

###############################################################################@
### rcvaTraits ----

ClinVar_rcvaTraits <- d$traits %>%
   select(rcvaId, id, type) %>%
   rename("t.id"="id", "traitType"="type") %>%
   mutate(
      rcvaId=as.integer(rcvaId),
      t.id=as.integer(t.id)
   ) %>%
   distinct()
d$traits <- NULL

###############################################################################@
### traits ----

ClinVar_traits <- d$traitNames %>%
   filter(trait.name.type=="Preferred") %>%
   select(trait.id, trait.name) %>% 
   distinct()
toAdd <- unique(setdiff(d$traitNames$trait.id, ClinVar_traits$trait.id))
if(length(toAdd) > 0){
   ClinVar_traits <- ClinVar_traits %>%
      bind_rows(
         d$traitNames %>%
            filter(trait.id %in% toAdd & !duplicated(trait.id)) %>%
            select(trait.id, trait.name)
      )
}
ClinVar_traits <- ClinVar_traits%>%
   rename("id"="trait.id", "name"="trait.name") %>%
   mutate(id=as.integer(id)) %>%
   filter(!duplicated(id))

###############################################################################@
### traitNames ----

ClinVar_traitNames <- d$traitNames %>%
   select(trait.id, trait.name, trait.name.type) %>%
   rename("t.id"="trait.id", "name"="trait.name", "type"="trait.name.type") %>%
   mutate(t.id=as.integer(t.id)) %>% 
   distinct()
d$traitNames <- NULL

###############################################################################@
### traitCref ----

ClinVar_traitCref <- d$traitXRef %>%
   select(trait.id, ID, DB, Type) %>%
   rename("t.id"="trait.id", "id"="ID", "db"="DB", "type"="Type") %>%
   mutate(t.id=as.integer(t.id)) %>%
   distinct()
d$traitXRef <- NULL
## Cleaning trait cross references
dbCleanTable <- read_tsv(
   here("scripts/DB-ID-Cleaning-Table.txt"),
   col_types="cccl"
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
### variants, rcvaVariant, varNames ----

##
ClinVar_variants <- d$measures %>%
   select(id=ID, type=Type) %>%
   mutate(id=as.integer(id)) %>%
   distinct()
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
ClinVar_rcvaVariant <- d$measures %>%
   select(id=ID, rcvaId) %>%
   rename("varId"="id") %>%
   mutate(
      varId=as.integer(varId),
      rcvaId=as.integer(rcvaId)
   ) %>%
   distinct()
d$measures <- NULL
##
ClinVar_varNames <- d$measureNames %>%
   select(measureId, name, type) %>%
   rename("varId"="measureId") %>%
   mutate(varId=as.integer(varId)) %>%
   distinct()
d$measureNames <- NULL
##
ClinVar_variants <- ClinVar_variants %>%
   left_join(
      ClinVar_varNames %>%
         arrange(desc(type)) %>%
         filter(!duplicated(varId)) %>%
         rename("id"="varId") %>%
         select(id, name),
      by="id"
   ) %>% 
   distinct()

###############################################################################@
### varEntrez, entrezNames ----

ClinVar_varEntrez <- d$measureRelationships %>%
   filter(DB=="Gene") %>%
   select(measureId, ID, type) %>%
   rename("varId"="measureId", "entrez"="ID") %>%
   mutate(
      varId=as.integer(varId),
      entrez=as.integer(entrez)
   ) %>%
   distinct()
##
ClinVar_entrezNames <- d$measureRelationships %>%
   as_tibble()%>%
   filter(DB=="Gene") %>%
   select(ID, name, symbol) %>%
   rename("entrez"="ID") %>%
   mutate(
      entrez=as.integer(entrez)
   ) %>%
   distinct(entrez, .keep_all=TRUE)
d$measureRelationships <- NULL

###############################################################################@
### varCytoLoc ----

ClinVar_varCytoLoc <- d$measureCytogeneticLocations %>%
   select(measureId, cytogenicLocation) %>%
   rename("varId"="measureId", "location"="cytogenicLocation") %>%
   mutate(varId=as.integer(varId)) %>%
   distinct()
d$measureCytogeneticLocations <- NULL

###############################################################################@
### varSeqLoc ----

ClinVar_varSeqLoc <- d$measureSequenceLocations %>%
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
      stop=as.integer(stop),
      variantLength=as.integer(variantLength)
   ) %>%
   distinct()
d$measureSequenceLocations <- NULL

###############################################################################@
### varXRef ----

ClinVar_varXRef <- d$measureXRef %>%
   select(measureId, ID, DB, Type) %>%
   rename("varId"="measureId", "id"="ID", "db"="DB", "type"="Type") %>%
   mutate(varId=as.integer(varId)) %>%
   distinct()
d$measureXRef <- NULL

###############################################################################@
### varAttributes ----

ClinVar_varAttributes <- d$measureAttributes %>%
   select(
      measureId, Type, integerValue, Change, value
   ) %>%
   rename("varId"="measureId") %>%
   mutate(
      varId=as.integer(varId),
      integerValue=as.integer(integerValue),
      value=ifelse(value=="", NA,value)
   ) %>% 
   distinct()
d$measureAttributes <- NULL
##
message(Sys.time())
message("... Done\n")

###############################################################################@
## Source information ----
###############################################################################@

ClinVar_sourceFiles <- read_tsv(
   file.path(sdir, "ARCHIVES/ARCHIVES.txt"),
   col_types="ccTl"
) %>%
   filter(inUse) %>% 
   select(url, current) %>% 
   mutate(current=as.Date(current))

###############################################################################@
## Check tables ----
###############################################################################@

confront_data(dm, lapply(names(dm), function(x) get(x)) %>% setNames(names(dm)))

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


