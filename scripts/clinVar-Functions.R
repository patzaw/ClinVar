############################
# lrbind <- function(l, item, mc.cores=1){
#   unique(do.call(rbind, mclapply(
#     l,
#     function(x) x[[item]],
#     mc.cores=mc.cores
#   )))
# }

############################
lrbind <- function(l, item){
    litem <- lapply(
        l,
        function(x) x[[item]]
    )
    litem <- litem[which(unlist(lapply(litem, function(x) !is.null(x))))]
    if(length(litem)==0){
        return(c())
    }
    cn <- colnames(litem[[1]])
    f <- cn[1]
    toRet <- unlist(lapply(litem, function(x) x[,f]))
    for(f in cn[-1]){
        toRet <- data.frame(toRet, unlist(lapply(litem, function(x) x[,f])), stringsAsFactors=F)
    }
    colnames(toRet) <- cn
    return(unique(toRet))
}

############################
readClinVar <- function(file, n=-1L){
    cvRaw <- readLines(file, n=n)
    encoding <- sub(
        "\".*$",
        "",
        sub("^[<][?]xml .*encoding=\"", "", cvRaw[1])
    )
    starts <- grep("<ClinVarSet .*?>", cvRaw)
    ends <- grep("<[/]ClinVarSet>", cvRaw)
    starts <- starts[1:length(ends)]
    cvList <- apply(
        cbind(starts, ends), 1,
        function(x) sub(
            "[^>]*$",
            "",
            sub(
                "^[^<]*",
                "",
                paste(cvRaw[x[1]:x[2]], collapse="\n")
            )
        )
    )
    attr(cvList, "encoding") <- encoding
    return(cvList)
}

############################
parseCvs <- function(node){
    
    ## Preparing the list to return
    toRet <- list()
    
    ## Just take current records
    if(xmlValue(node[["RecordStatus"]]) != "current"){
        invisible(NULL)
    }
    
    ## ClinVarSet
    cvsId <- xmlAttrs(node)["ID"]
    cvsTitle <- xmlValue(node[["Title"]])
    toRet$ClinVarSet <- data.frame(
        id=cvsId,
        title=cvsTitle,
        stringsAsFactors=F
    )
    rownames(toRet$ClinVarSet) <- c()
    
    ## ReferenceClinVarAssertion
    if("ReferenceClinVarAssertion" %in% names(node)){
        rcvaOutput <- parseRcva(
            node[["ReferenceClinVarAssertion"]],
            cvsId
        )
    }else{
        rcvaOutput <- NULL
    }
    toRet$ReferenceClinVarAssertion <- rcvaOutput$ReferenceClinVarAssertion
    toRet$rcvaAttributes <- rcvaOutput$rcvaAttributes
    toRet$observedIn <- rcvaOutput$observedIn
    
    # toRet$measureSet <- rcvaOutput$measureSet
    for(ln in names(rcvaOutput$measureSet)){
        toRet[[ln]] <- rcvaOutput$measureSet[[ln]]
    }
    
    # toRet$traitSet <- rcvaOutput$traitSet
    for(ln in names(rcvaOutput$traitSet)){
        toRet[[ln]] <- rcvaOutput$traitSet[[ln]]
    }
    
    ## ClinVarAssertion
    cvaOutput <- parseCva(node, cvsId)
    toRet$ClinVarAssertions <- cvaOutput$cva
    toRet$submitters <- cvaOutput$submitters
    toRet$cvaObservedIn <- cvaOutput$observedIn
    
    ##
    return(toRet)
    
}

############################
parseRcva <- function(node, cvsId){
    
    ## Preparing the list to return
    toRet <- list()
    
    ## ReferenceClinVarAssertion
    rcvId <- xmlAttrs(node)["ID"]
    rcvAcc <- xmlAttrs(node[["ClinVarAccession"]])["Acc"]
    assertion <- xmlAttrs(node[["Assertion"]])["Type"]
    if("ClinicalSignificance" %in% names(node)){
        rcvCSrev <- xmlValue(node[["ClinicalSignificance"]][["ReviewStatus"]])
        rcvCSdes <- xmlValue(node[["ClinicalSignificance"]][["Description"]])
        rcvCSexp <- xmlValue(node[["ClinicalSignificance"]][["Explanation"]])
    }else{
        rcvCSrev <- rcvCSdes <- rcvCSexp <- NA
    }
    toRet$ReferenceClinVarAssertion <- data.frame(
        cvs=cvsId,
        id=rcvId,
        accession=rcvAcc,
        assertion=assertion,
        reviewStatus=rcvCSrev,
        clinicalSignificance=rcvCSdes,
        explanation=rcvCSexp,
        stringsAsFactors=F
    )
    rownames(toRet$ReferenceClinVarAssertion) <- c()
    
    ## ReferenceClinVarAssertion AttributeSets
    rcvaAttributes <- getAttrSet(node)
    if(!is.null(rcvaAttributes)){
        rownames(rcvaAttributes) <- c()
        rcvaAttributes$rcva.id <- rcvId
    }
    toRet$rcvaAttributes <- rcvaAttributes
    
    ## Observed In
    toRet$observedIn <- parseObservedIn(node, rcvId)
    
    ## MeasureSet
    toRet$measureSet <- parseMeasureSet(node, rcvId)
    
    ## TraitSet
    toRet$traitSet <- parseTraitSet(node, rcvId)
    
    ##
    return(toRet)
    
}

############################
getAttrSet <- function(node, attToGet=c("Type", "integerValue")){
    toTake <- which(names(node)=="AttributeSet")
    if(length(toTake)==0){
        return(NULL)
    }else{
        toRet <- unique(do.call(rbind, lapply(
            node[toTake],
            function(x){
                toRet <- as.data.frame(
                    t(xmlAttrs(x[["Attribute"]])[attToGet]),
                    stringsAsFactors=F
                )
                colnames(toRet) <- attToGet
                toRet$value <- xmlValue(x[["Attribute"]])
                return(toRet)
            }
        )))
        rownames(toRet) <- c()
        return(toRet)
    }
}

############################
parseObservedIn <- function(node, rcvaId){
    toTake <- which(names(node)=="ObservedIn")
    if(length(toTake)==0){
        return(NULL)
    }else{
        toRet <- unique(do.call(rbind, lapply(
            node[toTake],
            function(x){
                return(data.frame(
                    origin=xmlValue(x[["Sample"]][["Origin"]]),
                    taxonomyId=xmlAttrs(x[["Sample"]][["Species"]])["TaxonomyId"],
                    species=xmlValue(x[["Sample"]][["Species"]]),
                    affectedStatus=xmlValue(x[["Sample"]][["AffectedStatus"]]),
                    numberTested=xmlValue(x[["Sample"]][["NumberTested"]]),
                    # sampleDescription=xmlValue(x[["Sample"]][["SampleDescription"]][["Description"]]),
                    ## Information for further exploration
                    nbObsInTraitSet=sum(names(x)=="TraitSet"),
                    nbSampleTraitSet=sum(names(x[["Sample"]])=="TraitSet"),
                    stringsAsFactors=F
                ))
            }
        )))
        toRet$rcvaId <- rcvaId
        rownames(toRet) <- c()
        return(toRet)
    }
}

############################
parseCvaObservedIn <- function(node, cvaId){
    toTake <- which(names(node)=="ObservedIn")
    if(length(toTake)==0){
        return(NULL)
    }else{
        toRet <- unique(do.call(rbind, lapply(
            node[toTake],
            function(x){
                return(data.frame(
                    origin=xmlValue(x[["Sample"]][["Origin"]]),
                    species=xmlValue(x[["Sample"]][["Species"]]),
                    affectedStatus=xmlValue(x[["Sample"]][["AffectedStatus"]]),
                    stringsAsFactors=F
                ))
            }
        )))
        toRet$cvaId <- cvaId
        rownames(toRet) <- c()
        return(toRet)
    }
}

############################
parseMeasureSet <- function(node, rcvaId){
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
            ## Measure identification
            measure <- as.data.frame(
                t(xmlAttrs(x)[c("Type", "ID")]),
                stringsAsFactors=F
            )
            colnames(measure) <- c("type", "id")
            toRet$measure <- measure
            ## Names of measure
            toTake <- which(names(x)=="Name")
            if(length(toTake)==0){
                measureNames <- NULL
            }else{
                measureNames <- unique(do.call(rbind, lapply(
                    x[toTake],
                    function(y){
                        return(data.frame(
                            name=xmlValue(y[["ElementValue"]]),
                            type=xmlAttrs(y[["ElementValue"]])["Type"],
                            stringsAsFactors=F
                        ))
                    }
                )))
                measureNames$measureId <- measure$id
            }
            toRet$measureNames <- measureNames
            ## Attributes
            toRet$measureAttributes <- getAttrSet(x, attToGet=c("Type", "integerValue", "Change"))
            if(!is.null(toRet$measureAttributes)){
                toRet$measureAttributes$measureId <- measure$id
            }
            ## Cytogenetic locations
            toTake <- which(names(x)=="CytogeneticLocation")
            if(length(toTake)==0){
                cytogeneticLocations <- NULL
            }else{
                cytogeneticLocations <- unique(do.call(rbind, lapply(
                    x[toTake],
                    function(y){
                        return(data.frame(
                            cytogenicLocation=xmlValue(y),
                            stringsAsFactors=F
                        ))
                    }
                )))
                cytogeneticLocations$measureId <- measure$id
            }
            toRet$cytogeneticLocations <- cytogeneticLocations
            ## Sequence locations
            sequenceLocations <- getSeqLoc(x)
            if(!is.null(sequenceLocations)){
                sequenceLocations$measureId <- measure$id
            }
            toRet$sequenceLocations <- sequenceLocations
            ## XREF
            xref <- getXRef(x)
            if(!is.null(xref)){
                xref$measureId <- measure$id
            }
            toRet$xref <- xref
            ## MeasureRelationships
            toRet$measureRelationships <- getMeasureRelationship(x, measure$id)
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
    rownames(toRet$measures) <- c()
    toRet$measures$rcvaId <- rcvaId
    
    toRet$measureNames <- do.call(rbind, lapply(
        measureSet,
        function(x) x$measureNames
    ))
    rownames(toRet$measureNames) <- c()
    
    toRet$measureAttributes <- do.call(rbind, lapply(
        measureSet,
        function(x) x$measureAttributes
    ))
    rownames(toRet$measureAttributes) <- c()
    
    toRet$measureCytogeneticLocations <- do.call(rbind, lapply(
        measureSet,
        function(x) x$cytogeneticLocations
    ))
    rownames(toRet$measureCytogeneticLocations) <- c()
    
    toRet$measureSequenceLocations <- do.call(rbind, lapply(
        measureSet,
        function(x) x$sequenceLocations
    ))
    rownames(toRet$measureSequenceLocations) <- c()
    
    toRet$measureXRef <- do.call(rbind, lapply(
        measureSet,
        function(x) x$xref
    ))
    rownames(toRet$measureXRef) <- c()
    
    toRet$measureRelationships <- do.call(rbind, lapply(
        measureSet,
        function(x) x$measureRelationships
    ))
    rownames(toRet$measureRelationships) <- c()
    
    ##
    return(toRet)
}

############################
getSeqLoc <- function(
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
                toRet <- as.data.frame(
                    t(xmlAttrs(x)[attToGet]),
                    stringsAsFactors=F
                )
                colnames(toRet) <- attToGet
                return(toRet)
            }
        )))
        rownames(toRet) <- c()
        return(toRet)
    }
}

############################
getXRef <- function(node){
    attToGet <- c("DB", "ID", "Type")
    toTake <- which(names(node)=="XRef")
    if(length(toTake)==0){
        return(NULL)
    }else{
        toRet <- unique(do.call(rbind, lapply(
            node[toTake],
            function(x){
                toRet <- as.data.frame(
                    t(xmlAttrs(x)[attToGet]),
                    stringsAsFactors=F
                )
                colnames(toRet) <- attToGet
                return(toRet)
            }
        )))
        rownames(toRet) <- c()
        return(toRet)
    }
}

############################
getMeasureRelationship <- function(node, measureId){
    toTake <- which(names(node)=="MeasureRelationship")
    if(length(toTake)==0){
        return(NULL)
    }else{
        toRet <- unique(do.call(rbind, lapply(
            node[toTake],
            function(x){
                toRet <- data.frame(
                    id="1",
                    type=xmlAttrs(x)["Type"],
                    stringsAsFactors=F
                )
                ##
                toTake <- which(names(x)=="Name")
                if(length(toTake)==0){
                    toAdd <- data.frame(
                        id="1",
                        name=NA,
                        name.type=NA
                    )
                }else{
                    toAdd <- unique(do.call(rbind, lapply(
                        x[toTake],
                        function(y){
                            return(data.frame(
                                id="1",
                                name=xmlValue(y[["ElementValue"]]),
                                name.type=xmlAttrs(y[["ElementValue"]])["Type"],
                                stringsAsFactors=F
                            ))
                        }
                    )))
                }
                toRet <- merge(toRet, toAdd, by="id", all=T)
                ##
                toTake <- which(names(x)=="Symbol")
                if(length(toTake)==0){
                    toAdd <- data.frame(
                        id="1",
                        symbol=NA,
                        symbol.type=NA
                    )
                }else{
                    toAdd <- unique(do.call(rbind, lapply(
                        x[toTake],
                        function(y){
                            return(data.frame(
                                id="1",
                                symbol=xmlValue(y[["ElementValue"]]),
                                symbol.type=xmlAttrs(y[["ElementValue"]])["Type"],
                                stringsAsFactors=F
                            ))
                        }
                    )))
                }
                toRet <- merge(toRet, toAdd, by="id", all=T)
                toAdd <- getXRef(x)
                colnames(toAdd) <- c("xref.db", "xref.id", "xref.type")
                toAdd$id <- "1"
                toRet <- merge(toRet, toAdd, by="id", all=T)
                rownames(toRet) <- c()
                toRet <- toRet[, setdiff(colnames(toRet), "id")]
                return(toRet)
            }
        )))
        rownames(toRet) <- c()
        toRet$measureId <- measureId
        return(toRet)
    }
}

############################
parseTraitSet <- function(node, rcvaId){
    ts <- node[["TraitSet"]]
    toTake <- which(names(ts)=="Trait")
    traitSet <- lapply(
        ts[toTake],
        function(x){
            toRet <- list()
            ## Trait identification
            trait <- as.data.frame(
                t(xmlAttrs(x)[c("Type", "ID")]),
                stringsAsFactors=F
            )
            colnames(trait) <- c("type", "id")
            toRet$trait <- trait
            ## Direct xref
            xref <- getXRef(x)
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
                        trait.name.type=xmlAttrs(y[["ElementValue"]])["Type"],
                        trait.name=xmlValue(y[["ElementValue"]]),
                        stringsAsFactors=F
                    )
                    toRet$xref <- getXRef(y)
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
                getXRef
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
    rownames(toRet$traits) <- c()
    toRet$traits$rcvaId <- rcvaId
    
    toRet$traitXRef <- do.call(rbind, lapply(
        traitSet,
        function(x) x$xref
    ))
    rownames(toRet$traitXRef) <- c()
    
    toRet$traitNames <- do.call(rbind, lapply(
        traitSet,
        function(x) x$traitNames
    ))
    rownames(toRet$traitNames) <- c()
    
    ##
    return(toRet)
    
}

############################
parseCva <- function(node, cvsId){
    toTake <- which(names(node)=="ClinVarAssertion")
    toRet <- lapply(
        node[toTake],
        function(x){
            cva <- data.frame(
                cvs=cvsId,
                id=xmlAttrs(x)["ID"],
                accession=xmlAttrs(x[["ClinVarAccession"]])["Acc"],
                clinicalSignificance=xmlValue(x[["ClinicalSignificance"]][["Description"]]),
                stringsAsFactors=F
            )
            submitters <- data.frame(
                submitter=xmlAttrs(x[["ClinVarSubmissionID"]])["submitter"],
                cvaId=cva$id,
                primary=T,
                stringsAsFactors=F
            )
            if("AdditionalSubmitters" %in% names(x)){
                toAdd <- unlist(lapply(
                    x[["AdditionalSubmitters"]][1:length(names(x[["AdditionalSubmitters"]]))],
                    function(y){
                        return(xmlAttrs(y)["SubmitterName"])
                    }
                ))
                toAdd <- data.frame(
                    submitter=toAdd,
                    cvaId=cva$id,
                    primary=F,
                    stringsAsFactors=F
                )
                submitters=rbind(submitters, toAdd)
            }
            observedIn <- parseCvaObservedIn(x, cva$id)
            toRet <- list(cva=cva, submitters=submitters, observedIn=observedIn)
            return(toRet)
        }
    )
    cva <- do.call(rbind, lapply(
        toRet,
        function(x) x$cva
    ))
    submitters <- do.call(rbind, lapply(
        toRet,
        function(x) x$submitters
    ))
    observedIn <- do.call(rbind, lapply(
        toRet,
        function(x) x$observedIn
    ))
    return(list(
        cva=cva,
        submitters=submitters,
        observedIn=observedIn
    ))
}
