library(here)

source(here("scripts/supporting-functions.R"))
plan(multisession, workers=4)

xmlf <- here("sources/ClinVar-head-10E6.xml.gz")

system.time({
   con <- file(xmlf)
   a <- extract_nsn_xml_elements(con=con, n=-1, tag="ClinVarSet")
   b <- future_lapply(a, function(x) as_list(read_xml(x)))
})

system.time({
   d <- process_xml_elements(file=xmlf, tag="ClinVarSet", by=5*10^5)
})

