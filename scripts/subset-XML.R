library(here)

f <- here("sources/ClinVarFullRelease.xml.gz")
n <- 10^6


h <- readLines(f, n=n)
ll  <- tail(grep("<[/]ClinVarSet>", h), 1)
htw <- c(
   h[1:ll],
   "</ReleaseSet>"
)
of <- gzfile(here("sources/ClinVar-head-10E6.xml.gz"))
writeLines(htw, of)
