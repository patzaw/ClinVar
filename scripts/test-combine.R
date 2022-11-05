odir <- "scripts/now"
saved <- list.files(odir, pattern="^subcvs-[[:digit:]]+[.]rds$")
savedi <- as.numeric(sub("subcvs-", "", sub( "[.]rds", "", saved)))
saved <- saved[order(savedi)]
savedi <- savedi[order(savedi)]
library(parallel)
a <- do.call(c, mclapply(file.path(odir, saved), readRDS, mc.cores=15))
