BiocManager::install(version = "3.23", ask = FALSE)

BiocManager::install(c(
  "optparse",
  "GenomicRanges",
  "GenomicFeatures",
  "rtracklayer",
  "dplyr",
  "stringr", 
  "data.table",
  "parallel",
  "TxDb.Hsapiens.UCSC.hg38.knownGene"
), ask = FALSE)
