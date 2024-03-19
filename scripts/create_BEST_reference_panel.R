major.GT.file <- "BEST/filtered_non-swga_AF_monoclonals.DEploid.major.GT.txt.gz"
major.GT <- read.delim(major.GT.file, na.strings = "./.", check.names = FALSE)

# simplify the genotype format
# assume no mixed infections present
major.GT[major.GT == "1/1"] <- 1
major.GT[major.GT == "0/0"] <- 0

panel <- major.GT

CHROM.POS <- strsplit(rownames(panel), split = ":", fixed = TRUE)
CHROM <- vapply(CHROM.POS, FUN = "[", FUN.VALUE = character(1), 1)
POS <- vapply(CHROM.POS, FUN = "[", FUN.VALUE = character(1), 2)

panel.file <- "BEST/filtered_non-swga_AF_monoclonals.DEploid.panel.txt.gz"
panel.handle <- gzfile(panel.file, open = "wb")
write.table(
  cbind(CHROM, POS, panel),
  panel.handle,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)
close(panel.handle)
