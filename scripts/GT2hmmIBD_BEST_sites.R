#!/usr/bin/env Rscript

# Format for genotype file:
# tab-delimited text file, with one single nucleotide polymorphism
# (SNP) per line. The first two columns are the chromosome and
# position, followed by one sample per column. A header line,
# giving the sample names, is required.
#
# Genotypes are coded by number:
# -1 for missing data, 0 for the first allele, 1 for the second, etc.
#
# SNPs and indels (if you trust them) can thus be
# treated on an equal footing. The variants must be in
# chromosome and position order, and can have between
# two and eight alleles (more, if you feel like changing
# max_allele in the code).

.add.index <- function(x, index.col, drop.col = NULL) {
  if (is.null(drop.col))
    drop.col <- index.col
  
  index <- do.call(paste, c(x[, index.col], sep = ":"))
  x <- x[, !(colnames(x) %in% drop.col), drop = FALSE]
  rownames(x) <- index
  
  as.matrix(x)
}

# assumes unphased data
# index.col order matters!
# returns a matrix
read.GT <- function(file,
                    index.col = NULL,
                    drop.col = NULL) {
  GT <- read.delim(file, na.strings = "./.", check.names = FALSE)
  
  if (!is.null(index.col)) {
    .add.index(GT, index.col, drop.col = drop.col)
  } else {
    as.matrix(GT)
  }
}


construct.hmmIBD <- function(GT, simplify = TRUE) {
  CHROM.POS <- strsplit(rownames(GT), split = ":", fixed = TRUE)
  POS <- vapply(CHROM.POS, FUN = "[", FUN.VALUE = character(1), 2)
  POS <- as.numeric(POS)
  
  # hmmIBD does not like non-integer chromosomes
  CHROM <- vapply(CHROM.POS, FUN = "[", FUN.VALUE = character(1), 1)
  CHROM <-
    vapply(
      strsplit(CHROM, split = "_", fixed = TRUE),
      FUN = "[",
      FUN.VALUE = character(1),
      2
    )
  # remove leading zeroes
  CHROM <-
    sub(
      pattern = "^0*",
      replacement = "",
      x = CHROM,
      perl = TRUE
    )
  CHROM <- as.numeric(CHROM)
  
  rownames(GT) <- NULL
  
  if (simplify) {
    # assumes that there are no mixed infections
    GT[GT == "1/1"] <- 1
    GT[GT == "0/0"] <- 0
  }
  
  cbind(CHROM, POS, GT)
}


write.hmmIBD <- function(hmmIBD, file)
  write.table(
    hmmIBD,
    file = file,
    quote = FALSE,
    sep = "\t",
    na = "-1",
    row.names = FALSE
  )


index.col <- c("CHROM", "POS")

samples.file <-
  "run_DEploid/filtered_non-swga_AF.DEploid.polyclonals.txt"

samples <- scan(file = samples.file,
                what = character(),
                quiet = TRUE)
    
panel.file <-
  "BEST/filtered_non-swga_AF_monoclonals.DEploid.panel.txt.gz"
panel <- read.GT(panel.file, index.col = index.col)
hmmIBDs <- construct.hmmIBD(panel, simplify = FALSE)

for (sample in samples) {
  message("Processing ", sample)

  GT.file <-
    paste0("BEST/",
           sample,
           ".final.hap")
  GT <- read.GT(GT.file, index.col = index.col)
  GT <- as.data.frame(GT)

  prop.file <-
    paste0("BEST/",
           sample,
           ".ibd.prop")

  if (file_test("-f", prop.file)) {
    prop <- read.table(prop.file)
    available.hap <- tail(prop, n = 1)
    used.hap <- available.hap[available.hap >= 0.01]
    hap.order <- order(used.hap, decreasing = TRUE)

    GT <- GT[, hap.order, drop = FALSE]
    colnames(GT) <- paste0(sample, "-", seq_along(GT))

  } else {
    hap.order <- 1

    GT <- GT[, hap.order, drop = FALSE]
    colnames(GT) <- sample
  }

  hmmIBD <- construct.hmmIBD(GT, simplify = FALSE)
  hmmIBDs <- merge(hmmIBDs, hmmIBD, by = index.col, all = TRUE, sort = FALSE)
}

hmmIBDs <- hmmIBDs[do.call(order, hmmIBDs[, index.col]), ]
hmmIBDs[, "POS"] <- as.integer(hmmIBDs[, "POS"])
hmmIBDs.file <-
  "BEST/filtered_non-swga_AF.DEploid.hmmIBD.txt"
write.hmmIBD(hmmIBDs, hmmIBDs.file)
