#!/usr/bin/env Rscript

.add.index <- function(x, index.col, drop.col = NULL) {
  if (is.null(drop.col))
    drop.col <- index.col
  
  index <- do.call(paste, c(x[, index.col], sep = ":"))
  x <- x[, !(colnames(x) %in% drop.col)]
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

# index.col order matters!
# returns a matrix
read.FORMAT <- function(file,
                        index.col = NULL,
                        drop.col = NULL) {
  FORMAT <- read.delim(file, na.strings = ".", check.names = FALSE)
  
  if (!is.null(index.col)) {
    .add.index(FORMAT, index.col, drop.col = drop.col)
  } else {
    as.matrix(FORMAT)
  }
}


calculate.PLAF <- function(GT, REF.AD, ALT.AD) {
  if (any(rownames(GT) != rownames(REF.AD)) ||
      any(rownames(REF.AD) != rownames(ALT.AD)))
    stop("Different row names between GT, REF.AD, and ALT.AD!")
  
  if (any(colnames(GT) != colnames(REF.AD)) ||
      any(colnames(REF.AD) != colnames(ALT.AD)))
    stop("Different column names between GT, REF.AD, and ALT.AD!")
  
  DP <- REF.AD + ALT.AD
  
  is.missing <- is.na(GT)
  
  ALT.AD[is.missing] <- NA
  DP[is.missing] <- NA
  
  rowSums(ALT.AD) / rowSums(DP)
}


write.PLAF <- function(PLAF, file) {
  CHROM.POS <- strsplit(names(PLAF), split = ":", fixed = TRUE)
  CHROM <- sapply(CHROM.POS, "[", 1)
  POS <- sapply(CHROM.POS, "[", 2)
  
  write.table(
    data.frame(CHROM, POS, PLAF),
    file,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )
}


write.PLAF.sites <- function(PLAF, file) {
  CHROM.POS <- strsplit(names(PLAF), split = ":", fixed = TRUE)
  CHROM <- sapply(CHROM.POS, "[", 1)
  POS <- sapply(CHROM.POS, "[", 2)
  
  write.table(
    cbind(CHROM, POS),
    file,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
  )
}


helper <-
  function(GT,
           REF.AD,
           ALT.AD,
           PLAF.file.prefix,
           samples = NULL) {
    if (!is.null(samples)) {
      GT <- GT[, samples]
      REF.AD <- REF.AD[, samples]
      ALT.AD <- ALT.AD[, samples]
    }
    
    message(paste0("Stats for: ", PLAF.file.prefix))
    PLAF <- calculate.PLAF(GT, REF.AD, ALT.AD)
    
    message(paste0("Number of sites left: ", length(PLAF)))
    PLAF <- PLAF[!is.na(PLAF)]
    message(paste0("Number of sites left with no missingness: ", length(PLAF)))
    
    
    # drop monomorphic sites
    PLAF <- PLAF[PLAF > 0]
    PLAF <- PLAF[PLAF < 1]
    
    message(paste0("Number of polymorphic sites: ", length(PLAF)))
    
    PLAF.file <- paste0(PLAF.file.prefix, ".PLAF.txt")
    write.PLAF(PLAF, PLAF.file)
    
    PLAF.sites.file <- paste0(PLAF.file.prefix, ".PLAF.sites.txt")
    write.PLAF.sites(PLAF, PLAF.sites.file)
  }


#AF.workflow <- function() {
#  index.col <- c("CHROM", "POS")
#  
#  GT.file <- commandArgs(trailingOnly = TRUE)
#  REF.AD.file <- sub("GT", "AD.0", GT.file, fixed = TRUE)
#  ALT.AD.file <- sub("GT", "AD.1", GT.file, fixed = TRUE)
#
#  GT <- read.GT(GT.file, index.col = index.col)
#  REF.AD <- read.FORMAT(REF.AD.file, index.col = index.col)
#  ALT.AD <- read.FORMAT(ALT.AD.file, index.col = index.col)
#  
#  infix <- sub(".*/", "", GT.file, perl = TRUE)
#  infix <- sub(".GT.txt.gz", "", infix, fixed = TRUE)
#  
#  PLAF.file.prefix <- paste0("DEploid_input/", infix)
  #PLAF.file.prefix <- paste0(dirname(GT.file), "/", infix)
#
#    
#  helper(GT, REF.AD, ALT.AD, PLAF.file.prefix)
#}

#AF.workflow()
AF.workflow <- function(GT.file, REF.AD.file, ALT.AD.file) {
  index.col <- c("CHROM", "POS")
  
  GT <- read.GT(GT.file, index.col = index.col)
  REF.AD <- read.FORMAT(REF.AD.file, index.col = index.col)
  ALT.AD <- read.FORMAT(ALT.AD.file, index.col = index.col)
  
  infix <- sub(".*/", "", GT.file, perl = TRUE)
  infix <- sub(".GT.txt.gz", "", infix, fixed = TRUE)
  
  PLAF.file.prefix <- paste0("DEploid_input/", infix)
  
  helper(GT, REF.AD, ALT.AD, PLAF.file.prefix)
}

AF.workflow("merged.bi.filt.GT.miss0.4.snps.GT.txt.gz", "merged.bi.filt.GT.miss0.4.snps.AD.0.txt.gz", "merged.bi.filt.GT.miss0.4.snps.AD.1.txt.gz")
