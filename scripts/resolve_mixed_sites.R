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


# assumes biallelic, unphased data
resolve.mixed <- function(GT, REF.AD, ALT.AD) {
  if (any(rownames(GT) != rownames(REF.AD)) ||
      any(rownames(REF.AD) != rownames(ALT.AD)))
    stop("Different row names between GT, REF.AD, and ALT.AD!")
  
  if (any(colnames(GT) != colnames(REF.AD)) ||
      any(colnames(REF.AD) != colnames(ALT.AD)))
    stop("Different column names between GT, REF.AD, and ALT.AD!")
  
  
  # remember homozygous GTs
  GT.missing <- is.na(GT)
  GT.ref <- GT == "0/0"
  GT.alt <- GT == "1/1"
  
  
  # the highest AD is the most likely GT
  alt.major <- REF.AD < ALT.AD
  
  GT[alt.major] <- "1/1"
  GT[!alt.major] <- "0/0"
  
  
  # restore homozygous GTs
  GT[GT.missing] <- NA
  GT[GT.ref] <- "0/0"
  GT[GT.alt] <- "1/1"
  
  
  GT
}


save.GT <- function(GT, file, compress = TRUE) {
  if (compress)
    file <- gzfile(file, open = "wb")
  
  write.table(
    GT,
    file = file,
    quote = FALSE,
    sep = "\t",
    na = "./."
  )
  
  if (compress)
    close(file)
}


index.col <- c("CHROM", "POS")

GT.file <- "DEploid_input/merged.bi.filt.GT.miss0.4.snps.GT.txt.gz"
GT <- read.GT(GT.file, index.col = index.col)

# forgot about phasing
GT <- sub("|", "/", GT, fixed = TRUE)

REF.AD.file <- "DEploid_input/merged.bi.filt.GT.miss0.4.snps.AD.0.txt.gz"
REF.AD <- read.FORMAT(REF.AD.file, index.col = index.col)
  
ALT.AD.file <- "DEploid_input/merged.bi.filt.GT.miss0.4.snps.AD.1.txt.gz"
ALT.AD <- read.FORMAT(ALT.AD.file, index.col = index.col)

# subset to DEploid sites and monoclonals
sites <-
  read.delim("DEploid_input/merged.bi.filt.GT.miss0.4.snps.PLAF.sites.txt",
             header = FALSE,
             col.names = c("CHROM", "POS"))
sites <- do.call(paste, c(sites, sep = ":"))

monoclonals <-
  scan(file = "run_DEploid/filtered_non-swga_AF.DEploid.monoclonals.txt",
       what = character(),
       quiet = TRUE)
GT <- GT[sites, monoclonals]
REF.AD <- REF.AD[sites, monoclonals]
ALT.AD <- ALT.AD[sites, monoclonals]

major.GT <- resolve.mixed(GT, REF.AD, ALT.AD)
  
major.GT.file <- "BEST/filtered_non-swga_AF_monoclonals.DEploid.major.GT.txt.gz"
save.GT(major.GT, major.GT.file)
