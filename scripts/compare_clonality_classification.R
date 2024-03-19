fws <- read.delim("metadata.tsv")

Fws.clonality <-
  data.frame(sample = fws[, "sample"], Fws_0.90 = fws[, "Fws"] > 0.90)

props <- read.table("run_DEploid/filtered_non-swga_AF.DEploid.props")

DEploid.clonality <-
  data.frame(sample = props[, 1],
             DEploid_0.99 = apply(props[, 2:length(props)] > 0.99, 1, any))


meta <- merge(Fws.clonality, DEploid.clonality)
meta <- meta[match(fws[, "sample"], meta[, "sample"]), ]

print("Any clonality classification difference between Fws and DEploid?")
print(table(meta[, "Fws_0.90"], meta[, "DEploid_0.99"]))

print("Any samples that are polyclonal by DEploid classification?")
print(meta[meta[, "Fws_0.90"] == TRUE &
             meta[, "DEploid_0.99"] == FALSE, "sample"])

print("Any samples that are polyclonal by Fws classification?")
print(meta[meta[, "Fws_0.90"] == FALSE &
             meta[, "DEploid_0.99"] == TRUE, "sample"])

# mark as polyclonals if either Fws or DEploid classified as polyclonal
polyclonals.index <- !meta[, "Fws_0.90"] | !meta[, "DEploid_0.99"]
polyclonals <- meta[polyclonals.index, "sample"]
monoclonals <- meta[!polyclonals.index, "sample"]


write.table(
  polyclonals,
  file = "run_DEploid/filtered_non-swga_AF.DEploid.polyclonals.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
write.table(
  monoclonals,
  file = "run_DEploid/filtered_non-swga_AF.DEploid.monoclonals.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
