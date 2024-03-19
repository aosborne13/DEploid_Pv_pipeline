library(ape)

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


GT.as.distance <-
  function(GT,
           method = "percentage",
           samples.in.column = TRUE) {
    if (samples.in.column)
      GT <- t(GT)
    
    dist.gene(GT, method = method)
  }


.get.labels <-
  function(distance,
           metadata,
           label.col,
           sample.col = "sample") {
    sample.order <- match(metadata[, sample.col], labels(distance))
    actual.length <- length(na.omit(sample.order))
    
    metadata <- metadata[order(sample.order), ]
    metadata <- metadata[1:actual.length, ]
    
    labels <- metadata[, label.col]
    labels[is.na(labels)] <- NaN
    
    factor(labels)
  }


.generate.label.palette <-
  function(labels, palette = "Tableau 10") {
    label.names <- levels(labels)
    label.palette <-
      palette.colors(length(label.names), palette = "Tableau 10")
    names(label.palette) <- label.names
    
    label.palette
  }


.generate.label.colours <-
  function(labels, label.palette)
    label.palette[labels]


plot.PCoA <-
  function(file,
           distance,
           metadata,
           label.col,
           sample.col = "sample",
           label.palette = NULL,
           legend.title = NA) {
    labels <- .get.labels(distance, metadata, label.col)
    if (is.null(label.palette))
      label.palette <- .generate.label.palette(labels)
    label.colours <- .generate.label.colours(labels, label.palette)
    
    PCoA <- pcoa(distance)
    vectors <- PCoA[["vectors"]]
    
    # expected variance under the broken stick model
    explained.variance <- PCoA[["values"]][, "Broken_stick"] * 100
    
    png(
      file,
      width = 7 + 2,
      height = 7,
      units = "in",
      res = 300,
      type = "cairo"
    )
    
    # environment setting to enable plotting legends outside of the plot
    par(mar = c(5, 4, 4, 11) + 0.1)
    
    plot(
      vectors[, 1],
      vectors[, 2],
      xlab = paste0("Coordinate 1 (", round(explained.variance[1], digits = 2), "%)"),
      ylab = paste0("Coordinate 2 (", round(explained.variance[2], digits = 2), "%)"),
      col = label.colours
    )
    
    legend(
      "bottomright",
      legend = levels(labels),
      fill = label.palette,
      title = legend.title,
      inset = c(-0.3, 0),
      xpd = TRUE
    )
    
    dev.off()
  }


.compute.parents.colours <-
  function(edge.colours, edges, mixed.colour) {
    # determine what the parents' colour would be based on
    # the colours passed on from the children
    tapply(edge.colours[, 1], edges[, 1], function(x) {
      # retrieve children's colours
      unique.colours <- na.omit(unique(x))
      
      # if the children's colours are mixed
      if (mixed.colour %in% unique.colours)
        return(mixed.colour)
      
      unique.length <- length(unique.colours)
      
      # if no information about the children's colours
      # we'll return back to this parent later
      if (unique.length == 0)
        return(NA)
      
      # if the children only have one colour
      if (unique.length == 1)
        return(unique.colours)
      
      # if the children have multiple colours
      mixed.colour
    })
  }


.update.parents.colours <-
  function(edge.colours, edges, mixed.colour) {
    parents.colours <-
      .compute.parents.colours(edge.colours, edges, mixed.colour)
    
    parent.as.child <- match(edges[, 2], names(parents.colours))
    
    # if the parents are children of other parents,
    # carry the colour of parents as the child's colour
    edge.colours[, 2] <-
      mapply(function(x, y) {
        if (is.na(y)) {
          return(x)
        } else {
          return(y)
        }
      }, x = edge.colours[, 2],
      y = parents.colours[parent.as.child])
    
    edge.colours[, 1] <- edge.colours[, 2]
    
    edge.colours
  }


.colour.NJ <- function(NJ, tip.colours) {
  edges <- NJ[["edge"]]
  
  # copy available tip colours as edge colours
  edge.colours <- tip.colours[edges]
  
  # copy colours on the tip to its direct parent
  dim(edge.colours) <- dim(edges)
  edge.colours[, 1] <- edge.colours[, 2]
  
  
  # trace back parent colours until they meet different colours
  mixed.colour <- "#E5E4E2"
  while (sum(is.na(edge.colours))) {
    edge.colours <-
      .update.parents.colours(edge.colours, edges, mixed.colour)
  }
  
  
  # finalise edge colours
  .update.parents.colours(edge.colours, edges, mixed.colour)
}


.plot.NJ <-
  function(file,
           NJ,
           labels = NULL,
           tip.color = par("col"),
           edge.color = NULL,
           tip.palette = NULL,
           legend.title = NA,
           rooted = TRUE,
           show.tip.label = TRUE) {
    if (rooted) {
      type <- "fan"
    } else {
      type <- "unrooted"
    }
    
    png(
      file,
      width = 13 + 1,
      height = 13,
      units = "in",
      pointsize = 9,
      res = 300,
      type = "cairo"
    )
    
    # environment setting to enable plotting legends outside of the plot
    par(mar = c(1, 1, 1, 12) + 0.1)
    
    plot(
      NJ,
      type = type,
      underscore = TRUE,
      lab4ut = "axial",
      show.tip.label = show.tip.label,
      edge.color = edge.color,
      tip.color = tip.color
    )
    
    if (!is.null(tip.palette))
      legend(
        "bottomright",
        legend = levels(labels),
        fill = tip.palette,
        cex = 1.5,
        title = legend.title,
        inset = c(-0.1, 0),
        xpd = TRUE
      )
    
    dev.off()
  }


plot.NJ <-
  function(file,
           distance,
           metadata,
           label.col,
           sample.col = "sample",
           tip.palette = NULL,
           legend.title = NA,
           rooted = TRUE,
           show.tip.label = TRUE,
           NJ = NULL) {
    labels <-
      .get.labels(distance, metadata, label.col, sample.col = sample.col)
    if (is.null(tip.palette))
      tip.palette <- .generate.label.palette(labels)
    tip.color <- .generate.label.colours(labels, tip.palette)
    
    if (is.null(NJ)) {
      NJ <- nj(distance)
    } else {
      message("NJ is supplied, distance will not be used")
    }
    
    edge.color <- .colour.NJ(NJ, tip.color)
    
    .plot.NJ(
      file,
      NJ,
      labels = labels,
      tip.color = tip.color,
      edge.color = edge.color,
      tip.palette = tip.palette,
      legend.title = legend.title,
      rooted = rooted,
      show.tip.label = show.tip.label
    )
  }


plot.NJs <-
  function(output.prefix,
           distance,
           metadata,
           label.col,
           tip.palette,
           legend.title) {
    NJ <- nj(distance)

    output.name <-
      paste0(output.prefix, "_NJ_unrooted_unlabelled.png")
    plot.NJ(
      output.name,
      distance,
      metadata,
      label.col,
      tip.palette = tip.palette,
      legend.title = legend.title,
      rooted = FALSE,
      show.tip.label = FALSE,
      NJ = NJ
    )
    
    output.name <-
      paste0(output.prefix, "_NJ_unrooted_labelled.png")
    plot.NJ(
      output.name,
      distance,
      metadata,
      label.col,
      tip.palette = tip.palette,
      legend.title = legend.title,
      rooted = FALSE,
      show.tip.label = TRUE,
      NJ = NJ
    )
    
    output.name <-
      paste0(output.prefix, "_NJ_rooted_unlabelled.png")
    plot.NJ(
      output.name,
      distance,
      metadata,
      label.col,
      tip.palette = tip.palette,
      legend.title = legend.title,
      rooted = TRUE,
      show.tip.label = FALSE,
      NJ = NJ
    )
    
    output.name <-
      paste0(output.prefix, "_NJ_rooted_labelled.png")
    plot.NJ(
      output.name,
      distance,
      metadata,
      label.col,
      tip.palette = tip.palette,
      legend.title = legend.title,
      rooted = TRUE,
      show.tip.label = TRUE,
      NJ = NJ
    )
  }


metadata <-
  read.delim("metadata.tsv",
             check.names = FALSE,
             encoding = "UTF-8")
  
major.monoclonal.GT.file <-
  "BEST/filtered_non-swga_AF_monoclonals.DEploid.major.GT.txt.gz"
major.monoclonal.GT <- read.GT(major.monoclonal.GT.file)


# give dist.gene less sites to process
major.monoclonal.GT <- major.monoclonal.GT[complete.cases(major.monoclonal.GT), ]
major.monoclonal.GT <- as.matrix(major.monoclonal.GT)

message("Now converting GT to distance...")
major.monoclonal.distance <- GT.as.distance(major.monoclonal.GT)
  



# label by treatment
major.monoclonal.legend.title <- "treatment"
major.monoclonal.label.col <- "treatment"

message("Now plotting PCoA...")
major.monoclonal.PCoA.file <-
  "BEST/filtered_non-swga_AF_monoclonals_DEploid_treatment_PCoA.png"
plot.PCoA(
  major.monoclonal.PCoA.file,
  major.monoclonal.distance,
  metadata,
  major.monoclonal.label.col,
  legend.title = major.monoclonal.legend.title
)

message("Now plotting NJs...")
major.monoclonal.output.prefix <-
  "BEST/filtered_non-swga_AF_monoclonals_DEploid_treatment"
plot.NJs(
  major.monoclonal.output.prefix,
  major.monoclonal.distance,
  metadata,
  major.monoclonal.label.col,
  NULL,
  major.monoclonal.legend.title
)

# label by release
major.monoclonal.legend.title <- "external_id"
major.monoclonal.label.col <- "external_id"

message("Now plotting PCoA...")
major.monoclonal.PCoA.file <-
  "BEST/filtered_non-swga_AF_monoclonals_DEploid_external_id_PCoA.png"
plot.PCoA(
  major.monoclonal.PCoA.file,
  major.monoclonal.distance,
  metadata,
  major.monoclonal.label.col,
  legend.title = major.monoclonal.legend.title
)

message("Now plotting NJs...")
major.monoclonal.output.prefix <-
  "BEST/filtered_non-swga_AF_monoclonals_DEploid_external_id"
plot.NJs(
  major.monoclonal.output.prefix,
  major.monoclonal.distance,
  metadata,
  major.monoclonal.label.col,
  NULL,
  major.monoclonal.legend.title
)
