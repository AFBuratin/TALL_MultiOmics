cimDiablo = function (object, color = NULL, color.Y, color.blocks, comp = NULL, 
          margins = c(2, 15), legend.position = "topright", transpose = T, 
          row.names = TRUE, col.names = TRUE, size.legend = 1.5, ...) 
{
  if (!is(object, "block.splsda")) 
    stop("cimDiablo is only available for 'block.splsda' objects")
  if (length(object$X) <= 1) 
    stop("This function is only available when there are more than 3 blocks")
  ncomp = min(object$ncomp)
  if (is.null(comp)) {
    comp = 1:ncomp
  }
  if (length(comp) > 1) {
    comp = unique(comp)
    if (!is.numeric(comp) || any(comp < 1)) 
      stop("invalid vector for 'comp'.", call. = FALSE)
    if (any(comp > ncomp)) 
      stop("the elements of 'comp' must be <= ", ncomp, 
           ".", call. = FALSE)
  }
  if (length(comp) == 1) {
    if (is.null(comp) || !is.numeric(comp) || comp <= 0 || 
        comp > ncomp) 
      stop("invalid value for 'comp'.", call. = FALSE)
    comp = c(comp, comp)
  }
  comp = round(comp)
  if (missing(color.Y)) {
    color.Y = color.mixo(1:nlevels(object$Y))
  }
  else {
    if (length(color.Y) != nlevels(object$Y)) 
      stop("'color.Y' needs to be of length ", nlevels(object$Y))
  }
  if (missing(color.blocks)) {
    color.blocks = brewer.pal(n = 12, name = "Paired")[seq(2, 
                                                           12, by = 2)]
  }
  else {
    if (length(color.blocks) != length(object$X)) 
      stop("'color.blocks' needs to be of length ", length(object$X))
  }
  X = object$X
  Y = object$Y
  indY = object$indY
  object$variates = c(object$variates[-indY], object$variates[indY])
  object$loadings = c(object$loadings[-indY], object$loadings[indY])
  object$loadings = lapply(object$loadings, function(x) {
    x[, comp, drop = FALSE]
  })
  keepA = lapply(object$loadings, function(i) apply(abs(i), 
                                                    1, sum) > 0)
  XDatList = mapply(function(x, y) {
    x[, y]
  }, x = X, y = keepA[-length(keepA)], SIMPLIFY = FALSE)
  XDat = do.call(cbind, XDatList)
  XDat[which(XDat > 2)] = 2
  XDat[which(XDat < -2)] = -2
  VarLabels = factor(rep(names(X), lapply(keepA[-length(keepA)], 
                                          sum)), levels = names(X))
  opar = par()[!names(par()) %in% c("cin", "cra", "csi", "cxy", 
                                    "din", "page")]
  par(mfrow = c(1, 1))
  res <- cim(XDat, transpose = transpose, color = color, row.names = row.names, 
             col.names = col.names, col.sideColors = color.blocks[as.numeric(VarLabels)], 
             row.sideColors = color.Y[as.numeric(Y)], margins = margins, clust.method = c("ward.D", "ward.D"),
             dist.method = c("correlation","correlation"), keysize = c(1,1), row.cex = 0.3, save = "png", name.save = "varSel",
             ...)
  if (!transpose) {
    legend(legend.position, c("Rows", c(levels(Y)[order(levels(Y))], 
                                        "", "Columns", names(X))), col = c(1, color.Y, 1, 
                                                                           1, color.blocks[1:nlevels(VarLabels)][match(levels(VarLabels), 
                                                                                                                       names(X))]), pch = c(NA, rep(19, nlevels(Y)), 
                                                                                                                                            NA, NA, rep(19, nlevels(VarLabels))), bty = "n", 
           cex = size.legend, text.font = c(2, rep(1, nlevels(Y)), 
                                            NA, 0.5, rep(0.5, nlevels(VarLabels))))
  }
  else {
    legend(legend.position, c("Rows", names(X), "", "Columns", 
                              c(levels(Y)[order(levels(Y))])), col = c(1, color.blocks[1:nlevels(VarLabels)][match(levels(VarLabels), 
                                                                                                                   names(X))], 1, 1, color.Y), pch = c(NA, rep(19, nlevels(VarLabels)), 
                                                                                                                                                       NA, NA, rep(19, nlevels(Y))), bty = "n", cex = size.legend, 
           text.font = c(2, rep(0.5, nlevels(VarLabels)), NA, 
                         0.5, rep(0.5, nlevels(Y))))
  }
  par(opar)
  return(invisible(res))
}

circosPlot = function (object=diablo.plot, comp = 1:min(object$ncomp), cutoff, color.Y, 
          color.blocks, color.cor, var.names = NULL, showIntraLinks = FALSE, 
          line = FALSE, size.legend = 0.8, ncol.legend = 1, size.variables = 0.25, 
          size.labels = 1, legend = TRUE, linkWidth = 1, ...) 
{
  Features = Exp = Dataset = Mean = linkColors = chrom = po = NULL
  figSize = 800
  segmentWidth = 25
  linePlotWidth = 90
  if (!is(object, "block.splsda")) 
    stop("circosPlot is only available for 'block.splsda' objects")
  if (length(object$X) < 2) 
    stop("This function is only available when there are more than 3 blocks\n    (2 in object$X + an outcome object$Y)")
  if (missing(cutoff)) 
    stop("'cutoff' is missing", call. = FALSE)
  if (missing(color.Y)) {
    color.Y = color.mixo(1:nlevels(object$Y))
  }
  else {
    if (length(color.Y) != nlevels(object$Y)) 
      stop("'color.Y' must be of length ", nlevels(object$Y))
  }
  if (missing(color.blocks)) {
    color.blocks = brewer.pal(n = 12, name = "Paired")
    if (length(object$X) > 6) {
      color.blocks <- colorRampPalette(color.blocks)(2 * 
                                                       length(object$X))
    }
  }
  else {
    if (length(color.blocks) != length(object$X)) 
      stop("'color.blocks' must be of length ", length(object$X))
    color.blocks.adj = adjustcolor(color.blocks, alpha.f = 0.5)
    color.blocks = c(rbind(color.blocks, color.blocks.adj))
  }
  if (missing(color.cor)) {
    color.cor = c(colors()[134], colors()[128])
  }
  else {
    if (length(color.cor) != 2) 
      stop("'color.cor' must be of length 2")
  }
  X = object$X
  Y = object$Y
  indY = object$indY
  object$variates = c(object$variates[-indY], object$variates[indY])
  object$loadings = c(object$loadings[-indY], object$loadings[indY])
  object$ncomp = c(object$ncomp[-indY], object$ncomp[indY])
  sample.X = lapply(object$loadings[-length(object$loadings)], 
                    function(x) {
                      1:nrow(x)
                    })
  if (is.null(var.names)) {
    var.names.list = unlist(sapply(object$loadings[-length(object$loadings)], 
                                   rownames))
  }
  else if (is.list(var.names)) {
    if (length(var.names) != length(object$loadings[-length(object$loadings)])) 
      stop.message("var.names", sample.X)
    if (sum(sapply(1:length(var.names), function(x) {
      length(var.names[[x]]) == length(sample.X[[x]])
    })) != length(var.names)) 
      stop.message("var.names", sample.X)
    var.names.list = var.names
  }
  else {
    stop.message("var.names", sample.X)
  }
  if (any(comp > min(object$ncomp))) {
    warning("Limitation to ", min(object$ncomp), " components, as determined by min(object$ncomp)")
    comp[which(comp > min(object$ncomp))] = min(object$ncomp)
  }
  comp = unique(sort(comp))
  invalid.linkWidth <- FALSE
  if (mode(linkWidth) == "numeric") {
    if (length(linkWidth) == 1) 
      linkWidth <- rep(linkWidth, 2)
    else if (length(linkWidth) != 2) 
      invalid.linkWidth <- TRUE
    linkWidth <- sort(linkWidth)
  }
  else {
    invalid.linkWidth <- TRUE
  }
  if (invalid.linkWidth) 
    stop("'linkWidth' must be a numeric of length 2 (or 1) specifying ", 
         "the range of widths used for link lines based on similarity measures.", 
         call. = FALSE)
  keepA = lapply(object$loadings, function(i) apply(abs(i)[, 
                                                           comp, drop = FALSE], 1, sum) > 0)
  cord = mapply(function(x, y, keep) {
    cor(x[, keep], y[, comp], use = "pairwise")
  }, x = object$X, y = object$variates[-length(object$variates)], 
  keep = keepA[-length(keepA)], SIMPLIFY = FALSE)
  simMatList = vector("list", length(X))
  for (i in 1:length(cord)) {
    for (j in 1:length(cord)) {
      simMatList[[i]][[j]] = cord[[i]] %*% t(cord[[j]])
    }
  }
  simMat = do.call(rbind, lapply(simMatList, function(i) do.call(cbind, 
                                                                 i)))
  Xdat = as.data.frame(do.call(cbind, X)[, colnames(simMat)])
  AvgFeatExp0 <- mutate(.data = Xdat, Y = Y)
  AvgFeatExp0 <- gather(data = AvgFeatExp0, Features, Exp, 
                        -Y)
  AvgFeatExp0 <- group_by(.data = AvgFeatExp0, Y, Features)
  AvgFeatExp0 <- summarise(.data = AvgFeatExp0, Mean = mean(Exp, 
                                                            na.rm = TRUE), SD = sd(Exp, na.rm = TRUE))
  AvgFeatExp0$Dataset <- factor(rep(names(X), unlist(lapply(cord, 
                                                            nrow))), levels = names(X))[match(AvgFeatExp0$Features, 
                                                                                              colnames(Xdat))]
  featExp <- group_by(.data = AvgFeatExp0, Dataset, Y)
  featExp <- arrange(.data = featExp, Mean)
  chr = genChr(featExp, color.blocks = color.blocks)
  Xblocks <- unique(paste0("chr", names(object$X)))
  chr$block <- factor(chr$chrom, levels = Xblocks, ordered = TRUE)
  chr <- chr[order(chr$block), ]
  chr$block <- NULL
  chr.names = unique(chr$chrom)
  db = segAnglePo(chr, seg = chr.names)
  db = data.frame(db)
  links = genLinks(chr, simMat, threshold = cutoff)
  if (nrow(links) < 1) 
    warning("Choose a lower correlation threshold to highlight\n    links between datasets")
  circleR = (figSize/2) - segmentWidth - linePlotWidth
  linksR = circleR - segmentWidth
  linePlotR = circleR + segmentWidth
  chrLabelsR = (figSize/2)
  ind.match = match(chr$name, unlist(sapply(object$loadings[-length(object$loadings)], 
                                            rownames)))
  chr$name.user = unlist(var.names.list)[ind.match]
  opar1 = par("mar")
  par(mar = c(2, 2, 2, 2))
  plot(c(1, figSize), c(1, figSize), type = "n", axes = FALSE, 
       xlab = "", ylab = "", main = "")
  drawIdeogram(R = circleR, cir = db, W = segmentWidth, show.band.labels = TRUE, 
               show.chr.labels = TRUE, chr.labels.R = chrLabelsR, chrData = chr, 
               size.variables = size.variables, size.labels = size.labels, 
               color.blocks = color.blocks, line = line, ...)
  if (nrow(links) > 0) 
    drawLinks(R = linksR, cir = db, mapping = links, col = linkColors, 
              lineWidth = linkWidth, drawIntraChr = showIntraLinks, 
              color.cor = color.cor)
  cTypes = levels(Y)
  lineCols = color.Y
  if (line == TRUE) {
    for (i in 1:length(chr.names)) {
      seg.name = gsub("chr", "", chr.names[i])
      expr = subset(featExp, featExp$Dataset == seg.name)
      expr = dcast(expr, formula = Features ~ Y, value.var = "Mean")
      expr = merge(expr, chr, by.x = "Features", by.y = "name")
      expr$po = (as.numeric(expr$chromStart) + as.numeric(expr$chromEnd))/2
      expr = rename(expr, seg.name = chrom, seg.po = po)
      cOrder = c(c(grep("seg.name", colnames(expr)), grep("seg.po", 
                                                          colnames(expr))), c(1:length(cTypes) + 1))
      expr = expr[, cOrder]
      subChr = subset(db, db$seg.name == chr.names[i])
      drawLinePlot(R = linePlotR, cir = subChr, W = linePlotWidth, 
                   lineWidth = 1, mapping = expr, col = lineCols, 
                   scale = FALSE)
    }
  }
  opar = par("xpd")
  par(xpd = TRUE)
  if (legend == TRUE) {
    legend(x = 5, y = (circleR/4), title = "Correlations", 
           c("Positive Correlation", "Negative Correlation"), 
           col = color.cor, pch = 19, cex = size.legend, bty = "n")
    if (line == TRUE) 
      legend(x = figSize - (circleR/3), y = (circleR/3), 
             title = "Expression", legend = levels(Y), col = lineCols, 
             pch = 19, cex = size.legend, bty = "n", ncol = ncol.legend)
    legend(x = figSize - (circleR/2), y = figSize, title = "Correlation cut-off", 
           legend = paste("r", cutoff, sep = "="), col = "black", 
           cex = size.legend, bty = "n")
    legend(x = -circleR/4, y = figSize, legend = paste("Comp", 
                                                       paste(comp, collapse = "-")), col = "black", cex = size.legend, 
           bty = "n")
  }
  par(xpd = opar, mar = opar1)
  return(invisible(simMat))
}