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
