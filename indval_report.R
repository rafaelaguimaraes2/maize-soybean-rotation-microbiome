indval_report <- function(spe, iva, alpha = 0.05, p.adj.method = "BH") {
  # Correction of the p-values for multiple testing
  pval.adj <- p.adjust(iva$pval, method = p.adj.method)
  # Table of the significant indicator species
  gr <- iva$maxcls[pval.adj <= alpha]
  iv <- iva$indcls[pval.adj <= alpha]
  pv <- iva$pval[pval.adj <= alpha]
  fr <- apply(spe > 0, 2, sum)[pval.adj <= alpha]
  fidg <- data.frame(
    group = dimnames(iva$indval)[[2]][gr],
    indval = iv,
    pvalue = pv,
    freq = fr
  )
  fidg <- fidg[order(fidg$group,-fidg$indval), ]
  fidg
}
