#' create manhatan plot
#'
#' @param gwa
#' @param adj_method adjustment method, see p.adjust
#' @param thresh_p significant threshold (default 0.05)
#' @param chr [char] filter to show only some chromosomes (show all if NA)
#'
#' @return plotly graph
manPlot <- function(gwa, adj_method, thresh_p = 0.05, chr = NA){
  logger <- logger$new("r-manPlot()")

  # P-Values adjustment
  logger$log("Adjust p-values ...")
  gwa$p_adj <- p.adjust(gwa$p,
                        method = adj_method,
                        n = nrow(gwa))
  maxSinifP <- max(gwa[gwa$p_adj < thresh_p, "p"], thresh_p/nrow(gwa)/100)
  minUnSinifP <- min(gwa[gwa$p_adj > thresh_p, "p"], 1)
  logger$log("Adjust p-values DONE")

  # threshold adjustment
  logger$log("Adjust p-values threshold...")
  thresh_pAdj <- uniroot(
    function(log10p){
      p <- 10^(log10p)
      p.adjust(p,
               method = adj_method,
               n = nrow(gwa)) - thresh_p
    },
    c(log10(maxSinifP), log10(minUnSinifP))
  )
  thresh_pAdj <- 10^(thresh_pAdj$root)
  logger$log("Adjust p-values threshold DONE")

  # filter according to "chr"
  if (!identical(chr, NA)) {
    gwa <- gwa[as.character(gwa$chr) %in% chr,]
  }

  # get significant SNP
  significantSNP <- gwa[gwa$p_adj <= thresh_p, "id"]

  p <- manhattanly(
    data.frame(CHR = gwa$chr,
               BP = gwa$pos,
               SNP = gwa$id,
               P = gwa$p),
    snp = "SNP",
    highlight = significantSNP,
    genomewideline = -log10(thresh_pAdj),
    suggestiveline = FALSE
  )
  p
}
