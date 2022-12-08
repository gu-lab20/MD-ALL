#' Title get_RNAseqCNV_plot
#'
#' @param RNAseqCNV_out
#' @param adjust
#' @param estimate_lab
#' @param scale_cols
#'
#' @return
#' @export
#'
#' @examples
get_RNAseqCNV_plot=function(RNAseqCNV_out,
                            adjust=TRUE,
                            estimate_lab = TRUE,
                            scale_cols = scaleCols,
                            dpRatioChromEdge = dpRatioChrEdge){
  #get parameters
  sample_name="TestSample"
  gender=as.character(RNAseqCNV_out$df_cnv_out$gender)[1]
  count_ns=RNAseqCNV_out$count_ns
  feat_tab_alt=RNAseqCNV_out$feat_tab_alt
  smpSNPdata=RNAseqCNV_out$smpSNPdata

  #adjust the gene expression according to the estimation of which chromosomes are diploid
  if (adjust == TRUE) {
    count_ns <-  adjust_dipl(feat_tab_alt, count_ns)
  }

  #calculate box plots
  box_wdt <- get_box_wdt(count_ns = count_ns, scaleCols = scale_cols)

  #adjust y axis limits
  ylim <- adjust_ylim(box_wdt = box_wdt, ylim = c(-0.4, 0.4))

  count_ns_final <- prep_expr(count_ns = count_ns, dpRatioChrEdge = dpRatioChromEdge, ylim = ylim)

  # filter low weighted genes for clearer visualization
  #count_ns_final <- filter_expr(count_ns_final = count_ns_final, cutoff = 0.6)


  gg_exp <- plot_exp(count_ns_final = count_ns_final, box_wdt = box_wdt, sample_name = sample_name, ylim = ylim,
                     estimate = estimate_lab, feat_tab_alt = feat_tab_alt, gender = gender)

  gg_snv <- plot_snv(smpSNPdata, sample_name = sample_name, estimate = estimate_lab)

  fig <- arrange_plots(gg_exp = gg_exp, gg_snv = gg_snv)

  fig
}
