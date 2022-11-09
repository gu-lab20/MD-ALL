#' gg_dimPlot_RNAseqtsne
#'
#' @param df_plot
#' @param cols
#' @param axis_by
#'
#' @return
#' @export
#'
#' @examples
gg_dimPlot_RNAseqtsne=function(df_plot,cols=subtypeCol,axis_by=10){

  names(df_plot)[grepl("FeatureN",names(df_plot))]

  N=unique(df_plot$FeatureN)
  i=unique(df_plot$perplexityN)
  # N=gsub("FeatureN","",word(names(df_plot)[2],-2,sep="_"))
  # i=gsub("PerplexityN","",word(names(df_plot)[2],-1,sep="_"))
  x_lab=paste0("tSNE dimension 1\nTop variable gene N = ", N, "; Perplexity=", i)
  y_lab="tSNE dimension 2"

  names(df_plot)[1:3]=c("feature","X","Y")

  col_in=get_cols_cat(df_plot$feature,cols_in=cols,cols_in_default)

  ggplot() +     xlab(x_lab) +    ylab(y_lab)+    theme_bw() +
    geom_point(data=df_plot, aes(X, Y, color=feature)) +
    scale_color_manual(values = col_in) +
    scale_x_continuous(breaks = round(seq(floor(min(df_plot$X)), ceiling(max(df_plot$X)), by = axis_by),1)) +
    scale_y_continuous(breaks = round(seq(floor(min(df_plot$Y)), ceiling(max(df_plot$Y)), by = axis_by),1)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.text=element_text(size=15),
          axis.title=element_text(size=15),
          legend.text = element_text( size=15),
          legend.title = element_text( size=15,face="bold")
    )
}
