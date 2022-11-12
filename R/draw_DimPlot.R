#' draw_DimPlot
#'
#' @param obj_in
#' @param group.by
#' @param cols
#' @param axis_by
#' @param reduction
#' @param highlightLevel
#' @param sizeHighlight
#' @param plot_title
#'
#' @return
#' @export
#'
#' @examples
draw_DimPlot=function(obj_in,group.by="diag",cols=subtypeCol(),axis_by=10,reduction="tsne",
                      plot_title=NULL,
                      highlightLevel=NULL,sizeHighlight=2){

  # if(!reduction %in% c("tsne","umap")){stop("Reduction not defined")}

  df_plot=get_embeding_feature(obj_in = obj_in,feature = group.by,reduction = reduction)

  if(reduction=="tsne"){
    if(is.null(plot_title)){plot_title="tSNE"}
    x_lab=paste0("tSNE_1\nGene N = ", unique(df_plot$FeatureN), "; Perplexity N=", unique(df_plot$perplexityN))

    p1=gg_dimPlot(df_in=df_plot,x="tSNE_1",y="tSNE_2",var_col="diag",cols=cols,
                  size=1,axis_by=axis_by,
                  highlightLevel=highlightLevel,sizeHighlight=sizeHighlight,
                  plot_title=plot_title,x_lab=x_lab,y_lab=NULL,legend_title="Group")  }

  if(grepl("umap",reduction)){
    if(is.null(plot_title)){plot_title="uMAP"}
    x_lab=paste0("uMAP_1\nGene N = ", unique(df_plot$FeatureN), "; Neighbor N=", unique(df_plot$n_neighbors))

    p1=gg_dimPlot(df_in=df_plot,x="uMAP_1",y="uMAP_2",var_col="diag",cols=cols,
                  size=1,axis_by=axis_by,
                  highlightLevel=highlightLevel,sizeHighlight=sizeHighlight,
                  plot_title=plot_title,x_lab=x_lab,y_lab=NULL,legend_title="Group")
  }
  p1
}
