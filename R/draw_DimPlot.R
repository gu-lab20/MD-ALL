#' draw_DimPlot
#'
#' @param obj_in
#' @param group.by
#' @param cols
#' @param axis_by
#' @param reduction
#'
#' @return
#' @export
#'
#' @examples
draw_DimPlot=function(obj_in,group.by="diag",cols=subtypeCol,axis_by=10,reduction="tsne"){

  df_plot=get_embeding_feature(obj_in = obj_in,feature = group.by,reduction = reduction)

  if(reduction=="tsne"){
    gg_dimPlot_RNAseqtsne(df_plot,cols = cols,axis_by = axis_by)
  }

}
