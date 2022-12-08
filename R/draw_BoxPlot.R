#' draw_BoxPlot
#'
#' @param obj_in
#' @param group.by
#' @param features
#' @param cols
#' @param assay_name_in
#' @param plot_title
#' @param x_lab
#' @param y_lab
#' @param highlightLevel
#' @param sizeHighlight
#'
#' @return
#' @export
#'
#' @examples
draw_BoxPlot=function(obj_in,group.by="diag",features=c("CRLF2"),
                      cols=subtypeCol(),
                      assay_name_in="vst",
                      plot_title="Box Plot",x_lab="Subtype",y_lab=NULL,
                      useGeneName=T,df_geneName=info_gtf_hg38,
                      highlightLevel=NULL,sizeHighlight=1.5){

  if(is.null(y_lab)){y_lab=features}

  if(!useGeneName){
    df_feature=get_features_df(obj_in,assay_name_in=assay_name_in,features=c(features,group.by))
    gg_boxPlot(df_feature,var_value = features,var_group = group.by,plot_title = "BoxPlot",x_lab = x_lab,y_lab = y_lab,highlightLevel = highlightLevel,cols = cols)
  }

  if(useGeneName){
    geneId=df_geneName$gene_id[df_geneName$gene_name %in% features][1]
    df_feature=get_features_df(obj_in,assay_name_in=assay_name_in,features=c(geneId,group.by))
    gg_boxPlot(df_feature,var_value = geneId,var_group = group.by,plot_title = "BoxPlot",x_lab = x_lab,y_lab = y_lab,highlightLevel = highlightLevel,cols = cols)
  }
}
