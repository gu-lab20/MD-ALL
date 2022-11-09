#' get_embeding_feature
#'
#' @param obj_in
#' @param assay_name_in
#' @param features
#' @param reduction
#'
#' @return
#' @export
#'
#' @examples
get_embeding_feature=function(obj_in,assay_name_in="vst",features="cluster_Phenograph_pca1",reduction="tsne"){
  if(!all(row.names(colData(obj_in$SE))==row.names(obj_in[reduction]))){stop("Sample names in colData and reduction not match")}

  df_feature=get_features_df(obj_in = obj_in,assay_name_in = assay_name_in,features = features)

  df_reduction=obj_in[[reduction]]
  df_out=bind_cols(df_feature,df_reduction)

  df_out
}
