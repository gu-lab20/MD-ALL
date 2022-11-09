#' get_features_df
#'
#' @param obj_in
#' @param assay_name_in
#' @param features
#'
#' @return
#' @export
#'
#' @examples
get_features_df=function(obj_in,assay_name_in="vst",features="cluster_Phenograph_pca1"){
  features1=features[features %in% names(colData(obj_in$SE))]
  if(length(features1) >=1){
    df_feature1=as.data.frame(colData(obj_in$SE)[features1])
    df_feature=df_feature1
  }

  features2=features[features %in% rownames(obj_in$SE)]
  if(length(features2) >=1){
    df_feature2=as.data.frame(t(assays(obj_in$SE[features2,])[[assay_name_in]]))
    df_feature=df_feature2
  }

  if(length(features1) >=1 & length(features2) >=1){
    if(!all(row.names(df_feature1)==row.names(df_feature2))){stop("df feature rownames not match")}
    df_feature=bind_cols(df_feature1,df_feature2)
  }

  df_feature
}
