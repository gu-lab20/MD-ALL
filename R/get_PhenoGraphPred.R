#' get_PhenoGraphPred
#'
#' @param obj_in
#' @param panelName
#' @param SampleLevel
#' @param type
#'
#' @return
#' @export
#'
#' @examples
get_PhenoGraphPred=function(obj_in,panelName="PhenographPred",SampleLevel="TestSample",type="value"){
  df_phenograph=obj_in[[panelName]]
  df_out=df_phenograph[df_phenograph$COH_sample %in% c(SampleLevel),]

  value_out=paste0("PhenoGraph Clustering Labeled Subtype: ",df_out$diag_pred,
                   " (FeatureN=",df_out$FeatureN,"; NeighborN=",df_out$top_neighborN,")")
  # if(type=="value"){out=value_out}
  # if(tolower(type) %in% c("df","dataframe")){out=df_out}
  return(list(df=df_out,value=value_out))
}
