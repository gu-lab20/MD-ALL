#' run_umap
#'
#' @param obj_in
#' @param n_neighbors
#' @param dims
#' @param variable_n
#' @param out_label
#' @param feature_panel
#'
#' @return
#' @export
#'
#' @examples
run_umap=function(obj_in,n_neighbors=30,dims=2,variable_n=800,out_label="umap",feature_panel="variable_genes",min_dist=0.3){
  #get features
  variable_genes_in=obj_in[[feature_panel]][1:min(variable_n,length((obj_in[[feature_panel]])))]

  #get running data
  indata=assays(obj_in$SE[variable_genes_in,])[["vst"]]
  indata1=t(indata)
  cat("Running uMAP: Feature N=", dim(indata1)[2], "; Sample N=", dim(indata1)[1],"; n_neighbors=", n_neighbors,  "\n")

  #running umap
  set.seed(10)
  fit_umap=umap(indata1,n_neighbors=n_neighbors,random_state=123,n_components=dims,min_dist=min_dist)

  #get output
  df_umap=as.data.frame(fit_umap$layout)
  names(df_umap)=c("uMAP_1","uMAP_2","uMAP_3")[1:ncol(df_umap)]
  df_umap$n_neighbors=n_neighbors
  df_umap$FeatureN=dim(indata1)[2]
  df_umap$COH_sample=row.names(df_umap)

  #summary output
  obj_in[[out_label]]=df_umap
  obj_in[[paste0(out_label,"_fit")]]=fit_umap

  obj_in
}
