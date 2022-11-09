#' run_tsne
#'
#' @param obj_in
#' @param perplexity_one
#' @param dims
#' @param variable_n
#' @param out_label
#' @param feature_panel
#'
#' @return
#' @export
#'
#' @examples
run_tsne=function(obj_in,perplexity_one=30,dims=2,variable_n=800,out_label="tsne",feature_panel="variable_genes"){
  variable_genes_in=obj_in[[feature_panel]][1:min(variable_n,length((obj_in[[feature_panel]])))]

  indata=assays(obj_in$SE[variable_genes_in,])[["vst"]]
  indata1=t(indata)
  cat("Run tSNE: Used Feature N=", dim(indata1)[2], "; Used Sample N=", dim(indata1)[1],"; Perplexity=", perplexity_one,  "\n")

  set.seed(10) # Sets seed for reproducibility
  tsne_out = Rtsne(indata1, dims = dims, perplexity = perplexity_one,
                   theta = 0.5,
                   max_iter = 5000,
                   check_duplicates = F,
                   partial_pca=T,
                   num_threads = 4)

  #summary output
  df_id=data.frame(COH_sample=colnames(indata))

  df_tsne=as.data.frame(tsne_out$Y)
  names(df_tsne)=c("tSNE_1","tSNE_2","tSNE_3")[1:ncol(df_tsne)]

  df_tsne$perplexityN=perplexity_one
  df_tsne$FeatureN=dim(indata1)[2]

  # names(df_tsne)=paste0(names(df_tsne),"_FeatureN",dim(indata1)[2],"_PerplexityN",perplexity_one)

  row.names(df_tsne)=colnames(indata)

  obj_in[[out_label]]=df_tsne
  obj_in
}
