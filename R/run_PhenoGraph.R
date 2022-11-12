#' run_PhenoGraph
#'
#' @param obj_in
#' @param feature_panel
#' @param variable_n
#' @param neighbor_k
#' @param ratio_cutoff
#' @param out_label
#'
#' @return
#' @export
#'
#' @examples
run_PhenoGraph=function(obj_in,
                        feature_panel="variable_genes",
                        variable_n = variable_n,
                        neighbor_k=neighbor_k,
                        ratio_cutoff=0.5,
                        out_label="PhenographPred"){

  variable_genes_in=obj_in[[feature_panel]][1:min(variable_n,length((obj_in[[feature_panel]])))]
  indata=assays(obj_in$SE[variable_genes_in,])[["vst"]]

  indata1=t(indata)
  cat("Run Phenograph: Used Feature N=", dim(indata1)[2], "; Used Sample N=", dim(indata1)[1],"; Neighbor_k=", neighbor_k,  "\n")

  df_diag_match=data.frame(COH_sample=row.names(indata1),diag=obj_in$SE$diag,stringsAsFactors = F) %>% mutate(obs=1:n())

  set.seed(10)
  PG_out=Rphenograph(indata1, neighbor_k)

  df_cluster=data.frame(
    obs=as.numeric(names(membership(PG_out[[2]]))),
    cluster=as.vector(membership(PG_out[[2]])),
    stringsAsFactors = F
  )

  df_cluster_detail=df_diag_match %>% left_join(df_cluster) %>%
    group_by(diag, cluster) %>% mutate(N_diagCluster=n()) %>%
    group_by(cluster) %>%  mutate(N_cluster=n()) %>%
    mutate(ratio = round(N_diagCluster/N_cluster, digits = 2)) %>%
    arrange(cluster, desc(N_cluster))

  df_cluster_diag=df_cluster_detail %>% select(diag,cluster,ratio,N_cluster,N_diagCluster,) %>% distinct() %>%
    filter(ratio>ratio_cutoff) %>% group_by(diag) %>% arrange(desc(N_diagCluster)) %>%
    mutate(diag_pred=ifelse(!is.na(cluster),diag,"NoClusterAsigned"),
           diagCluster_index=1:n(),
           diagCluster_index_max=n(),
           diag_pred_granular=ifelse(diagCluster_index_max==1,diag_pred,paste0(diag_pred,"_",diagCluster_index)),
           diag_pred_freq=paste0(diag_pred,"(",N_diagCluster,";",ratio*100,"%)"),
           diag_pred_granular_freq=paste0(diag_pred_granular,"(",N_diagCluster,";",ratio*100,"%)"),
    ) %>% select(diag,cluster,diag_pred,diag_pred_freq,diag_pred_granular,diag_pred_granular_freq)

  df_cluster_diag$diag=NULL

  df_diag_pred=df_cluster_detail %>%
    left_join(df_cluster_diag %>% mutate(ratio=NULL)) %>%
    mutate(obs=NULL,FeatureN=variable_n,top_neighborN=neighbor_k)

  df_diag_pred$diag_pred=ifelse(is.na(df_diag_pred$cluster),"NoClusterAsigned",df_diag_pred$diag_pred)
  df_diag_pred$diag_pred_freq=ifelse(is.na(df_diag_pred$cluster),"NoClusterAsigned",df_diag_pred$diag_pred_freq)
  df_diag_pred$diag_pred_granular=ifelse(is.na(df_diag_pred$cluster),"NoClusterAsigned",df_diag_pred$diag_pred_granular)
  df_diag_pred$diag_pred_granular_freq=ifelse(is.na(df_diag_pred$cluster),"NoClusterAsigned",df_diag_pred$diag_pred_granular_freq)

  row.names(df_diag_pred)=df_diag_pred$COH_sample

  obj_in[[out_label]]=df_diag_pred
  obj_in
}
