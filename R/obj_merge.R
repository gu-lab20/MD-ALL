#' obj_merge
#'
#' @param obj_in , object, object of reference
#' @param file_obj , string, file of object of reference
#' @param df_in , dataframe, count dataframe for merge
#'
#' @return
#' @export obj_merge
#'
#' @examples
#' obj_x=obj_merge(file_obj = "../tsne/234_featureCounts/obj_234_featureCountsMini.rds",df_in = df_count)
obj_merge=function(obj_in=NULL,file_obj=NULL,df_in,assay_name_in="vst"){
  if(!is.null(file_obj)){
    obj_in=readRDS(file_obj)
  }
  #get obj count data
  count_ref=assays(obj_in$SE)[[assay_name_in]]
  count_ref$feature=row.names(count_ref)

  #get merged count data

  if(!any(grepl("feature",names(df_in)))){df_in$feature=row.names(df_in)}

  features_overlap=df_in$feature[df_in$feature %in% row.names(count_ref)]
  count_merge=df_in %>% left_join(count_ref) %>%
    filter(feature %in% features_overlap)
  row.names(count_merge)=count_merge$feature
  count_merge$feature=NULL

  #get info data
  df_info_ref=as.data.frame(colData(obj_in$SE))
  df_info_in=data.frame(
    id=names(df_in)[2:ncol(df_in)],
    diag="TestSample",
    library="Unknown"
  )

  row.names(df_info_in)=names(df_in)[2:ncol(df_in)]

  df_info_all=bind_rows(df_info_in,df_info_ref)[c("diag")]

  #get output obj
  if(all(names(count_merge)==row.names(df_info_all))){
    matrix_list=list(counts=count_merge)
    names(matrix_list)=assay_name_in
    SE=SummarizedExperiment(assays=matrix_list, colData = df_info_all)
    obj_out=obj_in
    obj_out$SE=SE
  }

  if(!all(names(count_merge)==row.names(df_info_all))){
    stop("Names not match")
  }
  obj_out
}
