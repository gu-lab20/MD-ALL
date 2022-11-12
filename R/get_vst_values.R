#' get_vst_values
#'
#' @param obj_in , object, object of reference
#' @param file_obj , string, file of object of reference
#' @param df_count , dataframe, count dataframe for merge
#'
#' @return
#' @export get_vst_values
#'
#' @examples
get_vst_values=function(obj_in=NULL,file_obj=NULL,df_count){

  obj_x=obj_merge(obj_in=obj_in,file_obj = file_obj,df_in = df_count,assay_name_in = "counts")

  obj_x=run_vst(obj_in = obj_x)

  df_vst=as.data.frame(assays(obj_x$SE)[["vst"]])
  df_vst$feature=row.names(df_vst)

  df_vst_out=df_vst[c("feature",names(df_count)[2:ncol(df_count)])] %>% arrange(feature)

  df_vst_out
}

