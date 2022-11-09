#' run_vst
#'
#' @param obj_in , object, object of reference
#' @param assay_name_in , string, assay name used for running vst
#' @param assay_name_out , string, assay names used for storing the vst values
#' @param var_design , string, variable in coldata used in vst running formula design. will not affecting the vst values
#'
#' @return
#' @export run_vst
#'
#' @examples
#' obj_x=run_vst(obj_in = obj_x)
run_vst=function(obj_in,assay_name_in="counts",assay_name_out="vst",var_design="diag"){

  matrix_=as.matrix(assays(obj_in$SE)[[assay_name_in]])
  df_colData=as.data.frame(as.matrix(colData(obj_in$SE)))


  formula_in=formula(paste0("~",paste0(var_design,collapse = "+")))

  obj_deseq2 = DESeq2::DESeqDataSetFromMatrix(countData = matrix_,
                                              colData = df_colData,
                                              design = formula_in)
  cat("Running vst\n")
  obj_deseq2=DESeq2::varianceStabilizingTransformation(obj_deseq2,blind = T)
  df_vst=round(as.matrix(assay(obj_deseq2)),6)

  assays(obj_in$SE)[[assay_name_out]]=df_vst
  obj_in
}
