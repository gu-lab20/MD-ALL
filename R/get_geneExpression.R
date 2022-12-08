#' get_geneExpression
#'
#' @param df_vst
#' @param genes
#'
#' @return
#' @export
#'
#' @examples
get_geneExpression=function(df_vst,genes){
  names(df_vst)=c("feature","TestSample")
  geneId=info_gtf_hg38$gene_id[info_gtf_hg38$gene_name %in% genes]
  out=df_vst[df_vst$feature %in% geneId,] %>%
    left_join(info_gtf_hg38 %>%
                transmute(feature=gene_id,Gene=gene_name)) %>%
    mutate(Expression = TestSample) %>%
    select(Gene,Expression)
  out
}
