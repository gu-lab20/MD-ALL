#' get_BALL_mutation
#'
#' @param file_vcf
#'
#' @return
#' @export
#'
#' @examples
get_BALL_mutation=function(file_vcf){
  df_snv=vcf_to_snv(vcf_file = file_vcf,outfmt = 2) %>%
    mutate(mutation=paste0(chr,start,ref,var))

  df_snv1=BALL.mutation %>%
    mutate(mutation=paste0(chrom,start,refallele,varallele)) %>%
    filter(mutation %in% df_snv$mutation) %>%
    left_join(df_snv %>% select(mutation,depth,maf)) %>%
    mutate(mutation=F)

  out_text_BALLmutation=NA
  out_text_SubtypeDefiningMutation=NA

  if(nrow(df_snv1) > 0){
    out_text_BALLmutation=paste0(df_snv1$gene,": ",gsub("p.","",df_snv1$aaPos))

    df_snv1_1=df_snv1 %>% filter(!is.na(subtypesignature))
    if(nrow(df_snv1_1) > 0){
      out_text_SubtypeDefiningMutation=paste0(df_snv1_1$gene,": ",gsub("p.","",df_snv1_1$aaPos))
    }

  }

  return(list(BALL_snv=df_snv1,out_text_BALLmutation=out_text_BALLmutation,out_text_SubtypeDefiningMutation=out_text_SubtypeDefiningMutation))
}
