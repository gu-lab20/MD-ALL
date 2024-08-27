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
    mutate(mutationLabel=paste0(gene,":",gsub("p.","",aaPos),"(",sprintf("%.2f",maf),")"))

  out_text_BALLmutation=NA
  out_text_SubtypeDefiningMutation=NA
  subtype_mutation=NA

  if(nrow(df_snv1) > 0){
    out_text_BALLmutation=paste0(df_snv1$mutationLabel,collapse = "\n")

    df_snv1_1=df_snv1 %>% filter(!(is.na(subtypesignature) | subtypesignature==""))
    if(nrow(df_snv1_1) > 0){
      out_text_SubtypeDefiningMutation=paste0(df_snv1_1$mutationLabel,collapse = "\n")
      subtype_mutation=df_snv1_1$subtypesignature
    }
  }
  return(list(BALL_snv=df_snv1,out_text_BALLmutation=out_text_BALLmutation,
              out_text_SubtypeDefiningMutation=out_text_SubtypeDefiningMutation,
              subtype_mutation=subtype_mutation))
}














