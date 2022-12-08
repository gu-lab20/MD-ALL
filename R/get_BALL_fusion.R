#' get_BALL_fusion
#'
#' @param file_fusion
#' @param type
#'
#' @return
#' @export
#'
#' @examples
get_BALL_fusion=function(file_fusion,type){

  #read fusion file
  if(tolower(type) %in% c("fustioncatcher","fc")){
    df_fusion=as.data.frame(data.table::fread(input = file_fusion))
    df_fusion=df_fusion[c(1,2,5,6)]
    names(df_fusion)=c("gene1","gene2","feature1","feature2")
    df_fusion=df_fusion %>% group_by(gene1,gene2) %>%
      arrange(desc(feature1),desc(feature2)) %>%
      slice_head(n=1) %>%
      ungroup()
  }

  if(tolower(type) %in% c("cicero","c")){
    df_fusion=as.data.frame(data.table::fread(input = file_fusion))
    df_fusion=df_fusion[c(2,7,13,14)]
    names(df_fusion)=c("gene1","gene2","feature1","feature2")
    df_fusion=df_fusion %>%
      filter(!(is.na(gene1) | is.na(gene2))) %>% group_by(gene1,gene2) %>%
      slice_head(n=1) %>%
      ungroup()
  }

  df_fusion$gene1=gsub("@","",df_fusion$gene1);df_fusion$gene2=gsub("@","",df_fusion$gene2);

  #get fusions in file
  df_fusion1=bind_rows(
    df_fusion %>%
      mutate(fusion=paste0(gene1,"::",gene2),
             fusion_ordered=fusion) %>%
      select(fusion,fusion_ordered,feature1,feature2) %>%
      tidyr::separate_rows(fusion_ordered, sep = "::") %>% distinct(),

    df_fusion %>%
      distinct() %>%
      group_by(gene1,gene2) %>%
      mutate(fusion=paste0(gene1,"::",gene2),
             fusion_ordered=paste0(sort(unlist(strsplit(fusion,"::"))),collapse = "::")) %>%
      ungroup() %>%
      select(fusion_ordered,fusion,feature1,feature2)
  )

  #merge with info
  fusion_merge=BALL.fusion %>% left_join(df_fusion1) %>%
    filter(!is.na(fusion)) %>%
    select(fusion,feature1,feature2,fusionFeature,Subtype,note,subtypesignificance)

  names(fusion_merge)=c('FusionInFile',"feature1",'feature2','FusionFeature',"RelatedSubtype","SubtypeSignificance","ConfidenceToRelatedSubtype")

  if(tolower(type) %in% c("fustioncatcher","fc")){names(fusion_merge)[2:3]=c("Spanning_pairs","Spanning_unique_reads")}
  if(tolower(type) %in% c("cicero","c")){names(fusion_merge)[2:3]=c("readsA","readsB")}

  #get output
  fusion_merge

}


# file_fusion="tests/CRLF2NonPhLike.fusioncatcher"
# file_fusion="tests/CRLF2NonPhLike.cicero"
# head(df_fusion)




