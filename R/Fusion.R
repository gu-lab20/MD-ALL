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

  if(!file.exists(file_fusion)){fusion_merge=data.frame()}
  else if (length(readLines(file_fusion))<=1){fusion_merge=data.frame()}
  else {
    #read fusion file
    if(tolower(type) %in% c("fustioncatcher","fc")){
      df_fusion=as.data.frame(data.table::fread(input = file_fusion))
      df_fusion=df_fusion[c(1,2,5,6)]
      names(df_fusion)=c("gene1","gene2","feature1","feature2")
      df_fusion=df_fusion %>% group_by(gene1,gene2) %>%
        arrange(desc(feature2),desc(feature1)) %>%
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
      select(fusion,feature1,feature2,fusionFeature,Subtype,note,subtypesignificance) %>%
      mutate(subtype_info=paste0(Subtype,"@",note,",",subtypesignificance))

    if("CRLF2 fusion" %in% fusion_merge$note){
      fusion_merge=fusion_merge[!fusion_merge$note %in% c("adjacent to CRLF2"),]
    }

    fusion_merge=fusion_merge %>% group_by(fusion,feature1,feature2) %>%
      mutate(subtype_info=paste0(sort(subtype_info),collapse = ";")) %>%
      ungroup() %>%
      select(fusion,feature1,feature2,Subtype,subtype_info) %>%
      distinct() %>%
      arrange(desc(feature2),desc(feature1))


    # names(fusion_merge)=c('FusionInFile',"feature1",'feature2','FusionFeature',"RelatedSubtype","SubtypeSignificance","ConfidenceToRelatedSubtype")

    names(fusion_merge)=c("FusionInFile","feature1","feature2","PossibleSubtype","SubtypeInfo")

    if(tolower(type) %in% c("fustioncatcher","fc")){names(fusion_merge)[2:3]=c("Spanning_pairs","Spanning_unique_reads")}
    if(tolower(type) %in% c("cicero","c")){names(fusion_merge)[2:3]=c("readsA","readsB")}

    fusion_merge$method=type

    #get output
  }

  fusion_merge
}

