#' get_BALL_fusion
#'
#' @param file_fusion
#' @param type
#'
#' @return
#' @export
#'
#' @examples
get_BALL_fusion=function(file_fusion,type,cutoff=2){

  if(!file.exists(file_fusion)){fusion_merge=data.frame()}
  else if (length(readLines(file_fusion))<=1){fusion_merge=data.frame()}
  else {
    #read fusion file
    if(tolower(type) %in% c("fustioncatcher","fc")){
      df_fusion=as.data.frame(data.table::fread(input = file_fusion))
      df_fusion=df_fusion[c(1,2,5,6)]
      names(df_fusion)=c("gene1","gene2","feature1","feature2")
      df_fusion=df_fusion %>%
        filter(feature1>0 & feature2>0) %>%
        arrange(desc(feature2),desc(feature1)) %>%
        filter(!(is.na(gene1) | is.na(gene2))) %>%
        group_by(gene1,gene2) %>%
        slice_head(n=1) %>%
        ungroup()
    }

    if(tolower(type) %in% c("cicero","c")){
      df_fusion=as.data.frame(data.table::fread(input = file_fusion))
      df_fusion=df_fusion[c(2,7,13,14)]
      names(df_fusion)=c("gene1","gene2","feature1","feature2")
      df_fusion=df_fusion %>%
        filter(feature1>0 & feature2>0) %>%
        arrange(desc(feature2),desc(feature1)) %>%
        filter(!(is.na(gene1) | is.na(gene2))) %>%
        group_by(gene1,gene2) %>%
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
    # x=BALL.fusion
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

    if(any((fusion_merge$feature1>=cutoff & fusion_merge$feature2 >=cutoff) | fusion_merge$feature1+fusion_merge$feature2 >= 5)){
      # names(fusion_merge)=c('FusionInFile',"feature1",'feature2','FusionFeature',"RelatedSubtype","SubtypeSignificance","ConfidenceToRelatedSubtype")

      fusion_merge=fusion_merge %>% filter((feature1>=cutoff & feature2 >=cutoff) | feature1+feature2>=5)

      names(fusion_merge)=c("FusionInFile","feature1","feature2","PossibleSubtype","SubtypeInfo")
      if(tolower(type) %in% c("fustioncatcher","fc")){names(fusion_merge)[2:3]=c("Spanning_pairs","Spanning_unique_reads")}
      if(tolower(type) %in% c("cicero","c")){names(fusion_merge)[2:3]=c("readsA","readsB")}
      fusion_merge$method=type

    } else {
      fusion_merge=data.frame()
    }

    #get output
  }

  fusion_merge
}

# cat raw_files/fusion/COH005448_D1.cicero |cut -f1,2,3,4,5,7,8,9,10,13,14,27,6,11|grep ETV6
#cat raw_files/fusion/COH002003_D1.fusioncatcher |cut -f1,2,5,6,9,10,15,16|grep KMT2A

