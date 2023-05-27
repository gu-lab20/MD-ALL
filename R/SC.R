#' read_sc_file
#'
#' @param dfPath
#'
#' @return
#' @export
#'
#' @examples
read_sc_file=function(dfPath){
  count_sc=read.table(dfPath,stringsAsFactors = F,sep = "\t",header = T)
  row.names(count_sc)=unlist(count_sc[,1])
  count_sc=count_sc[,-1]
  count_sc
}

#' scCount2umap
#'
#' @param count_matrix
#'
#' @return
#' @export
#'
#' @examples
scCount2umap=function(count_matrix){
  obj_=CreateSeuratObject(count_matrix)
  obj_=NormalizeData(obj_,verbose = F)
  obj_=FindVariableFeatures(obj_,nfeatures = 2000,verbose = F)
  obj_=ScaleData(obj_,verbose = F)
  obj_=RunPCA(obj_,npcs = 50,verbose = F)
  obj_=RunUMAP(obj_,dims = 1:50,verbose = F)
  obj_
}

#' run_singleR
#'
#' @param inObj
#' @param ref_singleR
#' @param var_group
#'
#' @return
#' @export
#'
#' @examples
run_singleR=function(inObj,ref_singleR,var_group){
  for_singleR=GetAssayData(inObj)
  # dim(for_singleR)

  # print(dim(for_singleR))

  value_group=unlist(ref_singleR@colData[var_group])

  #cell level annotation
  # cat("Cell level annotation...\n")
  out_singleR = SingleR::SingleR(test = for_singleR, ref = ref_singleR,labels = value_group)
  #get output
  df_singleR_out=data.frame(
    barcode=row.names(inObj@meta.data),
    labels=out_singleR$labels,
    pruned.labels=out_singleR$pruned.labels,
    stringsAsFactors = F)

  return(list(obj_singleR=out_singleR,df_singleR_out=df_singleR_out))
}

#' get_singleR_label
#'
#' @param singleR_out
#' @param obj_in
#' @param reduction_method
#' @param levelsIn
#' @param label_levelsIn
#'
#' @return
#' @export
#'
#' @examples
get_singleR_label=function(singleR_out,
                          obj_in=NULL,reduction_method="umap",
                          levelsIn=c("Pro-B1","Pro-B2","Pre-B1","Pre-B2"),
                          label_levelsIn="Pro/Pre-B"){

  #get delta
  df_delta=data.frame(
    barcode=singleR_out$obj_singleR@rownames,
    deltaFromMedian=SingleR::getDeltaFromMedian(singleR_out$obj_singleR),
    stringsAsFactors = F
  )

  if(!is.null(obj_in)){
    #get embeding
    df_embed = as.data.frame(Embeddings(obj_in[[reduction_method]]))
    df_embed$barcode=row.names(df_embed)

    #get label
    df_label=data.frame(barcode=row.names(obj_in@meta.data),stringsAsFactors = F) %>%
      left_join(df_embed) %>%
      left_join(singleR_out$df_singleR_out) %>% mutate(pruned.labels=NULL) %>%
      left_join(df_delta)
  }

  if(is.null(obj_in)){
    df_label=data.frame(barcode=singleR_out$obj_singleR@rownames,stringsAsFactors = F) %>%
      left_join(singleR_out$df_singleR_out) %>% mutate(pruned.labels=NULL) %>%
      left_join(df_delta)
  }

  #get label percentage
  df_label_per=df_label %>% select(labels) %>% mutate(n=n()) %>% group_by(labels) %>%
    mutate(n_celltype=n(),
           percentage=round(n_celltype/n,4),
           label=paste0(labels,"\n(",sprintf("%.1f",percentage*100),"%)"),
           label_high=ifelse(labels %in% levelsIn,label_levelsIn,"Others")
    ) %>%
    select(labels,percentage,label,label_high) %>%
    distinct() %>%
    arrange(desc(percentage)) %>% ungroup()

  list(df_label=df_label,df_label_per=df_label_per)
}


#' get_singleR_topCells
#'
#' @param singleR_out
#' @param topPercentileOfDeltaToKeep
#' @param remove_valueLabelNonConsistent
#'
#' @return
#' @export
#'
#' @examples
get_singleR_topCells=function(singleR_out,topPercentileOfDeltaToKeep=0.9,remove_valueLabelNonConsistent=T){
  #get label using max score
  df_score=as.data.frame(singleR_out$obj_singleR$scores)
  df_score$barcode=singleR_out$obj_singleR@rownames
  df_score1=df_score %>% reshape2::melt(id.vars=c("barcode")) %>%
    group_by(barcode) %>% filter(value==max(value))
  names(df_score1)[2:3]=c("labelsFromMaxScore","maxScore")

  #get top cells using delta
  cutoff_delta_precentile=1-topPercentileOfDeltaToKeep
  df_delta=get_singleR_label(singleR_out)$df_label %>%
    group_by(labels) %>%
    mutate(
      delta_cutoff=quantile(deltaFromMedian,cutoff_delta_precentile),
      label_top=ifelse(deltaFromMedian>delta_cutoff,as.character(labels),NA)
    ) %>% ungroup()

  #get output
  df_out=df_score1 %>% left_join(df_delta %>% select(barcode,label_top))

  if(remove_valueLabelNonConsistent){
    df_out=df_out %>%
      mutate(label_top=ifelse(label_top==labelsFromMaxScore,label_top,NA))
  }

  df_out
}

#' get_singleR_heatmapScore
#'
#' @param singleR_out
#'
#' @return
#' @export
#'
#' @examples
get_singleR_heatmapScore=function(singleR_out){
  df_label=get_singleR_label(singleR_out = singleR_out)$df_label %>%
    left_join(get_singleR_topCells(singleR_out,topPercentileOfDeltaToKeep=1,remove_valueLabelNonConsistent=F)) %>%
    filter(labels==labelsFromMaxScore) %>% mutate(n=n()) %>%
    group_by(labels) %>% mutate(n_g=n(),diag_per=round(n_g/n,4)) %>%
    arrange(desc(n_g),desc(deltaFromMedian)) %>% ungroup() %>% mutate(obs=1:n())

  df_label_unique=df_label %>% select(labels,n_g,diag_per) %>% arrange(desc(n_g)) %>% distinct() %>%
    mutate(label_=paste0(letters[1:n()],".",labels),
           per=paste0(sprintf("%.1f",diag_per*100),"%"),
           label_per=paste0(labels," (",per,")"))

  names(df_label_unique)[1]="variable"

  matrix_score=bind_cols(
    data.frame(
      barcode=singleR_out$obj_singleR@rownames,
      stringsAsFactors = F
    ),
    as.data.frame(singleR_out$obj_singleR$scores)
  )
  row.names(matrix_score)=matrix_score$barcode

  matrix_score_=df_label %>% left_join(matrix_score) %>% arrange(obs)

  scale_0to1=function(x){
    x=as.numeric(x)
    (x-min(x))/(max(x)-min(x))}

  df_score_=matrix_score_[c("barcode","obs","labels",df_label_unique$variable)] %>%
    reshape2::melt(id.vars=c("barcode","obs","labels")) %>% left_join(df_label_unique) %>%
    group_by(barcode) %>%
    mutate(value1=scale_0to1(value)) %>% ungroup()

  df_score_
}

#' get_SC_subtypes
#'
#' @param count_matrix
#' @param levelsIn
#' @param SE_celltype
#' @param SE_BALL
#'
#' @return
#' @export
#'
#' @examples
get_SC_subtypes=function(count_matrix,levelsIn=c("Pro-B1","Pro-B2","Pre-B1","Pre-B2"),SE_celltype=SE_celltype,SE_BALL=SE_BALL){
  #get cols   ----------------------------------------------------------------------------------------------------------------------------------------------------
  #celltype
  df_col_celltype=bind_rows(
    as.data.frame(colData(SE_celltype)) %>% select(celltype,col) %>% distinct(),
    data.frame(celltype=c("Others","NotForBALLsubtyping"),col=c("powderblue","grey90"))
  )

  cols_celltype_raw=df_col_celltype$col
  names(cols_celltype_raw)=df_col_celltype$celltype

  #BALL subtype
  df_col_BALL=bind_rows(
    as.data.frame(colData(SE_BALL)) %>% select(diag,col) %>% distinct(),
    data.frame(diag=c("Others"),col=c("powderblue"))
  )

  cols_BALL_raw=df_col_BALL$col
  names(cols_BALL_raw)=df_col_BALL$diag


  #prepare single cell object  ----------------------------------------------------------------------------------------------------------------------------------------------------
  message('\nRunning single cell normalization, PCA and UMAP ...')
  obj_sc=scCount2umap(count_matrix)


  #celltype annotation ----------------------------------------------------------------------------------------------------------------------------------------------------
  message("\nRunning cell type annotation ...")
  singleR_celltype=run_singleR(inObj = obj_sc,ref_singleR = SE_celltype,var_group = "celltype")
  df_celltype=get_singleR_label(obj_in = obj_sc,singleR_out = singleR_celltype,levelsIn = levelsIn)$df_label

  #umap for all celltypes ----------------------------------------------------------------------------------------------------------------------------------------------------
  message("\nGet UMAP for  cell types ...")
  obj_sc$celltype=df_celltype$labels

  umap_celltype=DimPlot(object = obj_sc,group.by = "celltype",order = T,pt.size = 0.8,
                        cols = cols_celltype_raw[names(cols_celltype_raw) %in% obj_sc$celltype])+
    labs(title = "Cell types")+ guides(fill = guide_legend(title.position = "top",title.hjust = 0.5,title = "Cell types")) +
    theme(plot.margin = margin(t = 1,r = 1,b = 1,l = 1, "cm"))

  #barplot for cell types ----------------------------------------------------------------------------------------------------------------------------------------------------
  message("\nGet barplot for cell types ...")

  df_celltype_per=get_singleR_label(obj_in = obj_sc,singleR_out = singleR_celltype)$df_label_per

  barplot_celltype=
    gg_barplot(df_plot  =df_celltype_per ,x = "label_high",y = "percentage",group.by = "labels",bar_width = 0.5,
               cols = cols_celltype_raw[names(cols_celltype_raw) %in% df_celltype_per$labels],
               x_lab = "Cell type group",y_lab="Percentages",title = "Bar plot for cell types") +
    theme(legend.position = "right") +
    guides(fill = guide_legend(title.position = "top",title.hjust = 0.5,title = "Cell types"))+
    theme(plot.margin = margin(t = 2,r = 4,b = 2,l = 4, "cm"),
          plot.title = element_text(face = "bold",hjust = 0.5,size = 17))

  #select top cells ----------------------------------------------------------------------------------------------------------------------------------------------------
  message("\nGet top cells ...")
  df_celltype_top=get_singleR_topCells(singleR_out = singleR_celltype,topPercentileOfDeltaToKeep = 0.9)
  df_celltype_top1=df_celltype_top %>% filter(label_top %in% levelsIn)

  if(nrow(df_celltype_top1)<=1){
    message(paste0("Cells numbers for BALL subtyping=",nrow(df_celltype_top1),". Not Enough!"))
    list_out=list()

    p_out=plot_grid(plot_grid(umap_celltype,NULL,NULL,ncol  = 1,align = "hv"),plot_grid(barplot_celltype,NULL,ncol = 1),ncol=2)
  }

  if(nrow(df_celltype_top1)>1){

    #umap for celltypes of top cells ----------------------------------------------------------------------------------------------------------------------------------------------------
    message("\nGet UMAP of cell types for top cells ...")

    obj_sc$celltype_inBALLsubtyping=ifelse(colnames(x=obj_sc) %in% df_celltype_top1$barcode,obj_sc$celltype,"NotForBALLsubtyping")

    umap_celltypeForBALL=DimPlot(object = obj_sc,group.by = "celltype_inBALLsubtyping",order = T,pt.size = 0.8,
                                 cols = cols_celltype_raw[names(cols_celltype_raw) %in% obj_sc$celltype_inBALLsubtyping]) +
      labs(title = "Celltypes for BALL subtyping") +
      guides(fill = guide_legend(title.position = "top",title.hjust = 0.5,title = "Cell types"))+
      theme(plot.margin = margin(t = 1,r = 1,b = 1,l = 1, "cm"))

    #BALL annotation ----------------------------------------------------------------------------------------------------------------------------------------------------
    message("\nRunning BALL subtyping on single cell level ...")
    obj_subset=subset(obj_sc,cells=df_celltype_top1$barcode)
    row.names(SE_BALL)=rowData(SE_BALL)$gene_name

    singleR_BALL=run_singleR(inObj = obj_subset,ref_singleR = SE_BALL,var_group = "diag")

    df_BALL=get_singleR_topCells(singleR_out = singleR_BALL,topPercentileOfDeltaToKeep = 1)

    #umap of BALL subtypes  ----------------------------------------------------------------------------------------------------------------------------------------------------
    message("\nGet UMAP for BALL subtypes ...")
    obj_subset1=subset(obj_subset,cells=(df_BALL %>% filter(!is.na(label_top)))$barcode)

    df_BALL1=data.frame(barcode=row.names(obj_subset1@meta.data),stringsAsFactors = F) %>% left_join(df_BALL)

    obj_subset1$diag=df_BALL1$label_top

    umap_BALL=DimPlot(object = obj_subset1,group.by = "diag",order = T,pt.size = 0.8,
                      cols = cols_BALL_raw[names(cols_BALL_raw) %in% obj_subset1$diag]) +
      labs(title = "BALL subtyping")+
      xlim(min(df_celltype$UMAP_1), max(df_celltype$UMAP_1))+
      ylim(min(df_celltype$UMAP_2), max(df_celltype$UMAP_2)) +
      guides(fill = guide_legend(title.position = "top",title.hjust = 0.5,title = "BALL subtypes"))+
      theme(plot.margin = margin(t = 1,r = 1,b = 1,l = 1, "cm"))

    #heatmap of BALL subtypes  ----------------------------------------------------------------------------------------------------------------------------------------------------
    message("\nGet heatmap of SingleR scores for BALL subtypes ...")

    df_heatmap_score=get_singleR_heatmapScore(singleR_BALL) %>%
      mutate(label_new=paste0(variable,"\nN=",n_g," (",per,")"))

    df_label_per1=df_heatmap_score %>% select(obs,labels) %>% mutate(obs_max=max(obs)) %>% distinct() %>% mutate(group=1) %>% mutate(n=n()) %>% group_by(labels) %>%
      mutate(n_g=n(),
             `Cell number`=n_g,
             per=round(n_g/n,4),
             Percentage=paste0(sprintf("%.1f",per*100),"%"),
             `Subtypes, N (%)`=paste0(labels,", ",n_g," (",Percentage,")")
      ) %>%
      select(labels,obs_max,n_g,`Cell number`,Percentage,`Subtypes, N (%)`) %>% distinct()

    p_heatmap=gg_heatmap(df_in = df_heatmap_score,x = "obs",y = "label_",var_col = "value1",reverse_y = T,y_tick_label_var = "label_new",
                         title = "Scores (column scaled to 0-1)") +
      theme(legend.box.spacing = unit(0, "pt"),axis.text.y.right = element_text(margin = margin(l=-0.5,unit = "cm"))) +
      labs(title = "Heatmap for singleR scores")+
      theme(plot.margin = margin(t = 1,r = 1,b = 1,l = 1, "cm"),
            plot.title = element_text(face = "bold",hjust = 0.5,size = 17))

    p_out=plot_grid(plot_grid(umap_celltype,umap_celltypeForBALL,umap_BALL,ncol  = 1,align = "hv"),plot_grid(barplot_celltype,p_heatmap,ncol = 1),ncol=2)

    #output report  ----------------------------------------------------------------------------------------------------------------------------------------------------
    message('\nGet output report ...')

  }
  message("\nDone single cell BALL subtyping!")
  p_out
}
























































































