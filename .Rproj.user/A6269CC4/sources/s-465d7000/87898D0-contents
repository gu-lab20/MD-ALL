#' server_singleSample
#'
#' @param input
#' @param output
#' @param session
#'
#' @return
#' @export
#'
#' @examples
server_singleSample=function(input, output, session){
  showNotification(id = "showme_notif","Please upload the input file(s) to begin.",duration = 15,type = "message")

  #Count input  --------------------------------------------------------------------------------------
  observe({
    shinyjs::toggleState("run_all", condition = !(is.null(input$fileCount) & is.null(input$fileVcf)))
  })

  # output$fileCountName=reactive({req(input$fileCount);paste0("Uploaded file name: ",input$fileCount$name)})

  #GEP -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  output$body_GEP=renderUI({req(input$fileCount);body_GEP()})
  output$para_umap=renderUI({req(input$fileCount);para_umap()})
  output$para_phenograph=renderUI({req(input$fileCount);para_phenograph()})

  df_count=reactive({read_input(input$fileCount$datapath)})

  output$countTop_table=renderTable({req(df_count());head(df_count(),4)})
  output$countTail_table=renderTable({req(df_count());tail(df_count(),4)})

  #run vst and check first 4 lines
  df_vst=eventReactive({input$run_all},{
    shinyjs::disable("run_all")
    running_norm=showNotification("Running Normalization ...",duration = NULL,closeButton = F,type = "warning")
    on.exit(removeNotification(running_norm),add=T)
    get_vst_values(obj_in = obj_234_HTSeq,df_count = df_count())
  }
  )

  output$df_vstTop=renderTable({req(df_vst());head(df_vst(),4)})
  output$df_vstTail=renderTable({req(df_vst());tail(df_vst(),4)})

  output$downloadvst=renderUI({req(df_vst());downloadButton("df_vst.tsv", "Download normalized data")})
  output$df_vst.tsv = downloadHandler(
    filename = function(){"df_vst.tsv"},
    content = function(file){write.table(df_vst(),file,col.names = T,row.names = F,quote = F,na = "",sep="\t")}
  )

  GeneExpressionTabSum=reactive({get_geneExpression(df_vst = df_vst(),genes = c("CDX2","CRLF2","DUX4","HLF","NUTM1","MEGF10"))})

  #Imputation  ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  df_vst_i=reactive({
    shinyjs::disable("run_all")
    req(df_vst())
    running_phenograph=showNotification("Running Imputation ...",duration = NULL,closeButton = F,type = "warning")
    on.exit(removeNotification(running_phenograph),add=T)

    f_imputation(obj_ref = obj_234_HTSeq,df_in = df_vst())
  })

  #Get obj for umap and GEP   ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  obj_=eventReactive(input$run_all,{
    req(df_vst_i())
    obj_merge(obj_in = obj_1821,df_in = df_vst_i(),assay_name_in = "vst")
  })

  #GEP prediction --------------------------------------------------------------------------------------------------------------------------------------------------------------------
  phenographGeneN=reactive(as.numeric(unlist(strsplit(input$phenographGeneN,split = ","))))
  output$phenographGeneN=renderText(phenographGeneN())

  df_out_phenograph=eventReactive(input$run_all,{
    shinyjs::disable("run_all")
    req(df_vst_i())
    running_phenograph=showNotification("Running PhenoGraph Prediction ...",duration = NULL,closeButton = F,type = "warning")
    on.exit(removeNotification(running_phenograph),add=T)

    get_PhenoGraphPreds(obj_in = obj_(),feature_panel = "keyFeatures",SampleLevel = "TestSample",
                        neighbor_k = 10,
                        variable_n_list=phenographGeneN()
    )
  })

  df_out_svm=eventReactive(input$run_all,{
    running_svm=showNotification("Running SVM Prediction ...",duration = NULL,closeButton = F,type = "warning")
    on.exit(removeNotification(running_svm),add=T)
    get_SVMPreds(models_svm,df_in = df_vst_i())
  })

  df_pred=reactive(
    bind_rows(df_out_phenograph(),df_out_svm()) %>% mutate(N=sprintf("%04d",featureN))
  )

  p_predHeatmap=reactive(
    gg_tilePlot(df_in = df_pred(),x = "N",y = "method",var_col = "pred",x_tick_label_var = "featureN",title = "Prediction Heatmap")
  )

  output$predHeatmap=renderPlot(p_predHeatmap(),height = 200, width = 500)

  output$predHeatmap.pdf = downloadHandler(
    filename = function(){"predHeatmap.pdf"},
    content = function(file){ggsave(file,plot=p_predHeatmap(),width = 6,height = 2)}
  )

  observeEvent(df_pred(),enable("run_all"))


  #run umap --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  obj_umap=eventReactive(input$run_all,{
    shinyjs::disable("run_all")
    run_umap(obj_in = obj_(),out_label = "umap",
             n_neighbors = input$umapNeighbor,
             variable_n = input$umapGeneN,
             feature_panel = "keyFeatures")
  })

  p_umap=eventReactive(input$run_all,{
    running_umap=showNotification("Running UMAP ...",duration = NULL,closeButton = F,type = "warning")
    on.exit(removeNotification(running_umap),add=T)

    draw_DimPlot(obj_umap(),group.by = "diag_raw",reduction = "umap",highlightLevel = "TestSample")+
      theme(legend.position = "bottom")
  })


  output$umap=renderPlot(p_umap(),height = 500, width = 500)

  output$umap.pdf = downloadHandler(
    filename = function(){"umap.pdf"},
    content = function(file){ggsave(file,plot=p_umap(),width = 11,height = 9)}
  )
  observeEvent(p_umap(),enable("run_all"))


  #get box Plot ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  obj_boxplot=reactive({req(df_vst());
    obj_merge(obj_in = obj_1821,df_in = df_vst())
  })

  p_boxPlot=reactive({req(obj_boxplot());
    shinyjs::disable("run_all")
    draw_BoxPlot(obj_in = obj_boxplot(),group.by = "diag_raw1",
                 plot_title = "Box Plot of Gene Expression",
                 features = input$BoxPlotGeneName,highlightLevel = "TestSample")
  })

  output$BoxPlotGene=renderPlot({req(p_boxPlot())
    p_boxPlot()
  },
  height = 450, width = 700)

  observeEvent(obj_boxplot(),enable("run_all"))
  observeEvent(obj_boxplot(),
               showNotification("Done GEP Prediction.",duration = 3,type = "warning")
  )


  #RNAseqCNV --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  output$para_RNAseqCNV=renderUI({req(input$fileVcf);para_RNAseqCNV()})
  output$body_CNV=renderUI({req(input$fileVcf);body_CNV()})


  output$CNVparameters=reactive(
    paste0("Minimum Read Count of Gene: ",input$CNVminReadCnt,"; Mininum Depth of SNV: ",
           input$CNVminDepth,"; Min MAF:",input$CNVmafMin, "; Max MAF: ",input$CNVmafMax)
  )

  RNAseqCNV_out=eventReactive(input$run_all,{
    req(input$fileVcf)
    req(obj_boxplot())
    shinyjs::disable("run_all")

    running_cnv=showNotification("Running RNAseqCNV ...",duration = NULL,closeButton = F,type="warning")
    on.exit(removeNotification(running_cnv),add=T)

    run_RNAseqCNV(df_count = df_count(),
                  snv_file = input$fileVcf$datapath,
                  genome_version = "hg38",
                  minReadCnt = input$CNVminReadCnt,
                  minDepth = input$CNVminDepth,
                  mafRange = c(input$CNVmafMin,input$CNVmafMax))
  })

  chrom_n=reactive({
    req(!is.null(RNAseqCNV_out()))
    RNAseqCNV_out()$df_cnv_out$chrom_n})

  CNV_label=reactive({
    req(!is.null(RNAseqCNV_out()))
    paste0(RNAseqCNV_out()$df_cnv_out$gender,";\n",RNAseqCNV_out()$df_cnv_out$chrom_n,";",RNAseqCNV_out()$df_cnv_out$alterations)
  })

  output$cnvtext=renderTable(RNAseqCNV_out()$df_cnv_out %>% mutate(status=NULL,chrom_n=as.character(chrom_n)))

  p_cnv=eventReactive(input$run_all,{
    shinyjs::disable("run_all")

    running_cnv=showNotification("Running RNAseqCNV ...",duration = NULL,closeButton = F,type="warning")
    on.exit(removeNotification(running_cnv),add=T)

    get_RNAseqCNV_plot(RNAseqCNV_out=RNAseqCNV_out())
  })

  output$cnvplot=renderPlot(p_cnv())

  output$RNAseqCNV.pdf = downloadHandler(
    filename = function(){"RNAseqCNV.pdf"},
    content = function(file){ggsave(file,plot=p_cnv(),width = 20,height = 10)}
  )

  observeEvent(p_cnv(),enable("run_all"))

  observeEvent(p_cnv(),
               showNotification("Done RNAseqCNV Analysis.",duration = 3,type = "warning")
  )

  #Mutation  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  output$body_mutation=renderUI({
    req(input$fileVcf)
    body_mutation()
  })

  out_mutation=reactive({
    shinyjs::disable("run_all")
    req(input$fileVcf)
    getting_fc=showNotification("Extracting Mutations from VCF file ...",duration = NULL,closeButton = F,type = "warning")
    on.exit(removeNotification(getting_fc),add=T)
    get_BALL_mutation(file_vcf = input$fileVcf$datapath)
  })

  output$mutation_all=renderText(paste(out_mutation()$out_text_BALLmutation,collapse = ", "))
  output$mutation_BALLsubtypeDefining=renderText(paste(out_mutation()$out_text_SubtypeDefiningMutation,collapse = ", "))

  MutationTabSum=reactive({
    data.frame(SubtypeDefiningMutation=
                 paste0(out_mutation()$BALL_snv$gene,":",gsub("p.","",out_mutation()$BALL_snv$aaPos))[!is.na(out_mutation()$BALL_snv$subtypesignature)],
               RelatedSubtype=out_mutation()$BALL_snv$subtypesignature[!is.na(out_mutation()$BALL_snv$subtypesignature)],
               stringsAsFactors = F
    )
  })

  observeEvent(out_mutation(),enable("run_all"))


  #Fusion  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  output$body_fusion=renderUI({
    req(!(is.null(input$fileFusioncatcher) & is.null(input$fileCicero)))
    body_fusion()
  })

  #fusioncatcher
  fusionTableFusioncatcher=reactive({
    # shinyjs::disable("run_all")
    req(input$fileFusioncatcher)
    getting_fc=showNotification("Extracting fusions from Fusioncathcer output ...",duration = NULL,closeButton = F,type = "warning")
    on.exit(removeNotification(getting_fc),add=T)
    get_BALL_fusion(input$fileFusioncatcher$datapath,type = "fc")
  })
  output$fusionTableFusioncatcher=renderTable({fusionTableFusioncatcher() %>%
      mutate(PossibleSubtype=NULL,SubtypeInfo=NULL,method=NULL) %>% distinct()})
  observeEvent(fusionTableFusioncatcher(),enable("run_all"))

  #cicero
  fusionTableCicero=reactive({
    # shinyjs::disable("run_all")
    req(input$fileCicero)
    getting_c=showNotification("Extracting fusions from Cicero output ...",duration = NULL,closeButton = F,type = "warning")
    on.exit(removeNotification(getting_c),add=T)
    get_BALL_fusion(input$fileCicero$datapath,type = "c")
  })
  output$fusionTableCicero=renderTable({fusionTableCicero() %>%
      mutate(PossibleSubtype=NULL,SubtypeInfo=NULL,method=NULL) %>% distinct()})

  observeEvent(fusionTableFusioncatcher(),enable("run_all"))
  observeEvent(fusionTableCicero(),enable("run_all"))

  #Summarize   --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  output$body_summarise=renderUI({
    req(
      (!is.null(p_umap()) & !is.null(RNAseqCNV_out())) & (!is.null(fusionTableFusioncatcher())) & (!is.null(fusionTableCicero()))
    )
    body_summarise()
  })

  df_sum=reactive({
    req(!is.null(df_out_svm))
    get_subtype_final(
      id="TestSample",
      df_feateure_exp = GeneExpressionTabSum(),
      df_out_phenograph = df_out_phenograph(),df_out_svm = df_out_svm(),
      out_mutation = out_mutation(),
      chrom_n = chrom_n(),CNV_label = CNV_label(),
      fusion_fc = fusionTableFusioncatcher(),fusion_c = fusionTableCicero())

  })

  tab_geneticAlt=reactive({
    req(df_sum())
    data.frame(FusionCatcher=df_sum()$fusion_fc,Cicero=df_sum()$fusion_c,Mutation=df_sum()$Mutation_Sub_def,RNAseqCNV=CNV_label(),stringsAsFactors = F)
  })
  tab_GEP=reactive({
    req(df_sum())
    data.frame(PhenoGraph=df_sum()$subtype_phenograph,SVM=df_sum()$subtype_svm,stringsAsFactors = F)
  })
  tab_sum=reactive({
    req(df_sum())
    data.frame(Subtype=df_sum()$subtype_final,stringsAsFactors = F)
  })

  output$tab_geneticAlt=renderTable({tab_geneticAlt()})
  output$tab_GEP=renderTable({tab_GEP()})
  output$tab_sum=renderTable({tab_sum()})

  observeEvent(tab_geneticAlt(),enable("run_all"))
  observeEvent(tab_GEP(),enable("run_all"))
  observeEvent(tab_sum(),enable("run_all"))
}


#' server_batch
#'
#' @param input
#' @param output
#' @param session
#'
#' @return
#' @export
#'
#' @examples
server_batch=function(input, output, session){
  observe({
    shinyjs::toggleState("run_batch", condition = !(is.null(input$fileList)))
  })

  #input  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  output$ListingFileExample=renderUI({
    req(input$showListingFileExample)
    if (input$showListingFileExample){
      tableOutput("listingFileExample")
    }
  })

  output$listingFileExample=renderTable(file_listing_example)

  df_listing=reactive({
    req(input$fileList)

    df_=as.data.frame(vroom::vroom(input$fileList$datapath))
    names(df_)=tolower(names(df_))

    vars_in_listing=c('id','count','vcf','fusioncatcher','cicero')

    missingValueTest=any(is.na(as.matrix(df_)))
    missingValueTest1=any(trimws(as.character(as.matrix(df_)))=="")
    varNameTest=!all(vars_in_listing %in% names(df_))

    shinyFeedback::feedbackWarning("fileList",varNameTest | missingValueTest1 | missingValueTest,
                                   ifelse(varNameTest,paste0("Missing field: ",paste0(vars_in_listing[!vars_in_listing %in% names(df_)],collapse = ", ")),
                                          ifelse(missingValueTest,"Missing values in listing file, please check",
                                                 ifelse(missingValueTest1,"Missing values in listing file, please check"," "))))

    req(!varNameTest)
    req(!missingValueTest)
    req(!missingValueTest1)

    df_
  })

  output$listingCheck=renderText({
    req(input$fileList)

    df_=as.data.frame(vroom::vroom(input$fileList$datapath))
    names(df_)=tolower(names(df_))

    vars_in_listing=c('id','count','vcf','fusioncatcher','cicero')

    missingValueTest=any(is.na(as.matrix(df_)))
    missingValueTest1=any(trimws(as.character(as.matrix(df_)))=="")
    varNameTest=!all(vars_in_listing %in% names(df_))

    ifelse(varNameTest,paste0("Missing field: ",paste0(vars_in_listing[!vars_in_listing %in% names(df_)],collapse = " "),", please check"),
           ifelse(missingValueTest,"Missing values in listing file, please check",
                  ifelse(missingValueTest1,"Missing values in listing file, please check"," ")))
  })

  output$recordsInListingFile=renderText({
    paste0("Number of records in listing file: ",nrow(df_listing()))
  })
  output$TopListingDf=renderTable({req(df_listing());head(df_listing())})

  output$TopListing=renderUI({req(df_listing());TopListing()})

  output$para_batch=renderUI({req(df_listing());para_batch()})

  output$run_batch=renderUI({req(df_listing());run_batch()})

  #run Analysis  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  phenographGeneN=reactive(as.numeric(unlist(strsplit(input$phenographGeneN,split = ","))))

  outs_multiple=eventReactive(input$run_batch,{
    shinyjs::disable("run_batch")
    req(df_listing())

    running_batch=showNotification("Running batch-mode Analysis ...",duration = NULL,closeButton = F,type = "warning")
    on.exit(removeNotification(running_batch),add=T)

    run_multiple_samples(input$fileList$datapath,
                         featureN_PG = phenographGeneN(),minReadCnt = input$CNVminReadCnt,minDepth = input$CNVminDepth,mafmin = input$CNVmafMin,mafmax = input$CNVmafMax)
  }
  )

  observeEvent(outs_multiple(),disable("run_batch"))


  #download summary ----
  df_sums_out=reactive({
    req(outs_multiple())
    outs_multiple()$df_sums %>% select(sample_id,fusion_fc,fusion_c,Mutation_BALL,Mutation_Sub_def,RNAseqCNV_ChromN,RNAseqCNV_label,
                                       subtype_phenograph,subtype_phenograph_label,subtype_svm,subtype_svm_label,
                                       subtype_final)
  })

  output$out_multiple_sum.tsv = downloadHandler(
    filename = function(){"out_multiple_sum.tsv"},
    content = function(file){write.table(df_sums_out(),file,col.names = T,row.names = F,quote = F,na = "",sep="\t")}
  )

  observeEvent(df_sums_out(),disable("run_batch"))


  #output  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  observeEvent(df_listing(),{
    updateSelectInput(inputId = "selectedSample",choices = df_listing()$id)
  })

  df_sum_batch=reactive({
    req(outs_multiple())
    outs_multiple()$df_sums[outs_multiple()$df_sums$sample_id==input$selectedSample,]
  })

  CNV_label=reactive({
    req(!is.null(outs_multiple()))

    RNAseqCNV_out=outs_multiple()$RNAseqCNV_outs[[input$selectedSample]]

    paste0(RNAseqCNV_out$df_cnv_out$gender,";\n",RNAseqCNV_out$df_cnv_out$chrom_n,";",RNAseqCNV_out$df_cnv_out$alterations)
  })

  tab_geneticAlt_batch=reactive({
    # req(df_sum_batch())
    data.frame(FusionCatcher=df_sum_batch()$fusion_fc,Cicero=df_sum_batch()$fusion_c,
               Mutation=df_sum_batch()$Mutation_Sub_def,RNAseqCNV=CNV_label(),stringsAsFactors = F)
  })
  tab_GEP_batch=reactive({
    # req(df_sum_batch())
    data.frame(PhenoGraph=df_sum_batch()$subtype_phenograph,SVM=df_sum_batch()$subtype_svm,stringsAsFactors = F)
  })
  tab_sum_batch=reactive({
    # req(df_sum_batch())
    data.frame(Subtype=df_sum_batch()$subtype_final,stringsAsFactors = F)
  })

  output$tab_geneticAlt_batch=renderTable({tab_geneticAlt_batch()})
  output$tab_GEP_batch=renderTable({tab_GEP_batch()})
  output$tab_sum_batch=renderTable({tab_sum_batch()})

  #heatmap ----
  df_pred_batch=reactive({
    req(outs_multiple())
    bind_rows(outs_multiple()$df_outs_phenograph[[input$selectedSample]],outs_multiple()$df_outs_svm[[input$selectedSample]]) %>%
      mutate(N=sprintf("%04d",featureN))
  }
  )

  p_predHeatmap_batch=reactive(
    gg_tilePlot(df_in = df_pred_batch(),x = "N",y = "method",var_col = "pred",x_tick_label_var = "featureN",title = "Prediction Heatmap")
  )

  output$predHeatmap_batch=renderPlot(p_predHeatmap_batch(),height = 200, width = 500)

  output$predHeatmap_batch.pdf = downloadHandler(
    filename = function(){"predHeatmap_batch.pdf"},
    content = function(file){ggsave(file,plot=p_predHeatmap_batch(),width = 6,height = 2)}
  )

  #umap ----
  umaps_batch=reactive({
    req(outs_multiple())
    running_umap=showNotification("Running UMAP ...",duration = NULL,closeButton = F,type = "warning")
    on.exit(removeNotification(running_umap),add=T)

    run_multiple_umap(outs_multiple(),n_neighbors = input$umapNeighbor,variable_n = input$umapGeneN)
  })

  umap_batch=reactive({
    req(umaps_batch())
    umaps_batch()[[input$selectedSample]]
  })

  output$umap_batch=renderPlot(umap_batch(),height = 500, width = 500)

  output$umap_batch.pdf = downloadHandler(
    filename = function(){"umap_batch.pdf"},
    content = function(file){ggsave(file,plot=umap_batch(),width = 11,height = 9)}
  )

  observeEvent(umaps_batch(),disable("run_batch"))


  #Box plot ----
  obj_boxplot_batch=reactive({
    req(outs_multiple())
    obj_merge(obj_in = obj_1821,df_in = outs_multiple()$df_vsts[c("feature",input$selectedSample)],assay_name_in = "vst")
  })

  BoxPlotGene_batch=reactive({
    req(obj_boxplot_batch())
    req(input$BoxPlotGeneName_batch %in% info_gtf_hg38$gene_name[info_gtf_hg38$gene_id %in% rownames(obj_boxplot_batch()$SE)])
    draw_BoxPlot(obj_in = obj_boxplot_batch(),group.by = "diag_raw1",
                 plot_title = "Box Plot of Gene Expression",
                 features = input$BoxPlotGeneName_batch,highlightLevel = "TestSample")
  })

  output$BoxPlotGene_batch=renderPlot({
    req(BoxPlotGene_batch())
    BoxPlotGene_batch()
  },
  height = 450, width = 700)
  observeEvent(obj_boxplot_batch(),disable("run_batch"))

  #RNAseqCNV ----
  p_cnv_batch=reactive({
    req(outs_multiple())
    get_RNAseqCNV_plot(RNAseqCNV_out=outs_multiple()$RNAseqCNV_outs[[input$selectedSample]],sample_name = input$selectedSample)
  })

  output$cnvplot_batch=renderPlot(p_cnv_batch())

  output$RNAseqCNV_batch.pdf = downloadHandler(
    filename = function(){"RNAseqCNV_batch.pdf"},
    content = function(file){ggsave(file,plot=p_cnv_batch(),width = 20,height = 10)}
  )
  observeEvent(p_cnv_batch(),disable("run_batch"))

  observeEvent(p_cnv_batch(),
               showNotification("Done batch-mode Analysis.",duration = 5,type = "warning")
  )
}


server_countOnly=function(input, output, session){
  observe({
    shinyjs::toggleState("run_countOnly", condition = !(is.null(input$fileCountMatrix)))
  })
  #input  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  output$para_countOnly=renderUI({req(input$fileCountMatrix);para_countOnly()})
  output$run_countOnly=renderUI({req(input$fileCountMatrix);run_countOnly()})

  #run Analysis  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  phenographGeneN=reactive(as.numeric(unlist(strsplit(input$phenographGeneN,split = ","))))

  outs_countMatrix=eventReactive(input$run_countOnly,{
    shinyjs::disable("run_countOnly")
    req(input$fileCountMatrix)

    running_count=showNotification("Running Analysis for Count Matrix...",duration = NULL,closeButton = F,type = "warning")
    on.exit(removeNotification(running_count),add=T)

    run_countMatrix(file_countMatrix = input$fileCountMatrix$datapath,featureN_PG = phenographGeneN())
  }
  )

  observeEvent(outs_countMatrix(),disable("run_countOnly"))


  #download summary ----
  df_sums_countMatrix=reactive({
    req(outs_countMatrix())
    outs_countMatrix()$df_sums %>% select(sample_id,subtype_phenograph,subtype_phenograph_label,subtype_svm,subtype_svm_label)
  })

  output$out_countMatrix_sum.tsv = downloadHandler(
    filename = function(){"out_countMatrix_sum.tsv"},
    content = function(file){write.table(df_sums_countMatrix(),file,col.names = T,row.names = F,quote = F,na = "",sep="\t")}
  )

  observeEvent(df_sums_countMatrix(),disable("run_countOnly"))


  #output  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  observeEvent(df_sums_countMatrix(),{
    updateSelectInput(inputId = "selectedSample_countMatrix",choices = df_sums_countMatrix()$sample_id)
  })

  df_sum_countMatrix=reactive({
    req(df_sums_countMatrix())
    df_sums_countMatrix()[df_sums_countMatrix()$sample_id==input$selectedSample_countMatrix,]
  })


  tab_GEP_countMatrix=reactive({
    # req(df_sum_batch())
    data.frame(PhenoGraph=df_sum_countMatrix()$subtype_phenograph,SVM=df_sum_countMatrix()$subtype_svm,stringsAsFactors = F)
  })

  output$tab_GEP_countMatrix=renderTable({tab_GEP_countMatrix()})

  #heatmap ----
  df_pred_countMatrix=reactive({
    req(outs_countMatrix())
    bind_rows(outs_countMatrix()$df_outs_phenograph[[input$selectedSample_countMatrix]],outs_countMatrix()$df_outs_svm[[input$selectedSample_countMatrix]]) %>%
      mutate(N=sprintf("%04d",featureN))
  }
  )

  p_predHeatmap_countMatrix=reactive(
    gg_tilePlot(df_in = df_pred_countMatrix(),x = "N",y = "method",var_col = "pred",x_tick_label_var = "featureN",title = "Prediction Heatmap")
  )

  output$predHeatmap_countMatrix=renderPlot(p_predHeatmap_countMatrix(),height = 200, width = 500)

  output$predHeatmap_countMatrix.pdf = downloadHandler(
    filename = function(){"predHeatmap_countMatrix.pdf"},
    content = function(file){ggsave(file,plot=p_predHeatmap_countMatrix(),width = 6,height = 2)}
  )

  #umap ----
  umaps_countMatrix=reactive({
    req(outs_countMatrix())
    running_umap=showNotification("Running UMAP ...",duration = NULL,closeButton = F,type = "warning")
    on.exit(removeNotification(running_umap),add=T)

    run_multiple_umap(outs_multiple = outs_countMatrix(),n_neighbors = input$umapNeighbor,variable_n = input$umapGeneN)
  })

  umap_countMatrix=reactive({
    req(umaps_countMatrix())
    umaps_countMatrix()[[input$selectedSample_countMatrix]]
  })

  output$umap_countMatrix=renderPlot(umap_countMatrix(),height = 500, width = 500)

  output$umap_countMatrix.pdf = downloadHandler(
    filename = function(){"umap_countMatrix.pdf"},
    content = function(file){ggsave(file,plot=umap_countMatrix(),width = 11,height = 9)}
  )

  observeEvent(umaps_countMatrix(),disable("run_countOnly"))

  #Box plot ----
  obj_boxplot_countMatrix=reactive({
    req(outs_countMatrix())
    req(input$selectedSample_countMatrix %in% names(outs_countMatrix()$df_vsts))
    obj_merge(obj_in = obj_1821,df_in = outs_countMatrix()$df_vsts[c("feature",input$selectedSample_countMatrix)],assay_name_in = "vst")
  })


  BoxPlotGene_countMatrix=reactive({
    req(obj_boxplot_countMatrix())
    # req(input$BoxPlotGeneName_countMatrix %in% info_gtf_hg38$gene_name[info_gtf_hg38$gene_id %in% rownames(obj_boxplot_countMatrix()$SE)])
    draw_BoxPlot(obj_in = obj_boxplot_countMatrix(),group.by = "diag_raw1",
                 plot_title = "Box Plot of Gene Expression",
                 features = input$BoxPlotGeneName_countMatrix,highlightLevel = "TestSample")
  })

  output$BoxPlotGene_countMatrix=renderPlot({
    req(BoxPlotGene_countMatrix())
    BoxPlotGene_countMatrix()
  },
  height = 450, width = 700)
  observeEvent(obj_boxplot_countMatrix(),disable("run_countOnly"))





}



#' Title
#'
#' @param input
#' @param output
#' @param session
#'
#' @return
#' @export
#'
#' @examples
server_SC=function(input, output, session){
  #SC input   --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  observe({
    shinyjs::toggleState("run_sc", condition = !(is.null(input$fileSinglecell)))
  })

  count_sc=reactive({
    req(input$fileSinglecell)
    read_sc_file(input$fileSinglecell$datapath)
  })

  #SC output   --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  output$body_SCoutput=renderUI({
    req(input$fileSinglecell)
    body_SCoutput()
  })

  sc_report=eventReactive(input$run_sc,{
    shinyjs::disable("run_sc")
    running_sc=showNotification("Running B-ALL Subtyping for Single-cell Data ...",duration = NULL,closeButton = F,type = "warning")
    on.exit(removeNotification(running_sc),add=T)

    get_SC_subtypes(count_matrix = count_sc(),SE_celltype = SE_celltype,SE_BALL = SE_BALL)
  })

  output$scReport=renderPlot(sc_report(),width = 1400,height = 1300)

  output$downloadSCreport=renderUI({
    req(sc_report())
    downloadButton("Singlecell_BALL_subtyping.pdf","Download SC report")
  })

  output$Singlecell_BALL_subtyping.pdf = downloadHandler(
    filename = function(){"Singlecell_BALL_subtyping.pdf"},
    content = function(file){ggsave(file,plot=sc_report(),width = 24,height = 25)}
  )

  observeEvent(sc_report(),enable("run_sc"))
  observeEvent(sc_report(),
               showNotification("Done B-ALL Subtyping for Single-cell Data.",duration = 5,type = "warning")
  )
}


#' shinyApp_server
#'
#' @param input
#' @param output
#' @param session
#'
#' @return
#' @export
#'
#' @examples
shinyApp_server=function(input, output, session) {
  # increase upload limit to 500MB
  options(shiny.maxRequestSize = 500*1024^2)

  showNotification(id = "welcome_notif","Welcome to MD-ALL!",duration = 10,type = "message")

  # output$bulkRNAseq_selection=renderUI({
  #   radioButtons("bulkRNAseq_analysisType","Analysis Type:",selected = character(0),
  #                choices = c("Single Sample","Multiple Samples","Count Matrix Only"),inline = T)
  # })

  #Define analysis type -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  output$tabBody_bulkRNA=renderUI({
    if (input$bulkRNAseq_analysisType=="Single Sample"){
      tabBody_bulkRNA_singleSample()
    } else if (input$bulkRNAseq_analysisType=="Multiple Samples") {
      tabBody_bulkRNA_batch()
    } else if (input$bulkRNAseq_analysisType=="Count Matrix Only") {
      tabBody_bulkRNA_countONly()
    }
    # tabBox(width = NULL,
    #   title = "Bluk RNAseq",
    #   tabPanel("Single Sample",     tabBody_bulkRNA_singleSample()),
    #   tabPanel("Multiple Samples",  tabBody_bulkRNA_batch()),
    #   tabPanel("Count Matrix Only", tabBody_bulkRNA_countONly()),
    # )
  })

  observeEvent(input$run_all,disable("bulkRNAseq_analysisType"))
  observeEvent(input$run_batch,disable("bulkRNAseq_analysisType"))
  observeEvent(input$run_countOnly,disable("bulkRNAseq_analysisType"))


  server_singleSample(input, output, session)
  server_batch(input, output, session)
  server_countOnly(input, output, session)
  server_SC(input, output, session)
}

#' run_shiny_MDALL
#'
#' @return
#' @export
#'
#' @examples
run_shiny_MDALL=function(){
  # addResourcePath("image","img/")

  req_librarys=c("dplyr","stringr",
                 "DESeq2","SummarizedExperiment","Rphenograph",
                 "Seurat","SingleR",
                 "ggplot2","ggrepel","cowplot","umap",
                 "shiny","shinyjs","shinydashboard","shinyFeedback","MDALL")

  for(name in req_librarys){suppressMessages(suppressWarnings(library(name,quietly = T,verbose=F,character.only = T)))}

  shinyApp(ui=ui_dashboard,server=shinyApp_server)
}













