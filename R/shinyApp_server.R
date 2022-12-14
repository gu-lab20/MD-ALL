#' shinyApp_server
#'
#' @param input
#' @param output
#'
#' @return
#' @export
#'
#' @examples
shinyApp_server=function(input, output, session) {
  # increase upload limit to 500MB
  options(shiny.maxRequestSize = 500*1024^2)
  #Gene expression  --------------------------------------------------------------------------------------
  #Print first 4 lines of count
  output$fileCountName=reactive(input$fileCount$datapath)
  df_count=reactive(read_input(input$fileCount$datapath))

  output$countTop_table=renderTable({
    if(is.null(input$fileCount)){return(NULL)}
    head(df_count(),4)
  })

  output$countTail_table=renderTable({
    if(is.null(input$fileCount)){return(NULL)}
    tail(df_count(),4)
  })

  #run vst and check first 4 lines
  df_vst=eventReactive(input$run_vst,
                       get_vst_values(obj_in = obj_234_HTSeq,
                                      df_count = df_count())
  )

  output$df_vstTop=renderTable({
    if(is.null(df_vst())){return(NULL)}
    head(df_vst(),4)
  })

  output$df_vstTail=renderTable({
    if(is.null(df_vst())){return(NULL)}
    tail(df_vst(),4)
  })

  GeneExpressionTabSum=reactive({if(is.null(df_vst())){return(NULL)};get_geneExpression(df_vst = df_vst(),genes = c("CDX2","CRLF2","DUX4"))})

  #get box Plot ----

  obj_boxplot=eventReactive(input$run_vst,
                            obj_merge(obj_in = obj_234_HTSeq,df_in = df_vst()))

  p_boxPlot=reactive(draw_BoxPlot(obj_in = obj_boxplot(),group.by = "diag",features = input$BoxPlotGeneName,highlightLevel = "TestSample"))

  output$BoxPlotGene=renderPlot({if(is.null(df_vst())){return(NULL)};p_boxPlot()},height = 650, width = 900)

  #GEP --------------------------------------------------------------------------------------

  #run tsne -----------------------
  obj_=eventReactive(input$run_vst,
    obj_merge(obj_in = obj_2042_HTSeq,df_in = df_vst(),assay_name_in = "vst"))

  obj_tsne=eventReactive(input$run_tsne,
                     run_tsne(obj_in = obj_(),out_label = "tsne",perplexity_one = input$tsneperplexity,variable_n = input$tsneGeneN,feature_panel = "boruta_genes")
  )
  p1=eventReactive(input$run_tsne,
                   draw_DimPlot(obj_tsne(),group.by = "diag",reduction = "tsne",highlightLevel = "TestSample")
  )
  output$tsne=renderPlot(p1(),height = 700, width = 850)

  output$tsne.pdf = downloadHandler(
    filename = function(){"tsne.pdf"},
    content = function(file){ggsave(file,plot=p1(),width = 11,height = 9)}
  )

  df_tsne=reactive({
    if(is.null(obj_tsne())){NULL}
    get_embeding_feature(obj_in = obj_tsne(),features = c("diag"))
  })

  df_knn_tsne=reactive({
    if(is.null(df_tsne())){NULL}
    data.frame(Method="tSNE",
               GEP_prediction=knn_pred_one(indata = df_tsne(),var_y = "diag", var_x = c("tSNE_1","tSNE_2"))[which(df_tsne()$diag=="TestSample")],
               stringsAsFactors = F)
  })

  #run umap -----------------------
  obj_umap=eventReactive(input$run_umap,
                     run_umap(obj_in = obj_(),out_label = "umap",n_neighbors = input$umapNeighbor,variable_n = input$umapGeneN,feature_panel = "boruta_genes")
  )
  p2=eventReactive(input$run_umap,
                   draw_DimPlot(obj_umap(),group.by = "diag",reduction = "umap",highlightLevel = "TestSample")
  )
  output$umap=renderPlot(p2(),height = 700, width = 850)

  output$umap.pdf = downloadHandler(
    filename = function(){"umap.pdf"},
    content = function(file){ggsave(file,plot=p2(),width = 11,height = 9)}
  )

  df_umap=reactive({
    if(is.null(obj_umap())){NULL}
    get_embeding_feature(obj_in = obj_umap(),features = c("diag"),reduction = "umap")
  })

  df_knn_umap=reactive({
    if(is.null(df_umap())){NULL}
    data.frame(Method="uMAP",
               GEP_prediction=knn_pred_one(indata = df_umap(),var_y = "diag", var_x = c("uMAP_1","uMAP_2"))[which(df_umap()$diag=="TestSample")],
               stringsAsFactors = F)
  })

  #run phenograph -----------------
  obj_pheno=eventReactive(input$run_Phenograph,
                            run_PhenoGraph(obj_in = obj_(),feature_panel = "boruta_genes",variable_n = input$phenographGeneN,neighbor_k = input$phenographpNeighbor)
  )

  out_phenograph=eventReactive(input$run_Phenograph,
                               get_PhenoGraphPred(obj_pheno())
  )

  output$phenograph=renderText(out_phenograph()$value)

  df_phenograph=reactive({
    if(is.null(out_phenograph())){NULL}
    data.frame(Method="PhenoGraph Clustering",
               GEP_prediction=out_phenograph()$df$diag_pred,
               stringsAsFactors = F)
  })

  GEPpredictionTabSum=reactive({
    if(is.null(df_knn_tsne()) & is.null(df_knn_umap()) & is.null(df_phenograph())){NULL}
    bind_rows(df_phenograph(),df_knn_tsne(),df_knn_tsne())
  })

  #RNAseqCNV --------------------------------------------------------------------------------------
  output$CNVparameters=reactive(
    paste0("Minimum Read Count of Gene: ",input$CNVminReadCnt,"; Mininum Depth of SNV: ",input$CNVminDepth,"; Min MAF:",input$CNVmafMin, "; Max MAF: ",input$CNVmafMax)
  )

  RNAseqCNV_out=eventReactive(input$runRnaseqcnv,
                              run_RNAseqCNV(df_count = df_count(),
                                            snv_file = input$fileVcf$datapath,
                                            genome_version = "hg38",
                                            minReadCnt = input$CNVminReadCnt,
                                            minDepth = input$CNVminDepth,
                                            mafRange = c(input$CNVmafMin,input$CNVmafMax)),


  )

  output$cnvtext=renderTable(RNAseqCNV_out()$df_cnv_out)

  p_cnv=eventReactive(input$runRnaseqcnv,
                      get_RNAseqCNV_plot(RNAseqCNV_out=RNAseqCNV_out())
  )

  output$cnvplot=renderPlot(p_cnv())

  output$RNAseqCNV.pdf = downloadHandler(
    filename = function(){"RNAseqCNV.pdf"},
    content = function(file){ggsave(file,plot=p_cnv(),width = 20,height = 10)}
  )

  CNVTabSum=reactive({
    data.frame(
      PredictedChromosomeNumber=RNAseqCNV_out()$df_cnv_out$chrom_n,
      Chr21Alteration=ifelse(grepl("21",RNAseqCNV_out()$df_cnv_out$alterations),"Yes","No"),
      stringsAsFactors = F
    )
  })
  #Mutation  --------------------------------------------------------------------------------------
  out_mutation=eventReactive(input$getBALLmutation,
                             get_BALL_mutation(file_vcf = input$fileVcf$datapath)
  )

  output$mutation_all=renderText(paste(out_mutation()$out_text_BALLmutation,collapse = ", "))
  output$mutation_BALLsubtypeDefining=renderText(paste(out_mutation()$out_text_SubtypeDefiningMutation,collapse = ", "))

  MutationTabSum=reactive({
      data.frame(SubtypeDefiningMutation=
                   paste0(out_mutation()$BALL_snv$gene,":",gsub("p.","",out_mutation()$BALL_snv$aaPos))[!is.na(out_mutation()$BALL_snv$subtypesignature)],
                 RelatedSubtype=out_mutation()$BALL_snv$subtypesignature[!is.na(out_mutation()$BALL_snv$subtypesignature)],
                 stringsAsFactors = F
      )

  })

  #Fusion  --------------------------------------------------------------------------------------

  fusionTableFusioncatcher=reactive({if(is.null(input$fileFusioncatcher)){return(NULL)};get_BALL_fusion(input$fileFusioncatcher$datapath,type = "fc")})
  output$fusionTableFusioncatcher=renderTable({if(is.null(fusionTableFusioncatcher())){return(NULL)};fusionTableFusioncatcher()})

  fusionTableCicero=reactive({if(is.null(input$fileCicero)){return(NULL)};get_BALL_fusion(input$fileCicero$datapath,type = "c")})
  output$fusionTableCicero=renderTable({if(is.null(fusionTableCicero())){return(NULL)};fusionTableCicero()})

  FusionTabSum=reactive({
    if(is.null(fusionTableFusioncatcher()) & is.null(fusionTableCicero())){NULL}
    if((!is.null(fusionTableFusioncatcher())) | (!is.null(fusionTableCicero()))){
      bind_rows(fusionTableFusioncatcher()[c(1,4,5,6,7)],fusionTableCicero()[c(1,4,5,6,7)]) %>% distinct()}
  })

  #Summarize   --------------------------------------------------------------------------------------
  output$GeneExpressionTabSum=renderTable({if(is.null(GeneExpressionTabSum())){return(NULL)};GeneExpressionTabSum()})

  output$GEPpredictionTabSum=renderTable({if(is.null(GEPpredictionTabSum())){return(NULL)};GEPpredictionTabSum()})

  output$CNVTabSum=renderTable({if(is.null(CNVTabSum())){return(NULL)};CNVTabSum()})

  output$MutationTabSum=renderTable({if(is.null(MutationTabSum())){return(NULL)};MutationTabSum()})

  output$FusionTabSum=renderTable({if(is.null(FusionTabSum())){return(NULL)};FusionTabSum()})

}


















