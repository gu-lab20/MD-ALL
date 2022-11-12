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

  # # Input file validation
  # output$out_HTSeq1 <- renderTable((
  #   if(is.null(input$HTSeq1)){return()}
  #   else{
  #     input$HTSeq1
  #   }
  # ))

  #Print first 4 lines of file
  df_count=reactive(read_input(input$fileinput$datapath))

  output$countTop_table=renderTable({
    if(is.null(input$fileinput)){return(NULL)}
    head(df_count(),4)
  })

  output$countTtail_table=renderTable({
    if(is.null(input$fileinput)){return(NULL)}
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

  #run tsne -----------------------
  obj_=eventReactive(input$run_vst,
    obj_merge(obj_in = obj_2042_HTSeq,df_in = df_vst(),assay_name_in = "vst"))

  obj_tsne=eventReactive(input$run_tsne,
                     run_tsne(obj_in = obj_(),out_label = "tsne",feature_panel = "boruta_genes")
  )
  p1=eventReactive(input$run_tsne,
                   draw_DimPlot(obj_tsne(),group.by = "diag",reduction = "tsne",highlightLevel = "TestSample")
  )
  output$tsne=renderPlot(p1(),height = 700, width = 750)

  output$tsne.pdf = downloadHandler(
    filename = function(){"tsne.pdf"},
    content = function(file){ggsave(file,plot=p1(),width = 11,height = 9)}
  )

  #run umap -----------------------
  obj_umap=eventReactive(input$run_umap,
                     run_umap(obj_in = obj_(),out_label = "umap",n_neighbors = 30,variable_n = 800,feature_panel = "boruta_genes")
  )
  p2=eventReactive(input$run_umap,
                   draw_DimPlot(obj_umap(),group.by = "diag",reduction = "umap",highlightLevel = "TestSample")
  )
  output$umap=renderPlot(p2(),height = 700, width = 750)

  output$umap.pdf = downloadHandler(
    filename = function(){"umap.pdf"},
    content = function(file){ggsave(file,plot=p2(),width = 11,height = 9)}
  )

  #run phenograph -----------------
  obj_pheno=eventReactive(input$run_Phenograph,
                            run_PhenoGraph(obj_in = obj_(),feature_panel = "boruta_genes",variable_n = 800,neighbor_k = 30)
  )

  out_phenograph=eventReactive(input$run_Phenograph,
                               get_PhenoGraphPred(obj_pheno())
  )

  output$phenograph=renderText(out_phenograph())





}
