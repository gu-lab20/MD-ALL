library(stringr)
library(dplyr)
library(DESeq2)
library(shiny)
library(Rtsne)
library(caret)
library(ggrepel)
library(shinydashboard)
library(shinythemes)

gc(rm(list=ls()))



body <- dashboardBody(
  fluidRow(
    tabBox(side = "right", height = "300px", width = 4,
           title = tagList(shiny::icon("file-upload"), "File Upload"),
      fileInput("HTSeq1", "Upload HTSeq File", multiple = TRUE, accept = ".HTSeq", buttonLabel = "Select Count File..."),
      p(),tags$h5("Don't Forget to Check", p(), "File Status and VST Normalization"),
      ),
    tabBox(
      side = "right", height = "300px",width = 8,
      selected = "HTSeq Uploaded",
      title = tagList(shiny::icon("binoculars"), "Count Data Summary"),
      tabPanel("HTSeq Uploaded", tableOutput("out_HTSeq1")),
      tabPanel("Check the File", p(),strong("Top 5 rows:"), tableOutput("countTop_table")),
      tabPanel("Check Normalization", tags$h6("*Wait until 'top 5 rows' Generated"),p(),strong("Top 5 rows:"), tableOutput("df_vst_Table")),
      tabPanel("Download",actionButton("run_batch","Run Batch Analysis"), p(),downloadButton("downloadfile", "Download"))        
    )
  ),
  fluidRow(
    tabBox(side = "right", height = "800px", width = 12,
      title = tagList(shiny::icon("glasses"), "Visualization"),
      selected = "tSNE plot",
      tabPanel("Box plot", selectizeInput("gene","Type the gene name",choices=NULL,multiple=FALSE), 
               actionButton("run_boxplot", "Run Box Plot"),textOutput("box plot"),plotOutput("plot_boxplot")
      ),
      tabPanel("neighbor heatmap", radioButtons("FeatureN_interval","Select top MAD gene interval",choices = list(50,100,200), selected = 100), 
               # radioButtons("neighbors","Select Number of Nearest Neighbors",choices = list(5,10,20), selected = 10),
               actionButton("run_heatmap", "Run Heatmap"),textOutput("heatmap"),plotOutput("plot_heatmap")
      ),
      tabPanel("tSNE plot", sliderInput("tsne_perplexity","Select Perplexity",min=5, max=95, value=30, step=5), 
               actionButton("run_tSNE", "Run tSNE"),textOutput("tsne"),plotOutput("plot_tsne", click="plot_click"),
               hr(),p(textOutput("info")),p(tableOutput("plot_info"))
      )
    )
  ),
  fluidRow(
    tabBox(side = "right", height = "400px", width = 12,
           title = tagList(shiny::icon("calculator"), "Prediction"),
      # tabPanel("KNN Prediction",actionButton("run_KNN","run KNN Prediction"),tableOutput("KNN")),
      tabPanel("Phenograph Prediction", actionButton("run_phenograph","run Phenograph Prediction"),tableOutput("phenograph"))
    )
  )
)





shinyApp(
  ui = dashboardPage(skin = "purple",
    dashboardHeader(title = "GEP Prediction", titleWidth = 300),
    dashboardSidebar(disable=TRUE),
    body),
  
  server = function(input, output, session){
    # increase upload limit to 500MB
    options(shiny.maxRequestSize = 500*1024^2)
    
    # Input file validation
    output$out_HTSeq1 <- renderTable((
      if(is.null(input$HTSeq1)){return()}
      else{
        input$HTSeq1
      }
    ))
    
    #Print first 4 lines of file
    df_count=reactive(read_tsv(input$HTSeq1$datapath))
    output$countTop_table=renderTable({
      if(is.null(input$HTSeq1)){return(NULL)}
        head(df_count(),4)
    })
    
    #run vst and check first 4 lines
    df_vst=reactive(runVST_FromFile(input$HTSeq1$datapath))
    output$df_vst_Table=renderTable({
      if(is.null(input$HTSeq1)){return(NULL)}
      head(df_vst(),4)
      })
    
    # run tsne
      tsne_perplexity_value=reactive(input$tsne_perplexity)
      output$tsne=renderText(paste("tSNE plot, top MAD gene N=800, perplexity=",tsne_perplexity_value()))
      plot_out_tsne=eventReactive(input$run_tSNE,
                               draw_tsne_FromVst(df_vst(),diag_1409,for_tsne_1409_,top_feature_list_1409,tsne_perplexity_value()))
      output$plot_tsne=renderPlot(plot_out_tsne())
    
    
      tsne_table=eventReactive(input$run_tSNE,
                             tsne_df_FromVst(df_vst(),diag_1409,for_tsne_1409_,top_feature_list_1409,tsne_perplexity_value()))
    # tsne_perplexity_value=reactive(input$tsne_perplexity)
    # output$tsne=renderText(paste("tSNE plot, top MAD gene N=800, perplexity=",tsne_perplexity_value()))
    # 
    # plot_out_tsne=eventReactive(input$run_tSNE,
    #                             draw_tsne_MAD(tsne_table(),diag_1409,diagCol,800,tsne_perplexity_value(),10))
    # output$plot_tsne=renderPlot(plot_out_tsne())
    # data_plot=reactive(tsne_table() %>% left_join(diag_1409) %>% na.omit())
    # diag_col=reactive(diagCol[names(diagCol) %in% diag_1409$diag])
    # tsne_plot=eventReactive(input$run_tSNE,
    #                    ggplot() + 
    #                      xlab(paste0("tSNE dimension 1\ntop gene num = ", 800, "; perplexity=", tsne_perplexity_value())) +    ylab("tSNE dimension 2")+
    #                      theme_bw() +
    #                      
    #                      geom_point(data=subset(data_plot(), diag != "_Prediction" ), aes(X, Y, color=diag), size=1)+
    #                      geom_point(data=subset(data_plot(), diag == "_Prediction" ), aes(X, Y, color=diag), size=3, alpha=0.75)+
    #                      geom_text_repel(data=subset(data_plot(), diag == "_Prediction" ), aes(X, Y, label=COH_sample), fontface = "bold")+
    #                      scale_color_manual(values = diag_col()) +
    #                      
    #                      scale_x_continuous(breaks = round(seq(floor(min(data_plot()$X)), ceiling(max(data_plot()$X)), by = 10),1)) +
    #                      scale_y_continuous(breaks = round(seq(floor(min(data_plot()$Y)), ceiling(max(data_plot()$Y)), by = 10),1)) +
    #                      
    #                      theme(axis.title=element_text(size=10),
    #                            #axis.text=element_blank(),
    #                            #axis.ticks=element_blank(),
    #                            legend.text = element_text(size=10),
    #                            legend.title = element_text(size=10,face="bold")
    #                      ) +
    #                      scale_fill_discrete(name="Experimental\nCondition")+
    #                      guides(color = guide_legend(override.aes = list(size = 1)))
    #                    
    # )
    # output$plot_tsne=renderPlot(tsne_plot())
    
    # run heatmap
    FeatureN_interval_=reactive(input$FeatureN_interval)
    val=reactiveValues()
    observe(val$FeatureN_interval_value<-as.numeric(input$FeatureN_interval))
    heatmap_label=renderText(paste("neighbor heatmap, top MAD gene N= 100-2000, interval=", FeatureN_interval_(),", neighbors= 20"))
    output$heatmap=renderText(heatmap_label())
    plot_out_heatmap=eventReactive(input$run_heatmap,
                                   draw_neighbor_heatmap_FromVst(df_vst(),for_tsne_1409_,diag_1409,top_feature_list_1409,
                                                                 seq(100,2000,val$FeatureN_interval_value))
    )
    output$plot_heatmap=renderPlot(plot_out_heatmap())
    
    # Phenograph prediction
    phenograph_=eventReactive(input$run_phenograph,
                              Phenograph_pred_out(df_vst(),for_tsne_1409_,diag_1409,top_feature_list_1409,FeatureN_list = seq(500,2000,500),top_neighborN_list = c(5,10,20))
    )
    output$phenograph=renderTable(phenograph_())
    
    # KNN prediction
    # COH_sample=reactive(substr(basename(input$HTSeq1$name),1,12))
    # # output$COH=renderText(COH_sample())
    # 
    # KNN_=eventReactive(input$run_KNN,
    #                            KNN_pred_out(df_vst(),COH_sample(),for_tsne_1409_,diag_1409,top_feature_list_1409,seq(200,2000,500),c(20,30))
    #                            )
    # output$KNN=renderTable(KNN_())

    # check tsne dot info
    tsne_info=reactive(left_join(tsne_table(),diag_1409, by="COH_sample"))
    tsne_info_=reactive(dplyr::select(tsne_info(),"COH_sample","X","Y","diag.x","fusion"))
    output$plot_info=renderTable({
      req(input$plot_click)
      nearPoints(tsne_info_(), input$plot_click, threshold = 10, maxpoints=1)
    })
    

    # Box plot
    updateSelectizeInput(session, 'gene', choices = gene_list, server = TRUE)
    gene_sep=reactive(input$gene)
    # gene_input=reactive(input$gene)
    # gene_sep=eventReactive(input$gene,word(gene_input(),1,sep=" "))
    box_plot=eventReactive(input$run_boxplot,
                           draw_boxplot(for_tsne_1409_,df_vst(),gene_sep(),gene_ref_,diag_1409)
    )
    output$plot_boxplot=renderPlot(box_plot())


    # Download output
    output$downloadfile=downloadHandler(
      filename = function(){
        paste(input$HTSeq1, ".csv", sep = "")
      },
      content = function(file) {
        write.csv(phenograph_(), file, row.names = FALSE)
      }
    )
    # output$downloadfile <- downloadHandler(
    #   filename = "Downloads.zip",
    #   content = function(file){
    #     withProgress(message = "Writing Files to Disk. Please wait...", {
    #       temp <- setwd(tempdir())
    #       on.exit(setwd(temp))
    #       files <- c("phenograph.csv", "tsne.pdf")
    #       
    #       phenograph_ <- phenograph_() # Reading my reactive values in advance of launcing the future
    #       plot_out_tsne <- plot_out_tsne()
    #       
    #       future({
    #         write.csv(phenograph_, "phenograph.csv")
    #         pdf(plot_out_tsne, "tsne.pdf")
    #         zip(zipfile = file, files = files)
    #       }) 
    #     })
    #   }
    # )
    

    
    #Batch Analysis
    # plot_out_tsne1=eventReactive(input$run_batch,
    #                             draw_tsne_FromVst(df_vst(),diag_1409,for_tsne_1409_,top_feature_list_1409,tsne_perplexity_value()))
    # output$plot_tsne=renderPlot(plot_out_tsne1())
    # 
    # 
    # plot_out_heatmap1=eventReactive(input$run_batch,
    #                                draw_neighbor_heatmap_FromVst(df_vst(),for_tsne_1409_,diag_1409,top_feature_list_1409,
    #                                                              seq(100,2000,val$FeatureN_interval_value)
    #                                )
    # )
    # output$plot_heatmap=renderPlot(plot_out_heatmap1())
    # 
    # phenograph_1=eventReactive(input$run_batch,
    #                           Phenograph_pred_out(df_vst(),for_tsne_1409_,diag_1409,top_feature_list_1409,FeatureN_list = seq(500,2000,500),top_neighborN_list = c(5,10,20))
    # )
    # output$phenograph=renderTable(phenograph_1())
    # 
    # updateSelectizeInput(session, 'gene', choices = gene_list, server = TRUE)
    # box_plot1=eventReactive(input$run_batch,
    #                        draw_box_plot(df_vst(),diag_1409,for_tsne_1409_,gene_sep())
    #                        )
    # output$plot_boxplot=renderPlot(box_plot1())

    
    
  }
)