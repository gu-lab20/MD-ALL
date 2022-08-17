#set -------------
setwd("//isi-dcnl/user_data/zgu_grp/DryLab/GEP_prediction")
#load library --------
library(shiny)
#load functions -----------




ui <- fluidPage(
  #input file
  fileInput("file_count", "Upload count file",buttonLabel = "Select count file..."),
  textOutput("countTop_text"),
  tableOutput("countTop_table"),
  
  #tSNE
  actionButton("run_tSNE", "run tSNE"),
  textOutput("tsne"),
  plotOutput("plot")
  
)

server <- function(input, output, session) {
  tsne_label=reactive("tSNE plot, top MAD gene N=800, perplexity=20")
  output$tsne = renderText(tsne_label())

  output$countTop_text=renderText("Top 10 rows:")
  
  df_count=reactive(read_tsv(input$file_count$datapath))
  
  output$countTop_table=renderTable({
    if(is.null(input$file_count)){return(NULL)}
    head(df_count())
  }
  )
  plot_out=eventReactive(input$run_tSNE,
                         draw_tsne_FromVst(runVST_FromFile(input$file_count$datapath),diag_1449,vst_gene_1499,top_feature_list_1499))
  output$plot = renderPlot(plot_out())
}

shinyApp(ui, server)

# dataset <- reactive({
#   get(input$dataset, "package:ggplot2")
# })
# output$summmry <- renderPrint({
#   summary(dataset())
# })
# output$plot <- renderPlot({
#   plot(dataset)
# })
# 

# 
# 
# 
# # load file
# body <- dashboardBody(
#   fluidRow(
#     tabBox(
#       fileInput("HTSeq1", "Choose HTSeq File", multiple = TRUE, accept = ".HTSeq"),
#       fileInput("vst", "Choose vst Reference", accept = ".rdata")
#     ),
#   
#   
#   
#   
#     tabBox(
#       side = "right", height = "250px",
#       selected = "Tab1",
#       title = tagList(shiny::icon("binoculars"), "Upload Status"),
#       tabPanel("HTSeq Uploaded", tableOutput("out_HTSeq1")),
#       tabPanel("vst Reference Uploaded", tableOutput("out_vst")),
#       tabPanel("Tab3", textOutput("id"))
#   ),
# ),
#   fluidRow(
#     tabBox(
#       title = tagList(shiny::icon("asterisk"), "tabBox status"),
#       tabPanel("Tab1",
#                "Currently selected tab from first box:",
#                verbatimTextOutput("tabset1Selected")
#       ),
#       tabPanel("Tab2", textOutput("x1_head"))
#     )
#   )
# )
# 
# 
# 
# 
# 
# shinyApp(
#   ui = dashboardPage(
#     dashboardHeader(title = "GEP Prediction", titleWidth = 300),
#     dashboardSidebar(disable = TRUE),
#     body
#   ), 
#   
#   
#   
#   server = function(input, output){
#     # increase upload limit to 100MB
#     options(shiny.maxRequestSize = 100*1024^2)
#     # Input file validation
#     
#     output$out_HTSeq1 <- renderTable((
#       if(is.null(input$HTSeq1)){return()}
#       else{
#         input$HTSeq1
#       }
#     ))
#     output$out_vst <- renderTable((
#       if(is.null(input$vst)){return()}
#       else{
#         input$vst
#         }
#     ))
#     
#     
# 
#     id= function(){
#       word(basename(input$HTSeq1,1,sep="[.]")
#     }
#     output$id <- renderText(id)
# 
# 
#     
#     
#     # file= "test/COH000066_D1.DUX4patched.HTSeq"
#     # id=word(basename(file),1,sep="[.]")
#     # output$id <- renderText(id)
#     
#     
#     
#     
#     
#     # prepare count for reference + new sample-------------------
#      # count_one=reactive(read.table(input$HTSeq1,header = F,stringsAsFactors = F))
#      # names(count_one)=reactive(c("feature",id))
#      # count_panel_with_one=reactive(count_gene_234 %>% left_join(count_one))
#      # row.names(count_panel_with_one)=reactive(count_panel_with_one$feature)
#    #  
# 
# 
# 
#     # 
#     # #prepare diag with panel data
#     # if(!id %in% diag_234$COH_sample){
#     #   diag_in=rbind(
#     #     diag_234 %>% mutate(source="N2239") %>% select(COH_sample,diag,source),
#     #     data.frame(
#     #       COH_sample=id,
#     #       stringsAsFactors = F,
#     #       diag="test",
#     #       source="for_prediction"
#     #     )
#     #   )
#     # }
#     # 
#     # if(id %in% diag_234$COH_sample){
#     #   diag_in=diag_234
#     # }
#     # 
#     # diag_in$diag=as.factor(diag_in$diag)
#     # row.names(diag_in)=diag_in$COH_sample
#     # 
#     # #running vst------------------
#     # deseqobject=create_deseqObject_from_matrixdata(count_panel_with_one,diag_in,"COH_sample","diag")
#     # df_vst=run_vst(deseqobject)
#     # 
#     # #output data -----------
#     # df_vst_one=data.frame(
#     #   feature=row.names(df_vst),
#     #   vst=unlist(df_vst[id]),
#     #   stringsAsFactors = F
#     # )
#     # 
#     # names(df_vst_one)[2]=id   
#     
#         
#     
# 
#    #   
#    #   
#      
#      
#   }
# )
#      
# 
#      
# 
