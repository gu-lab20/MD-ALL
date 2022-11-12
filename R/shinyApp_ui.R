#' shinyApp_ui
#'
#' @return
#' @export
#'
#' @examples
shinyApp_ui=function(){dashboardPage(
  header(),
  siderbar(),
  body()
)
}


header = function(){
  dashboardHeader(title = "MD-ALL")
}

siderbar = function(){dashboardSidebar(disable = TRUE)}

body = function(){ dashboardBody(
  fluidRow(
    column(width = 3,

           box(width = NULL,status = "primary",
               title = tagList(shiny::icon("file-upload"), "File Upload"),
               fileInput("fileinput", "Upload Count File", multiple = TRUE, accept = ".HTSeq", buttonLabel = "Select Count File..."),
               textInput("sep", "Delimiter :", value = "\t"),checkboxInput("header", "Has Header? ", value = FALSE),
               actionButton("run_vst","Run Normalization"),downloadButton("downloadvst", "Download Normalization"),

               p(),
               actionButton("run_all","Run All")
           ),

           box(width = NULL,status = "warning",
               title="PhenoGraph",
               numericInput("phenographGeneN","GeneN",value = 800),numericInput("phenographpNeighbor","Neighbors",value = 30),
               p(),
               actionButton("run_Phenograph","Run PhenoGraph"),
           ),

           box(width = NULL,status = "success",
               title="tSNE",
               numericInput("tsneGeneN","GeneN",value = 800),numericInput("tsneperplexity","Perplexity",value = 30),
               p(),
               actionButton("run_tsne","Run tSNE"),

           ),

           box(width = NULL,status = "info",
               title="uMAP",
               numericInput("umapGeneN","GeneN",value = 800),numericInput("umapNeighbor","Neighbors",value = 30),
               p(),
               actionButton("run_umap","Run uMAP"),
           )


    ),


    column(width = 9,
           fluidRow(
             tabBox(width = 6,height = "300px",side = "right",
                    title = tagList(shiny::icon("binoculars"), "Data Checking"),
                    tabPanel("Check Normalization",
                             tags$h6("*Wait until 'top rows' Generated"),
                             p(),
                             strong("Top rows:"),
                             tableOutput("df_vstTop"),
                    ),

                    tabPanel("Check Raw data",
                             strong("Top rows:"),
                             tableOutput("countTop_table")
                    )

             ),

             box(width = 6,height = "300px", status = "warning",
                 title="PhenoGraph",
                 textOutput('phenograph')
             )
           ),
           fluidRow(
             box(width = 6,height = "800px",status = "success",
                 title = tagList(shiny::icon("glasses"), "tSNE"),
                 plotOutput("tsne",height = "700px",),
                 downloadButton("tsne.pdf", "Download")),


             box(width = 6,height = "800px",status = "info",
                 title = tagList(shiny::icon("glasses"), "uMAP"),

                 plotOutput("umap",height = "700px",),
                 downloadButton("umap.pdf", "Download"))
           )
    )
  )
)
}




