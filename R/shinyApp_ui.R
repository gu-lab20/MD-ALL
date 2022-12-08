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
  dashboardHeader(title = "MD-ALL: Molecular Diagnosis of Acute Lymphoblastic Leukemia")
}

siderbar = function(){dashboardSidebar(
  sidebarMenu(
    menuItem("Input Files", tabName = "inputSidebar", icon = icon("file"),
             fileInput("fileCount", "Upload Count File", multiple = TRUE,accept = "*.*", buttonLabel = "Select Count File..."),
             fileInput("fileVcf", "Upload VCF File", multiple = TRUE, buttonLabel = "Select VCF File..."),
             fileInput("fileFusioncatcher", "Upload FusionCatcher File", multiple = TRUE, buttonLabel = "Select FusionCatcher File..."),
             fileInput("fileCicero", "Upload Cicero File", multiple = TRUE, buttonLabel = "Select Cicero File...")
             # actionButton("run_all","Run All")
    ),
    menuItem("Gene Expression", tabName = "geneExpression", icon = icon("arrow-right")),
    menuItem("GEP-prediction", tabName = "GEP", icon = icon("palette")),
    menuItem("RNAseqCNV", tabName = "CNV", icon = icon("barcode")),
    menuItem("Gene mutation", tabName = "mutation", icon = icon("filter")),
    menuItem("Gene fusion", tabName = "fusion", icon = icon("bolt")),
    menuItem("Summarise", tabName = "summarise", icon = icon("layer-group"))
  )
)}

body_geneExpression=function(){fluidRow(
  fluidRow(
    box(width = 6,title=tagList(shiny::icon("binoculars"), "Raw Count Data"),status = "warning",height = "350px",
        textOutput('fileCountName'),
        fluidRow(
          column(width=6,strong("Top rows:"),tableOutput("countTop_table")),
          column(width=6,strong("Tail rows:"),tableOutput("countTail_table"))
        )
    ),
    box(width = 6,title="Normalized Data",status = "warning",height = "350px",
        actionButton("run_vst","Run Normalization"),
        downloadButton("downloadvst", "Download Normalization"),
        tags$h6("*Wait until 'top rows and tail rows' Generated"),
        fluidRow(
          column(width=6,strong("Top rows:"),tableOutput("df_vstTop")),
          column(width=6,strong("Tail rows:"),tableOutput("df_vstTail"))
        )
    ),

  ),
  fluidRow(
    box(width = 10,height = "900px",status = "info",
        title = tagList(shiny::icon("image"), "Box Plot of Normalized Gene Expression"),
        fluidRow(
          column(width=3,textInput("BoxPlotGeneName", "Input Gene Name", "CRLF2"))
        ),
        plotOutput("BoxPlotGene",height = "700px")
    )
  )
)
}

body_GEP=  function(){fluidRow(
  fluidRow(
    box(width = 6,height = "350px", status = "warning",title="PhenoGraph",
        fluidRow(
          column(width=2,numericInput("phenographGeneN","GeneN",value = 800)),
          column(width=2,numericInput("phenographpNeighbor","Neighbors",value = 30)),
          column(width=6,br(),actionButton("run_Phenograph","Run PhenoGraph"))
        ),
        textOutput('phenograph')
    )
  ),

  fluidRow(
    box(width = 6,height = "900px",status = "success",
        title = tagList(shiny::icon("image"), "tSNE"),
        fluidRow(
          column(width=2,numericInput("tsneGeneN","GeneN",value = 800)),
          column(width=2,numericInput("tsneperplexity","Perplexity",value = 30)),
          column(width=2,br(),actionButton("run_tsne","Run tSNE")),
          column(width=2,br(),downloadButton("tsne.pdf", "Download tSNE"))
        ),

        plotOutput("tsne",height = "700px")
    ),

    box(width = 6,height = "900px",status = "info",
        title = tagList(shiny::icon("image"), "uMAP"),
        fluidRow(
          column(width=2,numericInput("umapGeneN","GeneN",value = 800)),
          column(width=2,numericInput("umapNeighbor","Neighbors",value = 30)),
          column(width=2,br(),actionButton("run_umap","Run uMAP")),
          column(width=2,br(),downloadButton("umap.pdf", "Download uMAP"))
        ),
        plotOutput("umap",height = "700px")
    )
  )
)


}

body_CNV=function(){
  fluidRow(
           box(width = NULL,status = "warning",
               title="RNAseqCNV Parameters",

               fluidRow(
                 # column(width=1,radioButtons('CNVgenomeVesrion',label = "Genome Version",choices = c("hg38","hg19"))),
                 column(width=1,numericInput("CNVminReadCnt","Minimum Read Count of Gene",value = 3),numericInput("CNVminDepth","Mininum Depth of SNV",value = 20)),
                 column(width=3,box(title = "MAF range",numericInput("CNVmafMin","Min",value = 0.1),numericInput("CNVmafMax","Max",value = 0.9))),
                 column(width=1,),
                 column(width=1,),
               ),
               p(strong("Input Parameters: ")),textOutput('CNVparameters'),
               br(),
               actionButton("runRnaseqcnv","Run RNAseqCNV"),
           ),

           box(width = 12,status = "warning",
               title="RNAseqCNV results",
               tableOutput('cnvtext')
           ),

           box(width = 12,height = "800px",status = "success",
               title = tagList(shiny::icon("glasses"), "RNAseqCNV plot"),
               plotOutput("cnvplot",height = "800px",),
               downloadButton("RNAseqCNV.pdf", "Download CNV plot"))
  )
}

body_mutation=function(){
  fluidRow(
    column(width=4,
           box(width = NULL,status = "warning",title="B-ALL Mutations in Sample",
               actionButton("getBALLmutation","Get B-ALL Mutations"),
               br(),
               box(width = NULL,status = "success",height = "400px",title="All B-ALL Mutations",
                   textOutput("mutation_all")),
               box(width = NULL,status = "info",height = "200px",title="B-ALL Subtype Defining Mutations",
                   textOutput("mutation_BALLsubtypeDefining"))
               )
    )
  )
}

body_fusion=function(){
  fluidRow(
    column(width=8,
           box(width=12,status = "warning",title="B-ALL Fusions in Sample",
               br(),
               box(width = 12,status = "success",title="Fusions Detected by FusionCatcher",
                   tableOutput("fusionTableFusioncatcher")),
               box(width = 12,status = "info",title="Fusions Detected by Cicero",
                   tableOutput("fusionTableCicero"))
               )
           ),
  )
}

body_summarise=function(){fluidRow(
  box(width=12,status='warning',title="Summarise",
      br(),
      fluidRow(
        column(width = 5,
               box(width=12,height = "200px",title="Gene Expression",
                   br(),
                   tableOutput("GeneExpressionTabSum")),
               box(width=12,height = "200px",title="GEP prediction",
                   br(),
                   tableOutput("GEPpredictionTabSum")),
               box(width=12,height = "200px",title="CNV",
                   br(),
                   tableOutput("CNVTabSum")),
               box(width=12,height = "200px",title="Mutation",
                   br(),
                   tableOutput("MutationTabSum"))
        ),
        column(width = 7,
               box(width=12,title="Fusion",
                   br(),
                   tableOutput("FusionTabSum"))
        )
      ),
      fluidRow(
        box(width=12,title="Subtype Inference",
            br(),
            textOutput("subtypeInference"))
      )

  )
)}



body = function(){ dashboardBody(
  tabItems(
    tabItem(tabName = "geneExpression",
            body_geneExpression()
    ),
    tabItem(tabName = "GEP",
            body_GEP()
    ),
    tabItem(tabName = "CNV",
            body_CNV()
    ),
    tabItem(tabName = "mutation",
            body_mutation()
    ),
    tabItem(tabName = "fusion",
            body_fusion()
    ),
    tabItem(tabName = "summarise",
            body_summarise()
    )
  )
)
}




