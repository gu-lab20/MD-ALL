#Sidebar -----------------------------------------------------------------------------------------------------------------------------------------
#' sidebar of dashboard
#'
#' @return
#' @export
#'
#' @examples
sidebar = function(){dashboardSidebar(
  sidebarMenu(
    menuItem("Bulk RNAseq", tabName = "bulkRNAseq", icon = icon("bolt")),
    menuItem("Sing-cell RNAseq", tabName = "ScRNAseq", icon = icon("palette")),
    tags$footer(
      tags$p("Zunsong Hu@GuLab@COH",style="color:navy;"),
      # div(tags$img(src="image/Gu-lab-logo.jpg",width="200px")),
      style = "text-align:center; align: center;text-color:black;
                margin: 0px;bottom: 0;width: 100%;padding:10px;
                position: absolute;
                background:white;
                box-sizing:border-box;")
  )
)
}

# '#00adfc'

#body_inputBulkRNA -----------------------------------------------------------------------------------------------------------------------------------------
#' body_inputBulkRNA
#'
#' @return
#' @export
#'
#' @examples
body_inputBulkRNA=function(){  tabPanel("Input Data",
  fluidRow(
    # box(title = "BALL subtyping using Bulk RNAseq data",background = "aqua",width = 12,
    column(width = 4,
           box(title="Upload Count File",width = NULL,solidHeader = T,status = "primary",
               p(strong("Count File")),
               tags$h6("with the first column for ENSG gene IDs and second column for counts"),
               fileInput("fileCount", label = NULL, multiple = F,accept = "*.*", buttonLabel = "Select ..."),
               textOutput('fileCountName'),
               uiOutput("para_umap"),
               uiOutput("para_phenograph")
           )),

    column(width = 4,
           box(title = "Upload VCF File",width = NULL,solidHeader = T,status = "warning",
               p(strong("VCF File")),
               tags$h6("from variants calling software like GATK Haplotypercaller"),
               fileInput("fileVcf", label =NULL, multiple = F, buttonLabel = "Select ..."),
               uiOutput("para_RNAseqCNV")
           )),

    column(width = 4,
           box(title="Upload Fusion File",width = NULL,height = 300,solidHeader = T,status = "success",
               p(strong("FusionCatcher Output")),tags$h6("'final-list_candidate-fusion-genes.txt' from FusionCatcher output"),
               fileInput("fileFusioncatcher", label = NULL, multiple = F, buttonLabel = "Select ..."),

               p(strong("Cicero Output")),tags$h6("'final_fusions.txt' from Cicero output"),
               fileInput("fileCicero", label = NULL, multiple = F, buttonLabel = "Select ..."),
           ))
  ),

  fluidRow(
    column(width = 5),
    column(width = 2,actionButton("run_all","Run",class="btn-block btn-success")),
    column(width = 5)
  )
)

}

#' Title
#'
#' @return
#' @export
#'
#' @examples
para_umap=function(){box(title="UMAP Parameters",width=NULL,status = "warning",
                         # p(strong("Used Gene Number")),
                         # tags$h6("Recommand 1058"),
                         column(width = 6,
                                numericInput("umapGeneN",label = "Used Gene Number",value = 1058)),
                         # p(strong("Neighbors")),
                         # tags$h6("Recommand 10"),
                         column(width = 6,
                                numericInput("umapNeighbor",label = "Neighbors",value = 10))
)}

#' Title
#'
#' @return
#' @export
#'
#' @examples
para_phenograph=function(){box(title="PhenoGraph Parameters",width=NULL,status = "warning",
                               column(width = 6,
                               p(strong("Used Gene Number")),tags$h6("Choose from: 100,200,300,400,500,600,700,800,900,1000,1058. Separated by ,"),
                               textInput("phenographGeneN",label = NULL,value = "100")),
                               column(width = 6,
                               p("Input Gene Numbers:"),
                               textOutput("phenographGeneN"))
)}

#' Title
#'
#' @return
#' @export
#'
#' @examples
para_RNAseqCNV=function(){box(title="RNAseqCNV Parameters",width = NULL,status = "warning",
                              column(width = 6,numericInput("CNVminReadCnt","Minimum Read Count of Gene",value = 3)),
                              column(width = 6,numericInput("CNVminDepth","Mininum Depth of SNV",value = 20)),
                              column(width = 6,numericInput("CNVmafMin","Min of MAF",value = 0.1)),
                              column(width = 6,numericInput("CNVmafMax","Max of MAF",value = 0.85))
)}

#body_inputBulkRNA_batch -----------------------------------------------------------------------------------------------------------------------------------------
#' Title
#'
#' @return
#' @export
#'
#' @examples
body_inputBulkRNA_batch=function(){tabPanel(
  "Input Data",
  shinyFeedback::useShinyFeedback(),
  box(title="Upload Meta Table",width = NULL,solidHeader = T,status = "primary",
      p(strong("File of Meta Table")),
      tags$h6("Five columns are needed: id, count, VCF, fusioncatcher and cicero. Separated by tab (\\t)"),
      tags$h6("No missing value allowed"),
      checkboxInput("showListingFileExample","Show Meta Table Example"),
      uiOutput("ListingFileExample"),
      fileInput("fileList", label = NULL, multiple = F,accept = "*.*", buttonLabel = "Select File of Meta Table ...",width = "50%"),
      textOutput("listingCheck"),
      uiOutput("TopListing"),
  ),

  uiOutput("para_batch"),
  uiOutput("run_batch")
)
}

#' Title
#'
#' @return
#' @export
#'
#' @examples
TopListing=function(){
  box(title=NULL,width = NULL,status = "warning",
  textOutput("recordsInListingFile"),
  strong("Top rows:"),
  tableOutput("TopListingDf")
  )
}

#' Title
#'
#' @return
#' @export
#'
#' @examples
para_batch=function(){
  box(title="Parameters",width = NULL,solidHeader = T,status = "primary",
      column(width = 4,para_umap()),
      column(width = 4,para_phenograph()),
      column(width = 4,para_RNAseqCNV())
  )
}

#' Title
#'
#' @return
#' @export
#'
#' @examples
run_batch=function(){
  fluidRow(
    column(width = 5),
    column(width = 2,actionButton("run_batch","Run",class="btn-block btn-success")),
    column(width = 5)
  )
}

#body_outputBulkRNA_batch  -----------------------------------------------------------------------------------------------------------------------------------------

#' Title
#'
#' @return
#' @export
#'
#' @examples
body_outputBulkRNA_batch=function(){tabPanel(
  "Results",
  shinyFeedback::useShinyFeedback(),

  fluidRow(
    column(width=3,box(title = "Select Sample ID",width = NULL,solidHeader = T,status = "primary",
                       selectInput("selectedSample","Select Sample ID: ",choices = "sample"),
                       downloadButton("out_multiple_sum.tsv", "Download summary"),
    )),
    column(width=9,box(title = "Summary",width = NULL,solidHeader = T,status = "primary",
                       column(width = 6,box(width=NULL,title="Genetic Alteration",status = "info",
                                            tableOutput('tab_geneticAlt_batch'))),
                       column(width = 4,box(width=NULL,title="Gene Expression",status = "info",
                                            tableOutput("tab_GEP_batch"))),
                       column(width = 2,box(width=NULL,title = "Subtype Summary",status = "success",
                                            tableOutput("tab_sum_batch")))
    )),
  ),

  tabsetPanel(
    tabPanel("GEP",
             box(title = "GEP",width = NULL,solidHeader = T,status = "primary",
                 column(width = 5,
                        box(width = NULL,status = "success",collapsible = T,
                            title = tagList(shiny::icon("image"), "GEP Prediction Plots"),
                            downloadButton("predHeatmap_batch.pdf", "Download Heatmap"),
                            downloadButton("umap_batch.pdf", "Download UMAP"),
                            div(style = "text-align: center;",
                                plotOutput("predHeatmap_batch",height = "200px")),
                            div(style = "text-align: center;",
                                plotOutput("umap_batch",height = "500px"))
                        )
                 ),

                 column(width = 7,
                        box(width = NULL,status = "info",collapsible = T,
                            title = tagList(shiny::icon("image"), "Box Plot of Normalized Gene Expression"),
                            fluidRow(column(width=3,textInput("BoxPlotGeneName_batch", "Input Gene Name", "CRLF2"))),
                            div(style = "text-align: center;",
                                plotOutput("BoxPlotGene_batch",height = "660px"))
                        )
                 )
             ),),
    tabPanel("RNAseqCNV",
             box(title = "RNAseqCNV",width = NULL,solidHeader = T,status = "primary",
                 box(width = NULL,status = "success",
                     title = tagList(shiny::icon("glasses"), "RNAseqCNV plot"),
                     downloadButton("RNAseqCNV_batch.pdf", "Download CNV plot"),
                     div(style = "text-align: center;",
                         plotOutput("cnvplot_batch",height  = "700px"))
                 )
             )),
  )
)
}

#body_inputBulkRNA_countONly -----------------------------------------------------------------------------------------------------------------------------------------
#' Title
#'
#' @return
#' @export
#'
#' @examples
body_inputBulkRNA_countONly=function(){
  tabPanel(
    "Input Data",
    shinyFeedback::useShinyFeedback(),

    fluidRow(
      column(width = 6,
             box(title="Upload Count Matrix",width = NULL,solidHeader = T,status = "primary",
                 p(strong("File of Count Matrix")),
                 tags$h6("Rows for genes and columns for samples. ENSG gene IDs in first column."),
                 tags$h6("No missing value allowed."),

                 fileInput("fileCountMatrix", label = NULL, multiple = F,accept = "*.*", buttonLabel = "Select File of Count Matrix ...",width = "50%"),
                 textOutput("CountMatrixCheck"),
             )
             ),
      column(width = 6,
             uiOutput("para_countOnly")
             )
    ),

    uiOutput("run_countOnly")
  )
}

#' para_countOnly
#'
#' @return
#' @export
#'
#' @examples
para_countOnly=function(){
  box(title="Parameters",width = NULL,solidHeader = T,status = "primary",
      column(width = 12,para_umap()),
      column(width = 12,para_phenograph())
  )
}

#' run_countOnly
#'
#' @return
#' @export
#'
#' @examples
run_countOnly=function(){
  fluidRow(
    column(width = 5),
    column(width = 2,actionButton("run_countOnly","Run",class="btn-block btn-success")),
    column(width = 5)
  )
}

#body_outputBulkRNA_countONly  -----------------------------------------------------------------------------------------------------------------------------------------
#' body_outputBulkRNA_countONly
#'
#' @return
#' @export
#'
#' @examples
body_outputBulkRNA_countONly=function(){
  tabPanel(
    "Results",
    shinyFeedback::useShinyFeedback(),

    fluidRow(
      column(width=3,box(title = "Select Sample ID",width = NULL,solidHeader = T,status = "primary",
                         selectInput("selectedSample_countMatrix","Select Sample ID: ",choices = "sample"),
                         downloadButton("out_countMatrix_sum.tsv", "Download summary"),
      )),
      column(width=9,box(title = "Summary",width = NULL,solidHeader = T,status = "primary",
                         column(width = 12,box(width=NULL,title="Gene Expression",status = "info",
                                              tableOutput("tab_GEP_countMatrix")))

      )),
    ),

    box(title = "GEP",width = NULL,solidHeader = T,status = "primary",
        column(width = 5,
               box(width = NULL,status = "success",collapsible = T,
                   title = tagList(shiny::icon("image"), "GEP Prediction Plots"),
                   downloadButton("predHeatmap_countMatrix.pdf", "Download Heatmap"),
                   downloadButton("umap_countMatrix.pdf", "Download UMAP"),
                   div(style = "text-align: center;",
                       plotOutput("predHeatmap_countMatrix",height = "200px")),
                   div(style = "text-align: center;",
                       plotOutput("umap_countMatrix",height = "500px"))
               )),

        column(width = 7,
               box(width = NULL,status = "info",collapsible = T,
                   title = tagList(shiny::icon("image"), "Box Plot of Normalized Gene Expression"),
                   fluidRow(column(width=3,textInput("BoxPlotGeneName_countMatrix", "Input Gene Name", "CRLF2"))),
                   div(style = "text-align: center;",
                       plotOutput("BoxPlotGene_countMatrix",height = "660px")),
                      textOutput("testtab")
               ))
    )
  )

}


#body_GEP -----------------------------------------------------------------------------------------------------------------------------------------
#' body_GEP
#'
#' @return
#' @export
#'
#' @examples
body_GEP=function(){fluidRow(
  fluidRow(
    box(width = 6,status = "warning",collapsible = T,
        title=tagList(shiny::icon("binoculars"), "Raw Count Data"),
        fluidRow(
          column(width=6,strong("Top rows:"),tableOutput("countTop_table")),
          column(width=6,strong("Tail rows:"),tableOutput("countTail_table"))
        )
    ),
    box(width = 6,status = "warning",collapsible = T,
        title=tagList(shiny::icon("binoculars"), "Normalized Data"),
        # tags$h6("*Wait until 'top rows and tail rows' Generated"),
        fluidRow(
          column(width=6,strong("Top rows:"),tableOutput("df_vstTop")),
          column(width=6,strong("Tail rows:"),tableOutput("df_vstTail"))
        ),
        uiOutput("downloadvst")
    )
  ),

  fluidRow(
    column(width = 5,
           box(width = NULL,status = "success",collapsible = T,
               title = tagList(shiny::icon("image"), "GEP Prediction Plots"),
               downloadButton("predHeatmap.pdf", "Download Heatmap"),
               downloadButton("umap.pdf", "Download UMAP"),
               div(style = "text-align: center;",
                   plotOutput("predHeatmap",height = "200px")),
               div(style = "text-align: center;",
                   plotOutput("umap",height = "500px"))
           )
    ),

    column(width = 7,
           box(width = 12,status = "info",collapsible = T,
               title = tagList(shiny::icon("image"), "Box Plot of Normalized Gene Expression"),
               fluidRow(column(width=3,textInput("BoxPlotGeneName", "Input Gene Name", "CRLF2"))),
               div(style = "text-align: center;",
                   plotOutput("BoxPlotGene",height = "660px"))
           )
    )
  )
)
}

#body_RNAseqCNV -----------------------------------------------------------------------------------------------------------------------------------------
#' body_CNV
#'
#' @return
#' @export
#'
#' @examples
body_CNV=function(){fluidRow(
  box(width = 12,status = "warning",
      p(strong("Input Parameters: ")),
      textOutput('CNVparameters'),
  ),

  box(width = 12,status = "warning",
      title="RNAseqCNV results",
      tableOutput('cnvtext')
  ),

  box(width = 12,status = "success",
      title = tagList(shiny::icon("glasses"), "RNAseqCNV plot"),
      downloadButton("RNAseqCNV.pdf", "Download CNV plot"),
      div(style = "text-align: center;",
          plotOutput("cnvplot",height  = "700px"))
  )
)
}
#body_Genemutation -----------------------------------------------------------------------------------------------------------------------------------------
#' body_mutation
#'
#' @return
#' @export
#'
#' @examples
body_mutation=function(){fluidRow(
  column(width=6,
             box(width = NULL,status = "success",title="B-ALL Mutations",
                 textOutput("mutation_all")),
             box(width = NULL,status = "info",title="B-ALL Subtype Defining Mutations",
                 textOutput("mutation_BALLsubtypeDefining"))
  )
)
}
#body_Gene fusion -----------------------------------------------------------------------------------------------------------------------------------------
#' body_fusion
#'
#' @return
#' @export
#'
#' @examples
body_fusion=  function(){fluidRow(
  column(width=6,
             box(width = 12,status = "success",title="Fusions Detected by FusionCatcher",
                 tableOutput("fusionTableFusioncatcher")),
             box(width = 12,status = "info",title="Fusions Detected by Cicero",
                 tableOutput("fusionTableCicero"))
  )
)
}

#body_Summarise -----------------------------------------------------------------------------------------------------------------------------------------
#' body_summarise
#'
#' @return
#' @export
#'
#' @examples
body_summarise=function(){
      fluidRow(
        column(width = 5,box(width=NULL,title="Genetic Alteration",solidHeader = T,status = "info",
                             tableOutput('tab_geneticAlt'))),
        column(width = 3,box(width=NULL,title="Gene Expression",solidHeader = T,status = "info",
                             tableOutput("tab_GEP"))),
        column(width = 2,box(width=NULL,title = "Subtype Summary",solidHeader = T,status = "success",
                             tableOutput("tab_sum")))
      )
}

#Tab body bulkRNA single sample -----------------------------------------------------------------------------------------------------------------------------------------
#' Title
#'
#' @return
#' @export
#'
#' @examples
tabBody_bulkRNA_singleSample=function(){tabsetPanel(
body_inputBulkRNA(),
  tabPanel("GEP Prediction",uiOutput('body_GEP')),
  tabPanel("RNAseqCNV",uiOutput('body_CNV')),
  tabPanel("Gene mutation",uiOutput('body_mutation')),
  tabPanel("Gene fusion",uiOutput('body_fusion')),
  tabPanel("Summary",uiOutput('body_summarise'))
)}

#Tab body bulkRNA batch analysis -----------------------------------------------------------------------------------------------------------------------------------------
#' Title
#'
#' @return
#' @export
#'
#' @examples
tabBody_bulkRNA_batch=function(){tabsetPanel(
  body_inputBulkRNA_batch(),
  body_outputBulkRNA_batch()
)}

#Tab body bulkRNA count Matrix Only -----------------------------------------------------------------------------------------------------------------------------------------
#' Title
#'
#' @return
#' @export
#'
#' @examples
tabBody_bulkRNA_countONly=function(){tabsetPanel(
  body_inputBulkRNA_countONly(),
  body_outputBulkRNA_countONly()
)}

#body_inputSC -----------------------------------------------------------------------------------------------------------------------------------------
#' body_SCinput
#'
#' @return
#' @export
#'
#' @examples
body_SCinput=function(){fluidRow(
  # box(title = "BALL subtyping using single-cell RNAseq data",background = "teal",width = 12,
  fluidRow(column(width = 4,
                  box(title = "Upload single-cell count File",width = NULL,solidHeader = T,status = "primary",
                      p(strong("Single-cell Count File")),tags$h6("with rows for genes and columns for cells"),
                      fileInput("fileSinglecell",label =  NULL, multiple = F, buttonLabel = "Select ...")
                  ),
  ),
  column(8)),
  fluidRow(
    column(width = 1),
    column(width = 2,actionButton("run_sc","Run",class="btn-block btn-success")),
    column(width = 9)
  )
)
}

#body_outputSC -----------------------------------------------------------------------------------------------------------------------------------------
#' body_SCoutput
#'
#' @return
#' @export
#'
#' @examples
body_SCoutput=function(){fluidRow(
  h1("Single cell B-ALL subtyping report"),
  uiOutput("downloadSCreport"),
  plotOutput("scReport")
)
}

#dashboardBody -----------------------------------------------------------------------------------------------------------------------------------------
#' body of dashboard
#'
#' @return
#' @export
#'
#' @examples
body = function(){dashboardBody(
  tags$style("#shiny-notification-showme_notif {margin:20px;}"),
  useShinyjs(),
  tabItems(
    tabItem(tabName = "bulkRNAseq",
            # uiOutput("bulkRNAseq_selection"),
            radioButtons("bulkRNAseq_analysisType","Analysis Type:",selected = NULL,
                         choices = c("Single Sample","Multiple Samples","Count Matrix Only"),inline = T),
            uiOutput("tabBody_bulkRNA"),
    ),

    tabItem(tabName = "ScRNAseq",
            tabsetPanel(
              tabPanel("Input data",body_SCinput()),
              tabPanel("Single-cell BALL subtyping",uiOutput("body_SCoutput"))
            )
    )
  )
)
}

#ui_dashboard -----------------------------------------------------------------------------------------------------------------------------------------
#' ui_dashboard
#'
#' @return
#' @export
#'
#' @examples
ui_dashboard=function(){dashboardPage(
  dashboardHeader(title = "MD-ALL: Molecular Diagnosis of Acute Lymphoblastic Leukemia",
                  titleWidth = 600),
  sidebar(),
  body()
)
}






















