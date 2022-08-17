setwd("//isi-dcnl/user_data/zgu_grp/DryLab/GEP_prediction")

library(dplyr)
library(Rtsne)
library(caret)
library(ggrepel)
library(shiny)

ui <- fluidPage(
  fluidRow(
    column(3),
    column(4,tags$h1("First level"),
    tags$br(),
    tags$h6("Sixth level"),
    tags$hr(),
    tags$a(href= "http://www.google.com", "Google"),
    "when can I get the tSNE plot T.T")),
  fluidRow(
    column(3),column(4, tags$p("This is a",tags$strong("sad"),"story"),
    tags$p("another",tags$em("sad"),"story"),
    tags$code("cannot finish the task, sad"),
    tags$hr(),
    tags$img(height = 100,
        width = 100,
        src = "test_out/GEP_Prediction/COH000066_D1.tSNE_2D.png")))
)

server <- function(input, output) {}

shinyApp(ui = ui, server = server)

