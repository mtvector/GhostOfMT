library(shiny)
library(shinyFiles)
options(shiny.maxRequestSize=500*1024^2) 
# Define UI for slider demo application
shinyUI(pageWithSidebar(
  
  #  Application title
  headerPanel("Plot Favorite Genes Over Time"),
  
  # Sidebar with sliders that demonstrate various available options
  sidebarPanel(width=9,
               # file
               fileInput("filenames", label = "Data file input (support tab delimited file or csv).
                         Input should contain expression values - the first column contains gene names
                         and the first row (don't skip a space) contains standardized sample names. The
                         sample names should be structured as [TIMEPOINT]_[CONDITION1]_[CONDITIONS2],
                         for example: d1_mEpi3_100_NBCN. Need to open in browser to upload multiple files",multiple=T),
               
               # gene name input
               fileInput("genenames", label = "Gene list input. It should contains only one column, which provides gene names.
                         (support .csv,  .txt, .tab)")
  ),mainPanel=mainPanel(h4(textOutput("print0")))
))