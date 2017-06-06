print("Hello")
library(shiny)
library(shinyFiles)

shinyServer(function(input, output) {
  print(input)
})
