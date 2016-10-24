library(shiny)

shinyUI(fluidPage(
  titlePanel("Forecast probabilities"),
  fluidRow(
    column(width=6,
      sliderInput("nComp", "Number of components:", min=1, max=10, value=c(1, 9)),
      checkboxInput("delta6", "use an extra component for strains with a diameter of 6 mm", value=TRUE)),
    column(width=6,
      sliderInput("cbp", "Critical breakpoints:", min=7, max=40, value=c(17, 19)),
      textInput("sigma", "Standard deviation of technical error / mm:", value='1')
    )),
  fluidRow(
    column(width=6, tableOutput("mdlSmmry")),
    column(width=6, tableOutput("errorProbs"))
  ),
  fluidRow(
      plotOutput("masterPlot")
  )
))
