library(shiny)
source('../../Scripts/loadLibs.R')
source('../../Scripts/helper.R')

x <- as.double(as.matrix(read.table('tob.dat')))

shinyServer(function(input, output) {
  
  sigma <- reactive({  
    max(0, as.numeric(input$sigma))
  })
  
  mdl <- reactive({
    if (input$delta6) {
      nComp <- pmax(1, input$nComp - 1)
    } else {
      nComp <- input$nComp
    }
    getEligibleMdl(x, delta6=input$delta6, modelName='V', nComp=nComp[1]:nComp[2], sigma=sigma())
  })
  
  master <- reactive({
    return(masterPlot(x, mdl(), cbp=data.frame(S=input$cbp[1], R=input$cbp[2]), sigma()))
  })
  
  output$mdlSmmry <- renderTable({
    t <- getSummaryTable(mdl())
  }, include.rownames=FALSE)
  
  output$masterPlot <- renderPlot({
    grid.arrange(master()$plot)
  })
  
  output$errorProbs <- renderTable({
    formatErrorProbs(master()$errorProbs)
  })
})
