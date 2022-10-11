
library(shiny)

ui <- fluidPage(
  actionButton("name", "Press to go"),
  textInput("n", "Enter name of File", "noName"),
  textOutput("plot")
)

server <- function(input, output) {
  
  nameOfFile <- eventReactive(input$name,{
    input$n
  })
  
  output$plot <- renderText({
    nameOfFile()
  })
}

shinyApp(ui, server)