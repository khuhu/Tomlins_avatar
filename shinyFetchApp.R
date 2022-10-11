library(shiny)
library(DBI)
library(RMySQL)


#sqlFunction <- function(){
#  m <- dbDriver("MySQL")
#  con <- dbConnect(m,
#                   user='kevin',
#                   password='sat1840',
#                   host="localhost",
#                   dbname='AnnoDB')
#  dbdf <- dbGetQuery(con, paste0("SELECT * FROM ", tableName()[1]))
#  dbDisconnect(con)
#  dbdf
#}

ui <- fluidPage(
  headerPanel("Tomlins Lab Sequencing Data Pull"),
  sidebarLayout(
    sidebarPanel(fileInput("file", "Choose CSV File",accept = c("text/csv","text/comma-separated-values,text/plain",".csv")),
                 tags$hr()),
    mainPanel(
      tableOutput("contents"),
      tableOutput("sqlTable")
      )
  )
)


###might need to nest the RMySQL stuff into another function or something
server <- function(input, output, session) {
  data <- reactive({
    listOfSamples <- input$file
    read.csv(file=listOfSamples$datapath, stringsAsFactors = FALSE)
  })
  
  output$contents <- renderTable({
    data()
  })
  
  tableIdx <- reactive({read.table(file = "/home/kevhu/data/shinyAppTableIndex.txt", stringsAsFactors = FALSE, header = TRUE, sep = ",")})
  tableName <- reactive({
    if(is.null(data())){
      return()
    }
    tableIdx()$sqlTabs[which(tableIdx()$beds == unique(data()$Bed))]
    })
  
  #observe({
  #  if(!is.null(tableName())){
  #    df <- reactive({dbGetQuery(con,paste0("SELECT * FROM ", tableName()[1]))})
  #  }
  #})
  
  df <- reactive({
    if(is.null(tableName())){
      return()
    }
    else{
      #sqlFunction()
      m <- dbDriver("MySQL")
      con <- dbConnect(m,
                       user='kevin',
                       password='sat1840',
                       host="localhost",
                       dbname='AnnoDB')
      on.exit(dbDisconnect(con), add = TRUE)
      dbGetQuery(con,paste0("SELECT * FROM ", tableName()[1]))
    }
  })
  
  
  output$sqlTable <- renderTable({
    if(is.null(df())){
      return()
    }
    else{
      df()
    }
    })
  
  ###below tests if it can pull the table
  #output$sqlTable <- renderTable({
  #  if(is.null(tableName)){
  #    return()
  #  }
  #    m <- dbDriver("MySQL")
  #    con <- dbConnect(m,
  #                     user='kevin',
  #                     password='sat1840',
  #                     host="localhost",
  #                     dbname='AnnoDB')
  #    on.exit(dbDisconnect(con), add = TRUE)
  #    dummyVar <- paste0("'",data()$SampleName[1],"'")
  #    dummyVar2 <- paste0("'",data()$Barcode[1],"'")
  #    dbGetQuery(con,paste0("SELECT * FROM ", tableName()[1], " WHERE SampleName =", dummyVar," AND ","Barcode = ",dummyVar2," ;"))
  #})
  session$allowReconnect(TRUE)
  #on.exit(dbDisconnect(con), add = TRUE)
}

shinyApp(ui = ui, server = server)







#m <- dbDriver("MySQL")
#con <- dbConnect(m,
#                  user='kevin',
#                  password='sat1840',
#                  host="localhost",
#                  dbname='AnnoDB')
#on.exit(dbDisconnect(con), add = TRUE)
#finalTable <- NULL
#for(i in seq_along(data()$SampleName)){
#  dummyVar <- paste0("'",data()$SampleName[i],"'")
#  dummyVar2 <- paste0("'",data()$Barcode[i],"'")
#  res <- dbGetQuery(con,paste0("SELECT * FROM ", tableName()[1], " WHERE SampleName =", dummyVar," AND ","Barcode = ",dummyVar2," ;"))
#  finalTable <- rbind(finalTable, res)
#}

#tableIdx <- read.table(file = "/home/kevhu/data/shinyAppTableIndex.txt", stringsAsFactors = FALSE, header = TRUE, sep = ",")
#dbListConnections( dbDriver( drv = "MySQL"))
#lapply( dbListConnections( dbDriver( drv = "MySQL")), dbDisconnect)

