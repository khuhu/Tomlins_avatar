library(shiny)
library(hashmap)
library(xlsx)



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



server <- function(input, output, session) {
  data <- reactive({
    listOfSamples <- input$file
    read.csv(file=listOfSamples$datapath, stringsAsFactors = FALSE)
  })
  output$contents <- renderTable({
    data()
  })
  
  output$sqlTable <- renderTable({
    expTable()[10,]
  })
  
  expTable <- reactive({if(is.null(data())){
    return(NULL)
  }
    else{
      print("Starting Program")
      fullTable <- read.table("/mnt/DATA4/kevhu/fullRNACountTable2.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
      
      SampleName <- tryCatch(data()$SampleName, error = function(x) return(NULL))
      Barcode <- tryCatch(data()$Barcode, error = function(x) return(NULL))
      Bed <- tryCatch(data()$Bed, error = function(x) return(NULL))
      Report <- tryCatch(data()$Report, error = function(x) return(NULL))
  
      print("Creating Hash Table")
      
      key <- paste0(fullTable$SampleName,fullTable$Barcode, fullTable$ReportID)
      keyVals <- rep(0, length(key))
      hash <- hashmap(keys = key, values = keyVals)
      lookup <- paste0(SampleName, Barcode, Bed, Report)
      
      for(i in seq_along(lookup)){
        if(hash$has_key(lookup[i]) == TRUE){
          counter <- hash$find(lookup[i])
          hash$insert(lookup[i],c(counter + 1))
        }
      }
      
      listOfNames<- names(which(c(hash$data()) > 0))
      print("Look up done")
      
      exportTable <- fullTable[which(key %in% listOfNames),]
      
      #print(exportTable[1:10,])
      #return(exportTable)
      listOfBedsToSplit <- c(unique(exportTable$Bed))
      listOfSplitTables <- NULL
      
      
      for(i in seq_along(listOfBedsToSplit)){
        a <- NULL
        a <- exportTable[which(exportTable$Bed == listOfBedsToSplit[i]),]
        assign(paste0("RnaCountData", listOfBedsToSplit[i]), a)
        listOfSplitTables <- c(listOfSplitTables, paste0("RnaCountData", listOfBedsToSplit[i]))
      }
      
      #print(listOfSplitTables)
      
      listOfFinalExcelSheets <- NULL
      for(i in seq_along(listOfSplitTables)){
        sampCount <- unique(eval(as.name(listOfSplitTables[i]))$SampleName)
        print(sampCount)
        tmpTable <- NULL
        mappedReadsList <- NULL
        for(j in seq_along(sampCount)){
          tmpTable2 <- NULL
          tmpTable2 <- eval(as.name(listOfSplitTables[i]))[which(eval(as.name(listOfSplitTables[i]))$SampleName == sampCount[j]),]
          tmpTable2 <- tmpTable2[order(tmpTable2$AmpliconID),]
          mappedReadsList <- c(mappedReadsList, unique(tmpTable2$NumberOfMappedReads))
          if(j == 1){
            tmpTable <- cbind(tmpTable, tmpTable2[,c("AmpliconID")])
            tmpTable <- cbind(tmpTable, tmpTable2[,c("Gene")])
            tmpTable <- data.frame(tmpTable, stringsAsFactors = FALSE)
            tmpTable <- cbind(tmpTable, tmpTable2[,ncol(tmpTable2)])
            colnames(tmpTable)[1] <- c("AmpliconID")
            colnames(tmpTable)[2] <- c("Gene")
            colnames(tmpTable)[ncol(tmpTable)] <- sampCount[j]
          }
          else{
            tmpTable <- cbind(tmpTable, tmpTable2[,ncol(tmpTable2)])
            colnames(tmpTable)[ncol(tmpTable)] <- sampCount[j]
          }
        }
        
        ###need to add one more NA when I rerun the table to include the contig ID's
        mappedReadsList <- c(NA,NA, mappedReadsList)
        tmpTable <- rbind(mappedReadsList, tmpTable)
        assign(paste0("transformed",listOfBedsToSplit[i]), tmpTable)
        listOfFinalExcelSheets <- c(listOfFinalExcelSheets, paste0("transformed",listOfBedsToSplit[i]))
        print(eval(as.name(paste0("transformed",listOfBedsToSplit[i])))[1:10,])
      }
      for(i in seq_along(listOfFinalExcelSheets)){
        write.xlsx(x = eval(as.name(listOfFinalExcelSheets[i])), file = "/mnt/DATA4/kevhu/testServerAndi.xlsx",sheetName = as.character(listOfBedsToSplit[i]), row.names = FALSE,append = TRUE)
      }
      return(exportTable)
    }
  })
}
shinyApp(ui = ui, server = server)

