library(xlsx)
df <- read.xlsx("/mnt/DATA4/kevhu/PR_tissue_v1.panel_annotation.xlsx", sheetIndex = 1, sheetName = "Sheet1", stringsAsFactors = FALSE)
test.data.set <- read.xlsx("/mnt/DATA2/share/Users/Kevin/Dropbox/noName.Sri.xlsx", sheetIndex = 1, sheetName = "WG00196_02092016_Designed", stringsAsFactors = FALSE)


####need to change back to original format for shiny

key2 <- df$Amplicon_ID
keyVals2 <- df$Amplicon.Type
hash2 <- hashmap(key = key2, values = keyVals2)

lookup2 <- paste0(test.data.set$AmpliconID[2:nrow(test.data.set)])

row3 <- NULL
for(i in seq_along(lookup2)){
  if(hash2$has_key(lookup2[i]) == TRUE){
    toBind <- hash2$find(lookup2[i])
  }
  if(hash2$has_key(lookup2[i]) == FALSE){
    toBind <- "None"
  }
  row3 <- rbind(row3, toBind)
}

row3 <- rbind(NA, row3)





### format below
###
####
key2 <- paste0()
keyVals2 <- rep(0, length(key))
hash2 <- hashmap(keys = key, values = keyVals)

#Report <- tryCatch(a$ReportID, error = function(x) return(NULL))
lookup2 <- paste0(SampleName, Barcode, Bed, Report)

for(i in seq_along(lookup2)){
  if(hash$has_key(lookup2[i]) == TRUE){
    counter <- hash$find(lookup2[i])
    hash$insert(lookup2[i],c(counter + 1))
  }
}

listOfNames<- names(which(c(hash2$data()) > 0))