###uncomment if you don't have packages prior to this
#install.packages("RPostgreSQL")
#install.packages("jsonlite")
#install.packages("purrr")

library(jsonlite)
require("RPostgreSQL")

# create a connection
# save the password that we can "hide" it as best as we can by collapsing it
pw <- {
  "sat1840"
}

pw2 <- {
  "SAT_1840"
}

# loads the PostgreSQL driver
drv <- dbDriver("PostgreSQL")
# creates a connection to the postgres database
# note that "con" will be used later in each connection to the database
con <- dbConnect(drv, dbname = "umst_old",
                 host = "localhost", port = 5432,
                 user = "kevin", password = pw)


dbListTables(con)

###logging in using ALbert's credentials
con <- dbConnect(drv, dbname="iondb",
                 host="localhost", port= 5432,
                 user = "albert", password = pw2)


dbListTables(con)

res <- dbSendQuery(con, "SELECT * FROM rundb_pluginresult")

test <- dbFetch(res)

test2 <- test[,c("store")]
test2 <- data.frame(test2, stringsAsFactors = FALSE)
colnames(test2) <- c("store")
test2$store <- as.character(test2$store)
###result from grep shows us that the statistics is still there
grep("Uniformity of amplicon coverage", test2$store)
length(grep("Uniformity of amplicon coverage", test2$store))

###only 149/273 have the uniformity statistic

###tries to isolate the string of uniformity
test3 <- sub("*Uniformity of amplicon coverage", "", test2)
test3 <- strsplit(test2$store[5],"Uniformity")
test3

test2.1 <- test[,c("state","store")]
length(which(test2.1$state == "Completed")) #268 completed, so only 5 failed
test2.2 <- test2.1[which(test2.1$state == "Completed"),]
length(grep("Uniformity", test2.2$store))
length(test2.2$store)

###trying for just uniformity at this point
test4 <- sub(".*Uniformity of(.*?),", "", test2$store, perl = TRUE)
test4[3]

###test5 &6 lead to end product
test5 <- gsub(".*Uniformity of amplicon coverage(.*?),", "", test2$store)
test5[3]

test6 <- gsub("^(.*?)Uniformity(.*?):", "",test5)
test6[1:25]



###trying to read it in JSON parser

#mydf <- fromJSON(test2.2$store, simplifyDataFrame = TRUE)

#purrr::map(test2, jsonlite::fromJSON)
mydf <- jsonlite::stream_in(textConnection(gsub("\\n", "", test2.2$store)), flatten = TRUE)
mydf <- jsonlite::fromJSON(gsub("\\n", "", test2.2$store), flatten = TRUE)

###The flatten function works exceptionally well .... needs to be tested to see if it's being used correctly
mydf <- jsonlite::stream_in(textConnection(test2.2$store), flatten = TRUE)
mydf <- fromJSON(test2.2$store, flatten = TRUE)
mydf <- read_json("/mnt/DATA3/Yoda_data_dump/Auto_user_SATProton-143-10192016_SimpaMultifocalProstateCancer_OCPdevDNA_PrRNA_326_363/serialized_Auto_user_SATProton-143-10192016_SimpaMultifocalProstateCancer_OCPdevDNA_PrRNA_326.json",
                  simplifyVector = TRUE)
mydf <- jsonlite::stream_in(textConnection(dbFetch(res)$store), flatten = TRUE)

a <- lapply(paste0("[",test2.2$store,"]"), function(x) jsonlite::fromJSON(x))
#test
a <- lapply(paste0("[",test2.2$store,"]"), function(x) jsonlite::fromJSON(x)) 


#nope below doesn't help too much
b <- lapply(paste0("[",test2.2$store,"]"), function(x) jsonlite::stream_in(textConnection(x), simplifyVector = TRUE,flatten = TRUE)) 
b <- cbind(test, do.call(plyr::rbind.fill, lapply(paste0("[",test2.2$store,"]"), function(x) jsonlite::fromJSON(x))))
test <- a[[268]]
test <- bind_rows(a)


#test <- data.table::rbindlist(a, fill = TRUE)  
#test <- bind_rows(a)
dummy <- unlist(a)
grep("Uniformity", names(dummy))
dummy2 <- dummy[grep("Uniformity of amplicon coverage", names(dummy))]
unique(names(dummy2))

for(i in length(mydf$barcodes)){
  is.na(mydf$barcodes[i]$`Uniformity of base coverage`)
}
nrow(mydf)
colnames(mydf)
length(which(!is.na(mydf[35])))
mytestdf <- mydf[which(!is.na(mydf[35])),]


###testing tidyJSON
#install.packages("tidyjson")

library(tidyjson)
library(plyr)

a <- read_json(path = "/mnt/DATA5/Gandalf_data_keep/AC_PR_58_63_65_66_Good_115/plugin_out/coverageAnalysis_out/results.json", format = c("infer"))

attributes(a)
rm(a)



