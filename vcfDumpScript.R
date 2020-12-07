#vcfDir <- "/mnt/DATA4/kevhu/choLab/vcfs/test/"
#extensions <- c(".vcf.gz",".vcf.gz.tbi")
#directories <- c("/mnt/DATA3/Yoda_data_dump/Auto_user_SATProton-173-Cho_mouse_20171130_362_439/plugin_out/",
#                 "/mnt/DATA3/Yoda_data_dump/Auto_user_SATProton-174-Cho_mouse_20171213_363_441/plugin_out/");

###below is andi's stuff

vcfDir <- "/mnt/DATA4/kevhu/tmp2/";
extensions <- c(".vcf.gz",".vcf.gz.tbi")
setwd("/mnt/DATA3/Zeus_data_dump/AC17.2.1_Reanalysis_3_395/plugin_out/")

#!/usr/bin/env Rscript
.libPaths(c("/home/kevhu/R/x86_64-pc-linux-gnu-library/3.2",.libPaths()))
library(optparse)
library(stringr)

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="list of directories, should be full path with header 'Directories'", metavar="character"),
  make_option(c("-d", "--dir"), type="character", default="NULL", 
              help="path to output dir ", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

if (is.null(opt$dir)){
  print_help(opt_parser)
  stop("At least one argument must be supplied directory", call.=FALSE)
}

inputTab <- read.table(opt$input, header = TRUE, stringsAsFactors = FALSE)
directories  <- inputTab$Directories

setwd(directories[1])
locOfNames <-NULL
locOfNames <- c(locOfNames,system('sudo -kS find . -name *bc_summary.xls*  | grep -v "link"',
                input="sat1840", intern = TRUE))
locOfNames <- sub(".", directories, locOfNames)

library(jsonlite)
library(stringr)


for(i in seq_along(directories)){
  setwd(directories[i])
  listOfVcfs <- system('sudo -kS find . -name *TSVC_variants.vcf.gz*  | grep "IonXpress"',
                       input="sat1840", intern = TRUE)
  listOfVcfs <- listOfVcfs[order(listOfVcfs)]
  tableOfNames <- read.table(file = locOfNames, stringsAsFactors = FALSE, header = TRUE,sep = "\t")
  tableOfNames <- tableOfNames[-which(is.na(tableOfNames$Mean.Depth)), ]
  tableOfNames.justNames <- NULL
  ### onyl edited this below for Andi's stuff
  #for(j in seq_along(tableOfNames$Sample.Name)){
  #  tableOfNames.justNames <- c(tableOfNames.justNames, rep(tableOfNames$Sample.Name[j],2))
  #}
  for(j in seq_along(tableOfNames$Barcode.ID)){
    tableOfNames.justNames <- c(tableOfNames.justNames, rep(tableOfNames$Barcode.ID[j],2))
  }
  tableOfNames.justNames <- sapply(tableOfNames.justNames, FUN = function(x) str_replace_all(x," ",""))
  tableOfNames.justNames
  for(k in 1:length(listOfVcfs)){
    if(k %% 2 == 0){
      a <- extensions[2]
    }
    else{
      a <- extensions[1]
    }
    print(i)
    print(k)
    #system(paste0("ln -s  ", paste0(directories[i], sub("./","",listOfVcfs[k])), " ", vcfDir,tableOfNames.justNames[k],a))
    system(paste0("cp ", paste0(directories[i], sub("./","",listOfVcfs[k])), " ", vcfDir,"TSVC_variants_",tableOfNames.justNames[k],a))
  }
}


####testing purposes below
#paste0("cp ", listOfVcfs[1], " ", vcfDir,tableOfNames$Sample.Name[1],tableOfNames.justNames[1],extensions[1])
#system(paste0("cp ", listOfVcfs[1], " ", vcfDir,tableOfNames$Sample.Name[1],tableOfNames.justNames[1],extensions[1]))

