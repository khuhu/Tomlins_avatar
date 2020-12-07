#!/usr/bin/env Rscript

sampleVarCalls = read.table("/home/kevhu/data/allDataTest.txt",
                            header = TRUE,sep = '\t', quote = "", stringsAsFactors = FALSE)


test = apply(X = sampleVarCalls[,c("SAMPLE","BARCODE","BED","REPORTNUM","RunName")], MARGIN = 1, FUN = paste, collapse = "")
for(i in seq_along(sampleVarCalls$SAMPLE)){
  ##finding positional percentage
  posPct = sampleVarCalls[which((sampleVarCalls$Chr == sampleVarCalls$Chr[i]) & (sampleVarCalls$Start_position == sampleVarCalls$Start_position[i]) & (sampleVarCalls$End_position == sampleVarCalls$End_position[i]) & (sampleVarCalls$REF_1 == sampleVarCalls$REF_1[i])),]
  sampleVarCalls$new_TOT_OBS[i] = nrow(posPct)
  sampleVarCalls$new_TOT_NS[i] = length(unique(test))
  sampleVarCalls$new_TOT_POS_PCT[i] =  nrow(posPct)/length(unique(test))
  ##finding position variance
  varPCT = sampleVarCalls[which((sampleVarCalls$Chr == sampleVarCalls$Chr[i]) & (sampleVarCalls$Start_position == sampleVarCalls$Start_position[i]) & (sampleVarCalls$End_position == sampleVarCalls$End_position[i]) & (sampleVarCalls$REF_1 == sampleVarCalls$REF_1[i]) & (sampleVarCalls$ALT_1 == sampleVarCalls$ALT_1[i])),]
  sampleVarCalls$new_TOT_VAR_PCT[i] = nrow(varPCT)/length(unique(test))
}

write.table(x = sampleVarCalls ,file = "/home/kevhu/data/20170711_allData.tsv", sep = '\t', row.names = FALSE)
