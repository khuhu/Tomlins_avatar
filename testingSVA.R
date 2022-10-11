Sys.unsetenv("R_LIBS_USER")
Sys.setenv(R_LIBS_USER="/home/kevhu/R/x86_64-pc-linux-gnu-library/3.2/")

assign(".lib.loc", "/home/kevhu/R/x86_64-pc-linux-gnu-library/3.2/", envir = environment(.libPaths))


.libPaths()[1] 
source("https://bioconductor.org/biocLite.R")
biocLite("sva")
biocLite("pamr")
biocLite("limma")
biocLite("bladderbatch")

library(sva)
library(limma)
library(pamr)
library(bladderbatch)

data(bladderdata)
pheno = pData(bladderEset)
edata = exprs(bladderEset)
mod = model.matrix(~as.factor(cancer), data=pheno)
mod0 = model.matrix(~1,data=pheno)
n.sv = num.sv(edata,mod,method="leek")
n.sv
svobj = sva(edata,mod,mod0,n.sv=n.sv)
pValues = f.pvalue(edata,mod,mod0)
qValues = p.adjust(pValues,method="BH")

modSv = cbind(mod,svobj$sv)
mod0Sv = cbind(mod0,svobj$sv)
pValuesSv = f.pvalue(edata,modSv,mod0Sv)
qValuesSv = p.adjust(pValuesSv,method="BH")

fit = lmFit(edata,modSv)
contrast.matrix <- cbind("C1"=c(-1,1,0,rep(0,svobj$n.sv)),"C2"=c(0,-1,1,rep(0,svobj$n.sv)),"C3"=c(-1,0,1,rep(0,svobj$n.sv)))
fitContrasts = contrasts.fit(fit,contrast.matrix)

eb = eBayes(fitContrasts)
topTableF(eb, adjust="BH")


### testing for sva
dfDummy <- df;

for(x in allNames){
  dfDummy[[x]] <- dfResidRatio[[x]]/dfCorrectionRatio[[x]]
}

#testString <- "MYCLPool2"
#substr(testString, nchar(testString), nchar(testString))

for(i in seq_along(dfDummy$Gene)){
  dfDummy$PoolNumber[i] <- substr(dfDummy$Gene[i], nchar(dfDummy$Gene[i]), nchar(dfDummy$Gene[i]))
}
dfDummy$PoolNumber <- as.factor(dfDummy$PoolNumber)


###reshaping

#library(reshape2)
#dfDummy.amplicon <- melt(data = dfDummy, value.name = AmpliconId)


CNAmod <- model.matrix(~as.factor(PoolNumber), data = dfDummy)
CNAmod0 <- model.matrix(~1, data = dfDummy)

###quick fix for weird error
#tmplist <- NULL
#for(i in 1:nrow(CNAmod)){
#  tmplist[i] <- if(CNAmod[,2][i] == 1) 0 else 1
#}
#CNAmod <- cbind(CNAmod, tmplist)


dfDummy.values <- t(data.matrix(dfDummy[,allNames]))

CNA.num.sv <- num.sv(dfDummy.values, CNAmod, method = "leek")  
CNA.num.sv

str(edata)
str(mod)
str(CNAmod)
str(dfDummy.values)

### skipping enumerating the correct svs

CNA.svobj <- sva(dfDummy.values, mod = CNAmod, mod0 = CNAmod0)

CNA.svobj$n.sv
CNA.svobj$pprob.gam
CNA.svobj$pprob.b
