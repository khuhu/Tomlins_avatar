load(file = "/home/kevhu/data/segment.n190.allOCPv3.Robj")
segment.n190.allOCPv3 <-out
segment.n190.allOCPv3 <- segment.n190.allOCPv3[-1,]
rm(out)

load(file = "/home/kevhu/data/segment.n320.allCCP.Robj")
segment.n320.allCCP <- out
segment.n320.allCCP <- segment.n320.allCCP[-1,]
rm(out)

load(file = "/home/kevhu/data/segment.n633.allOCP1c.Robj")
segment.n633.allOCP1c <- out
segment.n633.allOCP1c <- segment.n633.allOCP1c[-1,]
rm(out)

load(file = "/home/kevhu/data/segment.n692.allOCPv2.Robj")
segment.n692.allOCPv2 <- out
segment.n692.allOCPv2 <- segment.n692.allOCPv2[-1,]
rm(out)


library(ggplot2)
library(gridExtra)

#OCPv3
allOCPv3.mean <- ggplot(data = segment.n190.allOCPv3, aes(seg.mean)) + geom_density() + ggtitle("OCPv3")
allOCPv3.sd <- ggplot(data = segment.n190.allOCPv3, aes(seg.sd)) + geom_density() + ggtitle("OCPv3")
allOCPv3.med <- ggplot(data = segment.n190.allOCPv3, aes(seg.median)) + geom_density() + ggtitle("OCPv3")

#CCP
allCCP.mean <- ggplot(data = segment.n320.allCCP, aes(seg.mean)) + geom_density() + ggtitle("CCP")
allCCP.sd <- ggplot(data = segment.n320.allCCP, aes(seg.sd)) + geom_density() + ggtitle("CCP")
allCCP.med <- ggplot(data = segment.n320.allCCP, aes(seg.median)) + geom_density() + ggtitle("CCP")

###This set has a weird kink on the mean and median

#OCP1c
allOCP1c.mean <- ggplot(data = segment.n633.allOCP1c, aes(seg.mean)) + geom_density() + ggtitle("OCPc1")
allOCP1c.sd <- ggplot(data = segment.n633.allOCP1c, aes(seg.sd)) + geom_density() + ggtitle("OCPc1")
allOCP1c.med <- ggplot(data = segment.n633.allOCP1c, aes(seg.median)) + geom_density() + ggtitle("OCPc1")

#OCPv2
allOCPv2.mean <- ggplot(segment.n692.allOCPv2, aes(seg.mean)) + geom_density() + ggtitle("OCPv2")
allOCPv2.sd <- ggplot(segment.n692.allOCPv2, aes(seg.sd)) + geom_density() + ggtitle("OCPv2")
allOCPv2.med <- ggplot(segment.n692.allOCPv2, aes(seg.median)) + geom_density() + ggtitle("OCPv2")


grid.arrange(allOCP1c.med,allOCPv2.med,allCCP.med,allOCPv3.med)
grid.arrange(allOCP1c.mean,allOCPv2.mean,allCCP.mean,allOCPv3.mean)
grid.arrange(allOCP1c.sd,allOCPv2.sd,allCCP.sd,allOCPv3.sd)

possArtifact <- segment.n633.allOCP1c[which((segment.n633.allOCP1c$seg.mean > 2.5) & (segment.n633.allOCP1c$seg.mean < 8)),]
plot(density(possArtifact$seg.mean))

###Thing to note is subset consists of subset of IDs -> not specific genes b/c the same amount of genes are uniqely represented



####looking at kind from segment.n320.allCCP.Robj

length(unique(possArtifact$ID))
length(unique(segment.n320.allCCP$ID))

###looking across samples for locations of the weird number of copy nubmers
possArtifact$tmp <- do.call(paste0,possArtifact[,c("chrom","loc.start","loc.end")])
segment.n190.allOCPv3$tmp <- do.call(paste0, segment.n190.allOCPv3[,c("chrom","loc.start","loc.end")])
segment.n633.allOCP1c$tmp <- do.call(paste0,segment.n633.allOCP1c[,c("chrom","loc.start","loc.end")])
segment.n692.allOCPv2$tmp <- do.call(paste0, segment.n692.allOCPv2[,c("chrom","loc.start","loc.end")])



#control
length(which(possArtifact$tmp %in% segment.n633.allOCP1c$tmp))
length(which(possArtifact$tmp %in% segment.n633.allOCP1c$tmp))/length(segment.n633.allOCP1c$tmp)
#test - these are also found - now I wonder how the distribution of these areas are like
length(which(possArtifact$tmp %in% segment.n190.allOCPv3$tmp))
length(which(possArtifact$tmp %in% segment.n190.allOCPv3$tmp))/length(segment.n190.allOCPv3$tmp)
length(which(possArtifact$tmp %in% segment.n692.allOCPv2$tmp))
length(which(possArtifact$tmp %in% segment.n692.allOCPv2$tmp))/length(segment.n692.allOCPv2$tmp)



OCPv3.subset <- segment.n190.allOCPv3[which(segment.n190.allOCPv3$tmp %in% possArtifact$tmp),]
ggplot(data = OCPv3.subset, aes("seg.mean")) + geom_density()
OCPv2.subset <- segment.n692.allOCPv2[which(segment.n692.allOCPv2$tmp %in% possArtifact$tmp),]
ggplot(data = OCPv2.subset, aes("seg.mean")) + geom_density()

###creating bed file to see dig into the 



###only about 21% if the total unique CNAs are being called there
length(unique(segment.n633.allOCP1c$tmp))
length(unique(segment.n190.allOCPv3$tmp))
length(unique(segment.n692.allOCPv2$tmp))
length(unique(possArtifact$tmp))


#






###possibly just a batch effects
length(grep(".1$", possArtifact$ID))
length(grep(".2$", possArtifact$ID))
length(grep(".3$", possArtifact$ID))

###almost twice as many first time runs as second and third time runs but ~6x as many for those high CNAs
length(grep(".1$", segment.n633.allOCP1c$ID))
length(grep(".2$", segment.n633.allOCP1c$ID))
length(grep(".3$", segment.n633.allOCP1c$ID))

possArtIds <- NULL
for(i in seq_along(possArtifact$ID)){
  possArtIds[i] <- substr(possArtifact$ID[i],1,nchar(possArtifact$ID[i])-2)
}

possArtIds <- unique(possArtIds)
possArtIds.2 <- NULL

for(i in seq_along(possArtIds)) {
  a <- NULL
  a <- grep(possArtIds[i], segment.n633.allOCP1c$ID)
  possArtIds.2[i] <- list(a)
}
 
#b/c above returns all regions by panels, if they're divisble by 3, then they were run 3 times
which(possArtIds.2[[]] %% 3 == 0)

counts <- NULL
for(i in 1:length(possArtIds.2)){
  counts[i] <- (lengths(possArtIds.2[i]) %% 3 == 0)
}

###above shows only 29/76 samples have 3 runs each
###next step is to check if the ones with .1 ending have the same CNA number as their 2 and 3 counterpart if that makes sense
###for the ones with all 3 there it'll be a little more difficult to handle



