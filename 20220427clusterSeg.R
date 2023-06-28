### per chromosome 
### do hclust - and from cutree can combine segments using weighted mean


library("PSCBS")
library(copynumber)
library(foreach)
library(doParallel)

clustSeg <- function(df, distVar =  0.2){
  res <- NULL
  for(i in unique(df$chromosome)){
    print(i)
    chrDf <- df[which(df$chromosome == i), ]
    treeGroups <- tryCatch(cutree(hclust(dist(chrDf$mean)), h = distVar),
                           error = function(x) return(NULL))
    if (is.null(treeGroups)) {
      res <- rbind(res, chrDf)
      next()
    }
    for (j in unique(treeGroups)) {
      print(j)
      tmpMeans <- chrDf$mean[which(treeGroups == j)]
      tmpWeights <- chrDf$nbrOfLoci[which(treeGroups == j)]/sum(chrDf$nbrOfLoci[which(treeGroups == j)])
      chrDf$mean[which(treeGroups == j)] <- weighted.mean(tmpMeans, tmpWeights)
    }
    res <- rbind(res, chrDf)
  }
  return(res)
}



# clustSeg(tmpSegRes, h = 0.15)
# 
# h <- 0.15
# tmpSegRes
# df <- tmpSegRes
# distVar = 0.2

clusterSegG <- function(df, distVar = 0.2){
  
  treeGroups <- tryCatch(cutree(hclust(dist(df$mean)), h = distVar),
                         error = function(x) return(NULL))
  if (is.null(treeGroups)) {
    next()
  }
  
  
  for (j in unique(treeGroups)) {
    tmpMeans <- df$mean[which(treeGroups == j)]
    tmpWeights <- df$nbrOfLoci[which(treeGroups == j)]/sum(df$nbrOfLoci[which(treeGroups == j)])
    df$mean[which(treeGroups == j)] <- weighted.mean(tmpMeans, tmpWeights)
  }
  
  return(df)
}


# tmpSegRes
# df <- tmpSegRes
# distVar = 0.2
# distVar2 = 0.5

clusterSegG_twostep <- function(df, distVar = 0.2, distVar2 = 0.5){
  
  df$mean <- 2^df$mean
  
  res <- NULL
  for(i in unique(df$chromosome)){
    print(i)
    chrDf <- df[which(df$chromosome == i), ]
    treeGroups <- tryCatch(cutree(hclust(dist(chrDf$mean)), h = distVar),
                           error = function(x) return(NULL))
    if (is.null(treeGroups)) {
      res <- rbind(res, chrDf)
      next()
    }
    for (j in unique(treeGroups)) {
      print(j)
      tmpMeans <- chrDf$mean[which(treeGroups == j)]
      tmpWeights <- chrDf$nbrOfLoci[which(treeGroups == j)]/sum(chrDf$nbrOfLoci[which(treeGroups == j)])
      chrDf$mean[which(treeGroups == j)] <- weighted.mean(tmpMeans, tmpWeights)
    }
    res <- rbind(res, chrDf)
  }
  
  ### does by chromosome then global
  
  df <- res
  
  treeGroups <- tryCatch(cutree(hclust(dist(df$mean)), h = distVar2),
                         error = function(x) return(NULL))
  if (is.null(treeGroups)) {
    next()
  }
  
  
  for (j in unique(treeGroups)) {
    tmpMeans <- df$mean[which(treeGroups == j)]
    tmpWeights <- df$nbrOfLoci[which(treeGroups == j)]/sum(df$nbrOfLoci[which(treeGroups == j)])
    df$mean[which(treeGroups == j)] <- weighted.mean(tmpMeans, tmpWeights)
  }
  
  df$mean <- log2(df$mean)
  
  return(df)
}



### altnerative plotLogRatio from noshad pipeline
### used to retrieve specific plody and purity segs

plotForPloidy <- function(p, P, mcnv, seg, w.fn, id, LRank) {
  # Clr <- .baseLogRatioVars(p, P, c(0,8), Msel=0, Csel=1)$Clr ## mapping lr to absolute CN
  f  <- modelFit(mcnv, seg, p, P, opts)
  pd <- plotDataKevin(mcnv$tile, mcnv$var, seg, f, gobj ,opts)
  pd$off
  plt  <- plotLogRatio(pd, p, P, sel.data = "segment", sel.col = "segment", dir = w.fn, sample = id)+ggtitle(sprintf("purity = %0.2f- ploidy = %0.2f", p, P))
  fldr  <-  sprintf("%s%s", w.fn, id)
  if (!file.exists(fldr)) {
    system(sprintf("mkdir %s", fldr))
  }
  ggsave(sprintf("%s/%s-%0.2f-%0.2f.png", fldr, LRank, p, P),plot = plt, device = "png", width = 500 ,height = 400, units = "mm")
}


plotLogRatio <- function(pd, purity=NULL, ploidy=NULL, sel.chr=NULL, sel.data="tile", sel.col="segment", lr.range=c(-4,4), C.range=c(0, 8), max.point=1000, dir = NULL, sample = NULL) {
  ## select data from chromosomes
  if (!is.null(sel.chr)) {
    cov <- pd$cov[chr %in% sel.chr]
    off <- pd$off[chr %in% sel.chr]
    seg <- pd$seg[chr %in% sel.chr]
  } else {
    cov <- pd$cov
    off <- pd$off
    seg <- pd$seg
  }
  ## subsample points
  cov <- cov[sample(nrow(cov))]
  cov <- cov[,head(.SD,max.point),by=seg]
  ## compute plot limits
  ymin <- lr.range[1]
  ymax <- lr.range[2]
  if (is.null(purity) || is.null(ploidy)) {
    Clr <- data.table(C=seq(lr.range[1], lr.range[2]), lr=seq(lr.range[1], lr.range[2]))       
  } else {
    vars <- .baseLogRatioVars(purity, ploidy, C.range)
    Clr <- vars$Clr
  }
  ## base plot (no data)
  plt <- .baseLogRatioPlot(off, Clr, ymin, ymax)
  ## data is tile
  if (sel.data=="tile") {
    if (sel.col=="segment") {
      plt <- plt +
        geom_point(aes(x=start.off, y=lr, col=factor(seg %% 3)), cov[(target)], size=0.75, alpha=1) +
        scale_color_manual(values=c("#6495ED", "#DD8080", "#CDE2B8"), guide=FALSE) +
        geom_point(aes(x=start.off, y=lr), cov[(!target)], size=0.75, alpha=1)
    } else if (sel.col=="sC") {
      plt <- plt +
        geom_point(aes(x=start.off, y=lr, col=sC), cov[(target)], size=0.75, alpha=1) +
        scale_color_gradient2(low="blue", mid="black", high="red", midpoint=2, limits=c(0,6), oob = scales::squish) +
        geom_point(aes(x=start.off, y=lr), cov[(!target)], size=0.75, alpha=1)
    } else if (sel.col=="C") {
      plt <- plt +
        geom_point(aes(x=start.off, y=lr, color=as.character(C)), cov[(target)], size=0.75, alpha=1) +
        scale_color_manual(values=STRING_COL, guide=FALSE) +
        geom_point(aes(x=start.off, y=lr), cov[(!target)], size=0.75, alpha=1)
    } else if (sel.col=="CK") {
      cov[C>1 & K==0, C:=-1]
      cov[,C:=as.character(C)]
      plt <- plt +
        geom_point(aes(x=start.off, y=lr, color=C, text=gene_names), cov, size=0.75, alpha=1) +
        scale_color_manual(values=STRING_COL, guide=FALSE)
    } else {
      plt <- plt + 
        geom_point(aes_string(x="start.off", y="lr", col=sel.col), cov, size=0.75) +
        scale_color_gradient(low="black", high="red", guide=FALSE)
    }
  }
  ## data is segment
  if (sel.data=="segment") {
    if (sel.col=="segment") {
      plt <- plt + 
        # geom_segment(aes(x=start.off, xend=end.off, y=lr, yend=lr, col=factor(seg %% 3)), seg, size=2) + 
        # # scale_color_manual(values=c("#6495ED", "#DD8080", "#CDE2B8"), guide=FALSE)
        geom_segment(aes(x=start.off, xend=end.off, y=lr, yend=lr), seg, size=2) + 
        scale_color_manual(values=c("#6495ED"), guide=FALSE)
    } else {
      plt <- plt +
        geom_segment(aes_string(x="start.off", xend="end.off", y="lr", yend="lr", col=lr-Clr[C==2]$lr), seg, size=2) +
        scale_color_gradient(low="black", high="red", guide=FALSE)
    }
  }
  # seg$sample <- paste(sample, purity, round(ploidy, 2), sep = "_")
  seg$sample <- paste(sample, sprintf("%0.2f", purity), sprintf("%0.2f", ploidy), sep = "_")
  write.table(seg, paste0(dir, "absoluteRes.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, append = TRUE)
  return(plt)    
}

### comment this out for lower or higher limit of tumor content

.opt.models <- function(data, seg, opts, purity) {
  purity <- round(purity, 2)
  opts$opt.p.lo       <- 0.4
  opts$opt.p.hi       <- 0.99
  opts$opt.grid.p.res <- 0.01
  grid <- .opt.grid(data, opts)
  grid$p0 <- round(grid$p0, 2)
  grid <- grid[p0 == purity]
  fine <- .opt.fine(data, grid, opts)
  # fine <- fine[p0 == purity]
  
  return(fine)
}


nameStripper <- function(df){
  require(stringr)
  df <- str_remove(df, "_MG_X.*")
  df <- str_remove(str_remove(df, "^X"), "_X.*")
  df <- tolower(str_remove_all(str_remove_all(str_remove(df, "X.*"), "_"), "\\."))
  df  <- str_remove(df, "o")
  df  <- str_replace_all(df, " ", "")
  return(df)
}



calcZscore <- function(segRes){
  sampZ1 <- NULL
  sampZ2 <- NULL
  z1_vector <- NULL
  z2_vector <- NULL
  p_vector1 <- NULL
  p_vector2 <- NULL
  nmean_vector <- NULL
  tmean_vector <- NULL
  q_vector1 <- NULL
  q_vector2 <- NULL
  
  # uniqIdx <- unique(match(str_replace_all(segRes$ID, "[[:punct:]]", ""), 
  #                         str_replace_all(str_replace_all(colnames(zscore_gc_oe_ratios), "[[:punct:]]", ""), " ", "")))
  
  uniqIdx <- unique(match(segRes$ID, colnames(zscore_gc_oe_ratios)))
  zScoreColnames <- colnames(zscore_gc_oe_ratios)[uniqIdx]
  for (i in seq_along(unique(segRes$ID))) {
    tmp_seg_res <- segRes[which(segRes$ID == unique(segRes$ID)[i]), ]
    tmp_seg_ranges <- GRanges(seqnames = Rle(tmp_seg_res$chrom), 
                              ranges = IRanges(start = tmp_seg_res$loc.start,
                                               end = tmp_seg_res$loc.end))
    print(unique(segRes$ID)[i])
    print(zScoreColnames[i])
    
    sampZ1 <- NULL
    sampZ2 <- NULL
    
    for (j in seq_along(tmp_seg_ranges)) {
      tmpOverlap <- findOverlaps(all_probes_grange, tmp_seg_ranges[j])
      tmpTest <- zscore_gc_oe_ratios[[zScoreColnames[i]]][queryHits(tmpOverlap)]
      
      if (!is.numeric(tmpTest)) {
        sampZ1 <- c(sampZ1, NA)
        sampZ2 <- c(sampZ2, NA)
        nmean_vector <- c(nmean_vector, NA)
        tmean_vector <- c(tmean_vector, NA)
        next()
      } else{
        tmpNormal <- zscore_gc_oe_ratios[queryHits(tmpOverlap), mouseNormal]
        tmpNormal2 <- unlist(apply(tmpNormal, 2, mean))
        
        z1_stat <- (mean(tmpTest) - mean(tmpNormal2))/sd(tmpNormal2)
        z2_stat <- (mean(tmpTest) - mean(tmpNormal2))/(sd(tmpNormal2) * mean(tmpTest))
        pVal_z1 <- 2.0*pnorm(-abs(z1_stat))
        pVal_z2 <- 2.0*pnorm(-abs(z2_stat))
        
        
        normal_seg_mean <- mean(tmpNormal2)
        tumor_seg_mean <- mean(tmpTest)
        sampZ1 <- c(sampZ1, pVal_z1)
        sampZ2 <- c(sampZ2, pVal_z2)
        z1_vector <- c(z1_vector, z1_stat)
        z2_vector <- c(z2_vector, z2_stat)
        p_vector1 <- c(p_vector1, pVal_z1)
        p_vector2 <- c(p_vector2, pVal_z2)
        nmean_vector <- c(nmean_vector, normal_seg_mean)
        tmean_vector <- c(tmean_vector, tumor_seg_mean)
      }
    }
    qVal_z1 <- (p.adjust(sampZ1, method = "BH"))
    qVal_z2 <- (p.adjust(sampZ2, method = "BH"))
    q_vector1 <- c(q_vector1, qVal_z1)
    q_vector2 <- c(q_vector2, qVal_z2)
    
  }
  return(list("z1_vector" = z1_vector, "z2_vector" = z2_vector,
              "p1_vector" = p_vector1, "p2_vector" = p_vector2,
              "q1_vector" = q_vector1, "q2_vector" = q_vector2,
              "nmean_vector"= nmean_vector, "tmean_vector" = tmean_vector))
}


calcZscore_log2 <- function(segRes){
  sampZ1 <- NULL
  sampZ2 <- NULL
  z1_vector <- NULL
  z2_vector <- NULL
  p_vector1 <- NULL
  p_vector2 <- NULL
  nmean_vector <- NULL
  tmean_vector <- NULL
  q_vector1 <- NULL
  q_vector2 <- NULL
  
  # uniqIdx <- unique(match(str_replace_all(segRes$ID, "[[:punct:]]", ""), 
  #                         str_replace_all(str_replace_all(colnames(zscore_gc_oe_ratios), "[[:punct:]]", ""), " ", "")))
  
  uniqIdx <- unique(match(segRes$ID, colnames(zscore_gc_oe_ratios)))
  zScoreColnames <- colnames(zscore_gc_oe_ratios)[uniqIdx]
  for (i in seq_along(unique(segRes$ID))) {
    tmp_seg_res <- segRes[which(segRes$ID == unique(segRes$ID)[i]), ]
    tmp_seg_ranges <- GRanges(seqnames = Rle(tmp_seg_res$chrom), 
                              ranges = IRanges(start = tmp_seg_res$loc.start,
                                               end = tmp_seg_res$loc.end))
    print(unique(segRes$ID)[i])
    print(zScoreColnames[i])
    
    sampZ1 <- NULL
    sampZ2 <- NULL
    
    for (j in seq_along(tmp_seg_ranges)) {
      tmpOverlap <- findOverlaps(all_probes_grange, tmp_seg_ranges[j])
      tmpTest <- zscore_gc_oe_ratios[[zScoreColnames[i]]][queryHits(tmpOverlap)]
      
      if (!is.numeric(tmpTest)) {
        sampZ1 <- c(sampZ1, NA)
        sampZ2 <- c(sampZ2, NA)
        nmean_vector <- c(nmean_vector, NA)
        tmean_vector <- c(tmean_vector, NA)
        next()
      } else{
        tmpNormal <- zscore_gc_oe_ratios[queryHits(tmpOverlap), mouseNormal]
        tmpNormal2 <- unlist(apply(tmpNormal, 2, mean))
        
        tmpTest <- log2(tmpTest)
        tmpNormal2 <- log2(tmpNormal2)
        
        z1_stat <- (mean(tmpTest) - mean(tmpNormal2))/sd(tmpNormal2)
        z2_stat <- (mean(tmpTest) - mean(tmpNormal2))/(sd(tmpNormal2) * mean(tmpTest))
        pVal_z1 <- 2.0*pnorm(-abs(z1_stat))
        pVal_z2 <- 2.0*pnorm(-abs(z2_stat))
        
        
        normal_seg_mean <- mean(tmpNormal2)
        tumor_seg_mean <- mean(tmpTest)
        sampZ1 <- c(sampZ1, pVal_z1)
        sampZ2 <- c(sampZ2, pVal_z2)
        z1_vector <- c(z1_vector, z1_stat)
        z2_vector <- c(z2_vector, z2_stat)
        p_vector1 <- c(p_vector1, pVal_z1)
        p_vector2 <- c(p_vector2, pVal_z2)
        nmean_vector <- c(nmean_vector, normal_seg_mean)
        tmean_vector <- c(tmean_vector, tumor_seg_mean)
      }
    }
    qVal_z1 <- (p.adjust(sampZ1, method = "BH"))
    qVal_z2 <- (p.adjust(sampZ2, method = "BH"))
    q_vector1 <- c(q_vector1, qVal_z1)
    q_vector2 <- c(q_vector2, qVal_z2)
    
  }
  return(list("z1_vector" = z1_vector, "z2_vector" = z2_vector,
              "p1_vector" = p_vector1, "p2_vector" = p_vector2,
              "q1_vector" = q_vector1, "q2_vector" = q_vector2,
              "nmean_vector"= nmean_vector, "tmean_vector" = tmean_vector))
}


segZscoresFilt <- function(segResults, zscores){
  
  # filts for significance, cn-change and length
  segmental_zscores <- cbind(segResults, "z-scores1" = zscores[["z1_vector"]], "z-scores2" = zscores[["z2_vector"]],
                             "p-val1" = zscores[["p1_vector"]], "p-val2" = zscores[["p2_vector"]],
                             "q-val1" = zscores[["q1_vector"]], "q-val2" = zscores[["q2_vector"]],
                             "normal_seg_mean" = zscores[["nmean_vector"]], "tumor_seg_mean" = zscores[["tmean_vector"]])
  return(segmental_zscores)
}
