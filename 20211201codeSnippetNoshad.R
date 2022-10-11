res <- foreach(i = 1:length(samples), .combine = rbind) %dopar% {
  print(i)
  lr     <- lrs[variable == samples[i], ]
  
  if (isEmpty(lr$variable)) {
    res <- NULL
    return(res)
  }
  
  seg    <- segs[ID == samples[i]]
  purity <- Pp[sample == samples[i], tc]
  sex    <- Pp$sex[i]
  ## make tiled genome object
  cnv      <- list()
  cnv$var  <- GRanges()
  cnv$tile <- GenomicRanges::GRanges(seqnames = lr$ChromNum, ranges = IRanges(start = lr$Start, end = lr$End),  AmpliconID = lr$AmpliconId, lr = lr$value)
  seg      <- GenomicRanges::GRanges(seqnames = seg$chrom,   ranges = IRanges(start = seg$loc.start, end = seg$loc.end),  AmpliconID = seg$ID, 
                                     lr.mean = seg$seg.mean, nlr = seg$num.mark)
  seqinfo(seg) <- gobj$seqi
  ## model cnv
  mcnv             <- cnv
  mcnv$tile$arm    <- seqnames(cnv$tile)
  mcnv$tile$target <- TRUE
  mcnv$tile$hq     <- TRUE
  mcnv$tile$gap    <- 0
  mcnv$stats <- .modelStats(mcnv$tile, mcnv$var)
  ## normal copy number
  if( sex == "male") {
    mcnv$tile$nC  <- ifelse(seqnames(mcnv$tile) %in% "chrX", 1, 2)
  } else {
    mcnv$tile$nC  <- 2
  }
  mcnv$tile$baf       <- NA_real_
  mcnv$tile$baf.depth <- NA_real_
  mcnv$tile$baf.n     <- NA_real_
  
  data <- .opt.data(mcnv, seg, opts)
  opt <- tryCatch(.opt.models(data, seg, opts, purity), 
                  error = function(e) NULL)
  
  if (is.null(opt)) {
    res <- NULL
    return(res)
  }
  
  # this is where the final likelihoods are so after sort, append to df
  opt <- opt[order(-L)]
  opt$Lrank <- 1:nrow(opt)
  
  lapply(1:nrow(opt), function(j) {
    purity <- opt$p0[j]
    ploidy <- opt$P0[j]
    LRank  <- str_pad(opt$Lrank[j], 2, pad = 0)
    w.fn   <- sprintf("/mnt/DATA5/tmp/kev/testNoshad3/")
    id     <- samples[i]
    plotForPloidy(purity, ploidy, mcnv, seg, w.fn, id, LRank)
  })
  
  tmpOpt <- cbind("Sample" = rep(samples[i], nrow(opt)), opt)
  return(tmpOpt)
}