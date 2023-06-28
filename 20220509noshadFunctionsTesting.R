

plotForPloidy(purity, ploidy, mcnv, seg, w.fn, id, LRank)

p <- purity
P <- ploidy
tmpSeg <- seg

plotForPloidy <- function(p, P, mcnv, seg, w.fn, id, LRank) {
  # Clr <- .baseLogRatioVars(p, P, c(0,8), Msel=0, Csel=1)$Clr ## mapping lr to absolute CN
  f  <- modelFit(mcnv, seg, p, P, opts)
  pd <- plotDataKevin(mcnv$tile, mcnv$var, seg, f, gobj ,opts)
  plt  <- plotLogRatio(pd, p, P, sel.data = "segment", sel.col = "segment", dir = w.fn, sample = id)+ggtitle(sprintf("purity = %0.2f- ploidy = %0.2f", p, P))
  fldr  <-  sprintf("%s%s", w.fn, id)
  if (!file.exists(fldr)) {
    system(sprintf("mkdir %s", fldr))
  }
  ggsave(sprintf("%s/%s-%0.2f-%0.2f.png", fldr, LRank, p, P),plot = plt, device = "png", width = 500 ,height = 400, units = "mm")
}

sel.chr=NULL
max.point=1000
C.range=c(0, 8)
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
  seg$sample <- paste(sample, purity, round(ploidy, 2), sep = "_")
  write.table(seg, paste0(dir, "absoluteRes.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, append = TRUE)
  return(plt)    
}



plotDataKevin(mcnv$tile, mcnv$var, seg, f, gobj ,opts)

tile <- mcnv$tile
var <- mcnv$var
seg <- tmpSeg

function(tile, var, seg, fit, gobj, opts) {
  snp <- var
  all.chr <- seqlevels(gobj$seqi)
  tmp <- cumsum(as.numeric(seqlengths(seqinfo(seg))))
  off <- c(0, head(tmp, -1))
  names(off) <- seqnames(seqinfo(seg))
  chr.off <- data.table(chr=names(off), chr.start=off+1, chr.end=tmp, chr.col=rep(c("A","B"), length.out=length(tmp)))
  
  cov.dt <- data.table(
    chr=as.character(seqnames(tile)),
    start=floor((start(tile)+end(tile))/2),
    end=floor((start(tile)+end(tile))/2),
    type="COV",
    seg=GenomicRanges::findOverlaps(tile, seg, select = "first"),
    as.data.table(mcols(tile))
  )
  cov.dt[,chr:=factor(chr, all.chr, ordered=TRUE)]
  cov.dt[,":="(
    start.off=start+off[chr],
    end.off=end+off[chr]
  )]
  
  if (length(snp)>0) {
    snp.dt <- data.table(
      chr=as.character(seqnames(snp)),
      start=start(snp),
      end=start(snp),
      type="BAF",
      idx=(1:length(var)),
      tile=findOverlaps(snp, tile, maxgap=opts$tile.shoulder-1, select="first"),
      seg=findOverlaps(snp, seg, maxgap=opts$tile.shoulder-1, select = "first"),
      as.data.table(mcols(snp)[,!(colnames(mcols(snp)) %in% c("REF", "ALT"))])
    )
  } else {
    snp.dt <- data.table(
      chr=character(0),
      start=integer(0),
      end=integer(0),
      type=character(0),
      idx=integer(0),
      tile=integer(0),
      seg=integer(0),
      t.AF=numeric(0)
    )
  }
  snp.dt[,chr:=factor(chr, all.chr, ordered=TRUE)]
  snp.dt[,":="(
    start.off=start+off[chr],
    end.off=end+off[chr]
  )]
  
  if (is.null(fit)) {
    fit <- data.table(
      seg=1:length(seg),
      C=NA_integer_,
      K=NA_integer_,
      lr=NA_real_,
      tL=NA_real_,
      aL=NA_real_,
      d=NA_real_,
      anom=NA_real_,
      mse=NA_real_,
      nlr=NA_real_,
      naf=NA_real_,
      len=NA_real_,
      sC=NA_real_
    )
  }
  seg.dt <- data.table(
    chr=as.character(seqnames(seg)),
    start=start(seg),
    end=end(seg),
    type="SEG",
    fit
  )
  seg.dt[,chr:=factor(chr, all.chr, ordered=TRUE)]
  seg.dt[,":="(
    start.off=start+off[chr],
    end.off=end+off[chr]
  )]
  setkey(seg.dt, seg)
  cov.dt$C <- seg.dt[J(cov.dt$seg)]$C
  snp.dt$C <- seg.dt[J(snp.dt$seg)]$C
  cov.dt$K <- seg.dt[J(cov.dt$seg)]$K
  snp.dt$K <- seg.dt[J(snp.dt$seg)]$K
  pd <- list(cov=cov.dt, snp=snp.dt, seg=seg.dt, off=chr.off)
  return(pd)
}