# Coordinate-based, resolution-agnostic, per-slice focal (queen-8), with a DEFINED duplicate rule
# (aggregate sub-units per cell -> representative composition -> focal -> broadcast). Replaces terra's
# arbitrary xyz-collapse. Shared routine for fit AND predict.
suppressMessages(library(data.table))

compute_focal_coord <- function(dt, classes, coord=c("X","Y"), slice="out_year",
                                res=NULL, complete_only=TRUE, agg=c("mean","wmean"), wcol=NULL) {
  agg<-match.arg(agg); X<-coord[1]; Y<-coord[2]; out<-list()
  for (sl in unique(dt[[slice]])) {
    d<-dt[get(slice)==sl]
    r <- if (is.null(res)) min(diff(sort(unique(d[[X]])))) else res      # infer grid spacing
    # 1. representative per-cell composition (aggregate duplicate sub-units)
    if (agg=="wmean" && !is.null(wcol)) cell<-d[,lapply(.SD,function(v) sum(v*get(wcol))/sum(get(wcol))),by=c(X,Y),.SDcols=classes]
    else                               cell<-d[,lapply(.SD,mean),by=c(X,Y),.SDcols=classes]
    setnames(cell,c(X,Y),c("cx","cy"))
    # 2. queen-8 neighbour lookups
    off<-CJ(dx=c(-r,0,r),dy=c(-r,0,r))[!(dx==0&dy==0)]
    nb<-cell[,.(nx=cx+off$dx, ny=cy+off$dy, ox=cx, oy=cy), by=.(cx,cy)]
    setkey(cell,cx,cy); j<-cell[nb, on=.(cx=nx, cy=ny)]                  # neighbour shares at (nx,ny); ox,oy = focal cell
    agg_j<-j[, c(list(npresent=sum(!is.na(get(classes[1])))),
                 lapply(.SD,function(v) sum(v,na.rm=TRUE))), by=.(ox,oy), .SDcols=classes]
    for (cc in classes) set(agg_j, j=paste0("focal_",cc),
        value = if (complete_only) fifelse(agg_j$npresent==8L, agg_j[[cc]]/8, 0) else agg_j[[cc]]/pmax(agg_j$npresent,1))
    setnames(agg_j,c("ox","oy"),c(X,Y))
    out[[as.character(sl)]]<-agg_j[,c(coord,paste0("focal_",classes)),with=FALSE][,(slice):=sl]
  }
  rbindlist(out,fill=TRUE)
}

if (sys.nframe()==0) {
  dp<-as.data.table(readRDS("output/dat_pixel_FULL_10km_2026-06-17.rds"))
  foc<-grep("^focal_",names(dp),value=TRUE); own<-intersect(sub("^focal_","",foc),names(dp)); foc<-paste0("focal_",own)
  dp[,ncell:=.N,by=.(X,Y,out_year)]
  rec<-compute_focal_coord(dp, classes=own, coord=c("X","Y"), slice="out_year", res=10000, complete_only=TRUE)
  m<-merge(dp[,c("X","Y","out_year","ncell",foc),with=FALSE], rec, by=c("X","Y","out_year"), suffixes=c(".t",".c"))
  D<-abs(as.matrix(m[,paste0(foc,".t"),with=FALSE]) - as.matrix(m[,paste0(foc,".c"),with=FALSE]))
  cat(sprintf("\nEQUIVALENCE coord vs terra (%d rows, %d classes):\n", nrow(m), length(own)))
  cat(sprintf("  ALL   : mean|d|=%.4f max|d|=%.3f\n", mean(D), max(D)))
  cl<-m$ncell==1
  cat(sprintf("  CLEAN (ncell==1, %d rows): mean|d|=%.4f max|d|=%.3f  <- tests window/edge equivalence\n", sum(cl), mean(D[cl,]), max(D[cl,])))
  cat(sprintf("  DUP   (ncell>1,  %d rows): mean|d|=%.4f max|d|=%.3f  <- terra arbitrary-collapse vs coord aggregate\n", sum(!cl), mean(D[!cl,]), max(D[!cl,])))
}
