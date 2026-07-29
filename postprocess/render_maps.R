#!/usr/bin/env Rscript
# Render observed|predicted share grid maps (EPSG:3035) for every class of one branch.
suppressMessages({library(data.table); library(ggplot2); library(scales)})
args <- commandArgs(trailingOnly=TRUE); BR <- args[1]
CACHE <- "output/report_work/cache"
MAPDIR <- "output/report_work/maps"
dir.create(file.path(MAPDIR,BR), showWarnings=FALSE, recursive=TRUE)
R <- readRDS(file.path(CACHE, paste0(BR,".rds")))
cats <- R$cats; J <- R$J
base <- data.table(x=R$coordX, y=R$coordY)
# sequential single-hue magnitude ramp (perceptually uniform, CVD-safe): light -> dark teal
seqcols <- c("#F2FbF8","#CDEBE3","#98D6C8","#5FBBA8","#2E9A86","#0F7A67","#075447")
theme_map <- theme_void(base_size=11) + theme(
  plot.title=element_text(size=11,hjust=.5,margin=margin(b=3)),
  legend.position="right", legend.key.height=unit(14,"pt"), legend.key.width=unit(8,"pt"),
  plot.background=element_rect(fill="white",colour=NA), panel.spacing=unit(4,"pt"),
  strip.text=element_text(size=10,face="bold"))
for (j in seq_len(J)) {
  cl <- cats[j]
  d <- rbind(
    data.table(base, val=R$obs[,j], panel="Observed"),
    data.table(base, val=R$pred[,j], panel="Predicted (model)")
  )
  mx <- max(d$val, na.rm=TRUE); if (!is.finite(mx) || mx<=0) mx <- 1e-6
  p <- ggplot(d, aes(x,y,fill=val)) + geom_tile(width=10000,height=10000) +
    facet_wrap(~panel, ncol=2) +
    scale_fill_gradientn(colours=seqcols, limits=c(0,mx), labels=percent_format(accuracy=1),
                         name="share", oob=squish) +
    coord_equal(expand=FALSE) + theme_map +
    ggtitle(sprintf("%s  —  %s", BR, cl))
  ggsave(file.path(MAPDIR, BR, sprintf("%02d_%s.png", j, gsub("[^A-Za-z0-9]+","_",cl))),
         p, width=7.2, height=3.5, dpi=96)
}
# dominant predicted class map (categorical) as an overview
dom_obs <- cats[max.col(R$obs, ties.method="first")]
dom_pred<- cats[max.col(R$pred, ties.method="first")]
dd <- rbind(data.table(base, cl=dom_obs, panel="Observed dominant"),
            data.table(base, cl=dom_pred, panel="Predicted dominant"))
pd <- ggplot(dd, aes(x,y,fill=cl)) + geom_tile(width=10000,height=10000) +
  facet_wrap(~panel, ncol=2) + coord_equal(expand=FALSE) + theme_map +
  theme(legend.position="bottom", legend.text=element_text(size=6), legend.key.size=unit(7,"pt")) +
  guides(fill=guide_legend(ncol=6,title=NULL)) + ggtitle(sprintf("%s — dominant class", BR))
ggsave(file.path(MAPDIR, BR, "00_dominant.png"), pd, width=7.6, height=4.6, dpi=96)
cat(sprintf("[%s] rendered %d class maps + dominant overview\n", BR, J))
