#!/usr/bin/env Rscript
# Driver-effect heatplot (EXACT recovered slopes): covariate x class, zero-sum log-odds effect.
suppressMessages({library(data.table); library(ggplot2)})
args <- commandArgs(trailingOnly=TRUE); BR <- args[1]
CACHE <- "output/report_work/cache"
MAPDIR <- "output/report_work/maps"
R <- readRDS(file.path(CACHE, paste0(BR,".rds")))
rot <- as.data.table(R$rot)
rot[, sig := sign(q025)==sign(q975)]
# tidy covariate labels; drop focal-lag autoregressive terms to a separate group for readability
rot[, grp := fifelse(grepl("^focal_", cov), "focal LU (t-1)",
              fifelse(grepl("^OC_|^ROO|^AWC|^VS_", cov), "soil", "driver"))]
# order covariates: drivers, soil, focal; classes by overall prevalence
cls_order <- R$cats[order(-colMeans(R$obs))]
cov_order <- rot[, .(m=mean(abs(median))), by=cov][order(match(gsub("focal_|OC_.*|ROO.*|AWC.*|VS_.*","",cov), NA), -m)]$cov
rot[, cov := factor(cov, levels=rev(unique(c(rot[grp=="driver"][order(-abs(median)),unique(cov)],
                                             rot[grp=="soil",unique(cov)], rot[grp=="focal LU (t-1)",unique(cov)]))))]
rot[, class := factor(class, levels=cls_order)]
lim <- quantile(abs(rot$median), .98)
p <- ggplot(rot, aes(class, cov, fill=pmax(pmin(median,lim),-lim))) +
  geom_tile(colour="white", linewidth=.3) +
  geom_point(data=rot[sig==TRUE], shape=21, size=.45, colour="#111", fill=NA, stroke=.25) +
  scale_fill_gradient2(low="#2166AC", mid="#F7F7F5", high="#B2182B", midpoint=0, name="effect\n(log-odds)") +
  labs(x=NULL, y=NULL, title=sprintf("%s — driver effects (exact recovered slopes; dot = 95%% CI excludes 0)", BR)) +
  theme_minimal(base_size=9) +
  theme(axis.text.x=element_text(angle=45,hjust=1,size=7), axis.text.y=element_text(size=6.5),
        panel.grid=element_blank(), plot.title=element_text(size=10),
        plot.background=element_rect(fill="white",colour=NA))
nh <- max(5, nrow(rot[,.N,by=cov])*0.13)
ggsave(file.path(MAPDIR, BR, "heatplot.png"), p, width=8.5, height=nh, dpi=100, limitsize=FALSE)
cat(sprintf("[%s] driver heatplot rendered (%d covs x %d classes)\n", BR, length(unique(rot$cov)), length(unique(rot$class))))
