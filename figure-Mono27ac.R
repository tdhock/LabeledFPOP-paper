data(Mono27ac, package="PeakSegDisk")
library(ggplot2)
library(data.table)
l <- function(chromStart, chromEnd, annotation){
  data.table(chrom="chr11", chromStart, chromEnd, annotation)
}

max.peaks <- 13L
fit <- PeakSegOptimal::PeakSegPDPAchrom(Mono27ac$coverage, max.peaks)
seg.dt <- data.table(fit$segments)
seg.dt[, mid := (chromStart+chromEnd)/2]
seg.dt[300000 < mid & mid < 320000 & status=="peak"]

new.labels <- rbind(
  l(180000, 190000, "noPeaks"),
  l(206000, 207000, "peakStart"),
  l(208000, 210000, "peakEnd"),
  l(300000, 308250, "peakStart"),
  l(308260, 320000, "peakEnd"))
all.labels <- rbind(Mono27ac$labels, new.labels)
model.segs <- seg.dt[peaks==max.peaks]
dot.dt <- model.segs[status == "peak"]
model.color <- "blue"
ann.colors <- c(
  noPeaks="#f6f4bf",
  peakStart="#ffafaf",
  peakEnd="#ff4c4c",
  peaks="#a445ee")
err.df <- PeakError::PeakErrorChrom(dot.dt, all.labels)
peak.y <- -1
(gg <- ggplot()+
   theme_bw()+
   scale_fill_manual(values=ann.colors)+
   geom_rect(aes(
     xmin=chromStart, xmax=chromEnd,
     ymin=-Inf, ymax=Inf,
     fill=annotation),
     color="grey",
     alpha=0.5,
     data=all.labels)+
   geom_rect(aes(
     xmin=chromStart, xmax=chromEnd,
     ymin=-Inf, ymax=Inf,
     linetype=status),
     color="black",
     fill=NA,
     size=1,
     data=err.df)+
   scale_linetype_manual(
     "error type",
     values=c(
       correct=0,
       "false negative"=3,
       "false positive"=1))+
   geom_step(aes(
     chromStart, count),
     color="grey50",
     data=Mono27ac$coverage)+
   geom_step(aes(
     chromStart, mean),
     color=model.color,
     data=model.segs)+
   geom_point(aes(
     chromEnd, peak.y),
     data=dot.dt,
     shape=1,
     color=model.color)+
   geom_segment(aes(
     chromStart, peak.y,
     xend=chromEnd, yend=peak.y),
     data=dot.dt,
     color=model.color,
     size=2))

(gg.zoom <- gg+coord_cartesian(xlim=c(200000, 300000)))
png("figure-Mono27ac-zoom.png", 10, 3, units="in", res=200)
print(gg.zoom)
dev.off()
png("figure-Mono27ac.png", 10, 3, units="in", res=200)
print(gg)
dev.off()
