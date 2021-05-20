works_with_R(
  "4.1.0",
  PeakSegOptimal="2018.5.25",
  PeakError="2017.6.19",
  "tdhock/FLOPART@ee8d599891f665efe356698babafa0c1d7cca0ba")

data(Mono27ac, package="FLOPART")
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
  l(100000, 200000, "noPeaks"),
  l(206000, 207000, "peakStart"),
  l(208000, 220000, "peakEnd"),
  l(300000, 308250, "peakStart"),
  l(308260, 320000, "peakEnd"))
all.labels <- rbind(Mono27ac$labels, new.labels)[order(chromStart)]
flopart <- FLOPART::FLOPART(Mono27ac[["coverage"]], all.labels, 1400)
FLOPART.segs <- flopart[["segments_dt"]]
FLOPART.peaks <- FLOPART.segs[status == "peak"]
model.color <- "blue"
peak.y <- -1
ann.colors <- c(
  noPeaks="#f6f4bf",
  peakStart="#ffafaf",
  peakEnd="#ff4c4c")
gg <- ggplot()+
  geom_rect(aes(
    xmin=chromStart, xmax=chromEnd,
    fill=annotation,
    ymin=-Inf, ymax=Inf),
    data=all.labels,
    color="grey",
    alpha=0.5)+
  theme_bw()+
  scale_fill_manual(values=ann.colors)+
  geom_step(aes(
    chromStart, count),
    color="grey50",
    data=Mono27ac[["coverage"]])+
  geom_step(aes(
    chromStart, mean),
    color=model.color,
    data=FLOPART.segs)+
  geom_point(aes(
    chromEnd, peak.y),
    color=model.color,
    shape=21,
    data=FLOPART.peaks)
gg+coord_cartesian(xlim=c(500000, 510000))

gg+coord_cartesian(xlim=c(200000, 250000))

gg+coord_cartesian(xlim=c(150000, 200000))

gg+coord_cartesian(xlim=c(300000, 330000))

seg.dt.list <- list(
  PeakSegOptimal=seg.dt[peaks==max.peaks],
  FLOPART=FLOPART.segs)
err.dt.list <- list()
model.segs.list <- list()
for(pkg in names(seg.dt.list)){
  pkg.segs <- seg.dt.list[[pkg]][, .(chromStart, chromEnd, mean, status)]
  pkg.peaks <- pkg.segs[status=="peak"]
  err.df <- PeakError::PeakErrorChrom(pkg.peaks, all.labels)
  err.dt.list[[pkg]] <- data.table(pkg, err.df)
  model.segs.list[[pkg]] <- data.table(pkg, pkg.segs)
}
err.dt <- do.call(rbind, err.dt.list)
model.segs <- do.call(rbind, model.segs.list)

model.peaks <- model.segs[status=="peak"]
(gg <- ggplot()+
   theme_bw()+
   theme(panel.spacing=grid::unit(0, "lines"))+
   scale_fill_manual(values=ann.colors)+
   facet_grid(pkg ~ ., labeller=label_both)+
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
     data=err.dt)+
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
     data=model.peaks,
     shape=1,
     color=model.color)+
   geom_segment(aes(
     chromStart, peak.y,
     xend=chromEnd, yend=peak.y),
     data=model.peaks,
     color=model.color,
     size=2))

gg.list <- list(
  out=gg,
  peakEnd=gg+coord_cartesian(xlim=c(200000, 250000)),
  noPeaks=gg+coord_cartesian(xlim=c(150000, 200000), ylim=c(-1, 5)))
for(g.name in names(gg.list)){
  gg.zoom <- gg.list[[g.name]]
  out.png <- sprintf("figure-Mono27ac-new-labels-%s.png", g.name)
  png(out.png, 10, 3, units="in", res=200)
  print(gg.zoom)
  dev.off()
}
