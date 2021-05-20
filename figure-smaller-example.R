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

lo.lim <- 145000
hi.lim <- 327000
some.coverage <- Mono27ac$coverage[lo.lim < chromEnd & chromEnd < hi.lim]
some.coverage[1, chromStart := lo.lim]
some.coverage[.N, chromEnd := hi.lim]

max.peaks <- 7L
fit <- PeakSegOptimal::PeakSegPDPAchrom(some.coverage, max.peaks)
seg.dt <- data.table(fit$segments)
new.labels <- rbind(
  l(180000, 200000, "noPeaks"),
  l(206000, 207000, "peakStart"),
  l(208000, 220000, "peakEnd"),
  l(300000, 308250, "peakStart"),
  l(308260, 320000, "peakEnd"))
all.labels <- rbind(Mono27ac[["labels"]], new.labels)[order(chromStart)]
some.labels <- all.labels[lo.lim <= chromStart & chromEnd < hi.lim][c(1,3,4,5)]
FLOPART::FLOPART_data(some.coverage, some.labels)
flopart <- FLOPART::FLOPART(some.coverage, some.labels, 1400)
FLOPART.segs <- flopart[["segments_dt"]]
FLOPART.peaks <- FLOPART.segs[status == "peak"]
model.color <- "blue"
peak.y <- -2
ann.colors <- c(
  noPeaks="orange",
  peakStart="#efafaf",
  peakEnd="#ff4c4c")
gg <- ggplot()+
  geom_rect(aes(
    xmin=chromStart, xmax=chromEnd,
    fill=annotation,
    ymin=-Inf, ymax=Inf),
    data=some.labels,
    color="grey",
    alpha=0.5)+
  theme_bw()+
  scale_fill_manual(values=ann.colors)+
  geom_step(aes(
    chromStart, count),
    color="grey50",
    data=some.coverage)+
  geom_step(aes(
    chromStart, mean),
    color=model.color,
    data=FLOPART.segs)+
  geom_point(aes(
    chromEnd, peak.y),
    color=model.color,
    shape=21,
    data=FLOPART.peaks)
gg

seg.dt.list <- list(
  FLOPART=FLOPART.segs)
for(show.peaks in 5:max.peaks){
  seg.dt.list[[paste(show.peaks, "peaks")]] <- seg.dt[peaks==show.peaks]
}
err.dt.list <- list()
model.segs.list <- list()
for(model in names(seg.dt.list)){
  pkg.segs <- seg.dt.list[[model]][, .(chromStart, chromEnd, mean, status)]
  pkg.peaks <- pkg.segs[status=="peak"]
  err.df <- PeakError::PeakErrorChrom(pkg.peaks, some.labels)
  err.dt.list[[model]] <- data.table(model, err.df)
  model.segs.list[[model]] <- data.table(model, pkg.segs)
}
err.dt <- do.call(rbind, err.dt.list)
model.segs <- do.call(rbind, model.segs.list)
model.peaks <- model.segs[status=="peak"]
gg <- ggplot()+
   theme_bw()+
   theme(panel.spacing=grid::unit(0, "lines"))+
   scale_fill_manual(values=ann.colors)+
   facet_grid(model ~ ., labeller=label_both)+
   geom_rect(aes(
     xmin=chromStart, xmax=chromEnd,
     ymin=-Inf, ymax=Inf,
     fill=annotation),
     color="grey",
     alpha=0.5,
     data=some.labels)+
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
     data=some.coverage)+
   geom_step(aes(
     chromStart, mean),
     color=model.color,
     data=model.segs)+
   geom_segment(aes(
     chromStart, mean,
     xend=chromEnd, yend=mean),
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
     size=2)+
  scale_y_continuous("Count of aligned DNA sequence reads")+
  scale_x_continuous("Position on chromosome (bases)")
out.png <- "figure-smaller-example.png"
png(out.png, 10, 5, units="in", res=200)
print(gg)
dev.off()
