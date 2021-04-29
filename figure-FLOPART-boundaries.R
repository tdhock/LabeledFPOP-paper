works_with_R(
  "4.0.4",
  PeakSegOptimal="2018.5.25",
  PeakError="2017.6.19",
  "tdhock/FLOPART@b14fd70274ad2313c0742003270d8d66da28456d")
library(ggplot2)
library(data.table)
l <- function(chromStart, chromEnd, annotation){
  data.table(chromStart, chromEnd, annotation)
}
d <- function(chromStart, chromEnd, count){
  data.table(chromStart, chromEnd, count)
}
cov.dt <- rbind(
  d(0, 10, 3),
  d(10, 20, 4),
  d(20, 30, 15),
  d(30, 40, 23),
  d(40, 50, 8),
  d(50, 60, 9),
  d(60, 70, 45),
  d(70, 80, 43),
  d(80, 90, 5),
  d(90, 100, 2))

(gg.cov <- ggplot()+
  geom_step(aes(
    chromStart+0.5, count),
    data=cov.dt))

lab.dt <- rbind(
  l(25, 42, "peakEnd"),
  l(43, 45, "noPeaks"),
  l(46, 48, "noPeaks"),
  l(48, 55, "noPeaks"),
  l(58, 80, "peakStart"))
ann.colors <- c(
  noPeaks="#f6f4bf",
  peakStart="#ffafaf",
  peakEnd="#ff4c4c",
  peaks="#a445ee")
gg.cov+
  theme_bw()+
  geom_rect(aes(
    xmin=chromStart+0.5,
    xmax=chromEnd+0.5,
    ymin=-Inf, ymax=Inf,
    fill=annotation),
    data=lab.dt,
    color="grey",
    alpha=0.5)+
  scale_fill_manual(values=ann.colors)

get.pos <- function(dt)dt[, c(chromStart, chromEnd)]
uniq.pos <- sort(unique(c(get.pos(lab.dt), get.pos(cov.dt))))
uniq.dt <- data.table(
  chromStart=uniq.pos[-length(uniq.pos)],
  chromEnd=uniq.pos[-1])

with.counts <- cov.dt[
  uniq.dt,
  .(chromStart=i.chromStart, chromEnd=i.chromEnd, count,
    index = seq_along(count), weight = i.chromEnd - i.chromStart),
  on=.(chromStart < chromEnd, chromEnd > chromStart)]
label_code <- FLOPART::get_label_code()
ann.colors <- c(
  noPeaks="#f6f4bf",
  peakStart="#ffafaf",
  peakEnd="#ff4c4c",
  peaks="#a445ee")
label.index.dt <- with.counts[
  lab.dt,
  .(start=x.index[1], end=x.index[.N],
    labelStart=i.chromStart,
    labelEnd=i.chromEnd,
    type=label_code[annotation], annotation),
  by=.EACHI,
  on=.(chromEnd > chromStart, chromStart < chromEnd)]

## should look the same.
ggplot()+
  geom_step(aes(
    chromStart+0.5, count),
    data=with.counts)+
  theme_bw()+
  geom_rect(aes(
    xmin=chromStart+0.5,
    xmax=chromEnd+0.5,
    ymin=-Inf, ymax=Inf,
    fill=annotation),
    data=label.index.dt,
    color="grey",
    alpha=0.5)+
  scale_fill_manual(values=ann.colors)

gg <- ggplot()+
  geom_point(aes(
    index, count),
    data=with.counts)+
  theme_bw()+
  geom_rect(aes(
    xmin=start-0.5,
    xmax=end+0.5,
    ymin=-Inf, ymax=Inf,
    fill=annotation),
    data=label.index.dt,
    color="grey",
    alpha=0.5)+
  scale_fill_manual(values=ann.colors)
png("figure-FLOPART-boundaries.png", width=4, height=2, res=200, units="in")
print(gg)
dev.off()
