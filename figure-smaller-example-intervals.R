library("ggplot2")
library(data.table)
data("Mono27ac.simple", package="FLOPART")
Mono27ac.simple
label.pen <- 1400
fit <- with(Mono27ac.simple, FLOPART::FLOPART(coverage, label, label.pen))
lapply(fit, head)

## Plot number of intervals stored in cost function, versus data
## point, for each cost status.
imat <- fit[["intervals_mat"]]
interval.dt <- data.table(
  intervals=as.integer(imat),
  status=c("peak", "background")[as.integer(col(imat))],
  data.i=as.integer(row(imat)))
ann.colors <- c(
  noPeaks="orange",
  peakStart="#efafaf",
  peakEnd="#ff4c4c")
gg <- ggplot()+
  scale_fill_manual("label", values=ann.colors)+
  geom_rect(aes(
    xmin=firstRow-0.5, xmax=lastRow+0.5,
    ymin=-Inf, ymax=Inf,
    fill=annotation),
    alpha=0.5,
    color="grey",
    data=fit[["label_dt"]])+
  geom_line(aes(
    data.i, intervals, color=status),
    size=1,
    data=interval.dt)
png("figure-smaller-example-intervals.png", width=20, height=3, units="in", res=200)
print(gg)
dev.off()
