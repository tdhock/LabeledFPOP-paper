library(data.table)

labeled.data <- readRDS("labeled-data.rds")

pid.chr <- "1.1"
sig.dt <- labeled.data$signals[pid.chr, on="pid.chr"]
lab.dt <- labeled.data$regions[pid.chr, on="pid.chr"]

change.colors <- c(`0breakpoints` = "#f6f4bf", `1breakpoint` = "#ff7d7d")
library(ggplot2)
ggplot()+
  scale_fill_manual(values=change.colors)+
  geom_rect(aes(
    xmin=min, xmax=max,
    ymin=-Inf, ymax=Inf,
    fill=annotation),
    alpha=0.5,
    data=lab.dt)+
  geom_point(aes(
    position, logratio),
    data=sig.dt)

ggplot()+
  scale_fill_manual(values=change.colors)+
  geom_rect(aes(
    xmin=first.i, xmax=last.i,
    ymin=-Inf, ymax=Inf,
    fill=annotation),
    alpha=0.5,
    data=lab.dt)+
  geom_point(aes(
    data.i, logratio),
    data=sig.dt)
