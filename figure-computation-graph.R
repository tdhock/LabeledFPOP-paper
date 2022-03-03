node.dt <- data.table::CJ(state=c(-1,1), data.num=1:15)
node.dt[, next.num := data.num+1]
edge.dt <- node.dt[node.dt, on=.(next.num=data.num), nomatch=0L]
edge.dt[, type := ifelse(i.state==state, "no change", "change")]
library(ggplot2)
gg <- ggplot()+
  geom_point(aes(
    data.num, state),
    shape=21,
    size=10,
    data=node.dt,
    fill="white",
    color="black")+
  geom_segment(aes(
    data.num, state,
    linetype=type,
    xend=next.num-0.1, yend=i.state+ifelse(i.state!=state, -i.state*0.1, 0)),
    arrow=grid::arrow(length=unit(0.1, "in"), type="closed"),
    data=edge.dt)
png("figure-computation-graph.png", width=10, height=2, units="in", res=200)
print(gg)
dev.off()
