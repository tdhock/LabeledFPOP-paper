if(!file.exists("signal.list.annotation.sets.RData")){
  download.file("http://members.cbio.mines-paristech.fr/~thocking/neuroblastoma/signal.list.annotation.sets.RData", "signal.list.annotation.sets.RData")
}
load("signal.list.annotation.sets.RData")

library(data.table)
ann.dt <- do.call(rbind, lapply(names(annotation.sets), function(set){
  data.table(set, annotation.sets[[set]])
}))
ann.dt[, pid.chr := paste0(profile.id, ".", chromosome)]
sig.nrow.vec <- sapply(signal.list, nrow)
sig.nrow.dt <- data.table(pid.chr=names(sig.nrow.vec), nrow=sig.nrow.vec)
counts.tall <- ann.dt[, .(
  labels=.N
), by=.(set, pid.chr, annotation)]
counts.wide <- dcast(
  counts.tall,
  set + pid.chr ~ annotation,
  value.var="labels")[sig.nrow.dt, on=.(pid.chr)]
one.of.each <- counts.wide[!is.na(`1breakpoint`) & !is.na(`0breakpoints`)]
one.of.each[, labels := `1breakpoint` + `0breakpoints`]

library(ggplot2)
ggplot()+
  geom_point(aes(
    labels, nrow),
    data=one.of.each)+
  scale_y_log10()

set.counts <- one.of.each[, .(sets=.N), by=pid.chr][order(sets)]
one.set <- set.counts[sets==1]
keep <- one.of.each[pid.chr %in% one.set$pid.chr]
ann.keep <- ann.dt[keep, on=.(set, pid.chr)]
ann.keep[, fold := rep(1:2, l=.N), by=.(pid.chr)]

sig.dt <- do.call(rbind, lapply(names(signal.list), function(pid.chr){
  data.table(pid.chr, signal.list[[pid.chr]])
}))
sig.keep <- sig.dt[keep, on=.(pid.chr)]
sig.keep[, change.after := floor(c(diff(position)/2, NA)+position)+0.5, by=pid.chr]
sig.keep[, data.i := 1:.N, by=pid.chr]
join.dt <- sig.keep[ann.keep, .(
  min,
  max,
  annotation,
  first.i=min(data.i),
  last.i=max(data.i)+1L
), by=.EACHI, on=.(pid.chr, change.after > min, change.after < max)]
setkey(join.dt, pid.chr, first.i)
join.dt[, next.first := c(first.i[-1], NA), by=pid.chr]

next.overlaps <- join.dt[next.first<last.i, .(pid.chr, annotation, first.i, last.i, next.first)]

ann.no.overlap <- join.dt[!pid.chr %in% next.overlaps$pid.chr, .(pid.chr, annotation, first.i, last.i, min, max)]

counts.tall <- ann.no.overlap[, .(
  labels=.N
), by=.(pid.chr, annotation)]
counts.wide <- dcast(
  counts.tall,
  pid.chr ~ annotation,
  value.var="labels")[sig.nrow.dt, on=.(pid.chr)]
one.of.each <- counts.wide[!is.na(`1breakpoint`) & !is.na(`0breakpoints`)]
one.of.each[, labels := `1breakpoint` + `0breakpoints`]
ggplot()+
  geom_point(aes(
    labels, nrow),
    data=one.of.each)+
  scale_y_log10()

sig.no.overlap <- sig.keep[pid.chr %in% ann.no.overlap$pid.chr, .(
  pid.chr, data.i, position, logratio)]

fwrite(sig.no.overlap, "labeled-signals.csv")
fwrite(ann.no.overlap, "labeled-regions.csv")

labeled.data <- list(
  signals=sig.no.overlap,
  regions=ann.no.overlap)
saveRDS(labeled.data, "labeled-data.rds", version=2)
