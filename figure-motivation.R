### Write down what package versions work with your R code, and
### attempt to download and load those packages. The first argument is
### the version of R that you used, e.g. "3.0.2" and then the rest of
### the arguments are package versions. For
### CRAN/Bioconductor/R-Forge/etc packages, write
### e.g. RColorBrewer="1.0.5" and if RColorBrewer is not installed
### then we use install.packages to get the most recent version, and
### warn if the installed version is not the indicated version. For
### GitHub packages, write "user/repo@commit"
### e.g. "tdhock/animint@f877163cd181f390de3ef9a38bb8bdd0396d08a4" and
### we use install_github to get it, if necessary.
works_with_R <- function(Rvers,..., local=TRUE){
  if(isTRUE(local)){
    local.lib <- file.path(getwd(), "library")
    dir.create(local.lib, showWarnings=FALSE, recursive=TRUE)
    .libPaths(local.lib)
  }
  pkg_ok_have <- function(pkg,ok,have){
    stopifnot(is.character(ok))
    if(!as.character(have) %in% ok){
      warning("works with ",pkg," version ",
              paste(ok,collapse=" or "),
              ", have ",have)
    }
  }
  pkg_ok_have("R",Rvers,getRversion())
  pkg.vers <- list(...)
  for(pkg.i in seq_along(pkg.vers)){
    vers <- pkg.vers[[pkg.i]]
    pkg <- if(is.null(names(pkg.vers))){
      ""
    }else{
      names(pkg.vers)[[pkg.i]]
    }
    if(pkg == ""){# Then it is from GitHub.
      ## suppressWarnings is quieter than quiet.
      if(!suppressWarnings(require(requireGitHub))){
        ## If requireGitHub is not available, then install it using
        ## devtools.
        if(!suppressWarnings(require(devtools))){
          install.packages("devtools")
          require(devtools)
        }
        install_github("tdhock/requireGitHub")
        require(requireGitHub)
      }
      requireGitHub(vers)
    }else{# it is from a CRAN-like repos.
      if(!suppressWarnings(require(pkg, character.only=TRUE))){
        install.packages(pkg)
      }
      pkg_ok_have(pkg, vers, packageVersion(pkg))
      library(pkg, character.only=TRUE)
    }
  }
}
works_with_R(
  "4.1.1",
  PeakSegOptimal="2018.5.25",
  PeakError="2021.7.1",
  ggplot2="3.3.5",
  data.table="1.14.0")

file.RData <- "Mono27ac.simple.RData"
if(!file.exists(file.RData)){
  u <- paste0(
    "https://github.com/tdhock/FLOPART/raw/main/data/",
    file.RData)
  download.file(u, file.RData)
}
(objs <- load(file.RData))
str(Mono27ac.simple)
max.peaks <- 7L
fit <- PeakSegOptimal::PeakSegPDPAchrom(
  Mono27ac.simple[["coverage"]], max.peaks)
seg.dt <- data.table(fit$segments)
model.color <- "blue"
peak.y <- -2
ann.colors <- c(
  noPeaks="orange",
  peakStart="#efafaf",
  peakEnd="#ff4c4c")
show.segs <- seg.dt[peaks >= 4]
model.peaks <- show.segs[status=="peak"]
err.dt <- model.peaks[, {
  PeakError::PeakErrorChrom(.SD, Mono27ac.simple[["label"]])
}, by=peaks]
gg <- ggplot()+
  ggtitle(
    "Optimal Poisson segmentation with up-down constraint graph (two nodes, two edges): one segment for each peak and background region")+
   theme_bw()+
   theme(panel.spacing=grid::unit(0, "lines"))+
   scale_fill_manual("label", values=ann.colors)+
   facet_grid(peaks ~ ., labeller=label_both)+
   geom_rect(aes(
     xmin=chromStart, xmax=chromEnd,
     ymin=-Inf, ymax=Inf,
     fill=annotation),
     color="grey",
     alpha=0.5,
     data=Mono27ac.simple[["label"]])+
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
     data=Mono27ac.simple[["coverage"]])+
   geom_step(aes(
     chromStart, mean),
     color=model.color,
     data=show.segs)+
   geom_segment(aes(
     chromStart, mean,
     xend=chromEnd, yend=mean),
     color=model.color,
     data=show.segs)+
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
out.png <- "figure-motivation.png"
png(out.png, 15, 4, units="in", res=200)
print(gg)
dev.off()
