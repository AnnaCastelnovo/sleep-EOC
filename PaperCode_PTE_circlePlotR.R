library(circlize)
library(R.matlab)
library(stringr)
library("readxl")

#add path
setwd("/.../CirclePlots/")
set.seed(14)
TFM <- readMat("TFM_reord.mat")
labelsh <- readMat("reord_labelHigh.mat")
labels <- readMat("reord_scouts2.mat")
label = labels$test
labelh = labelsh$reord.label
temp = list.files(pattern="*.mat")

for (i in 1:length(temp)) {
  assign('TFM', readMat(temp[i]))
  name = temp[i]
  
  tfm = TFM$reord.M
  dim(tfm)
  
  rownames(tfm) <- sapply( label, paste0)
  colnames(tfm) <- sapply( label, paste0)
  
  df <- data.frame(from = rep(rownames(tfm), ncol(tfm)),
                   to = rep(colnames(tfm), each = nrow(tfm)),
                   value = as.vector(tfm))
  
  
  mat = df
  nami = str_trim(substr(name[1], 1, (nchar(name)-4)))
  
  #replace with correct path
  s_name = paste("path/to/save/imgs",nami,".png")
  s_name = gsub(" ", "", s_name, fixed = TRUE)
  png(s_name, width = 2300, height = 2300, res=150) 

  nm = unique(unlist(rownames(tfm)))
  group = structure(substr(nm, nchar(nm), nchar(nm)), names = nm)
  
  circos.par(start.degree = 250)
  circos.par("canvas.xlim" = c(-1.3, 1.3), "canvas.ylim" = c(-1.3, 1.3))
  par(cex = 1.8)
  
  rainb = c(rainbow(length(nm)/2, start = 0.1, end =0.4), rainbow(length(nm)/2, start = 0.5, end =0.9))
  
  grid.col <- setNames(rainb, nm)

  chordDiagram(tfm, group = group, big.gap = 20, annotationTrack = "grid", preAllocateTracks = 1, directional = 1, grid.col = grid.col, direction.type = c("diffHeight", "arrows"), link.arr.type = "big.arrow",scale = TRUE)
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    
    circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.25), cex = 0.5)
  }, bg.border = NA)
  
  

  dev.off()
  
  circos.clear()
  
}




