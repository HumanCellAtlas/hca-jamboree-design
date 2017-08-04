rm(list=ls(all=TRUE))
library(ggplot2)
library(dplyr)
setwd("/home/jovyan/shared_scratch/group 5/task6/Seurat_pbmc/")

dirs <- list.dirs(getwd())
dirs <- dirs[3:length(dirs)]
dirs <- dirs[-2]
plotmat <- read.table(paste0(dirs[1],"/pbmc_purity.txt"),header=TRUE)
colname <- gsub("/home/jovyan/shared_scratch/group 5/task6/Seurat_pbmc/","",dirs[1])
plotmat <- cbind(plotmat,file=rep(colname,nrow(plotmat)))
for (i in 2:length(dirs)) {
  colname <- gsub("/home/jovyan/shared_scratch/group 5/task6/Seurat_pbmc/","",dirs[i])
  if (file.exists(paste0(dirs[i],"/pbmc_purity.txt"))) {
    purr <- read.table(paste0(dirs[i],"/pbmc_purity.txt"),header=TRUE)
    print(1)
    print(colnames(purr))
  } else if (file.exists(paste0(dirs[i],"/purity.txt"))) {
    purr <- read.table(paste0(dirs[i],"/purity.txt"),header=TRUE)
  } else {
    next
  }
  purr <- cbind(purr,file=rep(colname,nrow(purr)))
  plotmat <- rbind(plotmat,purr)
}
down <- c(1,1,0,1,0,1,0,1,0,1,0,1,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0)
sub <- c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1)
downsub <- c(1,0,1,0,1,0,1,0,1,0,1,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0)
files_down <- gsub("/home/jovyan/shared_scratch/group 5/task6/Seurat_pbmc/","",dirs[down==1])
files_sub <- gsub("/home/jovyan/shared_scratch/group 5/task6/Seurat_pbmc/","",dirs[sub==1])
files_downsub <- gsub("/home/jovyan/shared_scratch/group 5/task6/Seurat_pbmc/","",dirs[downsub==1])
subsamp_down <- subset(plotmat,file %in% files_down)
subsamp_sub <- subset(plotmat,file %in% files_sub)
subsamp_downsub <- subset(plotmat,file %in% files_downsub)

png("/home/jovyan/shared_scratch/group 5/task6/subsamp_down.png",width=1000,height=480)
ggplot(aes(y = purity, x = file, fill = factor(cell.cluster)), data = subsamp_down) + geom_boxplot()
dev.off()

subsamp_sub$file <- ifelse(subsamp_sub$file=="Sub10","Sub90",
                           ifelse(subsamp_sub$file=="Sub20","Sub80",
                                  ifelse(subsamp_sub$file=="Sub30","Sub70",
                                         ifelse(subsamp_sub$file=="Sub40","Sub60",
                                                ifelse(subsamp_sub$file=="Sub60","Sub40",
                                                       ifelse(subsamp_sub$file=="Sub70","Sub30",
                                                              ifelse(subsamp_sub$file=="Sub80","Sub20",
                                                                     ifelse(subsamp_sub$file=="Sub90","Sub10",
                                                                            ifelse(subsamp_sub$file=="Sub50","Sub50","All_cells")))))))))
png("/home/jovyan/shared_scratch/group 5/task6/subsamp_sub.png",width=1000,height=480)
ggplot(aes(y = purity, x = file, fill = factor(cell.cluster)), data = subsamp_sub) + geom_boxplot()
dev.off()

png("/home/jovyan/shared_scratch/group 5/task6/subsamp_subdown.png",width=1000,height=480)
ggplot(aes(y = purity, x = file, fill = factor(cell.cluster)), data = subsamp_downsub) + geom_boxplot()
dev.off()
