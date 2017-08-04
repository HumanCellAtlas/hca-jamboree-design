# define centers
library(FNN)
data.use <- GetCellEmbeddings(pbmc,reduction.type = "pca",dims.use = 1:10)
n.cells  <- nrow(x = data.use)
cluster.ids <- GetClusters(pbmc)
neighbors.to.consider <- 150
my.knn <- get.knn(
  data <- as.matrix(x = data.use),
  k = neighbors.to.consider
)
nn.ranked <- cbind(1:n.cells, my.knn$nn.index[, 1:(k.param-1)])
neighbors <- my.knn$nn.index
cellnames <- row.names(data.use)
purity <- lapply(1:nrow(neighbors), function(i) {
  curr.cell <- row.names(data.use)[i]
  cell.cluster <- cluster.ids[i,2]
  neighbor.clusters <- cluster.ids[neighbors[i,],2]
  purity <- mean(neighbor.clusters==cell.cluster)
  cbind(curr.cell,cell.cluster,purity)
})
out.mat <- do.call(rbind,purity)