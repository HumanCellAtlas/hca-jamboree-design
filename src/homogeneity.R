


homogeneity <- function(data, cluster){
  library(edgeR)
  cluster_id <- unique(cluster)
  result <- data.frame(mean_gini=unlist(sapply(cluster_id,
         function(cluster_id, data,cluster){
           mean(gini(t(data[,cluster==cluster_id])), na.rm = TRUE)},
         data,cluster)),
         cluster_id = cluster_id)

}
