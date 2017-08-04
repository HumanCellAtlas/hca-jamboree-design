

pbmc_all = read.table("pbmc8k_S1_merged_1.dge.txt", sep="\t", header=T, row.names = 1, stringsAsFactors = F)
total_counts = colSums(pbmc_all)



pbmc_1 = generateSubsampledMatrix(pbmc_all, 0.1, seed=1234)
total_counts = colSums(pbmc_1)

pbmc_2 = generateSubsampledMatrix(pbmc_all, 0.5, seed=1234)
total_counts = colSums(pbmc_2)

pbmc = CreateSeuratObject(project = "10X_pbmc", raw.data = pbmc_2)

pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
pbmc <- ScaleData(object = pbmc)

pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, 
               genes.print = 5)

length(pbmc@var.genes)
write.table(pbmc@var.genes, file="pbmc_0.5_var_genes.txt", sep="\t", quote=F, row.names=F, col.names=F)
write.table(pbmc@dr$pca@cell.embeddings, file="pbmc_0.5_pca_cells.txt", sep="\t", quote=F)
write.table(pbmc@dr$pca@gene.loadings, file="pbmc_0.5_pca_genes.txt", sep="\t", quote=F)

pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:10, 
                     resolution = 0.6, print.output = 0, save.SNN = TRUE)

pbmc <- RunTSNE(object = pbmc, dims.use = 1:10, do.fast = TRUE)


pdf("pbmc_0.5_pca_clusters.pdf")
data.plot = as.data.frame(cbind(pbmc@dr$pca@cell.embeddings, pbmc@ident))
ggplot(data.plot, aes(x = PC1, y = PC2)) + 
  geom_point(aes(colour = pbmc@ident),  size = 1, alpha = 0.8) + 
  scale_color_hue(l = 65) + 
  scale_size(range = c(1, 1)) +
  theme_bw() +
  theme(legend.title=element_blank(),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        axis.text=element_text(size=17, colour = "black"),
        axis.title=element_text(size=19, colour = "black"),
        legend.text=element_text(size=15, colour = "black"),
        legend.key = element_rect(colour = NA),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"),
        plot.margin = unit(c(0.2, 0, 0, 0.2), "lines"))+
  coord_fixed(ratio=0.7)
dev.off()



pdf("pbmc_0.5_tsne_clusters.pdf")
data.plot = as.data.frame(cbind(pbmc@dr$tsne@cell.embeddings, pbmc@ident))
ggplot(data.plot, aes(x = tSNE_1, y = tSNE_2)) + 
  geom_point(aes(colour = pbmc@ident),  size = 1, alpha = 0.8) + 
  scale_color_hue(l = 65) + 
  scale_size(range = c(1, 1)) +
  theme_bw() +
  theme(legend.title=element_blank(),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        axis.text=element_text(size=17, colour = "black"),
        axis.title=element_text(size=19, colour = "black"),
        legend.text=element_text(size=15, colour = "black"),
        legend.key = element_rect(colour = NA),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"),
        plot.margin = unit(c(0.2, 0, 0, 0), "lines"))+
  coord_fixed(ratio=0.8)
dev.off()


meta = read.table("pbmc8k.meta.txt", sep="\t", header=T)
m1 = match(pbmc@cell.names, rownames(meta))
new_meta = meta[m1,'types']
names(new_meta) = pbmc@cell.names
pbmc = AddMetaData(pbmc, new_meta, "type")
write.table(pbmc@meta.data, file="pbmc_0.5_meta.txt", quote=F, sep="\t")


pdf("tnse_type.pdf")
TSNEPlot(pbmc, group.by="type")
dev.off()


library(FNN)
data.use <- GetCellEmbeddings(pbmc,reduction.type = "pca",dims.use = 1:10)
n.cells  <- nrow(x = data.use)
cluster.ids <- GetClusters(pbmc)
neighbors.to.consider <- 150
my.knn <- get.knn(
  data <- as.matrix(x = data.use),
  k = neighbors.to.consider
)
k.param = neighbors.to.consider
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


