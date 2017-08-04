library(data.table);
library(ggplot2);
##library(cellrangerRkit);
##library(DESeq2);
library(rpca);
##library(rsvd);
library(Rtsne);
library(irlba);

diff.exp <- function(sle, ra) {
  sle.ranges <- seq(1,ncol(sle),100);
  sle.sum <- sapply(sle.ranges, function(i) {
    if(i+100>ncol(sle)) {
      rowSums(as.matrix(exprs(sle[,i:ncol(sle)])));
    } else {
      rowSums(as.matrix(exprs(sle[,i:(i+100)])));
    }
  })
  
  ra.ranges <- seq(1,ncol(ra),100);
  ra.sum <- sapply(ra.ranges, function(x) {
    if(x+100>ncol(ra)) {
      rowSums(as.matrix(exprs(ra[,x:ncol(ra)])));
    } else {
      rowSums(as.matrix(exprs(ra[,x:(x+100)])));
    }
  })
  
  ##  browser();
  colData <- data.frame(disease=c(rep("sle",ncol(sle.sum)),rep("ra",ncol(ra.sum))));
  
  dds <- DESeqDataSetFromMatrix(countData = cbind(sle.sum, ra.sum),
                                colData = colData,
                                design = ~disease)
  
  featureData <- fData(sle);
  (mcols(dds) <- DataFrame(mcols(dds), featureData))
  
  dds <- DESeq(dds);
  res <- results(dds);
  
  res2 <- data.frame(res, featureData$symbol)
  
  return(res2[order(res2$padj),])
}

pbmc8k <- fread("pbmc8k_S1_merged_1.dge.txt");
pbmc8k.rownames <- pbmc8k$V1;
pbmc8k <- pbmc8k[,-1,with=F];
rownames(pbmc8k) <- pbmc8k.rownames;

##pure.cell.types <- readRDS("all_pure_select_11types2.rds");
##load("cell.type.cor.classification.rda");
##pure.major.cell.types <- pure.cell.types;
##pure.major.cell.types$pure_id <- pure.cell.types$pure_id[c(2,7,8,9,10)];
##pure.major.cell.types$pure_avg <- pure.cell.types$pure_avg[c(2,7,8,9,10),];

pbmc8k.colsums <- colSums(pbmc8k);
pbmc8k.median.sum <- median(pbmc8k.colsums);
pbmc8k.norm <- log2(sweep(pbmc8k, 2, pbmc8k.median.sum/pbmc8k.colsums,"*")+1);
##total.norm.most.variable.genes <- total.norm[rownames(training.norm.most.variable.genes),];

##total.matched <- total[match(rownames(cell.type.classification$cell.types.avg), fData(total)$symbol), ]
##total.match <- cor(as.matrix(exprs(total.matched)),cell.type.classification$cell.types.avg, method='spearman')
##total.labels <- colnames(cell.type.classification$cell.type)[as.numeric(apply(total.match, 1, which.max))]

##genes <- c("CD3D", "MS4A1", "LYZ", "CD19", "FCGR3A",
##           "CD8A", "NKG7", "CD14", "CCR7","IL7R","MS4A7","GNLY","FCER1A","CST3", "CD79A","S100A8","LDHB","CCL5","GZMB",
##           "CLEC9A", "CD34", "CD1C", "CD27", "CD74", "FTL");

##genes.id <- match(genes, rownames(pbmc8k));

pbmc8k.matched <- pbmc8k.norm[match(rownames(cell.type.classification$cell.types.avg), rownames(pbmc8k.norm)),];
pbmc8k.match <- cor(pbmc8k.matched,cell.type.classification$cell.types.avg, method='spearman',use='complete.obs')
pbmc8k.labels <- colnames(cell.type.classification$cell.type)[as.numeric(apply(pbmc8k.match, 1, which.max))]

##pbmc8k.match <- cor(pbmc8k[genes.id,],
##                   t(pure.major.cell.types$pure_avg[,match(genes,pure.major.cell.types$pure_use_gene_names)]),
##                   method='spearman')
##pbmc8k.match[which(is.na(pbmc8k.match))] <- 0;
##pbmc8k.pure.major.cell.types <- pure.major.cell.types$pure_id[apply(pbmc8k.match,1,which.max)];

##genes.matched <- intersect(pure.cell.types$pure_use_gene_names,rownames(pbmc8k));
##pbmc8k.match.all <- cor(pbmc8k[match(genes.matched, rownames(pbmc8k))],
##                    t(pure.cell.types$pure_avg[,match(genes.matched, pure.cell.types$pure_use_gene_names)]),
##                    method='spearman')
##pbmc8k.pure.all.cell.types <- pure.cell.types$pure_id[apply(pbmc8k.match.all,1,which.max)];

use_genes <- which(rowSums(pbmc8k)!=0);

umi=colSums(pbmc8k);
ribo <- grep("^RP",rownames(pbmc8k));
ribo.perc <- colSums(pbmc8k[ribo,])/umi;

pbmc8k.regress <- t(lm(t(pbmc8k.norm[use_genes,])~umi+ribo.perc)$residual);

pbmc8k.pca <- irlba::irlba(pbmc8k.regress[order(apply(pbmc8k.regress,1,var),decreasing=T)[1:2000],], 20);
pbmc8k.tsne <- Rtsne(pbmc8k.pca$v);

pbmc8k.df <- data.frame(pbmc8k.tsne$Y,
                        umi=umi,
                        types=pbmc8k.labels,
                        types.all=pbmc8k.pure.all.cell.types);
colnames(pbmc8k.df) <- c("tsne1","tsne2","umi","types","all.types");

ggsave("pbmc.cell.pdf",ggplot(aes(tsne1, tsne2, color=types), data=pbmc8k.df)+geom_point()+theme_bw()+
         theme(panel.border = element_blank(), panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")));

write.table(pbmc8k.df, file="pbmc8k.meta.txt", row.names=T, col.names=T, quote=F, sep="\t");