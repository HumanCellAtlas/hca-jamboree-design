library(subSeq);
library(data.table);

pbmc <- fread("pbmc8k_S1_merged_1.dge.txt");
rownames <- pbmc$V1;
pbmc <- pbmc[,-1,with=F];

pbmc.10 <- pbmc[,-sample(colnames(pbmc),.1*ncol(pbmc)),with=F];
pbmc.20 <- pbmc[,-sample(colnames(pbmc),.2*ncol(pbmc)),with=F];
pbmc.30 <- pbmc[,-sample(colnames(pbmc),.3*ncol(pbmc)),with=F];
pbmc.40 <- pbmc[,-sample(colnames(pbmc),.4*ncol(pbmc)),with=F];
pbmc.50 <- pbmc[,-sample(colnames(pbmc),.5*ncol(pbmc)),with=F];
pbmc.60 <- pbmc[,-sample(colnames(pbmc),.6*ncol(pbmc)),with=F];
pbmc.70 <- pbmc[,-sample(colnames(pbmc),.7*ncol(pbmc)),with=F];
pbmc.80 <- pbmc[,-sample(colnames(pbmc),.8*ncol(pbmc)),with=F];
pbmc.90 <- pbmc[,-sample(colnames(pbmc),.9*ncol(pbmc)),with=F];

pbmc.10.subsample.90 <- generateSubsampledMatrix(pbmc.10, sum(pbmc.90)/sum(pbmc.10), seed=1234)
pbmc.20.subsample.90 <- generateSubsampledMatrix(pbmc.20, sum(pbmc.90)/sum(pbmc.20), seed=1234)
pbmc.30.subsample.90 <- generateSubsampledMatrix(pbmc.30, sum(pbmc.90)/sum(pbmc.30), seed=1234)
pbmc.40.subsample.90 <- generateSubsampledMatrix(pbmc.40, sum(pbmc.90)/sum(pbmc.40), seed=1234)
pbmc.50.subsample.90 <- generateSubsampledMatrix(pbmc.50, sum(pbmc.90)/sum(pbmc.50), seed=1234)
pbmc.60.subsample.90 <- generateSubsampledMatrix(pbmc.60, sum(pbmc.90)/sum(pbmc.60), seed=1234)
pbmc.70.subsample.90 <- generateSubsampledMatrix(pbmc.70, sum(pbmc.90)/sum(pbmc.70), seed=1234)
pbmc.80.subsample.90 <- generateSubsampledMatrix(pbmc.80, sum(pbmc.90)/sum(pbmc.80), seed=1234)

write.table(cbind(rownames,pbmc.10), file="pbmc_downsampled_cells_10%.txt",row.names=F, col.names=T,quote=F, sep="\t")
write.table(cbind(rownames,pbmc.20), file="pbmc_downsampled_cells_20%.txt",row.names=F, col.names=T,quote=F, sep="\t");
write.table(cbind(rownames,pbmc.30), file="pbmc_downsampled_cells_30%.txt",row.names=F, col.names=T,quote=F, sep="\t");
write.table(cbind(rownames,pbmc.40), file="pbmc_downsampled_cells_40%.txt",row.names=F, col.names=T,quote=F, sep="\t");
write.table(cbind(rownames,pbmc.50), file="pbmc_downsampled_cells_50%.txt",row.names=F, col.names=T,quote=F, sep="\t");
write.table(cbind(rownames,pbmc.60), file="pbmc_downsampled_cells_60%.txt",row.names=F, col.names=T,quote=F, sep="\t");
write.table(cbind(rownames,pbmc.70), file="pbmc_downsampled_cells_70%.txt",row.names=F, col.names=T,quote=F, sep="\t");
write.table(cbind(rownames,pbmc.80), file="pbmc_downsampled_cells_80%.txt",row.names=F, col.names=T,quote=F, sep="\t");
write.table(cbind(rownames,pbmc.90), file="pbmc_downsampled_cells_90%.txt",row.names=F, col.names=T,quote=F, sep="\t");


write.table(cbind(rownames,pbmc.10.subsample.90), file="pbmc_downsampled_cells_10%_subsampled_90%.txt", row.names=F, col.names=T, quote=F, sep="\t")
write.table(cbind(rownames,pbmc.20.subsample.90), file="pbmc_downsampled_cells_20%_subsampled_90%.txt",row.names=F, col.names=T, quote=F, sep="\t")
write.table(cbind(rownames,pbmc.30.subsample.90), file="pbmc_downsampled_cells_30%_subsampled_90%.txt",row.names=F, col.names=T, quote=F, sep="\t")
write.table(cbind(rownames,pbmc.40.subsample.90), file="pbmc_downsampled_cells_40%_subsampled_90%.txt",row.names=F, col.names=T, quote=F, sep="\t")
write.table(cbind(rownames,pbmc.50.subsample.90), file="pbmc_downsampled_cells_50%_subsampled_90%.txt",row.names=F, col.names=T, quote=F, sep="\t")
write.table(cbind(rownames,pbmc.60.subsample.90), file="pbmc_downsampled_cells_60%_subsampled_90%.txt",row.names=F, col.names=T, quote=F, sep="\t")
write.table(cbind(rownames,pbmc.70.subsample.90), file="pbmc_downsampled_cells_70%_subsampled_90%.txt",row.names=F, col.names=T, quote=F, sep="\t")
write.table(cbind(rownames,pbmc.80.subsample.90), file="pbmc_downsampled_cells_80%_subsampled_90%.txt",row.names=F, col.names=T, quote=F, sep="\t")

pbmc.subsample.10 <- generateSubsampledMatrix(pbmc, .1, seed=1234)
pbmc.subsample.20 <- generateSubsampledMatrix(pbmc, .2, seed=1234)
pbmc.subsample.30 <- generateSubsampledMatrix(pbmc, .3, seed=1234)
pbmc.subsample.40 <- generateSubsampledMatrix(pbmc, .4, seed=1234)
pbmc.subsample.50 <- generateSubsampledMatrix(pbmc, .5, seed=1234)
pbmc.subsample.60 <- generateSubsampledMatrix(pbmc, .6, seed=1234)
pbmc.subsample.70 <- generateSubsampledMatrix(pbmc, .7, seed=1234)
pbmc.subsample.80 <- generateSubsampledMatrix(pbmc, .8, seed=1234)
pbmc.subsample.90 <- generateSubsampledMatrix(pbmc, .9, seed=1234)

##write.table(cbind(rownames,pbmc.50), file="pbmc_downsampled_cells_50%_subsampled.txt")
write.table(cbind(rownames,pbmc.subsample.10), file="pbmc_downsampled_cells_10%_subsampled_90%.txt", row.names=F, col.names=T, quote=F, sep="\t")
write.table(cbind(rownames,pbmc.subsample.20), file="pbmc_downsampled_cells_20%_subsampled_90%.txt",row.names=F, col.names=T, quote=F, sep="\t")
write.table(cbind(rownames,pbmc.subsample.30), file="pbmc_downsampled_cells_30%_subsampled_90%.txt",row.names=F, col.names=T, quote=F, sep="\t")
write.table(cbind(rownames,pbmc.subsample.40), file="pbmc_downsampled_cells_40%_subsampled_90%.txt",row.names=F, col.names=T, quote=F, sep="\t")
write.table(cbind(rownames,pbmc.subsample.50), file="pbmc_downsampled_cells_50%_subsampled_90%.txt",row.names=F, col.names=T, quote=F, sep="\t")
write.table(cbind(rownames,pbmc.subsample.60), file="pbmc_downsampled_cells_60%_subsampled_90%.txt",row.names=F, col.names=T, quote=F, sep="\t")
write.table(cbind(rownames,pbmc.subsample.70), file="pbmc_downsampled_cells_70%_subsampled_90%.txt",row.names=F, col.names=T, quote=F, sep="\t")
write.table(cbind(rownames,pbmc.subsample.80), file="pbmc_downsampled_cells_80%_subsampled_90%.txt",row.names=F, col.names=T, quote=F, sep="\t")
write.table(cbind(rownames,pbmc.subsample.90), file="pbmc_downsampled_cells_80%_subsampled_90%.txt",row.names=F, col.names=T, quote=F, sep="\t")
