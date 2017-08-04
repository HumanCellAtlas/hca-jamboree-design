library(data.table);
library(DESeq2);

file <- "pbmc_downsampled_cells_10%.txt";
pbmc <- fread(file);

