#' @title Function to simulate multigroup single cell datasets 
#'
#' @description 
#'
#' @param no_clust Number of cell groups to simulate
#' @param no_cells Total number of cells to simulate
#' @param mean_expr Vector of mean expression values - one for each gene
#' @param var_expr Vector of variance values - in the same order as mean_expr
#' @param proportions Vector of percentages which dictates the size of each group relative to the total number of cell
#' @param cvRatio Vector of scaling factors to increase noise withtin each cell population
#' @param distances Vector of proprtions that dictates the number of genes differentially expressed relative to the total number of genes
#' @param log2FC Scaling factor for strength of differential expression
#' @param genes_up Vector of proporations indicating the relative number of genes up-regulated in each group
#' 
#' @return A matrix of simulated expression counts
#'
#' @examples
#'

# To do
# Add variance in library size
# Add covariance structure
# Add technical noise 

Simu_data <- function(
  no_clust = 2,
  no_cells = 1000,
  mean_expr = NULL,
  var_expr = NULL,
  proportions = NULL,
  cvRatio = NULL,
  distances = NULL,
  log2FC = NULL,
  genes_up = NULL,
  ...
){
  if(is.null(mean_expr) & is.null(var_expr)) stop("Please supply mean and variance estimates.")
  
  # Generate distances vector if none was submitted
  if(is.null(distances)){
    distances = runif(no_clust)
  }
  
  # Generate log2FC vector if none was submitted
  if(is.null(log2FC)){
    log2FC = rep(1, no_clust)
  }
  
  # Generate cvRatio vector if none was submitted
  if(is.null(cvRatio)){
    cvRatio = rep(1, no_clust)
  }
  
  # Generate genes_up vector if none was submitted
  if(is.null(genes_up)){
    genes_up = rep(0.5, no_cells)
  }
  
  # Fit the mean variance trend 
  spl <- smooth.spline(log(var_expr) ~ log(mean_expr))
  
  # Generate mean and variance expression vector per group
  means <- list()
  vars <- list()
  for(i in 1:no_clust){
   cur_means <- mean_expr
   cur_var <- var_expr
   sam <- sample(1:length(mean_expr), length(mean_expr)*distances[i])
   sam_up <- sample(1:length(sam), length(sam)*genes_up[i])
   cur_means[sam[sam_up]] <- mean_expr[sam[sam_up]]*log2FC[i]
   cur_means[sam[-sam_up]] <- mean_expr[sam[-sam_up]]/log2FC[i]
   cur_var[sam] <- exp(predict(object = spl, log(cur_means[sam]))$y)
   means[[i]] <- cur_means
   vars[[i]] <- cur_var
  }
  
  source("shared_scratch/group2/Task6/nbSimulations.R")
  
  # Generate the cell groups
  cells = simCountsNB(means = as.numeric(means[[1]]), 
                      variances = as.numeric(vars[[1]]), 
                      nSamples = no_cells*proportions[1],
                      cvRatio = cvRatio[1])
  for(i in 2:no_clust){
    cur_counts <- simCountsNB(means = means[[i]], 
                              variances = vars[[i]], 
                              nSamples = no_cells*proportions[i],
                              cvRatio = cvRatio[i])
    cells <- cbind(cells, cur_counts)
  }
  
  # Introduce variance in cell sizes
  # lib.size <- rlnorm(ncol(cells), mean = 1, sd = 0.6)
  # cells <- t(t(cells)/lib.size)
  
  return(cells)
}