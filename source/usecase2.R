#hypomethylated regions

compute_cell_type_scores <- function(mean_score_matrix, sd_score_matrix, invert = FALSE){
    
  if(invert){
    sign_for_score <- -1
    direction <- "higher"
    sign_for_label <- "+"
  } 
  else{
    sign_for_score <- 1
    direction <- "lower"
    sign_for_label <- "-"
  } 
  #the three metrices described above
  mean_score <- t(apply(sign_for_score * mean_score_matrix, 
                       1, rank, ties.method = "max") - 1)
  mean_1_sd <- t(apply(sign_for_score * (mean_score_matrix - sd_score_matrix), 
                                  1, rank, ties.method = "max") - 1)
  mean_2_sd <- t(apply(sign_for_score * (mean_score_matrix - (2 * sd_score_matrix)), 
                                  1, rank, ties.method = "max") - 1)
  
  #the worst rank
  worst_rank <- pmax(mean_score, mean_1_sd, mean_2_sd)
  
  #combine regions with ranks
  result <- list(mean_score, mean_1_sd, mean_2_sd, worst_rank)
  names(result) <- c(paste("number of cell types with score", direction, "than average"), 
                     paste("number of cell types with score", direction, "than (average", sign_for_label, "SD)"),
                     paste("number of cell types with score", direction, "than (average", sign_for_label, "2*SD)"),
                     "worst rank")
  return(result)
}

#generate a list of cell type signatures
generate_cell_type_signatures <- function(unique_biosources, regions, ranks, 
                                          min.num.of.regions = 500,
                                          max.num.of.regions = NULL){
  if(!is.null(max.num.of.regions)) 
    if(min.num.of.regions > max.num.of.regions) 
      stop("max.num.of.regions needs to be equal to or larger than min.num.of.regions.")
  
  foreach(biosource = unique_biosources,
          .final = function(x) setNames(x, unique_biosources)) %do% {
            regions_ranks <- cbind(regions, ranks)
            result <- regions_ranks %>% 
              select(CHROMOSOME, START, END, celltypes_scoring_better = UQ(biosource)) %>%   
              top_n(min.num.of.regions, -celltypes_scoring_better)  
            if(!is.null(max.num.of.regions)){
              #return a random sample of regions
              result <- result %>% sample_n(max.num.of.regions)
            } 
            return(result)
          }
}
