#hypomethylated regions

compute_cell_type_hypo_meth_scores <- function(mean_score_matrix, sd_score_matrix){
  
  #the three metrices described above
  mean_meth <- t(apply(mean_score_matrix, 
                       1, rank, ties.method = "max") - 1)
  mean_meth_minus_1_sd <- t(apply((mean_score_matrix - sd_score_matrix), 
                                  1, rank, ties.method = "max") - 1)
  mean_meth_minus_2_sd <- t(apply((mean_score_matrix - (2 * sd_score_matrix)), 
                                  1, rank, ties.method = "max") - 1)
  
  #the worst rank
  worst_rank <- pmax(mean_meth, mean_meth_minus_1_sd, mean_meth_minus_2_sd)
  
  #combine regions with ranks
  result <- list(mean_meth, mean_meth_minus_1_sd, mean_meth_minus_2_sd, worst_rank)
  names(result) <- c("number of cell types with lower average methylation", 
                     "number of cell types with lower (average methylation - SD)",
                     "number of cell types with lower (average methylation - 2*SD)",
                     "worst rank")
  return(result)
}

#hypermethylated regions 

convert_hypo_to_hyper <- function(hypo_meth_ranks){
  #reverse ranking for turning hypo into hypermethylation ranks
  hyper_meth_ranks <- lapply(hypo_meth_ranks, 
                             function(rank_matrix) { max(rank_matrix[1,]) - rank_matrix })
  
  names(hyper_meth_ranks) <- str_replace(names(hyper_meth_ranks), "lower", "higher")
  return(hyper_meth_ranks)
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
