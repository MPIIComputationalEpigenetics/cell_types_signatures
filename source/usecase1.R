usecase_1_request_score_matrices <- function(selected_experiments, column_name, 
                                              reference_genome = "GRCh38",
                                              experiments_meta)
{
  chromosomes <- grep("_", deepblue_extract_ids(deepblue_chromosomes(reference_genome)), invert = TRUE, value = TRUE)

  # We also need to define the name of the column we need in the WIG files
  experiments_columns <- deepblue_select_column(selected_experiments, column_name)
  
  request_ids <- foreach(chromosome = chromosomes,
                         .final = function(x) setNames(x, chromosomes)) %do% {
                           
                           #Select BLUEPRINT Ensembl regulatory build regions
                           #deepblue_list_annotations(genome = reference_genome) for all annotations
                           ensembl_reg_build <- deepblue_select_annotations(annotation_name = "Blueprint Ensembl Regulatory Build",
                                                                            genome = reference_genome,
                                                                            chromosome = chromosome)
                           
                           #alternatively if you want to use 1kb tiling 
                           #deepblue_tiling_regions(size = 1000, genome = reference_genome, chromosome = chromosome)
                           
                           # request the score matrices we want
                           experiments_mean_request_id <- 
                             deepblue_score_matrix(experiments_columns = experiments_columns,
                                                   aggregation_function = "mean",
                                                   aggregation_regions_id = ensembl_reg_build)
                           
                           experiments_sd_request_id <- 
                             deepblue_score_matrix(experiments_columns = experiments_columns,
                                                   aggregation_function = "sd",
                                                   aggregation_regions_id = ensembl_reg_build)
                           
                           #collect all request ids in pairs
                           list(experiments_mean_request_id, experiments_sd_request_id)
                         }
  return(request_ids)
}

check_status <- function(request_ids){
  table(unlist(lapply(deepblue_info(unlist(request_ids)), function(x) x$state)))
}

usecase_1_download_score_matrices <- function(request_ids,
                                             experiments_meta)
{
  #check status of all requests
  state <- check_status(request_ids)
  
  while(names(state) != c("done")){
    print(state)
    message("At least one request is not completed. Waiting...")
    Sys.sleep(10)
  }
  
  # once finished, download the score matrices and merge them
  download_and_combine_score_matrices <- function(request_ids, i){
    # NOTE bind_rows will match columns by name such that we don't have to worry about column order
    foreach(pair_of_request_ids = request_ids, .combine = bind_rows) %do%
    {
      request_id <- pair_of_request_ids[[i]]
      if(deepblue_info(request_id)$state != "done") 
        stop(paste(request_id, "is not reported as done. Please check its (error) state."))
      deepblue_download_request_data(request_id)
    }
  }
  experiments_mean_score_matrix <- download_and_combine_score_matrices(request_ids, 1)
  experiments_sd_score_matrix <- download_and_combine_score_matrices(request_ids, 2)
  
  # extract regions coordinates
  regions <- experiments_mean_score_matrix[,1:3]
  
  # convert the score matrix to a numeric matrix (omitting chr, start, end columns)
  experiments_mean_score_numeric <- as.matrix(experiments_mean_score_matrix[,4:ncol(experiments_mean_score_matrix)])
  experiments_sd_score_numeric <- as.matrix(experiments_sd_score_matrix[,4:ncol(experiments_sd_score_matrix)])
  
  # the column order is not guaranteed to be identical in these two matrices, thus we reorder the second one
  experiments_sd_score_numeric <- experiments_sd_score_numeric[,colnames(experiments_mean_score_numeric)]
  
  # group by cell type
  unique_biosources <- unique(experiments_meta$user_celltype)
  
  # multiple columns in the matrices correspond to the same cell type. here we merge them by mean, omitting NAs.
  aggregate_by_cell_type <- function(score_matrix){
    foreach(biosource = unique_biosources, 
            .combine = cbind,
            .final = function(x) { 
              colnames(x) <- unique_biosources; 
              return(x)}) %do% {
                experiments_with_current_biosource <- experiments_meta[which(experiments_meta$user_celltype == biosource), "name"]
                subset_of_matrix <- score_matrix[,experiments_with_current_biosource]
                if(is.null(ncol(subset_of_matrix))) return(subset_of_matrix)
                else return(rowMeans(subset_of_matrix, na.rm = TRUE))
              }
  }
  
  biosources_mean_score_matrix <- aggregate_by_cell_type(experiments_mean_score_numeric)
  
  #we replace NAs with zeros
  biosources_mean_score_matrix[is.na(biosources_mean_score_matrix)] <- 0
  
  #we identify regions without methylation signal in any of the cell types
  signal_regions <- which(rowSums(biosources_mean_score_matrix) > 0)
  
  #and keep only those
  biosources_mean_score_matrix <- biosources_mean_score_matrix[signal_regions,]
  
  #follow the same strategy for the SD matrix
  biosources_mean_sd_matrix <- aggregate_by_cell_type(experiments_sd_score_numeric)
  biosources_mean_sd_matrix[is.na(biosources_mean_sd_matrix)] <- 0
  biosources_mean_sd_matrix <- biosources_mean_sd_matrix[signal_regions,]
  
  #make sure to also filter the regions coordinates data table
  regions <- regions[signal_regions,]
  
  return(list(mean_df = experiments_mean_score_matrix, 
              mean_matrix = biosources_mean_score_matrix, 
              sd_df = experiments_sd_score_matrix,
              sd_matrix = biosources_mean_sd_matrix,
              regions = regions))
}

