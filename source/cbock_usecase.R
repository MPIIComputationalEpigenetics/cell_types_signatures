# Load all package dependencies. install if missing!
library(DeepBlueR)
library(ggplot2)
library(dplyr)
library(foreach)
library(xlsx)
library(matrixStats)
library(stringr)
library(data.table)
library(gplots)
library(RColorBrewer)
library(doParallel)
library(testthat)

#register parallel processing, enable for speeding up signatures computation
#num_of_cores <- 37 #adjust as needed
#cl <- parallel::makeCluster(37, outfile = "")
#registerDoParallel(cl)

# Goal: establish a comprehensive list of marker regions that epigenetically identify a cell type of interest (e.g. for the development of cell type specific biomarkers)

### Use Case 4: Obtaining region-specific DNA methylation data for a custom selection of cell types ###

# 1. Export a summary table that lists metadata for all samples that have sufficient DNA methylation data to be included in the analysis
# 2. The user manually adds a new column "custom_cell_type" to this table (in Excel), defining which biologically related samples should be aggregated into which cell types ("NA" for cell samples that are to be ignored)
# 3. The annotation table is imported into DeepBlue, and the use case continues with step 2 of Use Case 1.

# first we identify suitable biosource terms
biosource_blood <- deepblue_get_biosource_related("hematopoietic cell")

# we identify experimental files in BLUEPRINT for those biosource terms
selected_experiments <- deepblue_list_experiments(genome = "GRCh38",
                                                  biosource = deepblue_extract_names(biosource_blood),
                                                  project = "BLUEPRINT Epigenome",
                                                  epigenetic_mark = "DNA methylation")

# keep only 'call' files
selected_experiments <- selected_experiments[grepl(pattern = "call\\.", selected_experiments$name),]

# keep only CpG methylation calls
selected_experiments <- selected_experiments[grepl(pattern = "CPG_methylation", selected_experiments$name),]

# download meta data
experiments_meta <- deepblue_meta_data_to_table(deepblue_extract_ids(selected_experiments))
experiments_meta$user_celltype <- experiments_meta$biosource_name

# write to XLSX
write.xlsx(experiments_meta, file = "experiments_metadata.xlsx")

# wait for user to edit the file and add custom cell type names
readline("Press enter when you are finished editing the metadata file experiments_metadata.xlsx")

# read edited meta data
experiments_meta <- read.xlsx(file = "experiments_metadata.xlsx", sheetIndex = 1)

# plot cell type sample numbers
plot_data <- experiments_meta %>% group_by(user_celltype) %>% summarize(count = n())

ggplot(plot_data, aes(x = user_celltype, y = count)) + geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


### Use Case 1: Obtaining region-specific DNA methylation data for blood cell types ###
# 1. Within cell_type_class=blood, group all samples by cell_type
# 2. Select a region set of interest (e.g. a 1kb tiling of the genome or the Ensembl segmentation map)
# 3. Annotate each region with the mean and standard deviation (alternatively: median, 5th percentile, and 95th percentile) of the region's DNA methylation status in all samples of a given cell_type
# 4. Export the resulting [region] x [cell_type] table for manual filtering in R

# NOTE: not sure what you mean by grouping by cell type. grouping is not supported in DeepBlue and has to be done
# after downloading the data in R. please elaborate how we should group multiple samples, e.g. mean of the median and SD?

# We split the matrix generation by chromosome to make it more efficient.
# First we ask DeepBlue for the chromosome names it uses.
# NOTE: could be you also would want to remove the sex chromosomes here but we left them in

source("source/usecase1.R")

dna_meth_request_ids <- usecase_1_request_score_matrices(selected_experiments = selected_experiments,
                                                         column_name = "VALUE",
                                                         reference_genome = "GRCh38",
                                                         experiments_meta = experiments_meta)

check_status(dna_meth_request_ids)

dna_meth_data <- usecase_1_download_score_matrices(request_ids = dna_meth_request_ids,
                                                   experiments_meta = experiments_meta)

### Use Case 2: Automated selection of cell-type specific biomarkers (this could also be done in R) ###

# 1. Based on the table produced in Use Case 1, add for each cell_type the following calculated measures of cell-type specificity:
# - number of cell types in which the region's mean DNA methylation is lower than the region's DNA methylation in the selected cell_type
# - number of cell types in which the region's mean DNA methylation minus 1x the standard deviation is lower than the region's mean DNA methylation in the selected cell_type
# - number of cell types in which the region's mean DNA methylation minus 2x the standard deviation is lower than the region's mean DNA methylation in the selected cell_type
# 2. Rank the regions based on the worst rank of the three metrics (each one ranked individually and then taking the row-wise maximum)
# 3. Return the top-500 regions that are hypomethylated (lower methylation) in the selected cell_type based on the consensus ranking
# 4. Repeat steps 1 to 3 with a focus on hypermethylated regions (higher methylation) in the selected cell_type

source("source/usecase2.R")

hypo_meth_ranks <- compute_cell_type_scores(dna_meth_data$mean_matrix,
                                            dna_meth_data$sd_matrix)

hyper_meth_ranks <- compute_cell_type_scores(dna_meth_data$mean_matrix,
                                             dna_meth_data$sd_matrix,
                                             invert = TRUE)

#check if the rank computations are correct using a random region and cell type
source("source/testing.R")
testing(dna_meth_data$mean_matrix, dna_meth_data$sd_matrix, hypo_meth_ranks, hyper_meth_ranks)

#get unique biosources for extracting signatures
unique_biosources <- as.character(unique(experiments_meta$user_celltype))

#min.num.of.regions -> include at least 500 regions, include remaining regions with the same rank score
#max.num.of.regions -> if more than x regions have been selected draw a random subset of x
#choose small x for example to make heatmaps feasible.
hyper_signatures <- generate_cell_type_signatures(unique_biosources, dna_meth_data$regions,
                                                  hyper_meth_ranks$`worst rank`,
                                                  min.num.of.regions = 100,
                                                  max.num.of.regions = NULL)
hypo_signatures <- generate_cell_type_signatures(unique_biosources, dna_meth_data$regions,
                                                 hypo_meth_ranks$`worst rank`,
                                                 min.num.of.regions = 100,
                                                 max.num.of.regions = NULL)

#plot size of individual signatures and heatmaps of a subset of the signatures

# plot number of signature regions for each cell type.

hyper_signatures_df <- rbindlist(hyper_signatures, idcol = "cell_type")
hyper_signatures_df$status <- "hyper"

hypo_signatures_df <- rbindlist(hypo_signatures, idcol = "cell_type")
hypo_signatures_df$status <- "hypo"

ggplot(data = rbind(hypo_signatures_df, hyper_signatures_df),
       aes(x = cell_type)) +
  geom_bar(stat = "count") +
  theme_bw() +
  facet_wrap(~status, ncol = 1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# bonus: plot a heatmap with smaller signatures

hyper_signatures_heatmap <- generate_cell_type_signatures(unique_biosources, dna_meth_data$regions,
                                                          hyper_meth_ranks$`worst rank`,
                                                          min.num.of.regions = 100,
                                                          max.num.of.regions = 100)
hypo_signatures_heatmap <- generate_cell_type_signatures(unique_biosources, dna_meth_data$regions,
                                                         hypo_meth_ranks$`worst rank`,
                                                         min.num.of.regions = 100,
                                                         max.num.of.regions = 100)
source("source/heatmap.R")
plot_heatmap(signatures = hypo_signatures_heatmap,
             mean_df_matrix = dna_meth_data$mean_df,
             experiments_meta = experiments_meta,
             signature_score_numeric = signature_score_numeric,
             filename = "DNA_meth_hypo")

plot_heatmap(signatures = hyper_signatures_heatmap,
             mean_df_matrix = dna_meth_data$mean_df,
             experiments_meta = experiments_meta,
             signature_score_numeric = signature_score_numeric,
             filename = "DNA_meth_hyper")

# bonus: tell DeepBlue to perform Gene Ontology enrichment with these signatures and write results to disk
# warning: this will take quite a while thus commented out

source("source/go_enrichment.R")
#go_enrichment(hypo_signatures, unique_biosources, "DNA_meth_hypo_GO_enrichment.xlsx")
#go_enrichment(hyper_signatures, unique_biosources, "DNA_meth_hypo_GO_enrichment.xlsx")

### Use Case 3: Obtaining region-specific open chromatin data (based on signal intensity) ###

# -> this use case is the same as Use Case 1, but based on ChIP-seq H3K27ac signal intensity (highest peak size in the region)

# first we identify suitable biosource terms
biosource_blood <- deepblue_get_biosource_related("hematopoietic cell")

# we identify experimental files in BLUEPRINT for those biosource terms
selected_experiments_h3k27ac <- deepblue_list_experiments(genome = "GRCh38",
                                                          biosource = deepblue_extract_names(biosource_blood),
                                                          project = "BLUEPRINT Epigenome",
                                                          epigenetic_mark = "H3K27ac")

# keep only 'bedgraph' files
selected_experiments_h3k27ac <- selected_experiments_h3k27ac[grepl(pattern = "bedgraph", selected_experiments_h3k27ac$name),]

# download meta data
experiments_meta_h3k27ac <- deepblue_meta_data_to_table(deepblue_extract_ids(selected_experiments_h3k27ac))
experiments_meta_h3k27ac$user_celltype <- experiments_meta_h3k27ac$biosource_name

# write to XLSX
write.xlsx(experiments_meta_h3k27ac, file = "experiments_meta_data_h3k27ac.xlsx")

# wait for user to edit the file and add custom cell type names
readline("Press enter when you are finished editing the metadata file experiments_meta_data_h3k27ac.xlsx")

# read edited meta data
experiments_meta_h3k27ac <- read.xlsx(file = "experiments_meta_data_h3k27ac.xlsx", sheetIndex = 1)

# plot cell type sample numbers
plot_data <- experiments_meta_h3k27ac %>% group_by(user_celltype) %>% summarize(count = n())

ggplot(plot_data, aes(x = user_celltype, y = count)) + geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

h3k27ac_request_ids <- usecase_1_request_score_matrices(selected_experiments = selected_experiments_h3k27ac,
                                                        column_name = "VALUE",
                                                        reference_genome = "GRCh38",
                                                        experiments_meta = experiments_meta_h3k27ac)

check_status(h3k27ac_request_ids)

h3k27ac_data <- usecase_1_download_score_matrices(request_ids = h3k27ac_request_ids,
                                                  experiments_meta = experiments_meta_h3k27ac)

lower_h3k27ac_ranks <- compute_cell_type_scores(h3k27ac_data$mean_matrix,
                                                h3k27ac_data$sd_matrix)

higher_h3k27ac_ranks <- compute_cell_type_scores(h3k27ac_data$mean_matrix,
                                                 h3k27ac_data$sd_matrix,
                                                 invert = TRUE)

unique_biosources_h3k27ac <- as.character(unique(experiments_meta_h3k27ac$user_celltype))

#min.num.of.regions -> include at least 500 regions, include remaining regions with the same rank score
#max.num.of.regions -> if more than 500 regions have been selected draw a random subset of 500
#reduced to 100 regions for plotting a heatmap
h3k27ac_higher_signatures <- generate_cell_type_signatures(unique_biosources_h3k27ac, h3k27ac_data$regions,
                                                           higher_h3k27ac_ranks$`worst rank`,
                                                           min.num.of.regions = 100,
                                                           max.num.of.regions = 100)
h3k27ac_lower_signatures <- generate_cell_type_signatures(unique_biosources_h3k27ac, h3k27ac_data$regions,
                                                          lower_h3k27ac_ranks$`worst rank`,
                                                          min.num.of.regions = 100,
                                                          max.num.of.regions = 100)

#cleanup
if(exists("cl")) stopCluster(cl)