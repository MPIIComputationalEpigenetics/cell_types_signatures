go_enrichment <- function(signatures, unique_biosources, file){
  
  foreach(biosource = unique_biosources) %do% {
    signature_string <- paste(apply(hypo_signatures[[biosource]][,-4], #keep only chr start end
                                    1, paste, collapse = "\t", sep = ""), #concatenate to string
                              collapse = "\n")
    
    qid <- deepblue_input_regions(signature_string, genome = "GRCh38") #submit to DeepBlue to get a query id
    
    enrichment <- deepblue_enrich_regions_go_terms(qid, gene_model = "gencode v23") #perform enrichment
    
    enrichment_result <- deepblue_download_request_data(enrichment) #download results
    enrichment_result_df <- enrichment_result$enrichment$go_terms
    enrichment_result_df$p_adjust <- p.adjust(enrichment_result_df$p_value)
    enrichment_result_df <- enrichment_result_df %>% arrange(p_adjust) %>% filter(p_value < 0.05)
    
    #write to disk
    if(biosource == "effector memory CD8-positive, alpha-beta T cell, terminally differentiated")
      biosource <- "term. diff. eff. mem. CD8-positive, alpha-beta T cell"
    write.xlsx2(enrichment_result_df, file = file, sheetName = biosource, append = TRUE)
  }

}