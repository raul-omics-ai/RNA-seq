########## 17/02/2025 ##########

# ==================================================== #
# FULL AUTOMATIC RNA-SEQ ANÁLISIS: FROM QC TO PATHWAYS #
# ==================================================== #

# This function will be eprform all analysis of RNA-seq analysis, from the quality control step from the
# raw count matrix to pathway analysis, doing the differential expression analysis in between. This 
# scripts will be able to select the model to perform the DE analysis (DESeq2 or limma-voom).

# The main dataset to perform some test will be the lactation dataset from https://www.nature.com/articles/ncb3117#Sec25
#full_rnaseq_automate_analysis(count_matrix = raw_counts[-56], 
#                              metadata = metadata, 
#                              experimental_design_variable = "tissue", 
#                              where_to_save = "~/Documentos/", 
#                              title = "Test_V4_full_RNASEQ_pipeline",
#                              model = "DESeq2", 
#                              organism = "mmu", 
#                              pathway_analysis = F, 
#                              omic = "rna", 
#                              topology_analysis = T, 
#                              new_pathway_txt = T, 
#                              path_txt = "~/Documentos/09_scripts_R/kegg_mmu_pathways.txt") 

full_rnaseq_automate_analysis <- function(count_matrix, 
                                          metadata, 
                                          experimental_design_variable = NULL,
                                          min_counts = 10,
                                          where_to_save = NULL,
                                          title = "RNAseq_analysis",
                                          model = c("limma", "DESeq2"),
                                          organism = "hsa",
                                          gene_key_type = "ENSEMBL",
                                          pathway_analysis = TRUE,
                                          omic = c("rna", "microrna"),
                                          topology_analysis = TRUE,
                                          kegg_db_path = NULL, 
                                          new_pathway_txt = FALSE, 
                                          path_txt = NULL
                                          ){
  # ======================================================================== #
  # ==== BLOCK 0: Loading packages and setting up the working directory ==== 
  # ======================================================================== #
  source("~/Documentos/09_scripts_R/print_centered_note_v1.R")
  source("~/Documentos/09_scripts_R/Automate_Saving_ggplots.R")
  source("~/Documentos/09_scripts_R/create_sequential_dir.R")
  source("~/Documentos/09_scripts_R/QC_Count_Matrix_RNAseq_v6.R")
  source("~/Documentos/09_scripts_R/Differential_Expression_Analysis_Limma_voom_v4.R")
  source("~/Documentos/09_scripts_R/Differential_Expression_Analysis_EdgeR_QL_v3.R")
  source("~/Documentos/09_scripts_R/enrichment_analysis_v5.R")
  source("~/Documentos/09_scripts_R/spia_topology_based_enrichment_analysis_V2.R")
  
  print_centered_note(toupper("Setting up working directory "))
  where_to_save <- ifelse(is.null(where_to_save), getwd(), where_to_save)
  output_dir <- create_sequential_dir(path = where_to_save, name = title)
  
  wd <- getwd()
  
  # ====================================== #
  # ==== BLOCK 1: RAW COUNT MATRIX QC ==== 
  # ====================================== #
  print_centered_note("STARTING RAW MATRIX QC STEP")
  
  clean_count_mat <- qc_count_matrix_rnaseq(count_mat = count_matrix, 
                         metadata = metadata, 
                         variable_experimental_design = experimental_design_variable, 
                         where_to_save =  output_dir, 
                         title = "INITIAL_QC", 
                         min_counts = min_counts, 
                         organism = organism, 
                         gene_keytype = gene_key_type,
                         omic = omic)
  
  # ========================================================================== #
  # ==== BLOCK 2: DIFFERENTIAL EXPRESSION ANALYSIS WITH DESEQ2/LIMMA-VOOM ====
  # ========================================================================== #
  print_centered_note("STARING THE DIFERENTIAL EXPRESSION ANALYSIS")
  
  print("Loading clean count matrix")
  #clean_count_mat <- read.xlsx(file.path(output_dir,"01_INITIAL_QC", "QC_Report.xlsx"), sheet = "Filtered_Counts", rowNames = T)
  if(length(model) == 2){
    print_centered_note("PERFORMING THE DE ANALYSIS WITH BOTH MODELS")
    print_centered_note(paste0("DEA MODEL SELECTED: ", model[1]))
    dge_df <- limma_voom_dea(count_matrix = clean_count_mat, 
                   metadata = metadata, 
                   where_to_save = output_dir, 
                   title = paste0("02_DEA_", model[1]), 
                   VOOM_transformation = TRUE, 
                   variable_experimental_design = experimental_design_variable,
                   organism = organism, 
                   gene_keytype = gene_key_type)
    
    print_centered_note(paste0("DEA MODEL SELECTED: ", model[2]))
    dge_df <- edgeR_QL_dea(count_matrix = clean_count_mat, 
                 metadata = metadata, 
                 variable_experimental_design = experimental_design_variable, 
                 where_to_save = output_dir,
                 title = paste0("DEA_", model[2]), 
                 organism = organism, 
                 gene_keytype = gene_key_type,
                 omic = omic)
  }# Both models key
  if(length(model) == 1){
    if(model == "DESeq2"){
      print_centered_note(paste0("DEA MODEL SELECTED: ", model))
      dge_df <- edgeR_QL_dea(count_matrix = clean_count_mat, 
                   metadata = metadata, 
                   variable_experimental_design = experimental_design_variable, 
                   where_to_save = output_dir,
                   title = paste0("DEA_", model), 
                   organism = organism, 
                   gene_keytype = gene_key_type,
                   omic = omic)
    } # DESeq2 DEA key
    if(model == "limma"){
      print_centered_note(paste0("DEA MODEL SELECTED: ", model))
      dge_df <- limma_voom_dea(count_matrix = clean_count_mat, 
                     metadata = metadata, 
                     where_to_save = output_dir, 
                     title = paste0("02_DEA_", model), 
                     VOOM_transformation = TRUE, 
                     variable_experimental_design = experimental_design_variable,
                     organism = organism, 
                     gene_keytype = gene_key_type)
      
    } # limma-voom DEA key
  }# Only one model key
  
  # =================================== #
  # ==== BLOCK 3: PATHWAY ANALYSIS ====
  # =================================== #
  if(pathway_analysis){
    print_centered_note("PERFORMING PATHWAY ANALYSIS")
    print("Loading DE Gene List")
    #dge_df <- read.xlsx(file.path(output_dir, paste0("02_DEA_", model), "DEA_Results.xlsx"), sheet = 6)
    
    if(length(model) == 1){
      if(model == "DESeq2"){
        enrichment_function(dge_df = dge_df, 
                            adj_p_val_colname = "FDR", 
                            where_to_save = output_dir, 
                            title = "03_PATHWAY_ANALYSIS")
      }# Pathway analysis for DESeq2 model key
      if(model == "limma"){
        enrichment_function(dge_df = dge_df, 
                            where_to_save = output_dir, 
                            title = "03_PATHWAY_ANALYSIS")
        
      }# Pathway analysis for limma model key
    }# If key for pathway analysis selection method
  }# Pathway analysis if key
  
  # ===================================================== #
  # ==== BLOCK 4: TOPOLOGY BASED ENRICHMENT ANALYSIS ==== #
  # ===================================================== #
  if(topology_analysis){
    print_centered_note(toupper("Performing Topology Based Enrichment Analysis"))
    SPIA_automatic_enrichment(dea_dataframe = dge_df, pval_adjust_colname = "FDR",
                              keytype = "SYMBOL", gene_column = "Symbol", 
                              where_to_save = output_dir, 
                              title = "TOPOLOGY_BASED_ENRICHMENT_ANALYSIS", 
                              organism = organism, 
                              kegg_db_path = kegg_db_path, new_pathway_txt = new_pathway_txt, 
                              path_txt = path_txt
    )
  }# topology if key
  print_centered_note("END OF COMPLETE RNASEQ ANALYSIS")
}# Main function key

