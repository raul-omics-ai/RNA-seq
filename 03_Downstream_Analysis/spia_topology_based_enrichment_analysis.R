########## 27/05/2025 ##########
#' ============================================================================ #
#' ==== AUTOMATE IMPACT ANALYSIS FOR TOPOLOGICAL-BASED ENRICHMENT ANALYSIS ==== #
#' ============================================================================ #
#'
#' Performs automated Signaling Pathway Impact Analysis (SPIA) to identify
#' significantly impacted KEGG pathways using topology-aware enrichment analysis.
#' Includes visualization, pathway mapping, and similarity analysis between pathways.
#'
#' @param dea_dataframe A data frame containing differential expression analysis results.
#'                    Must contain logFC values and adjusted p-values.
#' @param pval_adjust_colname Character string specifying the column name for adjusted p-values.
#'                          Default is "FDR".
#' @param keytype Character string specifying the gene identifier type in the input data.
#'              Default is "ENTREZID".
#' @param gene_column Character string specifying the column name containing gene identifiers.
#'                  Default is "ENTREZID".
#' @param where_to_save Character string specifying the directory path for saving results.
#'                    If NULL, uses current working directory. Default is NULL.
#' @param title Character string used for naming output files and directories.
#'            Default is "Topology_Based_Enrichment_Analysis".
#' @param org Character string specifying the organism. Options: "mmu" (mouse) or "hsa" (human).
#'          Default is "mmu".
#' @param significance Numeric value specifying the significance threshold for adjusted p-values.
#'                   Default is 0.05.
#' @param kegg_db_path Character string specifying the path to pre-existing KEGG database CSV file.
#'                   Default is NULL.
#' @param new_pathway_txt Logical indicating whether to create a new KEGG pathway database from text file.
#'                      Default is FALSE.
#' @param path_txt Character string specifying the path to KEGG pathway text file when new_pathway_txt = TRUE.
#'               Default is NULL.
#'
#' @return Invisible NULL. The function generates:
#'   - SPIA results in Excel format
#'   - Two-way evidence plot
#'   - Pathview pathway visualizations
#'   - Jaccard index dendrogram for significant pathways
#'   - Processed data files in organized directory structure
#'
#' @details
#' This function performs the following steps:
#' \enumerate{
#'   \item Loads required packages and sets up directory structure
#'   \item Converts gene identifiers to ENTREZID if necessary
#'   \item Runs SPIA analysis using pathway topology
#'   \item Generates visualizations (evidence plot, pathway maps)
#'   \item Calculates Jaccard similarity between significant pathways
#'   \item Creates hierarchical clustering dendrogram of pathways
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage with mouse data
#' SPIA_automatic_enrichment(dea_dataframe = my_deg_results, 
#'                          org = "mmu")
#'
#' # Human data with custom parameters
#' SPIA_automatic_enrichment(dea_dataframe = human_deg_results,
#'                          org = "hsa",
#'                          significance = 0.01,
#'                          title = "My_SPIA_Analysis")
#' }
#'
#' @export
#' @importFrom openxlsx write.xlsx
#' @importFrom dplyr filter mutate select left_join group_by summarise
#' @importFrom tidyr separate_rows
#' @importFrom SPIA spia plotP
#' @importFrom pathview pathview
#' @importFrom vegan vegdist
#' @importFrom stats hclust dist
#' @importFrom AnnotationDbi mapIds
#' @importFrom org.Mm.eg.db org.Mm.eg.db
#' @importFrom org.Hs.eg.db org.Hs.eg.db
SPIA_automatic_enrichment <- function(dea_dataframe,
                                      pval_adjust_colname = "FDR",
                                      keytype = "ENTREZID",
                                      gene_column = "ENTREZID",
                                      where_to_save = NULL,
                                      title = "Topology_Based_Enrichment_Analysis",
                                      org = "mmu",
                                      significance = 0.05,
                                      kegg_db_path = NULL,
                                      new_pathway_txt = FALSE,
                                      path_txt = NULL){
  
  # ============================================================================== #
  # ==== STEP 0: LOADING PACKAGES, CUSTOM FUNCTIONS AND SET WORKING DIRECTORY ==== #
  # ============================================================================== #
  source("~/Documentos/09_scripts_R/print_centered_note_v1.R")
  source("~/Documentos/09_scripts_R/Automate_Saving_ggplots.R")
  source("~/Documentos/09_scripts_R/create_sequential_dir.R")
  source("~/Documentos/09_scripts_R/automate_saving_dataframes_xlsx_format.R")
  print_centered_note('LOADING PACKAGES ')
  list_of_packages = c("openxlsx", "edgeR", "dplyr", "AnnotationDbi", "clusterProfiler",
                       "patchwork", "ggplot2", "tidyr", "SPIA", "pathview", "vegan",
                       "ggdendro")
  
  new_packages = list_of_packages[!(list_of_packages %in% installed.packages())]
  if(length(new_packages) > 0){install.packages(new_packages)}
  
  invisible(lapply(list_of_packages, FUN = library, character.only = T))
  rm(list_of_packages, new_packages)
  print('Packages loaded successfully.')
  
  print_centered_note(toupper("Setting up working directory "))
  where_to_save <- ifelse(is.null(where_to_save), getwd(), where_to_save)
  output_dir <- create_sequential_dir(path = where_to_save, name = title)
  figures <- create_sequential_dir(path = output_dir, name = "Figures")
  dataprocessed <- create_sequential_dir(path = output_dir, name = "Report")
  wd <- getwd()
  
  # ============================================= #
  # ==== STEP 1: Change KEYTYPE to ENTREZID  ==== #
  # ============================================= #
  
  if(keytype != "ENTREZID"){
    print_centered_note(paste0("Transforming ", keytype, "keytype to ENTREZID"))
    if(org == "mmu"){
      print("Selected specie: Mus musculus")
      library(org.Mm.eg.db)
      org.db <- org.Mm.eg.db
    }
    if(org == "hsa"){
      print("Selected specie: Homo sapiens")
      library(org.Hs.eg.db)
      org.db <- org.Hs.eg.db
    }
    
    dea_dataframe$ENTREZID <- AnnotationDbi::mapIds(org.db, 
                                                 keys = dea_dataframe[[gene_column]],
                                                 column = "ENTREZID",
                                                 keytype = toupper(keytype))
  }# if for entrezid conversion
  
  # =========================================== #
  # ==== STEP 2: SETTING UP SPIA FUNCTION  ==== #
  # =========================================== #
  print_centered_note(toupper("Setting up SPIA function "))

  # 1. Creating a vector with the logFC for DE genes
  deg_df <-dea_dataframe[dea_dataframe[[pval_adjust_colname]] < significance,]
  deg_df<-deg_df[!duplicated(deg_df$ENTREZID),]
  
  dea_logfc = deg_df$logFC
  names(dea_logfc)<-as.vector(deg_df$ENTREZID)
  
  # 2. Creating a vector with the background set of all the genes represented on the platform.
  background_entrez <- dea_dataframe$ENTREZID
  
  # 3. Call SPIA function
  spia_result <- spia(de=dea_logfc, 
                      all=background_entrez, 
                      organism= org, 
                      nB = 2000)
  
  print("Saving SPIA dataframe")
  save_dataframe(spia_result, title = paste0("SPIA_Result_", title), 
                 folder = dataprocessed)
  
  # ================================ #
  # ==== STEP 4: Visualizations ==== #
  # ================================ #
  print_centered_note(toupper("Creating SPIA Visualizatons"))
  # 1. Two-way evidence plot
  print("1. Two-way evidence plot")
  spia_result2 <- spia_result
  spia_result2$ID <- spia_result2$Name 
  spia_result2 <- na.omit(spia_result2) # Removing NA values
  threshold <- ifelse(min(spia_result2$pGFdr) > .05, round(min(spia_result2$pGFdr), digits = 3) + 0.001, 0.05)
  
  png(filename = file.path(figures, paste0("01_SPIA_Two-way_Evidence_Plot_thr_", threshold, ".png")), 
      width = 3000, height = 2000, res = 300, units = "px")
  plotP(spia_result2, threshold=threshold)
  dev.off()
  
  # 2. Pathview
  print("2. Pathview KEGG Pathways")
  pathwvie_folder <- create_sequential_dir(figures, name = "Pathview_Pathways")
  setwd(pathwvie_folder)
  
  tmp = sapply(spia_result[spia_result$pGFdr < threshold, "ID"], 
               function(pid) pathview (gene.data = dea_logfc, pathway.id = pid, 
                                       species = org, same.layer = F))
  
  setwd(wd)
  
  # ========================================================== #
  # ==== STEP 5: Obtain all genes from each KEGG pathway  ==== #
  # ========================================================== #
  # 1. Get Kegg Database with all genes in the database
  if(new_pathway_txt){
    print_centered_note("Creating the new TXT for KEGG")
    source("~/Documentos/09_scripts_R/crear_kegg_database_V1.R")
    kegg_db <- procesar_kegg_pathways(path_txt = path_txt, 
                                            output_csv = file.path(dataprocessed, paste0("02_Kegg_db_", org, "_",Sys.Date(),".csv")),
                                            organismo = org)
    
    } else{
      print_centered_note("Loading CSV for KEGG pathways")
      kegg_db <- read.csv(kegg_db_path)
    }# ifelse kegg database key
  
  # 2. Merge SPIA pathways with database to get all database genes
  spia_result_merge <- spia_result %>%
    mutate(ID = paste0(org, ID)) %>%
    left_join(kegg_db,
              by=c("ID" = "ID"))
  save_dataframe(spia_result_merge, title = paste0("SPIA_Result_All_Genes", title), 
                 folder = dataprocessed)
  # ==================================================================== #
  # ==== STEP 6: Get the Jaccad Index from all significant pathways ==== #
  # ==================================================================== #
  print_centered_note(toupper("Calculating Jaccard Index"))
  
  spia_result_merge.short <- spia_result_merge %>%
    filter(pGFdr < threshold) %>%
    select(Name, ID, Pathway, SYMBOL)

  if(nrow(spia_result_merge.short) > 2){
    gene_lists <- spia_result_merge.short %>%
    tidyr::separate_rows(SYMBOL, sep = ",\\s*") %>%  # Divide por coma y espacio
    group_by(Pathway) %>%
    summarise(Genes = list(unique(SYMBOL)))     # Lista única de genes por ruta
  
  # Extraer todas las rutas y genes únicos
  all_genes <- unique(unlist(gene_lists$Genes))
  
  # Construir matriz vacía
  inc_matrix <- sapply(gene_lists$Genes, function(gl) all_genes %in% gl) %>%
    t()
  colnames(inc_matrix) <- all_genes
  if(org == "mmu"){
    rownames(inc_matrix) <- gsub(" - Mus musculus \\(house mouse\\)$", "",gene_lists$Pathway)
  }
  if(org == "hsa"){
    rownames(inc_matrix) <- gsub(" - Homosapiens \\(human\\)$", "",gene_lists$Pathway)
  }
  inc_matrix <- as.data.frame(inc_matrix)
  
  # Cálculo de la distancia de Jaccard dist devuelve un objeto de clase 'dist'
  jacc_dist <- vegan::vegdist(inc_matrix, method = "jaccard") # Esto es distancia no similaridad (1-jaccard index)
  jacc_mat  <- as.matrix(jacc_dist)                  # convertir a matriz cuadrada
  
  # Clustering
  print_centered_note(toupper("Saving Jaccard Index Dendrogram"))
  hc <- hclust(jacc_dist, method = "single")
  png(filename = file.path(figures,"02_Jaccard_Index_Single_Hierarchical_Clustering.png"), 
      res =300, width = 3250, height = 2500)
  print(plot(hc, main = "Dendrograma de rutas basado en Jaccard", xlab = "", sub = ""))
  dev.off()
  
  print_centered_note(toupper("End of the script"))
  } else {print_centered_note(toupper("End of the script"))} #ifelse statement key
  
}# Key of the function
