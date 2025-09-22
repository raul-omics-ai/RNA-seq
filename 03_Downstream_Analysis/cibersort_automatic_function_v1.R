########## 17/09/2025 ##########
# ================================================ #
# ==== AUTOMATIC DECONVOLUTION WITH CIBERSORT ==== #
# ================================================ #

#' Automatic CIBERSORT Deconvolution
#'
#' Performs immune deconvolution using the CIBERSORT algorithm with integrated visualizations.
#'
#' @param count_matrix Gene expression matrix with genes in rows and samples in columns.
#' @param meta Metadata with sample information. Must have the same samples as count_matrix.
#' @param signature_matrix Signature matrix for deconvolution. Must have the genes in rownames.
#' @param where_to_save Directory to save results (default: NA uses working directory).
#' @param title Analysis title (default: "CIBERSORT_deconvolution_Result").
#' @param specie Data species ("mouse" or "human").
#' @param transform2human If TRUE, converts mouse genes to human (default: TRUE).
#' @param grouping_variable Variable in meta for grouping in boxplots.
#' @param P Number of permutations for CIBERSORT (default: 100).
#' @param ordered If TRUE, sorts samples in heatmap (default: FALSE).
#' @param sorting_variable Variable(s) for sample sorting.
#'
#' @return List with deconvolution results and visualizations. Saves files to disk.
#'
#' @examples
#' cibersort_automatic_function(
#'   count_matrix = counts,
#'   meta = metadata,
#'   signature_matrix = sig_matrix,
#'   title = "MyAnalysis",
#'   grouping_variable = "condition"
#' )
cibersort_automatic_function <- function(count_matrix, 
                                        meta, 
                                        signature_matrix, 
                                        where_to_save = NA,
                                        title = "CIBERSORT_deconvolution_Result",
                                        specie = "mouse",
                                        transform2human = TRUE,
                                        grouping_variable = NULL,
                                        P = 100,
                                        ordered = FALSE,
                                        sorting_variable = NULL){
  # ============================================================================== #
  # ==== STEP 0: LOADING PACKAGES, CUSTOM FUNCTIONS AND SET WORKING DIRECTORY ==== 
  # ============================================================================== #
  source("~/Documentos/09_scripts_R/print_centered_note_v1.R")
  source("~/Documentos/09_scripts_R/Automate_Saving_ggplots.R")
  source("~/Documentos/09_scripts_R/create_sequential_dir.R")
  source("~/Documentos/09_scripts_R/automate_saving_dataframes_xlsx_format.R")
  source("~/Documentos/09_scripts_R/CIBERSORT.R")
  print_centered_note('LOADING PACKAGES ')
  list_of_packages = c("openxlsx", "dplyr", "patchwork", "ggplot2", "tidyr", "immunedeconv", 
                       "AnnotationDbi", "homologene", "pheatmap", "reshape2", "viridisLite", "viridis", 
                       "ggsignif", "patchwork", "rlang")
  
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
  
  print("Creating the report")
  report <- data.frame(matrix(nrow = 0, ncol = 1))
  colnames(report) <- "Report"
  report[nrow(report)+1, 1] <- paste0("# Genes in Orignial Count Matrix: ", nrow(count_matrix), " genes")
  report[nrow(report)+1, 1] <- paste0("# Unique SYMBOLS in Count Matrix: ", length(unique(rownames(count_matrix))), " genes")
  
  # =================================== #
  # ==== STEP 1: PREPARE THE INPUT ==== 
  # =================================== #
  print_centered_note(toupper("Checking the checkpoints "))
  # Check Point 1 ==== #
  # All samples are in both datasets? Are them sorted?
  if(!all(colnames(count_matrix) %in% rownames(meta))){
    stop(paste0("Please review if all samples are in both datasets"))
  }#if key for sorting metadata and count_matrix

  # Checkpoint 2 === #
  # All SYMBOLs in count matrix are unique?
  if(!length(unique(rownames(count_matrix))) == nrow(count_matrix)){
    # Dejo el código comentado por si acaso por si quiero automatizar la transformación de ENSEMBL a SYMBOL para ratón.
    #c1 <- count_matrix %>%
    #  tibble::rownames_to_column(var = "ENSEMBL") %>%
    #  mutate("SYMBOL" = mapIds(org.Mm.eg.db, keys = ENSEMBL, column = "SYMBOL", keytype = "ENSEMBL")) %>%
    #  relocate(ENSEMBL, SYMBOL) %>%
    #  tibble::column_to_rownames(var = "SYMBOL") %>%
    #  dplyr::select(-ENSEMBL)
    
    stop(paste0("There are ", nrow(count_matrix), " genes in the dataset for ", length(unique(rownames(count_matrix))), "\nReview count matrix to remove duplicates"))
  }# if key for checking samples in both datasets 
  
  # Checkpoint 3 === #
  # The signature matrix for deconvolution has rownames?
  tiene_rownames <- function(df) {
    rn <- rownames(df)
    
    # Caso 1: No hay rownames definidos
    if (is.null(rn)) return(FALSE)
    
    # Caso 2: Rownames son iguales a la secuencia por defecto (1, 2, 3, ...)
    if (identical(rn, as.character(seq_len(nrow(df))))) {
      return(FALSE)
    }
    
    # Caso 3: Rownames personalizados
    return(TRUE)
  }# function tiene_rownames key for check rownames 
  
  if (!tiene_rownames(signature_matrix)) {
    stop("The 'signature_matrix' dataframe does not have rownames (Gene Symbols).")
  }
  
  if (!tiene_rownames(count_matrix)) {
    stop("The 'count_matrix' dataframe does not have rownames (Gene Symbols).")
  }
  
  # Checkpoint 4
  # La variable está dentro de metadatos
  # Verificar que las variables de ordenamiento existan en meta
  if (!all(sorting_variable %in% colnames(meta))) {
    stop("Sorting variable doesn't found in metadata")
  }
  
  print_centered_note(toupper("Transforming mouse genes to human genes"))
  if(transform2human & specie == "mouse"){
    mouse2human_orthologues <- homologene::mouse2human(rownames(count_matrix)) 
    save_dataframe(mouse2human_orthologues, title = "Human_Orthologues_DataFrame", folder = dataprocessed)
  
    report[nrow(report)+1, 1] <- paste0("# Human Orthologues mapped: ", nrow(mouse2human_orthologues), " genes")
    report[nrow(report)+1, 1] <- paste0("# Unique Human Orthologues mapped: ", length(unique(mouse2human_orthologues$humanGene)), " genes")
    report[nrow(report)+1, 1] <- paste0("Duplicated Human Orthologues: ", 
                                        paste(mouse2human_orthologues$humanGene[duplicated(mouse2human_orthologues$humanGene)], 
                                              collapse = ", "))
    report[nrow(report)+1, 1] <- paste0("Duplicated Genes Management: Keep the gene with higher expression")
    
    orthologues_count_mat <- mouse2human_orthologues %>% 
      dplyr::select(-mouseID, - humanID) %>%
      left_join(count_matrix %>%
                  tibble::rownames_to_column("mouseGene"),
                by = c("mouseGene" = "mouseGene")) %>%
      dplyr::select(-mouseGene)
    
    orthologues_count_mat <- orthologues_count_mat %>%
      group_by(humanGene) %>%
      slice_max(order_by = rowSums(across(where(is.numeric))), n = 1) %>%
      ungroup() %>%
      tibble::column_to_rownames(var = "humanGene")
    
    save_dataframe(orthologues_count_mat, title = "Mouse2Human_count_matrix_clean", folder = dataprocessed, rowNames = T)
    
  }# key to transform mouse genes to human genes
  
  # ============================================== #
  # ==== STEP 2: DECONVOLUTION WITH CIBERSORT ==== 
  # ============================================== # 
  print_centered_note(toupper("Running Deconvolution"))
  res_ciber <- cibersort(sig_matrix = signature_matrix, 
                         mixture_file = orthologues_count_mat, 
                         perm = P, #The number of permutation
                         QN = F # Quantile normalization of input mixture, default = FALSE for RNA-Seq data
                         )
  print("Saving CIBERSORT result")
  save_dataframe(as.data.frame(res_ciber), title = "CIBERSORT_Result", folder = dataprocessed, rowNames = T)
  
  print("Saving report")
  save_dataframe(report, title = "Report", folder = dataprocessed)
  # ================================ #
  # ==== STEP 3: VISUALIZATIONS ==== 
  # ================================ #
  print_centered_note(toupper("Creating the visualizations"))
  # 1) Merge metadata with cibersort output
  plot_matrix <- reshape2::melt(as.matrix(res_ciber))
  colnames(plot_matrix) <- c("SampleID","CellType", "Score")
  
  plot_matrix <- merge(plot_matrix, meta %>% 
                     tibble::rownames_to_column("SampleID"), 
                   by = c("SampleID" = "SampleID"))
  
  print("Saving the visualization dataframe")
  save_dataframe(plot_matrix, title = "Metadata_Cibersort_Result_Merge", folder = dataprocessed)

  # 2) Heatmap
  print("Creating the heatmap")
  hmap_mat <- res_ciber %>%
    as.data.frame() %>%
    dplyr::select(-c("P-value", "Correlation", "RMSE"))
  
  if(ordered){
    samples_ordered <- meta %>%
      tibble::rownames_to_column(var = "SampleID") %>%
      arrange(across(all_of(sorting_variable))) %>%
      pull(SampleID)
      
    X_ordered <- t(hmap_mat)
    X_ordered <- X_ordered[, samples_ordered]
    
    annotation_col <- meta
    
    hmap <- pheatmap(X_ordered,
             color = viridis(100, option = "D"),
             annotation_col = annotation_col,
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             scale = "none")
    save_ggplot(hmap, title = "Immune_Deconvolution_Heatmap_Sorted", folder = figures, width = 6000, height = 3000)
    
  } else{
    annotation_col <- meta
    
    hmap <- pheatmap(t(hmap_mat),
                     color = viridis(100, option = "D"),  # "D" = viridis clásico
                     annotation_col = annotation_col,
                     cluster_rows = FALSE,
                     cluster_cols = TRUE,
                     scale = "none")
    
    save_ggplot(hmap, title = "Immune_Deconvolution_Heatmap_Clustered", folder = figures, width = 6000, height = 3000)
    
  } # if else key for heatmap ordered samples
  
  # 3) Stacked plot
  print("Creating the Stacked plot")
  staks_plots <- lapply(sorting_variable, function(var) {
    plot_matrix %>%
      dplyr::filter(!CellType %in% c("RMSE", "P-value", "Correlation")) %>%
      arrange(across(all_of(sorting_variable))) %>%
      ggplot(aes(x = SampleID, y = Score, fill = CellType)) +
      geom_bar(position = "stack", stat = "identity") +
      facet_wrap(vars(!!sym(var)), scales = "free_x") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      theme(legend.position = "bottom") +
      labs(y = "Proportion", x = NULL)
  })
  
  for(plot in seq(1:length(staks_plots))){
    title = paste0("Stacked_Plot_", sorting_variable[plot])
    save_ggplot(plot = staks_plots[[plot]], title = title, folder = figures, width = 6000, height = 3000)
  }
  
  # 4) Individual plots by
  print("Creating the Boxplots")
  make_boxplots <- function(data, grouping_variable) {
    boxplots <- lapply(unique(data$CellType), function(ct) {
      data %>%
        dplyr::filter(CellType == ct,
               !CellType %in% c("RMSE", "P-value", "Correlation")) %>%
        ggplot(aes(x = !!sym(grouping_variable), y = Score)) +
        geom_boxplot(aes(colour = !!sym(grouping_variable))) +
        theme(legend.position = "none") +
        labs(title = ct, y = "Score", x = NULL) +
        geom_signif(comparisons = list(c("C57", "rd16")),
                    map_signif_level = FALSE)
    })
    return(boxplots)
  }
  
  boxplots <- make_boxplots(plot_matrix, grouping_variable = grouping_variable)
  
  boxplot_folder <- create_sequential_dir(path = figures, name = "Individual_Boxplots")
  
  lapply(seq(1:length(unique(plot_matrix$CellType))), function(ct){
    cell_type <- unique(plot_matrix$CellType)[ct]
    title = paste0("Boxplot_CellType_", cell_type)
    
    save_ggplot(boxplots[[ct]], title = title, folder = boxplot_folder, width = 2000, height = 1500)
    
  })

  combined <- wrap_plots(boxplots, ncol = 5)  # ajusta ncol como quieras
  save_ggplot(combined, title = paste0("Boxplot_CellTypes_by_", grouping_variable), folder = figures,
              width = 6000, height = 4000)
  
  print_centered_note(toupper("End of the script"))
  # Al final de la función, retornar resultados útiles
  invisible(list(
    cibersort_results = res_ciber,
    plots = list(heatmap = hmap, boxplots = boxplots),
    report = report
  ))
}# Main function key
