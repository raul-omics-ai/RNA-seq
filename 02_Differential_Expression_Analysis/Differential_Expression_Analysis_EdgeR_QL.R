########## 28/11/2024 ##########

# EdgeR quasi-likelihood data analysis for RNA and microRNA data.

## I'm going to create a new pipeline for the differential expression analysis, this time instead using limma-voom, I'll use
# edgeR with the Quasi-likelihood pipeline. I read that limma-voom pipeline it's better for large datasets with many samples, that can
# manage variance much better, so is that reason why I want to create another pipeline.

# All pipeline was based in the article from Chen et. al 2020: https://bioconductor.org/packages/release/workflows/vignettes/RnaSeqGeneEdgeRQL/inst/doc/edgeRQL.html

# "The edgeR's quasi-likelihood pipeline uses negative binomial generalized linear model. This method provides stricter error rate control 
# than other negative binomial based pipelines including the tradicional edgeR pielines or DESeq2".
# "The analysis approach in this article can be applied to any RNA-seq study and it's especially appropriate for biological experiments
# with multiple treatment factors and with small numbers of biological replicates."

# NOTA PARA EL FUTURO: Creo que hay bastante ruido en mi anális porque hay bastantes genes identificados como predichos (5.370 son predichos)
# entonces creo que eliminando estos genes, podríamos ver mejor las diferencias a nivel transcriptomico entre los dos grupos.
# Aunque hay 5.370 genes, la verdad que cuando hago el QC, muchos de estos se iran, entonces la verdad que los que se quedan si lso
# buscas en google si que te salen, entonces no se que hacer...
edgeR_QL_dea <- function(count_matrix, 
                         metadata, 
                         variable_experimental_design = NULL,
                         where_to_save = NULL, 
                         title = "DEA_edgeR_QL", 
                         organism = 'hsa', 
                         gene_keytype = "ENSEMBL",
                         omic = c("rna", "microrna")){
  ########## loading packages and custom functions ##########
  source("~/Documentos/09_scripts_R/print_centered_note_v1.R")
  source("~/Documentos/09_scripts_R/create_sequential_dir.R")
  source("~/Documentos/09_scripts_R/Automate_Saving_ggplots.R")
  
  print_centered_note('LOADING PACKAGES ')
  list_of_packages = c("ggplot2", "edgeR","openxlsx","dplyr", "ggplotify", "RColorBrewer", "circlize","ggrepel",
                        'AnnotationDbi', "graphics", "DEGreport", "patchwork","pheatmap", "grid", "gridExtra")
  
  new_packages = list_of_packages[!(list_of_packages %in% installed.packages())]
  if(length(new_packages) > 0){install.packages(new_packages)}
  
  invisible(lapply(list_of_packages, FUN = library, character.only = T))
  rm(list_of_packages, new_packages)
  print('Packages loaded successfully.')
  
  ########## Configuración del directorio de salida ##########
  print("Setting up the directory to save results")
  where_to_save <- ifelse(is.null(where_to_save), getwd(), where_to_save)
  
  output_dir <- create_sequential_dir(path = where_to_save, name = title)

  print(paste0("All files will be saved in ", output_dir))
  
  ########## Creating the reports ##########
  print_centered_note("SETTING UP REPORT FILE ")
  wb_title = paste('DEA_REPORT',title, Sys.Date(), sep=' ')
  wb = createWorkbook(title = wb_title)
  addWorksheet(wb, sheetName = 'Overview')
  
  report = data.frame(matrix(nrow=0, ncol=1))
  colnames(report) = 'Report'
  
  for(i in variable_experimental_design){
    my_kkk = as.data.frame(table(metadata[i]))
    report[nrow(report)+1,] = paste(my_kkk[[2]],my_kkk[[1]], sep=' ', collapse = ',')
  }
  
  report[nrow(report)+1,] = paste0('# of ENSEMBL IDs: ', nrow(count_matrix), ' genes')
  rm(my_kkk)
  ########## step 0: sort rows and columns in metadata and count_matrix ##########
  print_centered_note("SORTING ROWS AND COLUMNS IN METADATA AND COUNT_MATRIX")
  id = match(colnames(count_matrix), rownames(metadata))
  
  metadata = metadata[id,]
  if(all(colnames(count_matrix) == rownames(metadata))){
    print("Rows and columns succesfully sorted. Keep going.")}else{stop("Please, Check if both datasets have the same samples.")}
  
  rm(id)
  
  ########## step 1: set the experimental design ##########
  print_centered_note("SETTING THE EXPERIMENTAL DESIGN")
  design_variable_idx <- which(colnames(metadata) == variable_experimental_design)
  
  if(!is.factor(metadata[[design_variable_idx]])){
    cat("\tSetting", variable_experimental_design, "as a factor.\n\tPlease check levels after analysis.\n")
    metadata[[design_variable_idx]] <- factor(metadata[[design_variable_idx]])
  }else{
    print('Setting the design matrix')
  }
  print("Creating the Experimental Design matrix")
  
  group = make.names(metadata[[design_variable_idx]])
  
  design = model.matrix(~0+group)
  colnames(design) = gsub("group", "", colnames(design))
  colnames(design) = make.names(colnames(design))
  rownames(design) = colnames(count_matrix)
  
  print('Saving experimental design in the report.')
  addWorksheet(wb, sheetName = 'Experimental_design')
  writeDataTable(wb, 'Experimental_design', x = as.data.frame(design), rowNames = T)
  
  ########## step 2: create contrast matrix ##########
  print_centered_note("CREATING CONTRAST MATRIX ")
  # Identificar los niveles de los grupos en el diseño experimental
  levels_group <- levels(metadata[[design_variable_idx]])
  
  # Asignar el primer nivel como "basal" y los siguientes como "tratamientos"
  basal_group <- make.names(levels_group[1])           # Primer nivel es basal
  treatment_groups <- make.names(levels_group[-1])     # Los otros niveles son tratamientos
  
  # Crear la matriz de contrastes para cada tratamiento contra basal
  contrasts <- sapply(treatment_groups, function(trt) paste0(trt, "-", basal_group))
  contrast_matrix <- makeContrasts(contrasts = contrasts, levels = design)
  
  print("Contrast matrix created:")
  print(contrast_matrix)
  
  addWorksheet(wb, sheetName = 'Contrast_matrix')
  writeData(wb, sheet = 'Contrast_matrix', x = as.data.frame(contrast_matrix), rowNames = T)
  ########## step 3: Differential expression analysis ##########
  print_centered_note("MAKING THE DIFFERENTIAL EXPRESSION ANALYSIS")
  # 3.1) crear el objeto DGEList con los metadatos
  d = DGEList(count_matrix)
  d$samples$group = group
  
  # 3.2) Fitting the dge object if the count matrix is from a rna o microRNA experiment
  omic <- tolower(omic)
  if(omic == "rna"){
    print(paste0("Selected omic is ", toupper(omic)))
    # Adding gene annotations
    if(organism == 'hsa'){
      library(org.Hs.eg.db)
      print("Selecting Homo sapiens annotations")
      orgdb <- org.Hs.eg.db
    }
    
    if(organism == "mmu"){
      library(org.Mm.eg.db)
      print("Selecting Mus musculus annotations")
      orgdb <- org.Mm.eg.db
    }
    
    d$genes <- data.frame("Symbol" = mapIds(orgdb, rownames(d), keytype=gene_keytype, column="SYMBOL"),
                          "Gene.Name" = mapIds(orgdb, rownames(d), keytype = gene_keytype, column = "GENENAME"))
    
  } # if else omic == "rna" key
  
  if(omic == "microrna"){
    print("Selected omic is microRNA")
    d$miRNA_family <- gsub("^(([^-]+-){2}[^-]+).*", "\\1", rownames(count_matrix)) 
  } #if else omic == "microrna" key
  
  # Añadir un punto para el control de los genes duplicados cuando se hayan mapeado las lecturas únicamente sobre los genes codificantes.
  writeData(wb, sheet = "Overview", x = report)
  
  # NOTE: In the original tutorial, they said: "Whatever the filtering rule, it should be independent of the information in the 
  # targets file. It should not make any reference to which RNA libraries belong to which group, because doing so would bias the 
  # subsequent differential expression analysis."
  
  # 3.3) Calculate size normalization:
  print("Calculating normalize factor for each library")
  d <- calcNormFactors(d)

  # A normalization factor below one indicates that a small number of high count genes are monopolizing the sequencing, causing 
  # the counts for other genes to be lower than would be usual given the library size. This is a sign that these samples contain a 
  # number of very highly upregulated genes.
  addWorksheet(wb, 'Library_size_per_sample')
  writeDataTable(wb, 'Library_size_per_sample', d$samples, rowNames = T)
  
  # 3.4) Dispersion calculation
  print("Calculating the Dispersion")
  d <- estimateDisp(d, design, robust=TRUE)
  print("Fitting the NB - General Lineal Model")
  fit <- glmQLFit(d, design, robust=TRUE)
  
  #### BCF Plot x ggplot2
  # BCV PLOT
  #Extraer datos para la dispersión biológica
  bcv_df <- data.frame(
    AveLogCPM = d$AveLogCPM,  
    BCV = sqrt(d$tagwise.dispersion),  # Usamos la dispersión por gen en lugar de la tendencia
    Dispersion_Type = "Tagwise Dispersion"
  )
  
  # Extraer la curva de tendencia de dispersión
  trend_df <- data.frame(
    AveLogCPM = d$AveLogCPM,  
    BCV = sqrt(d$trended.dispersion),
    Dispersion_Type = "Trended Dispersion"
  )
  common_dispersion <- data.frame(
    AveLogCPM = range(d$AveLogCPM),  # Para que la línea horizontal cubra todo el eje X
    BCV = sqrt(d$common.dispersion),
    Dispersion_Type = "Common Dispersion"
  )
  dispersion_df <- rbind(bcv_df, trend_df, common_dispersion)
  
  bcv_plot <- ggplot(dispersion_df, aes(x = AveLogCPM, y = BCV, color = Dispersion_Type)) +
    geom_point(data = subset(dispersion_df, Dispersion_Type == "Tagwise Dispersion"), alpha = 0.5) +  # Puntos de genes individuales
    geom_line(data = subset(dispersion_df, Dispersion_Type == "Trended Dispersion"), linewidth = 1) +  # Línea de tendencia
    geom_hline(data = common_dispersion, aes(yintercept = BCV, color = Dispersion_Type), size = 1) +  # Línea horizontal
    labs(
      title = paste0("Common Dispersion: ", round(d$common.dispersion, 3)),
      x = "Average Log CPM",
      y = "Biological Coefficient of Variation (BCV)",
      color = "Dispersion Type"  # Título de la leyenda
    ) +
    scale_color_manual(values = c(
      "Tagwise Dispersion" = "black", 
      "Trended Dispersion" = "blue",
      "Common Dispersion" = "red"
    )) +
    theme_minimal()+
    theme(legend.position = "bottom")
  
  rm(dispersion_df, bcv_df, trend_df, common_dispersion)
  #### plotQLDisp x ggplot2
  # Extraer datos de dispersión Cuasi-Likelihood
  ql_raw_df <- data.frame(AveLogCPM = fit$AveLogCPM,
                          QLRaw = sqrt(sqrt(fit$deviance/fit$df.residual.zeros)), # Esto es de color negro,
                          QLSqueezed = sqrt(sqrt(fit$var.post)), # Esto es de color rojo
                          Trend_Dips = sqrt(sqrt(fit$var.prior)) # Esto es de color azul
                          )
  
  Genewise_QL_Dispersion <- ggplot(ql_raw_df, aes(x = AveLogCPM)) +
    geom_point(aes(y = QLRaw, color = "Raw"), alpha = 0.5) +  # Puntos negros
    geom_point(aes(y = QLSqueezed, color = "Squeezed"), alpha = 0.5) +  # Puntos rojos
    geom_line(aes(y = Trend_Dips, color = "Trend"), size = 1) +  # Línea azul
    labs(
      x = "Average Log2 CPM",
      y = "Quarter-Root Mean Deviance",
      color = "Dispersion Type"  # Título de la leyenda
    ) +
    scale_color_manual(values = c(
      "Raw" = "black", 
      "Trend" = "blue",
      "Squeezed" = "red"
    )) +
    theme_minimal() +
    theme(legend.position = "bottom")  # Leyenda en la parte inferior
  
  rm(ql_raw_df)
  
  combined_plot <- bcv_plot + Genewise_QL_Dispersion + plot_layout(ncol = 2) +
    plot_annotation(
      tag_levels = 'A',
      tag_prefix = 'Fig. ',
      title = "NB GLM Fitted Parameters",
      subtitle = "Biological Coefficient Variation Plot (Fig. A) \nMean-Variance relationship after fitting GLM (Fig. B)"
    ) & 
    theme(plot.subtitle = element_text(size = 14, hjust = 0.5),
          plot.title = element_text(
            family = "lobster", 
            size = 18,
            hjust = 0.5,
            face = "bold", 
            color = "#2a475e"
          ),
          plot.tag = element_text(size = 15, face = "bold", family="lobster"))
  
  save_ggplot(plot = combined_plot, title = "BCV_Genewise_QL_Dispersion_Plots", 
              folder = output_dir, width = 3500, height = 2000)
  
  rm(combined_plot, bcv_plot, Genewise_QL_Dispersion)
  # NOTE: "this moderation (QL) reduces the uncertainty of the estimates and improves testing power"
  #### STEP 4: Testing differentially expression genes ####
  print_centered_note("PERFORMING DIFFERENTIAL EXPRESSION ANALYSIS")
  addWorksheet(wb, sheetName = paste0("DEA_summary"))
  dea_summary <- data.frame(
    Contrast = character(),
    Total_Genes = numeric(),
    PValue_Significant = numeric(),
    Adjusted_PValue_Significant = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:(length(levels_group) - 1)) {
    for (j in (i + 1):length(levels_group)) {
      group1 <- levels_group[i]
      group2 <- levels_group[j]
      contrast_name <- paste(group1, "_vs_", group2, sep="")
      dea_summary[nrow(dea_summary)+1, 1] = contrast_name
      print(paste0("Making contrast between ", contrast_name))
      addWorksheet(wb, sheetName = paste0("DEA_", contrast_name))
      contrast <- makeContrasts(contrasts = paste(group1, "-", group2, sep=""), levels=design)
      test <- glmQLFTest(fit, contrast=contrast)
      print(paste0("Saving results for contrasts ", contrast_name))
      res.tab <- topTags(test, n=Inf)$table
      
      # Calcular los valores de expresión media para cada grupo en el contraste
      group_means <- as.data.frame(
        sapply(c(group1, group2), function(group) {
          round(rowMeans(d$counts[, d$samples$group == group]),1)
        })
      )
      colnames(group_means) <- paste0("Mean_Counts_", c(group1, group2))
      
      group_means <- group_means[match(rownames(res.tab), rownames(group_means)),]
      
      res.tab <- cbind(res.tab, group_means)
      
      # Reordenar las columnas
      if(omic == "rna"){
        res.tab <- res.tab[, c(1,2,8,9,3,4,5,6,7)]
      }
      if(omic == "microrna"){
        res.tab <- res.tab[, c(1,2,6,7,3,4,5)]
      }
      
      # Resumen de genes significativos
      total_genes <- nrow(res.tab)
      significant_genes <- nrow(res.tab %>% dplyr::filter(PValue < 0.05 & abs(logFC) >= 1))
      significant_adj_genes <- nrow(res.tab %>% dplyr::filter(FDR < 0.05 & abs(logFC) >= 1))
      
      dea_summary <- rbind(
        dea_summary,
        data.frame(
          Contrast = contrast_name,
          Total_Genes = total_genes,
          PValue_Significant = significant_genes,
          Adjusted_PValue_Significant = significant_adj_genes
        )
      )
      
      writeDataTable(wb, sheet = paste0("DEA_", contrast_name), res.tab, rowNames = TRUE)
      
      #### visualizations
      print_centered_note("CREATING VISUALIZATIONS OF DE RESULTS")
      # P-value distribution
      print("1.Pvalue Distribution")
      p_val_distr <- degQC(counts = d$counts, groups = d$samples$group, pvalue = res.tab$PValue)
      
      save_ggplot(plot = p_val_distr, title = paste0("P_VALUES_Distribution_Along_Variance_Expression_", contrast_name),
                  folder = output_dir,  width = 2000, height = 2000)
      
      # Volcanoplot
      print("2.Volcano Plot")
      res.tab <- res.tab %>%
        mutate(
          threshold = case_when(
            logFC > 1 & PValue < 0.05 ~ "UP",
            logFC < -1 & PValue < 0.05 ~ "DOWN",
            TRUE ~ "NA"
          ),
          threshold = factor(threshold, levels = c("DOWN", "UP", "NA"))
        )
      
      if(omic == "rna"){
        res.tab <- res.tab %>%
          tibble::rownames_to_column("ENSEMBL")
        
        # Selección de los 10 genes más significativos
        top10_genes <- res.tab %>%
          dplyr::filter(PValue < 0.05) %>%
          dplyr::filter(!is.na(Symbol)) %>%
          arrange(desc(abs(logFC))) %>%
          slice_head(n = 10) %>%
          pull(Symbol)
        
        # Etiquetas para los genes seleccionados
        res.tab <- res.tab %>%
          mutate(labels_top10 = ifelse(Symbol %in% top10_genes, Symbol, NA))
        
      }
      if(omic == "microrna"){
        res.tab <- res.tab %>%
          tibble::rownames_to_column("microRNA")
        
        # Selección de los 10 microRNA más significativos
        top10_genes <- res.tab %>%
          dplyr::filter(PValue < 0.05) %>%
          arrange(desc(abs(logFC))) %>%
          slice_head(n = 10) %>%
          pull(microRNA)
        
        # Etiquetas para los genes seleccionados
        res.tab <- res.tab %>%
          mutate(labels_top10 = ifelse(microRNA %in% top10_genes, microRNA, NA))
        
      }
      
      
      # Paso 2: Crear el gráfico
      volcano_plot <- ggplot(res.tab, aes(x = logFC, y = -log10(PValue))) +
        geom_point(aes(color = threshold), size = 1) +
        geom_point(data = res.tab %>% dplyr::filter(!is.na(labels_top10)), color = "red", size = 3) +
        geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "blue") +  # Línea p-valor < 0.01
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") + # Línea p-valor < 0.05
        geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "darkgray") +  # Líneas logFC
        geom_text_repel(
          aes(label = labels_top10),
          size = 4,
          show.legend = FALSE,
          seed = 42,
          box.padding = 0.5,
          max.overlaps = Inf,
          segment.size = 0.5
        ) +
        scale_color_manual(values = c("DOWN" = "#0094a4", "UP" = "#FAA99D", "NA" = "grey")) +
        scale_y_continuous(
          breaks = c(0, -log10(0.05), -log10(0.01), seq(from = 2, to = ceiling(max(-log10(res.tab$PValue))))),
          labels = c("0", "**", "***", seq(from = 2, to = ceiling(max(-log10(res.tab$PValue)))))
        ) +
        labs(
          x = expression(log[2]("Fold Change")),
          y = expression(-log[10]("P-value"))
          #title = paste0("Volcano Plot ", title),
          #subtitle = contrast_name
        ) +
        theme_minimal(base_size = 14) +
        theme(
          legend.position = "none",
          axis.text = element_text(face = "bold"),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, face = "italic")
        )
      
      save_ggplot(plot = volcano_plot, 
                  title = paste0("Volcano_Plot_", title, "_", contrast_name), 
                  folder = output_dir,
                  width = 4000, height = 2000)
      
      
      # heatmap
      print("3.Heatmap")
      if(omic == "rna"){
        significant_genes <- res.tab %>%
          dplyr::filter(FDR < 0.05 & abs(logFC) > 1) %>%
          pull(ENSEMBL)
        
        cpm_vis = cpm(d[which(rownames(d) %in% significant_genes),], log=T)
        ensembl_to_symbol <- setNames(d$genes[which(rownames(d$genes) %in% significant_genes), "Symbol"], rownames(d$genes)[which(rownames(d$genes) %in% significant_genes)])  # ENSEMBL como nombres, SYMBOL como valores
        rownames(cpm_vis) <- ensembl_to_symbol[rownames(cpm_vis)]
      }
      if(omic == "microrna"){
        significant_genes <- res.tab %>%
          dplyr::filter(FDR < 0.05 & abs(logFC) > 1) %>%
          pull(microRNA)
        
        cpm_vis = cpm(d[which(rownames(d) %in% significant_genes),], log=T)
      }
      
      # Paso 3: Crear las anotaciones
      # `metadata` debe tener las muestras como filas y las anotaciones como columnas.
      # Asegúrate de que las filas coincidan con los nombres de las columnas de `expression_data`.
      
      # only if there are any FDR < 0.05 genes, the heatmap will be created.
      if(length(significant_genes) > 0){
        annotation_col <- metadata
        sig_genes_heatmap <- as.ggplot(
          pheatmap(
            mat = cpm_vis,
            scale = "row",                # Estandarizar por filas
            cluster_rows = TRUE,          # Agrupar genes
            cluster_cols = TRUE,          # Agrupar muestras
            annotation_col = annotation_col,
            #annotation_colors = annotation_colors,
            show_rownames = FALSE,        # Opcional: Ocultar nombres de genes si son demasiados
            show_colnames = TRUE,         # Mostrar nombres de muestras
            color = colorRampPalette(c("blue", "white", "red"))(50),  # Gradiente de colores
            #main = paste0("Heatmap: ", title, " - ", contrast_name)  # Título del gráfico
          )
        )
        
        print("4.Combined Plot")
        combined_plot <- sig_genes_heatmap + (volcano_plot / p_val_distr)
        combined_plot <- combined_plot +
          plot_annotation(
            tag_levels = 'A',
            tag_prefix = 'Fig. ',
            title = paste0('DE Genes from Contrast ', contrast_name),
            subtitle = "Heatmap of standarize gene expression (Fig. A) \n VolcanoPlot (Fig. B) \nP-Value Distribution (Fig. C)",
          ) & 
          theme(plot.subtitle = element_text(size = 14, hjust = 0.5),
                plot.title = element_text(
                  family = "lobster", 
                  size = 18,
                  hjust = 0.5,
                  face = "bold", 
                  color = "#2a475e"
                ),
                plot.tag = element_text(size = 15, face = "bold", family="lobster"))
        
        save_ggplot(plot = sig_genes_heatmap, title = "Heatmap_Sig_Genes",
                    folder = output_dir, width = 5000, height = 3000)
        
        save_ggplot(plot = combined_plot, title = "SignifGenes_Visualizations",
                    folder = output_dir, width = 5000, height = 3000)
      } # if statemente heatmap check
      
    } # j for-loop 
  } # i for-loop
  cat("\n")
  writeDataTable(wb, sheet = paste0("DEA_summary"), dea_summary)
  print("Saving Excel Worbook")
  saveWorkbook(wb, file = file.path(output_dir, "DEA_Results.xlsx"),overwrite = T)
  
  print_centered_note("END OF DE ANALYSIS WITH EdgeR QUASI-LIKELIHOOD ")
  return(res.tab)
} # Llave de la función
