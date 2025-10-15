########## 12/11/2024 ##########

#·············································································#
################# LIMMA-VOOM DIFFERENTIAL EXPRESSION ANALYSIS ################# 
#·············································································#

# In this script I'm going to create a function to perform a differnetial expression analysis with limma-voom pipeline based on my
# previous script with limma-voom to add more updates, such as the DGE, more visualizations, omit the QC because I've got another script,
# and add pathway analysis.
limma_voom_dea <- function(count_matrix, 
                           metadata, 
                           where_to_save = NULL, 
                           title = "DEA_limma_voom", 
                           VOOM_transformation = T,
                           variable_experimental_design = NULL, 
                           organism = 'hsa', 
                           gene_keytype = "ENSEMBL"){
  ########## loading packages ##########
  source("~/Documentos/09_scripts_R/print_centered_note_v1.R")
  print_centered_note('LOADING PACKAGES ')
  list_of_packages = c("ggplot2", "edgeR","openxlsx","dplyr", "RColorBrewer", "circlize","ggrepel",
                       'AnnotationDbi', "limma", "graphics", "grid", "gridExtra", "ggplotify", "DEGreport")
  
  new_packages = list_of_packages[!(list_of_packages %in% installed.packages())]
  if(length(new_packages) > 0){install.packages(new_packages)}
  
  invisible(lapply(list_of_packages, FUN = library, character.only = T))
  rm(list_of_packages, new_packages)
  print('Packages loaded successfully.')
  
  ########## Configuración del directorio de salida ##########
  # initial checkpoint #
  if(is.null(variable_experimental_design)){
    stop("Please, check if the variable_experimental_design argument was assigned propertly")
  }
  
  where_to_save <- ifelse(is.null(where_to_save), getwd(), where_to_save)
  output_dir <- file.path(where_to_save, title)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  ########## Creating the reports ##########
  wb_title = paste('DEA_REPORT_LIMMA', Sys.Date(), sep='_')
  wb = createWorkbook(title = wb_title)
  addWorksheet(wb, sheetName = 'Overview')
  
  report = data.frame(matrix(nrow=0, ncol=1))
  colnames(report) = 'Report'
  
  for(i in variable_experimental_design){
    my_kkk = as.data.frame(table(metadata[i]))
    report[nrow(report)+1,] = paste(my_kkk[[2]],my_kkk[[1]], sep=' ', collapse = ',')
  }
  
  report[nrow(report)+1,] = paste0('Number of genes: ', nrow(count_matrix), ' genes')
  writeData(wb, sheet = "Overview", x = report)
  
  rm(i, my_kkk)
  ########## step 0: sort rows and columns in metadata and count_matrix ##########
  print_centered_note("SORTING ROWS AND COLUMNS IN METADATA AND COUNT_MATRIX")
  
  id = match(colnames(count_matrix), rownames(metadata))
  
  metadata = metadata[id,]
  if(all(colnames(count_matrix) == rownames(metadata))){
    print("Rows and columns succesfully sorted. Keep going")}  else{
    stop("Please, Check if both datasets have the same samples.")
  }
  
  rm(id)
  
  ########## step 1: set the experimental design ##########
  print_centered_note("SETTING THE EXPERIMENTAL DESIGN")
  design_variable_idx <- which(colnames(metadata) == variable_experimental_design)
  if(!is.factor(metadata[[design_variable_idx]])){
    cat("Setting", variable_experimental_design, "as a factor.\nPlease check levels after analysis.\n")
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
  print_centered_note("SETTING UP THE CONTRAST MATRIX ")
  
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
  
  ########## 3.1) CREAR EL OBJETO DGEList CON LOS METADATOS ########## 
  print("Creating the DGEList Object")
  d = DGEList(count_matrix)
  d$samples$group = group
  
  ########## 3.2) FACTOR DE NORMALIZACIÓN ##########
  print("Calculating Normalization Factors by TMM method")
  d = calcNormFactors(d, method="TMM") ## Calcular los factores de normalización Trimmed-Mean of M-values method
  
  addWorksheet(wb, sheetName = "Library_overview")
  writeData(wb, sheet = "Library_overview", as.data.frame(d$samples), rowNames = T)
  ########## 3.3) VOOM TRANSFORMATOIN ##########
  if(VOOM_transformation){
    print('Voom transformation will be applied to RNA-seq data')

    y = voom(d, design, plot=F, save.plot = T)
    print("Voom transformation applied succesfully")
    
    ########## 3.4) AJUSTE DEL MODELO LINEAL ##########
    print('Fitting the linear model')
    vfit = lmFit(y, design)
    vfit = contrasts.fit(vfit, contrast_matrix)
    efit = eBayes(vfit)
    
    ########## VISUALIZATIONS
    print("Making Voom Transformation Visualizations")
    # 1.Mean-Variance Relationship 
    voom_df <- data.frame(
      x = y$voom.xy$x,
      y = y$voom.xy$y
    )
    voom_line_df <- data.frame(
      x = y$voom.line[[1]],  # Coordenadas x de la línea
      y = y$voom.line[[2]]   # Coordenadas y de la línea
    )
    
    p_voom <- ggplot() +
      geom_point(data = voom_df, aes(x = x, y = y), color = "blue", alpha = 0.5) +  # Puntos
      geom_line(data = voom_line_df, aes(x = x, y = y), color = "red", size = 1) + # Línea de tendencia
      labs(title = "Mean-Variance Before Fitting",x = y$voom.xy$xlab, y = y$voom.xy$ylab) +
      theme_minimal()
    
    # 2. Mean-Variance Trend
    mv_df <- data.frame(
      x = efit$Amean,
      y = sqrt(efit$sigma),
      hline = sqrt(sqrt(efit$s2.prior))
    )
  
    p_trend <- ggplot(mv_df, aes(x = x, y = y)) +
      geom_point(alpha = 0.5, color = "black") +
      geom_line(data = mv_df, aes(x = x, y = hline), color = "blue", size = 1) + # Línea de tendencia
      labs(title = "Mean-Variance After Fitting", x = "Average Log-Expression", y = "Sqrt(sigma)") +
      theme_minimal()
    
    #3.Effect of normalization
    norm_data <- as.data.frame(y$E)
    norm_data_long <- reshape2::melt(norm_data)
    
    nsamples = ncol(norm_data)
    sample_col <- hcl.colors(nsamples, "Temps")
    group_col = hcl.colors(length(unique(d$samples$group)),'Temps')
    
    p_norm <- ggplot(norm_data_long, aes(x = value, color = variable)) +
      geom_density(size = 1, alpha = 0.7) +
      scale_color_manual(values = hcl.colors(nsamples, "Temps")) +
      labs(title = "Voom Normalization", x = "log-CPM", y = "Density") +
      theme_minimal() + theme(legend.position = "none")
    
    
    print("Saving Voom Transformations Plot")
    final_plot <- p_voom + p_trend + p_norm + plot_layout(ncol = 3)
    png(filename = file.path(output_dir, "01_Voom_Transformations_Plots.png"), res = 300, width = 2500, height = 1000)
    print(final_plot)
    dev.off()
    rm(final_plot, p_voom, p_trend, p_norm, norm_data, norm_data_long, mv_df,voom_df, voom_line_df)
  }else{print("Voom transformation won't be applied to RNA-seq data")}
  
  ########## step 4: Crear expresión promedio de cada gen en cada grupo ##########
  print_centered_note("CALCULATING AVERAGE EXPRESSION PER GENE PER GROUP")
  
  # Calcular la expresión promedio por gen y grupo
  avg_expression = apply(y$E, 1, function(gene_counts) {
    tapply(gene_counts, group, mean)
  })
  avg_expression_df = as.data.frame(t(avg_expression))   # Convertir la matriz resultante a un formato adecuado para agregarlo al reporte

  colnames(avg_expression_df) = paste0("Norm_Avg_Expression_", colnames(avg_expression_df))
  avg_expression_df = avg_expression_df %>%
    tibble::rownames_to_column("Gene")
  
  tfit <- treat(vfit, lfc = 1)
  dt <- decideTests(tfit)
  summary_deg = data.frame(summary(dt), stringsAsFactors = F)
  summary_deg = as.data.frame(apply(summary_deg, 2, as.character))
  summary_deg[nrow(summary_deg)+1,] = c("Note: This are DE genes with significant adj p-value < 0.05 & abs(log2foldchange) >= 1 was applied", NA, NA)
  
  addWorksheet(wb, sheetName = "DEA_summary")
  writeData(wb, sheet = "DEA_summary", x = summary_deg)
  
  ########## step 4: Extracting DE genes & visualizations ##########
  print_centered_note("EXTRACTING DIFFERENTIALY EXPRESSED GENES ")
  for(i in 1:ncol(contrast_matrix)){
    comparison = colnames(contrast_matrix)[i]
    print(paste0("Extracting DE genes for comparisson: ", comparison))
    deg = topTreat(tfit, coef = i, n= Inf)
    deg = deg %>%
      tibble::rownames_to_column("Gene") %>%
      arrange(P.Value)
    
    # In the following part I would like to introduce a adicional dataframe separetly from excel sheets 
    # to save the information about multimapping genes, but this will be in the second version.
    if(organism == "hsa"){
      print('Transforming ENSEMBL codes in SYMBOL for Human')
      library(org.Hs.eg.db)
      org.db <- org.Hs.eg.db
    }
    if(organism == "mmu"){
      print('Transforming ENSEMBL codes in SYMBOL for Mouse')
      library(org.Mm.eg.db)
      org.db <- org.Mm.eg.db
    }
    
    symbol = AnnotationDbi::mapIds(org.db, keys = deg$Gene, column = c('SYMBOL'), keytype = gene_keytype,
                                     multiVals = 'first') 
    genename = AnnotationDbi::mapIds(org.db, keys = deg$Gene, column = c('GENENAME'), keytype = gene_keytype,
                                     multiVals = 'first') 
    deg$Symbol <- symbol
    deg$Genename <- genename
    
    deg = merge(deg, avg_expression_df, by = "Gene")
    deg = deg %>% arrange(P.Value)
    deg = deg[c(1,2,3,9,10,4,5,6,7,8)]
    
    ########## step 5: visualizations ##########
    print(paste0("Making Visualizations for comparison: ", comparison))
    
    # P-value distribution
    # Más adelante me gustaría añadir esta función de distribución de p-valores a
    #print("1.Pvalue Distribution").
    #p_val_distr <- degQC(counts = d$counts, groups = d$samples$group, pvalue = tfit$p.value)
    #png(filename = file.path(output_dir, paste0("04_P_VALUES_Distribution_Along_Variance_Expression_", contrast_name, ".png")), res = 300, width = 2000, height = 2000)
    #print(p_val_distr)
    #dev.off()
    ########## VolcanoPlot
    print("1.Volcano Plot")
    top10 = deg %>%
      dplyr::filter(adj.P.Val < 0.05) %>%
      arrange(desc(abs(logFC)))
    top10_genes <- head(top10[!is.na(top10$Symbol),"Symbol"], 10)
    
    deg$threshold = ifelse(deg$logFC > 1 & deg$P.Value < 0.05, 'UP',
                               ifelse(deg$logFC < -1 & deg$P.Value < 0.05, 'DOWN','NA'))
    deg$threshold = factor(deg$threshold, levels = c('DOWN', 'UP', 'NA'))
    
    deg$labels_top10 = ifelse(deg$Symbol %in% top10_genes, deg$Symbol, NA)

    volcano_plot = ggplot(deg, aes(x = logFC, y = -log10(P.Value), label= labels_top10)) +
      geom_point(aes(color=threshold),size=1) +
      geom_point(data = deg[-which(is.na(deg$labels_top10)),], color = "red", size=3)+
      geom_hline(yintercept = -log10(0.01), linetype = "dashed") +  # Línea de corte
      scale_color_manual(values=c('#0094a4',"#FAA99D", 'grey'))+
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +  # Línea de corte
      geom_vline(xintercept = -1, linetype = "dashed") +
      geom_vline(xintercept = +1, linetype = "dashed") +  # Línea central
      labs(x = expression(paste("log"[2],"(Fold Change)")), y = expression(paste("-log"[10],"(P-value)")))+
      geom_text_repel(data=deg[-which(is.na(deg$labels_top10)),], aes(label= labels_top10), 
                      size=4,
                      show.legend = F,
                      min.segment.length = 0, 
                      seed = 42, 
                      box.padding = 0.5, 
                      max.overlaps = Inf,
                      xlim = c(-Inf, Inf),
                      segment.size=1,
                      # Do not repel from top or bottom edges.
                      ylim = c(-Inf, NA))+
      scale_y_continuous(breaks=c(0,-log10(0.05),-log10(0.01), seq(from=2, to =round(1.5*max(-log10(deg$P.Value)),digits = 0))),
                         labels=c('0', '**', '***', as.character(seq(from=2, to =round(1.5*max(-log10(deg$P.Value)),digits = 0))))
      )+
      theme(legend.position = 'none', 
            axis.text.x = element_text(face="bold",size=14),
            axis.text.y = element_text(face="bold", size=14))
    volcano_plot <- volcano_plot + theme(legend.position = "none") # Opcional: oculta la leyenda si se desea
    
   ########## Heatmap
   print("2.Heatmap")
   de_genes <- deg %>%
     dplyr::filter(adj.P.Val < 0.05 & abs(logFC) > 1)
   
   if(!is.null(de_genes)){
     cpm_vis = cpm(d[which(rownames(d) %in% de_genes$Gene),], log=T)
     
     ensembl_to_symbol <- setNames(de_genes$Symbol, de_genes$Gene)  # ENSEMBL como nombres, SYMBOL como valores
     rownames(cpm_vis) <- ensembl_to_symbol[rownames(cpm_vis)]
     
     heatmap_colors <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100)
     
     #group_colors <- setNames(hcl.colors(length(unique(d$samples$group)), "Set2"), unique(d$samples$group))
     #annotation_colors <- list(Group = group_colors)  # Lista de colores para la anotación
     
     annotation_col <- metadata
     heatmap_grob <- as.ggplot(
       pheatmap(cpm_vis, 
                cluster_rows = TRUE, 
                scale = "row",
                cluster_cols = TRUE,  
                color = heatmap_colors, 
                show_rownames = FALSE, 
                show_colnames = TRUE,
                fontsize_row = 8,
                fontsize_col = 10,
                #main = paste("Heatmap of Differentially Expressed Genes -", comparison),
                annotation_col = annotation_col,  
                #annotation_colors = annotation_colors,  
                legend_labels = c("logCPM")
       )
     )
     print("4.Combined Plot")
     combined_plot <- heatmap_grob +volcano_plot
     combined_plot <- combined_plot + plot_annotation(
       title = paste0('DE Genes from Contrast ', comparison),
       tag_levels = 'A',
       tag_prefix = 'Fig. ',
       theme = theme(
         plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
         plot.tag = element_text(size = 15, face = "bold", family="lobster")))
     
     png(filename = file.path(output_dir, "02_SignifGenes_Visualizations.png"), 
         res = 300, width = 5000, height = 3000)
     print(combined_plot)
     dev.off()
   }
   
   
   print(paste0("Saving DE genes for comparisson ", comparison))
   addWorksheet(wb, sheetName = paste0("DEA_Result_",comparison))
   writeDataTable(wb, paste0("DEA_Result_",comparison), deg)
   } # For loop key
  
  print("Saving the Workbook")
  saveWorkbook(wb, file.path(output_dir, paste0(wb_title,".xlsx")), overwrite = T)
  print_centered_note("End of the script")
}
  
