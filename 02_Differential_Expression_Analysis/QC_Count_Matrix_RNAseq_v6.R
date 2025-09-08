########## 05/11/2024 ##########

#·····················································#
##### Quality Control of count matrix RNA-seq data ####
#·····················································#

# I'm going to create a function to make a QC of any count matrix of RNA-seq data, create a report of library size in each
# sample and some visualizations to diagnostic the data. Also, this function will return a low expression genes filtered matrix
# to start working in DE analysis or whatever.

# For filtering, I'm going to use the function filterByExpr() and an mandatory argument will be the experimental design.

########## prerequesites ##########
# 1) Count matrix: rows must be the genes in ENSEMBL ID and columns must be sample names.
# 2) Metadata: Rownames must be the sample names. All metainformation about samples must be characters or factors.

########## inputs ##########
# 1) count_mat -> A csv file with count matrix with rownames as ENSEMBL IDs and colnames with sample names.
# 2) metadata -> A csv file with metadata with rownames as samplefiles and different variable information for each samples.
#                 Variables with other information must be factors ordered with the treatment at first levels and reference at second.
# 3) variable_experimental_design -> A variable name to set the experimental design for filtering low expression genes.

qc_count_matrix_rnaseq <- function(count_mat, 
                                   metadata, 
                                   variable_experimental_design = NULL,
                                   where_to_save = NULL, 
                                   title = "RNAseq_analysis", 
                                   min_counts = 10, 
                                   organism = "hsa",
                                   gene_keytype = "ENSEMBL",
                                   omic = c("rna", "microrna")){
  ########## STEP 1: loading packages ##########
  source("~/Documentos/09_scripts_R/print_centered_note_v1.R")
  source("~/Documentos/09_scripts_R/Automate_Saving_ggplots.R")
  source("~/Documentos/09_scripts_R/create_sequential_dir.R")
  
  print_centered_note("LOADING PACKAGES ")
  
  list_of_packages = c('limma', 'edgeR', 'DESeq2',"DEGreport", "tidyr","scales",'ggplot2', 'cowplot', 
                       'BiocManager', 'openxlsx',"AnnotationDbi",'grDevices', 'graphics', "dplyr",
                       'RColorBrewer', 'RNAseqQC', 'ensembldb', 'purrr', "patchwork", "pheatmap",
                       "ggplotify", "ggrepel")
  
  new_packages = list_of_packages[!(list_of_packages %in% installed.packages())]
  if(length(new_packages) > 0){install.packages(new_packages)}
  
  invisible(lapply(list_of_packages, FUN = library, character.only = T))
  rm(list_of_packages, new_packages)
  print('Packages loaded successfully.')
  
  # Checkpoints 
  if(length(variable_experimental_design) == 0){
    stop("Set the experimetal design variable to make the design matrix")
  }
  
  # Configuración del directorio de salida
  where_to_save <- ifelse(is.null(where_to_save), getwd(), where_to_save)
  output_dir <- create_sequential_dir(where_to_save, name = title)

  ########## STEP 2: creating the report ##########
  print_centered_note("CREATING THE OVERALL REPORT")
  
  print('Starting creating the report')
  report_wb <- createWorkbook()
  addWorksheet(report_wb, 'Overview')
  addWorksheet(report_wb, "Experimental_design")
  addWorksheet(report_wb, 'Library_size_per_sample')
  
  report = data.frame(matrix(nrow=0, ncol=1))
  colnames(report) = 'Report'
  
  ########## STEP 3: preparing count_mat and metadata ##########
  # First check: Are rownames of metadata in the same order as colnames in count matrix?
  print_centered_note("PREPARING COUNT MATRIX AND METADATA")
  
  id = match(colnames(count_mat), rownames(metadata))
  
  metadata = metadata[id,]
  if(!all(colnames(count_mat) == rownames(metadata))){
    stop('Please, prepare both datasets before to continue with the QC analysis.')
  }else{
    print('Both dataset are sorted. Continue with the analysis.')
  }
  rm(id)
  # Second check: Are all values in count_mat integers?
  # Función para verificar y redondear números decimales en una matriz de cuentas
  check_and_round <- function(count_mat) {
    # Verificar si hay algún número decimal en la matriz
    tiene_decimales <- any(count_mat %% 1 != 0)
    
    if (tiene_decimales) {
      print("Decimal numbers were found in the count matrix. Rounding to integers...")
      count_mat <- round(count_mat, digits = 0)
    } else {
      print("The count matrix doesn't have decimal numbers. Continue as well.")
    }
    
    return(count_mat)
  }
  print('Checking if there are decimal numbers in count matrix')
  count_mat <- check_and_round(count_mat)
  
  ########## STEP 2: set the experimental design ##########
  # Description of the dataset
  chr_idx <- sapply(metadata, is.character)
  
  print("Creating a summary for each categorical variable")
  for(column in colnames(metadata)[chr_idx]){
    print(paste0("Creating the report for ", column))
    my_kkk <- as.data.frame(table(metadata[column]))
    report[nrow(report)+1,] = paste0(column,": ", paste(my_kkk[[2]],my_kkk[[1]], sep=' ', collapse = ', '))
  }
  rm(column, my_kkk)
  # First get the index of that variable_experimental_design 
  
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
  rownames(design) = colnames(count_mat)
  
  print('Saving experimental design in the report.')
  writeDataTable(report_wb, 'Experimental_design', x = as.data.frame(design), rowNames = T)
  
  ########## STEP 3: cleaning count matrix ##########
  omic <- tolower(omic)
  
  if(organism == "hsa"){
    print("Selected organism is Human")
    library(org.Hs.eg.db)
    org.db <- org.Hs.eg.db
    anotation_hub <- "AH89426"
  }
  
  if(organism == "mmu"){
    print("Selected organism is Mouse")
    library(org.Mm.eg.db)
    org.db <- org.Mm.eg.db
    annotation_hub <- "AH95775"
  }
  
  if(omic == "rna"){
    print_centered_note("REMOVING GENES WITHOUT SYMBOL AND DUPLICATES FOR RNA-SEQ DATA")
    print("Creating the annotations")
    annotations <- data.frame(gene_keytype = rownames(count_mat),
                              "SYMBOL" = mapIds(org.db, 
                                                keys = rownames(count_mat), 
                                                column = "SYMBOL",
                                                keytype = gene_keytype,
                                                multiVals = "first"),
                              "GENENAME" = mapIds(org.db,
                                                  keys = rownames(count_mat),
                                                  column = "GENENAME",
                                                  keytype = gene_keytype,
                                                  multiVals = "first"))
    colnames(annotations)[1] <- gene_keytype
    addWorksheet(report_wb, sheetName = "Gene Annotations")
    writeDataTable(report_wb, sheet = "Gene Annotations" , x = annotations, rowNames = F)
    count_mat <- count_mat %>%
      tibble::rownames_to_column(var = gene_keytype)
    
    report[nrow(report)+1,] = paste(paste0("# Total ", gene_keytype, " IDs:"), nrow(count_mat))
    
    count_mat <- merge(count_mat, annotations, by = gene_keytype)
    idx_na <- which(is.na(count_mat$SYMBOL))
    report[nrow(report)+1,] = paste("# Genes without SYMBOL:", length(idx_na))
    idx_predicted <- grep(pattern = "predicted", count_mat$GENENAME)
    report[nrow(report)+1,] = paste("# Predicted Genes:", length(idx_predicted))
    
    print("Filtering predicted and NA ENSEMBLIDs")
    count_mat <- count_mat %>%
      dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
      dplyr::filter(!grepl("predicted", GENENAME, ignore.case = TRUE))
    
    dup <- nrow(count_mat)
    
    # Filtrar duplicados por mayor expresión
    print("Filtering duplicates based on higher expression")
    count_mat <- count_mat %>%
      group_by(SYMBOL) %>%
      dplyr::slice_max(order_by = rowSums(across(where(is.numeric))), n = 1) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(gene_keytype) %>%
      as.data.frame()
    
    rownames(count_mat) <- count_mat[,gene_keytype, drop = T]
    dup2 <- dup - nrow(count_mat)
    report[nrow(report)+1,] = paste("# Duplicated Genes :", dup2)
    report[nrow(report)+1,] = paste("Final # Genes :", nrow(count_mat))
    report[nrow(report)+1,] = paste("Selection Protocol: Based on Higher expression in across all samples")
    
  }
  
  # NOTA PARA EL FUTURO: Estaría bien introducir una parte también de preprocesamiento de las
  # matrices de cuentas de microRNA, porque podemos darle una matriz de microRNA maduros para
  # que sumen las cuentas de los microRNAs maduros que procedan de distintos precursores
  
  ########## STEP 4: setting the DGE objects ##########
  print_centered_note("CREATING THE DGE OBJECTS ")
  character_cols <- sapply(count_mat, is.character)
  
  print("Creating the DGEList Object")
  dge_raw = DGEList(count_mat[,character_cols == FALSE],
                    group = group)
  
  omic <- tolower(omic)
  if(omic == "rna"){
    print(paste0("Selected omic is ", toupper(omic)))
    # Adding gene annotations
    if(organism == 'hsa'){
      library(org.Hs.eg.db)
      print("Selecting Homo sapiens annotations")
      orgdb <- org.Hs.eg.db
      anotation_hub <- "AH89426"
    }
    
    if(organism == "mmu"){
      library(org.Mm.eg.db)
      print("Selecting Mus musculus annotations")
      orgdb <- org.Mm.eg.db
      anotation_hub <- "AH95775"
    }
    
    dge_raw$genes <- data.frame("Symbol" = mapIds(orgdb, rownames(dge_raw), keytype=gene_keytype, column="SYMBOL"),
                          "Gene.Name" = mapIds(orgdb, rownames(dge_raw), keytype = gene_keytype, column = "GENENAME"))
    
  } # if else omic == "rna" key
  
  if(omic == "microrna"){
    print("Selected omic is microRNA")
    dge_raw$miRNA_family <- gsub("^(([^-]+-){2}[^-]+).*", "\\1", rownames(count_mat[,character_cols == FALSE])) 
  } #if else omic == "microrna" key
  
  if(omic == "rna"){
    if(gene_keytype == "ENSEMBL"){
      dds_raw <- tryCatch(
        {
          suppressWarnings({
            # Ejecutar make_dds y suprimir advertencias
            make_dds(count_mat[,character_cols == FALSE], metadata, ah_record = anotation_hub)
          })
        },
        error = function(e) {
          # Código que se ejecuta si ocurre un error
          message("Error en make_dds: ", e$message)
          NULL  # Devuelve NULL si ocurre un error
        }
      )
    }
  } # rna if key for make dds
  
  if(omic == "microrna"){
    dds_raw <- dds_raw <- DESeqDataSetFromMatrix(countData = count_mat[,character_cols == FALSE],
                                                 colData = metadata, design = design)
  }
  #### STEP 5: RAW COUNTS VISUALIZATIONS ####
  # computing mean and median valuo of library size of the hole experiment
  print('Computing library size of hole experiment')
  samplenames = colnames(dge_raw)
  L = mean(dge_raw$samples$lib.size) * 1e-6 ## Tamaño medio de la librería de todo el experimento
  M = median(dge_raw$samples$lib.size) * 1e-6
  
  # Adding this information to the report
  report[nrow(report)+1,] = paste('Mean library size (M):', round(L, digits = 2))
  report[nrow(report)+1,] = paste('Median library size (M):', round(M, digits = 2))
  
  # Calculating cutoffs just for visualization
  lcpm.cutoff = log2(10/M + 2/L)
  nsamples = ncol(dge_raw)
  sample_col <- hcl.colors(nsamples, "Temps")
  group_col = hcl.colors(length(unique(dge_raw$samples$group)),'Temps')
  lcpm_pre = cpm(dge_raw, log=T)
  
  # FILTERING GENES ####
  print_centered_note('FILTERING LOW EXPRESSION GENES ')
  # Remember that the filtering thershold should be independent of the experimental design
  keep.exprs = edgeR::filterByExpr(dge_raw, min.total.count = min_counts) 
  d = dge_raw[keep.exprs,, keep.lib.sizes=FALSE] # Keep.lib.size es para recalcular el tamaño de la librería tras la eliminación de los genes
  c = dds_raw[keep.exprs,]
  lcpm_post <- cpm(d, log=TRUE)
  
  print_centered_note("CALCULATE NORMALIZATION FACTORS AND NORMALIZING GENE EXPRESSION")
  d = calcNormFactors(d, method="TMM") ## Calcular los factores de normalización Trimmed-Mean of M-values method
  
  lcpm_post_norm = cpm(d, log=TRUE)
  # Adding this information to the report
  report[nrow(report)+1,] = paste('QC filter: FilterByExperiment', "(min", min_counts, "counts/gene)")
  report[nrow(report)+1,] = paste('Genes Prefilter:', dim(dge_raw)[1]) # Genes Pre Filtering
  report[nrow(report)+1,] = paste('Genes Removed:', dim(dge_raw)[1] - dim(d)[1]) # Genes Post Filtering
  report[nrow(report)+1,] = paste('Genes Postfilter:', dim(d)[1]) # Genes Post Filtering
  
  writeDataTable(report_wb, sheet = 'Overview', report, rowNames = F)
  writeDataTable(report_wb, 'Library_size_per_sample', d$samples, rowNames = T)
  
  # visualizations ####
  print_centered_note("CREATING PREFILTERING QC PLOTS ")
  
  # PREFILTERING PLOTS ####
  # 1. DENSITY PLOT
  print("1.Density Plot")
  densities <- lapply(1:ncol(lcpm_pre), function(i) {
    den <- density(lcpm_pre[, i])
    data.frame(x = den$x, y = den$y, sample = as.factor(i))
  })
  density_data <- bind_rows(densities)
  counts_density <- ggplot(density_data, aes(x = x, y = y, color = sample)) +
    geom_line(linewidth = 1.5) +  # Trazar las líneas
    scale_color_manual(values = sample_col) +  # Definir los colores
    geom_vline(xintercept = lcpm.cutoff, linetype = "dashed", color = "black") +  # Línea vertical
    labs(
      #title = "A. Raw data",
      x = "Log(CPM)",
      y = "Density"
    ) +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12)
    ) +
    ylim(0, 0.26) +# Limitar el eje Y
    theme(legend.position = "none")
  
  # 2.GENE DETECTION
  print("2.Gene Detection Plot")
  gene_detection <- plot_gene_detection(dds_raw) + 
    geom_vline(xintercept = 100, linetype = "dotted", linewidth = 1)
    
  # 3.GENE BIOTYPES
  if(omic == "rna"){
    print("3.Gene Biotype Plot")
    gene_biotypes <- plot_biotypes(dds_raw) + 
      theme(axis.text.x = element_text(face="bold", size=4, angle=45),
            legend.position = "bottom") + 
      xlab(NULL) +
      theme(axis.text.x = element_blank())
    # + scale_y_continuous(labels = label_number(scale = 1e-02))
  }
  
  # 4.MEAN_VARIANCE_PLOT
  if(omic == "rna"){vsd_raw <- vst(dds_raw)}
  if(omic == "microrna"){vsd_raw <- varianceStabilizingTransformation(dds_raw)}
  print("4.Mean-Variance Plot")
  mean_variance_plot <- mean_sd_plot(vsd_raw)+ theme(legend.position = "none") + ylab("Standard Deviation (SD)")
    
  # 7.NORMALIZATION BOXPLOT
  print("5.Normalization Plot")
  df_lcpm <- as.data.frame(lcpm_pre) %>%
    mutate(Gene = rownames(.)) %>%
    pivot_longer(cols = -Gene, names_to = "Sample", values_to = "LogCPM")
  
  unnormalized <- ggplot(df_lcpm, aes(x = Sample, y = LogCPM, fill = Sample)) +
    geom_boxplot(outlier.shape = NA) +  # No mostrar outliers
    scale_fill_manual(values = sample_col) +  # Colores personalizados
    labs(y = "Log(CPM)", x = "") +
    #theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1,  vjust = -0.025),  # Rotar etiquetas del eje X
      legend.position = "none"  # Ocultar la leyenda si no es necesaria
    )
  
  # 6.PREFILTER COMBINED PLOT
  print("6.Pre-Filtering Low Expression Genes Combined Plot")
  if(omic == "rna"){
    pre_filter_plots <- (counts_density | gene_detection) / gene_biotypes  / unnormalized / mean_variance_plot
  }
  if(omic == "microrna"){
    pre_filter_plots <- (counts_density | gene_detection) / unnormalized / mean_variance_plot
    
  }
  
  #### POSTFILTERING VISUALIZATIONS ####
  print_centered_note("CREATING POSTFILTERING QC PLTOS")
  
  # 1. DENSITY PLOT
  print("1.Density Plot")
  densities_post <- lapply(1:ncol(lcpm_post), function(i) {
    den <- density(lcpm_post[, i])
    data.frame(x = den$x, y = den$y, sample = as.factor(i))
  })
  density_data_post <- bind_rows(densities_post)
  counts_density_post <- ggplot(density_data_post, aes(x = x, y = y, color = sample)) +
    geom_line(size = 1.5) +  # Trazar las líneas
    scale_color_manual(values = sample_col) +  # Definir los colores
    geom_vline(xintercept = lcpm.cutoff, linetype = "dashed", color = "black") +  # Línea vertical
    labs(
      #title = "B. Clean data",
      x = "Log(CPM)",
      y = "Density"
    ) +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12)
    ) +
    ylim(0, 0.26) +# Limitar el eje Y
    theme(legend.position = "none")
  
  # 2. Gene Detection
  print("2.Gene Detection Plot")
  post_gene_detection <- plot_gene_detection(c)
  
  if(omic == "rna"){
    # 3. Gene Biotype
    print("3.Gene Biotype Plot")
    post_gene_biotypes <- plot_biotypes(c) + 
      theme(axis.text.x = element_text(face="bold", size=4, angle=45),
            legend.position = "bottom") + 
      xlab(NULL) +
      theme(axis.text.x = element_blank())+
      ylab(NULL)
    }
  
  # 4. Mean Variance Plot
  print("4.Mean-Variance Plot")
  if(omic == "rna"){vsd <- vst(c)}
  if(omic == "microrna"){vsd <- varianceStabilizingTransformation(c)}
  
  post_mean_variance_plot <- mean_sd_plot(vsd) + theme(legend.position = "none") + ylab(NULL)
  # 5. Normalization
  print("5.Normalization plot")
  df_lcpm_norm <- as.data.frame(lcpm_post_norm) %>%
    mutate(Gene = rownames(.)) %>%
    pivot_longer(cols = -Gene, names_to = "Sample", values_to = "LogCPM")
  
  normalized <- ggplot(df_lcpm_norm, aes(x = Sample, y = LogCPM, fill = Sample)) +
    geom_boxplot(outlier.shape = NA) +  # No mostrar outliers
    scale_fill_manual(values = sample_col) +  # Colores personalizados
    labs(y = "Log(CPM)", x = "") +
    #theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1,  vjust = -0.025),  # Rotar etiquetas del eje X
      legend.position = "none"  # Ocultar la leyenda si no es necesaria
    ) + 
    ylab(NULL)
  
  print("6.Post-Filtering Low Expression Genes Combined Plot")
  if(omic == "rna"){
    post_filter_plots <- (counts_density_post | post_gene_detection) / post_gene_biotypes /normalized / post_mean_variance_plot
  }
  if(omic == "microrna"){
    post_filter_plots <- (counts_density_post | post_gene_detection) / normalized / post_mean_variance_plot
    
  }
  
  comparison <- (pre_filter_plots | post_filter_plots) +
    plot_annotation(
      title = "Raw Count Matrix QC Plots",
      subtitle = "Pre-filtering (left) vs. Post-filtering (right)"
    ) & 
    theme(plot.subtitle = element_text(size = 14, hjust = 0.5),
          plot.title = element_text(
            family = "lobster", 
            size = 18,
            hjust = 0.5,
            face = "bold", 
            color = "#2a475e"
          ))
  
  print("7.Saving the QC Combined Plot")
  save_ggplot(plot = comparison, title = "QC_RAW_COUNTS_QC_Plots", folder = output_dir,
              width = 5000, height = 3000)
  
  cat("\n")
  print("Saving the Filtered_Counts")
  addWorksheet(report_wb, sheetName = 'Filtered_Counts')
  writeDataTable(report_wb, sheet = "Filtered_Counts", x = as.data.frame(d$counts), rowNames = T)
  
  print("Saving the Log(CPM)_Counts")
  addWorksheet(report_wb, sheetName = "LogCPM_Filtered_Counts")
  writeDataTable(report_wb, sheet = "LogCPM_Filtered_Counts", x = as.data.frame(lcpm_post_norm), rowNames = T)
  
  #### STEP 7: CORRELATION BETWEEN SAMPLES ####
  print_centered_note("COMPUTING CORRELATION BETWEEN SAMPLES")
  c$lib.size <- d$samples$lib.size
  c$norm.factors <- d$samples$norm.factors
  
  # Calculate all correlations 
  cts_corr <- d$counts %>% 
    # we use Spearman's correlation, a non-parametric metric based on ranks
    cor(method = "spearman")

  print("Saving Correlation Between Samples")
  addWorksheet(report_wb, sheetName = "Correlation Between Samples")
  writeDataTable(report_wb, 
                 sheet = "Correlation Between Samples", 
                 x = as.data.frame(cts_corr), 
                 rowNames = T)
  print("1.Correlation Pheatmap")
  # Crear una tabla de anotación sin categorizar variables continuas
  annotation_col <- as.data.frame(colData(c))
  
  # Definir colores continuos personalizados
  annotation_colors <- list(
    lib.size = colorRampPalette(c("white", "blue"))(100),  # Gradiente azul
    norm.factors = colorRampPalette(c("white", "red"))(100)  # Gradiente rojo
  )
  
  # Generar el heatmap con colores continuos
  correlation_samples <- as.ggplot(pheatmap(
    cts_corr,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    #color = colorRampPalette(c("blue", "white", "red"))(50)  # Escala de color del heatmap
  ))
  
  # Lista para almacenar los MDS plots
  mds_plots <- list()
  
  # Crear los MDS plots y almacenarlos en la lista
  print('2.Creating MDSplots')
  for (cat in colnames(metadata)) {
    cat(paste0("\t","MDSplot for ", cat, " variable", "\n"))
    
    main_title <- paste('MDSplot colored by', cat)
    col.group <- factor(metadata[[cat]])  # Convertir a factor
    
    if (nlevels(col.group) > 12) {
      print('MDSplot cannot continue')
      next
    }
    
    # Crear un vector de colores nombrado sin modificar los niveles del factor
    color_palette <- if (nlevels(col.group) == 2) {
      c('#ff1b6b', '#103783')
    } else {
      brewer.pal(nlevels(col.group), 'Dark2')
    }
    
    color_mapping <- setNames(color_palette, levels(col.group))  # Asignar colores a cada nivel
    
    # Crear MDS plot con ggplot2
    mds_data <- plotMDS(lcpm_post_norm, plot = FALSE)
    df_mds <- data.frame(Dim1 = mds_data$x, Dim2 = mds_data$y, 
                         Sample = colnames(lcpm_post_norm), 
                         Group = metadata[[cat]])
    
    p <- ggplot(df_mds, aes(x = Dim1, y = Dim2, color = Group, label = Sample, shape = Group)) +
      geom_point(size = 4) +
      #geom_text(vjust = -1) +
      geom_text_repel() +
      scale_color_manual(values = color_mapping, name = cat) +  # Usamos el mapeo correcto
      labs(title = main_title, x = "Dimension 1", y = "Dimension 2") #+theme_minimal()
    
    mds_plots[[cat]] <- p  # Guardar en la lista
    
    save_ggplot(plot = p, title = paste0("MDSplot_", cat), folder = output_dir, dpi = 300, wid = 3000, height = 2000)
    
  }# For loop key for MDSplots
  
  # Organizar todos los MDS plots en una columna
  mds_grid <- wrap_plots(mds_plots, ncol = 1)  # 1 columna (N filas)
  
  # Unir la figura final en dos columnas: pheatmap a la izquierda y MDS plots a la derecha
  final_plot <- (correlation_samples | mds_grid) + 
    plot_annotation(
    tag_levels = 'A',
    tag_prefix = 'Fig. ',
    title = "Group and Replicate Variability",
    subtitle = "Heatmap of Spearman Correlation between samples (right) \nMDSplot of similarity between samples (left)"
  ) & 
    theme(plot.subtitle = element_text(size = 18, hjust = 0.5),
          plot.title = element_text(
            family = "lobster", 
            size = 22,
            hjust = 0.5,
            face = "bold", 
            color = "#2a475e"
          ),
          plot.tag = element_text(size = 15, face = "bold", family="lobster"))
  
  print("3.Saving the combined plot")
  save_ggplot(plot = final_plot, title = "Heatmap_MDSplot", folder = output_dir,
              width = 10000, height = 6000)
  
  print("Saving QC Report Workbook")
  saveWorkbook(report_wb, file = file.path(output_dir,'QC_Report.xlsx'), overwrite = T)
  print_centered_note('END OF THE QUALITY CONTROL STEP')
  
  return(as.data.frame(d$counts))
}
