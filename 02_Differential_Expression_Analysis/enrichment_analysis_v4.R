########## 12/11/2024 ##########

# This scripts will be focus in the enrichment analysis with GO terms, KEGG, MSigDB, Wikipathways and GSEA analysis. I will explore also some
# visualizations.

# For enrichment analysis, only can take care of genes with SYMBOL, so all those ENSEMBL IDs without SYMBOL,  won't be taken into consideration
enrichment_function <- function(dge_df, 
                                symbol_colname = 'Symbol', 
                                logfc_colname = 'logFC', 
                                adj_p_val_colname = "adj.P.Val", 
                                where_to_save = NULL, 
                                title = "Enrichment_analysis",
                                organism = 'mmu',  
                                ontologies = c("BP", "CC", "MF"), 
                                kegg = TRUE, 
                                wikip = TRUE,
                                kegg_pathview = FALSE,
                                universe = TRUE){
  ########## loading packages and custom funcitons ##########
  source("~/Documentos/09_scripts_R/print_centered_note_v1.R")
  source("~/Documentos/09_scripts_R/Automate_Saving_ggplots.R")
  print_centered_note('LOADING PACKAGES ')
  list_of_packages = c("openxlsx", "limma", "edgeR", "dplyr", "AnnotationDbi", "clusterProfiler",
                       "KEGG.db", "enrichplot", "graphics", "patchwork", "ggnewscale", "fgsea", 
                       "GseaVis", "gridExtra", "pathview", "ggplot2", "patchwork", "ggradar", "tidyr")
  
  new_packages = list_of_packages[!(list_of_packages %in% installed.packages())]
  if(length(new_packages) > 0){install.packages(new_packages)}
  
  invisible(lapply(list_of_packages, FUN = library, character.only = T))
  rm(list_of_packages, new_packages)
  print('Packages loaded successfully.')
  
  ########## Configuración del report ##########
  print_centered_note("Setting Up Enrichment Report")
  wb <- createWorkbook()
  addWorksheet(wb, sheetName = "Overview")
  
  report <- data.frame(matrix(ncol=1,nrow = 0))
  colnames(report) = "Report"
  
  report[nrow(report)+1,1] <- paste0("Number of genes in the dataset: ", nrow(dge_df))
  
  ########## Configuración del directorio de salida ##########
  where_to_save <- ifelse(is.null(where_to_save), getwd(), where_to_save)
  output_dir <- file.path(where_to_save, title)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  ########## step 0:check if there are NAs ##########
  print("Removing those rows without SYMBOL")
  symbol_idx <- which(colnames(dge_df) == symbol_colname)
  logfc_idx <- which(colnames(dge_df) == logfc_colname)
  pval_adj_idx <- which(colnames(dge_df) == adj_p_val_colname)
  
  clean_df <- dge_df[!is.na(dge_df[[symbol_idx]]), ]
  report[nrow(report)+1,1] <- paste0("Number of genes with SYMBOL: ", nrow(clean_df))
  
  if(all(!is.na(clean_df[[symbol_idx]]))){
    print("All rows without SYMBOLS have been removed")
  }
  
  print("Filtering no significant genes & Sorting the dataframe by log2foldchange")
  clean_df <- clean_df[which(clean_df[[pval_adj_idx]] < 0.05 & abs(clean_df[[logfc_idx]]) >= 1),]
  clean_df <- clean_df[order(clean_df[[pval_adj_idx]]), ]
  report[nrow(report)+1,1] <- paste0("Number of DE genes: ",nrow(clean_df))
  if(nrow(clean_df) == 0){
    stop("There are any DE genes, so can't continue with the enrichment analysis")
  }
  
  #writeData(wb, sheet = "Overview", x = report)
  
  ########## step 1: selecting genes Enrichment analysis ##########
  print("Starting with the Enrichment analysis")
  deGenes <- clean_df[[symbol_idx]] 
  geneUniverse <- na.omit(dge_df[[symbol_idx]])
  
  # Genes must be annotated in ENTREZ ID
  dir.create(path = file.path(output_dir, "02_GSEA_plots"), recursive = T, showWarnings = F)
  
  if(organism == 'mmu'){
    library(org.Mm.eg.db)
    print("Converting SYMBOLS in ENTREZ ID in mouse")
    org.db = org.Mm.eg.db
  }
  if(organism == "hsa"){
    library(org.Hs.eg.db)
    print("Converting SYMBOLS in ENTREZ ID in human")
    org.db = org.Hs.eg.db
  }
  
  #### go gene presparation 
  entrez_genes <- AnnotationDbi::mapIds(org.db, keys = clean_df[[symbol_idx]], column = "ENTREZID", keytype = "SYMBOL")
  if(universe){
    entrez_universe <- unlist(AnnotationDbi::mapIds(org.db, keys = dge_df[[symbol_idx]], column = "ENTREZID", keytype = "SYMBOL"))
  } else{
    entrez_universe = NULL
  }
  
  #### gsea gene preparation 
  original_gene_list <- dge_df[[logfc_idx]] 
  names(original_gene_list) <- dge_df[[symbol_idx]] 
  gene_list = sort(original_gene_list, decreasing = TRUE) # sort the list in decreasing order (required for clusterProfiler)
  
  kegg_list <- dge_df[[logfc_idx]]
  names(kegg_list) <- AnnotationDbi::mapIds(x = org.db, keys = dge_df[[symbol_idx]], column = "ENTREZID", keytype = "SYMBOL") 
  kegg_list <- sort(kegg_list, decreasing = T)
  
  ##################### GRAPHS GSE GO x ChatGPT
  if (!dir.exists(file.path(output_dir, "01_Individual_plots"))) {
    dir.create(file.path(output_dir, "01_Individual_plots"))
  }
  
  enrichment_plots <- list()
  
  for (ontology in ontologies) {
    message(paste0("Running GO enrichment for ", ontology))
    
    # ---- Enrichment GO ----
    print_centered_note(paste0("GO ENRICHMENT ANALYSIS FOR ONTOLOGY ", ontology, " "))
    ge.go <- enrichGO(gene = entrez_genes,
                      universe = entrez_universe,
                      OrgDb         = org.db,
                      ont           = ontology,
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)
    ge.go@result$Description  <- ifelse(
      nchar(ge.go@result$Description) > 47,  # Si la longitud original es mayor a 47
      paste0(substr(ge.go@result$Description, 1, 47), "..."),  # Trunca y añade "..."
      ge.go@result$Description  # Si no, deja el texto igual
    )
    
    tab.ge.go <- as.data.frame(ge.go)
    tab.ge.go$Description <- ifelse(
      nchar(tab.ge.go$Description) > 47,  # Si la longitud original es mayor a 47
      paste0(substr(tab.ge.go$Description, 1, 47), "..."),  # Trunca y añade "..."
      tab.ge.go$Description  # Si no, deja el texto igual
    )
    
    # Guardar resultados en Excel
    addWorksheet(wb, sheetName = paste0("GO_Enrichment_", ontology))
    writeDataTable(wb, sheet = paste0("GO_Enrichment_", ontology), x = tab.ge.go)
    
    report[nrow(report)+1,1] <- paste0("Number of GO terms in ontology ",ontology, ": ",nrow(tab.ge.go))
    
    # --- Guardar los gráficos de GO individualmente ----
    ans.go2 <- pairwise_termsim(ge.go)
    
    print("1.Enriched Map")
    set.seed(123)
    enrichmap_go <- emapplot(ans.go2, layout="kk", group_category = T, group_legend = T, 
                             cex_label_group = 1.5, node_label = "group", nCluster = 5,
                             repel = T)
    
    save_ggplot(plot =  enrichmap_go, 
                title = paste0("GO_EnrichMapClustered_", ontology), 
                folder = file.path(output_dir, "01_Individual_plots"), width = 14, height = 10,
                units = "in")
    #ggsave(file.path(output_dir, "01_Individual_plots", paste0("01_GO_EnrichMapClustered_", ontology, ".png")), 
    #       enrichmap_go, dpi = 300, width = 14, height = 10)
    
    #print("1.Treeplot") # La función treeplot ha dejado de funcionar con ggplot >= 3.4
    #tree_go <- treeplot(ans.go2, 
    #                    showCategory = 10, 
    #                    fontsize = 8) + 
    #  theme(legend.position = "left")
    #
    #ggsave(file.path(output_dir, "01_Individual_plots", paste0("01_Treeplot_GO_", ontology, ".png")),
    #       tree_go, dpi = 300, width = 18, height = 14)
    
    print("2.Barplot")
    barplot_go <- barplot(ge.go, showCategory = 20) + 
      ggtitle(paste0("GO Terms Barplot - ", ontology)) + 
      theme(axis.text.y = element_text(size = 14))  # Aumenta el tamaño de las etiquetas del eje Y
    
    save_ggplot(plot = barplot_go, 
                title = paste0("GO_Barplot_", ontology), 
                folder = file.path(output_dir, "01_Individual_plots"), width = 10, height = 6, 
                units = "in")
    #ggsave(file.path(output_dir, "01_Individual_plots", paste0("02_GO_Barplot_", ontology, ".png")), 
    #        barplot_go, dpi = 300, width = 10, height = 6)
    
    print("3.Dotplot")
    dotplot_go <- dotplot(ge.go, showCategory = 20, orderBy = "p.adjust", decreasing = FALSE) + 
      ggtitle(paste0("GO Terms Dotplot - ", ontology)) + 
      theme(axis.text.y = element_text(size = 14))  # Aumenta el tamaño de las etiquetas del eje Y
    
    save_ggplot(plot = dotplot_go, 
                title = paste0("GO_Dotplot_", ontology), 
                folder = file.path(output_dir, "01_Individual_plots"), 
                width = 10, height = 6, 
                units = "in")
    #ggsave(file.path(output_dir, "01_Individual_plots", paste0("03_GO_Dotplot_", ontology, ".png")),
    #       dotplot_go, dpi = 300, width = 10, height = 6)
    
    print("4.Combined Plot")
    enriched_comb_go <- enrichmap_go + (barplot_go / dotplot_go) + 
      plot_annotation(tag_levels = 'A', 
                      tag_prefix = "Fig. ",
                      title = paste0("Summary of GO ", ontology, " Enriched Terms"),
                      subtitle = "EnrichMap of enriched terms (Fig. A) \n Barplot of # genes of enriched terms (Fig. B) \nDotplot of gene ratio of enriched terms (Fig. C)"
                      ) & theme(plot.tag = element_text(size = 16),
                                plot.subtitle = element_text(size = 18, hjust = 0.5),
                                plot.title = element_text(
                                  family = "lobster", 
                                  size = 20,
                                  hjust = 0.5,
                                  face = "bold", 
                                  color = "#2a475e"))
      enrichment_plots[[ontology]] <- enriched_comb_go
    
    # Guardar figura combinada GO
    save_ggplot(plot = enriched_comb_go, 
                  title = paste0("GO_", ontology, "_enrichment_plot"), 
                  folder = output_dir, 
                  width = 35, height = 20, 
                  units = "in")
    #ggsave(file.path(output_dir, paste0("01_GO_", ontology, "_enrichment_plot.png")), 
    #       enriched_comb_go, dpi = 300, width = 35, height = 20)
    
    # ---- GSEA Analysis ----
    print_centered_note(paste0("GSEA ONTOLOGY ", ontology," "))
    gse <- clusterProfiler::gseGO(geneList=gene_list, 
                                  ont = ontology, 
                                  keyType = "SYMBOL", 
                                  #nPerm = 10000, 
                                  minGSSize = 3, 
                                  eps = 0,
                                  maxGSSize = 800, 
                                  pvalueCutoff = 0.05,
                                  seed = T,
                                  verbose = TRUE, 
                                  OrgDb = org.db,
                                  by="fgsea")
    
    tab.gse.go <- as.data.frame(gse)
    if (nrow(tab.gse.go) > 0) {
      print(paste0("Saving GO GSEA termos for ", ontology))
      addWorksheet(wb, sheetName = paste0("GO_GSEA_", ontology))
      writeDataTable(wb, sheet = paste0("GO_GSEA_", ontology), x = tab.gse.go)
      
      report[nrow(report)+1,1] <- paste0("Number of GSEA terms in ontology ", ontology, ": ", nrow(tab.gse.go))
      
      # GSE dotplot
      print("1.GSE-Dotplot")
      dotplot_gse <- dotplot(gse, showCategory = 10, split = ".sign") + 
        facet_grid(. ~ .sign) + 
        theme(axis.text.y = element_text(size = 14),  # Aumenta tamaño etiquetas eje Y
              strip.text = element_text(size = 14))  # Aumenta tamaño etiquetas de facet
      
      save_ggplot(plot = dotplot_gse, 
                  title = paste0("GSEA_Dotplot_", ontology), 
                  folder = file.path(output_dir, "01_Individual_plots"), 
                  width = 10, height = 6, 
                  units = "in")
      #ggsave(file.path(output_dir, "01_Individual_plots", paste0("04_GSEA_Dotplot_", ontology, ".png")),
      #       dotplot_gse, dpi = 300, width = 10, height = 6)
      
      # GSE barplot
      print("2.GSE-Barplot")
      top_paths <- tab.gse.go[order(tab.gse.go$NES, decreasing = TRUE), ]
      top_up <- head(top_paths, 10)
      top_down <- tail(top_paths, 10)
      top_paths <- rbind(top_up, top_down)
      
      top_paths$Regulation <- ifelse(top_paths$NES > 0, "Upregulated", "Downregulated")
      
      top_paths_gg <- ggplot(top_paths, aes(reorder(Description, NES), NES)) +
        geom_col(aes(fill = Regulation)) +
        scale_fill_manual(values = c("Upregulated" = "#006e90", "Downregulated" = "#f18f01")) +
        coord_flip() +
        labs(x = "Pathway", y = "Normalized Enrichment Score") +
        theme(axis.text.y = element_text(size = 14))
      
      save_ggplot(plot = top_paths_gg, 
                  title = paste0("GSE_Barplot_", ontology), 
                  folder = file.path(output_dir, "01_Individual_plots"), 
                  width = 10, height = 6, 
                  units = "in")
      #ggsave(file.path(output_dir, "01_Individual_plots", 
      #                 paste0("05_GSE_Barplot_", ontology, ".png")),
      #       top_paths_gg, dpi = 300, width = 10, height = 6)
      
      #### volcanoplot
      print("3.GSE-VolcanoPlot")
      source("~/Documentos/09_scripts_R/volcano_gsea_v1.R")
      volcano_gsea <- volcanoGsea_v1(data = gse, pvalue.cutoff = 0.05)
      
      save_ggplot(plot = volcano_gsea, 
                  title = paste0("GSE_VolcanoPlot_", ontology), 
                  folder = file.path(output_dir, "01_Individual_plots"), 
                  width = 10, height = 6, 
                  units = "in")
      #ggsave(file.path(output_dir, "01_Individual_plots", paste0("06_GSE_VolcanoPlot_", ontology, ".png")), 
      #       volcano_gsea,dpi = 300, width = 10, height = 6)
      
      #### combined plot
      print("4.GSE-Combined Plot")
      enriched_comb_gse <- volcano_gsea / (dotplot_gse + top_paths_gg)  + 
        plot_annotation(tag_levels = 'A', 
                        tag_prefix = "Fig. ",
                        title = paste0("Gene Set Enrichment Analysis of ", ontology, " Terms"),
                        subtitle = "VolcanoPlot (Fig. A) \nDotplot of gene ratios of enriched terms (Fig. B) \nBarplot of NES of enriched terms (Fig. C)"
        ) & theme(plot.tag = element_text(size = 16),
                    plot.subtitle = element_text(size = 18, hjust = 0.5),
                    plot.title = element_text(
                      family = "lobster", 
                      size = 20,
                      hjust = 0.5,
                      face = "bold", 
                      color = "#2a475e"))
      enrichment_plots[[paste0(ontology, "_GSEA")]] <- enriched_comb_gse
      
      save_ggplot(plot = enriched_comb_gse, 
                  title = paste0("GSEA_GO_", ontology, "_enrichment_plot"), 
                  folder = output_dir, 
                  width = 20, height = 13, 
                  units = "in")
      #ggsave(file.path(output_dir, paste0("02_GSEA_GO_", ontology, "_enrichment_plot.png")), 
      #       enriched_comb_gse, dpi = 300, width = 20, height = 13)
      #### gsea plots
      print("Saving TOP 10 Upregulated and Downregulated pathways NES from GO-GSEA")
      for(path in 1:nrow(top_paths)){
        print(paste0("Saving GSEA plot for ", top_paths$Description[path], " GO term"))
        png(file.path(output_dir, paste0("02_GSEA_plots/GSEAPLOT_",top_paths$Description[path] ,"_GO_", ontology,  ".png")), res = 300, width = 4000, height = 2000)
        print(gseaplot2(gse, geneSetID = top_paths$ID[path], title = top_paths$Description[path], pvalue_table = T))
        dev.off()
      } # Foor loop key for generate GSEA plots
    } # If key for GSE visualizations
    else{print("Any pathway was recovered in GSEA")}
  } # For loop for ontologies 
  
  #### KEGG GSEA ####
  print_centered_note("PERFORMING GSEA ON KEGG")
  gse.kegg <- gseKEGG(geneList=kegg_list, 
                      keyType = "ncbi-geneid", 
                      nPerm = 10000, 
                      minGSSize = 3, 
                      maxGSSize = 800, 
                      pvalueCutoff = 0.05, 
                      verbose = TRUE, 
                      eps = 0,
                      organism = organism, 
                      pAdjustMethod = "BH",
                      use_internal_data = T)
  
  tab.kegg <- as.data.frame(gse.kegg)
  
  print("Saving KEGG GSEA resutls")
  addWorksheet(wb, sheetName = "KEGG_GSEA")
  writeDataTable(wb, sheet = "KEGG_GSEA", x = tab.kegg)
  report[nrow(report)+1,1] <- paste0("Number of KEGG terms: ", nrow(tab.kegg))
  
  if(nrow(tab.kegg) > 0){
    # KEGG Visualizations ----
    print_centered_note("CREATING KEGG VISUALIZATIONS ")
    
    # --- KEGG dotplot
    print("1.KEGG GSEA-Dotplot")
    KEGG_plot <- dotplotGsea(data = gse.kegg,topn = 10, pajust = 0.05)$plot
    save_ggplot(plot = KEGG_plot, 
                title = "KEGG_GSEA_Dotplot", 
                folder = output_dir, 
                width = 12, height = 6, 
                units = "in")
    #ggsave(filename = file.path(output_dir, "01_Individual_plots","08_KEGG_GSEA_Dotplot.png"), 
    #       plot = KEGG_plot,
    #       dpi = 300, width = 12, height = 6)
    
    print("2.KEGG GSEA-Volcanoplot")
    volcano_kegg <- volcanoGsea_v1(data = gse.kegg, pvalue.cutoff = 0.05)
    
    save_ggplot(plot = volcano_kegg, 
                title = "KEGG_GSE_VolcanoPlot", 
                folder = file.path(output_dir, "01_Individual_plots"), 
                width = 10, height = 6, 
                units = "in")
    #ggsave(file.path(output_dir, "01_Individual_plots", paste0("09_KEGG_GSE_VolcanoPlot.png")), 
    #       volcano_kegg,dpi = 300, width = 10, height = 6)
    
    print("3.KEGG GSEA Combined Plot")
    
    kegg_combined <- volcano_kegg /KEGG_plot +
      plot_annotation(tag_levels = 'A', 
                      tag_prefix = "Fig. ",
                      title = "Gene Set Enrichment Analysis of KEGG Terms",
                      subtitle = "VolcanoPlot (Fig. A) \nDotplot of gene ratios of enriched terms (Fig. B) \nBarplot of NES of enriched terms (Fig. C)"
                      ) & theme(plot.tag = element_text(size = 16),
                                plot.subtitle = element_text(size = 18, hjust = 0.5),
                                plot.title = element_text(
                                  family = "lobster", 
                                  size = 20,
                                  hjust = 0.5,
                                  face = "bold", 
                                  color = "#2a475e"))
    save_ggplot(plot = kegg_combined, 
                title = "KEGG_GSEA_CombinedPlot", 
                folder = output_dir, 
                width = 12, height = 10, 
                units = "in")
    #ggsave(file.path(output_dir, "03_KEGG_GSEA_CombinedPlot.png"), kegg_combined, 
    #       dpi = 300, width = 12, height = 10)
    
    if(kegg_pathview){
      dir.create(path = file.path(output_dir, "03_KEGG_paths"), recursive = T, showWarnings = F)
      wd <- getwd()
      setwd(file.path(output_dir, "03_KEGG_paths"))
      #### pathview
      for(pathway in 1:nrow(tab.kegg)){
        pathway_name <- tab.kegg$Description[pathway]
        print(paste0("Saving KEGG pathway ", pathway_name))
        pathview_plot <- pathview(gene.data = kegg_list, 
                                  pathway.id = tab.kegg$ID[pathway], 
                                  species = organism)
        
        png(filename = file.path(output_dir, "03_KEGG_paths", paste0("KEGG_", pathway_name, ".png", collapse = "_")), 
            res = 300, width = 3000, height = 3000)
        print(pathview_plot)
        dev.off()
      }# Pathview path plot generator key
      setwd(wd)
    }# check if pathview generatior plot key
  }# check key for kegg_gsea resutls
  
  # ---- Wikipathways GSEA ---- 
  print_centered_note("PERFORMING GSEA ON WIKIPATHWAYS")
  if(organism == "mmu"){
    wp.organism = "Mus musculus"
  }
  if(organism == "hsa"){
    wp.organism = "Homo sapiens"
  }
  gse.wp <- gseWP(kegg_list, organism = wp.organism, eps = 0)
  readable.wp <- setReadable(gse.wp, org.db, keyType="ENTREZID")
  
  tab.wp <- as.data.frame(readable.wp)
  print("Saving WK GSEA resutls")
  addWorksheet(wb, sheetName = "WP_GSEA")
  writeDataTable(wb, sheet = "WP_GSEA", x = tab.wp)
  
  report[nrow(report)+1,1] <- paste0("Number of WP terms: ", nrow(tab.wp))
  writeData(wb, sheet = "Overview", x = report)
  # ---- WikiPathways Visualizations ----
  if(nrow(tab.wp) > 0){
    print_centered_note("CREATING WP VISUALIZATIONS ")
    #### plot NES plot
    print("1.WP GSEA-Barplot")
    tab.wp <- tab.wp[order(tab.wp$NES, decreasing = TRUE), ]
    top_up <- head(tab.wp, 10)
    top_down <- tail(tab.wp, 10)
    
    top_paths <- rbind(top_up, top_down)
    
    top_paths$Regulation <- ifelse(top_paths$NES > 0, "Upregulated", "Downregulated")
    
    barplot_wp <- ggplot(top_paths, aes(reorder(Description, NES), NES)) +
      geom_col(aes(fill = Regulation)) +
      scale_fill_manual(values = c("Upregulated" = "#006e90", "Downregulated" = "#f18f01")) +
      coord_flip() +
      labs(x = "Pathway", y = "Normalized Enrichment Score")+ 
      theme(axis.text.y = element_text(size = 14))
    
    save_ggplot(plot = barplot_wp, 
                title = "WP_GSEA_Barplot", 
                folder = file.path(output_dir, "01_Individual_plots"), 
                width = 10, height = 6, 
                units = "in")
    #ggsave(filename =file.path(output_dir, "01_Individual_plots", "11_WP_GSEA_Barplot.png"), 
    #       barplot_wp,
    #       dpi = 300,  width = 10, height = 6)
    
    # WP dotplot
    print("2.WP GSEA-Dotplot")
    wp_dotplot <- dotplotGsea(data = gse.wp,topn = 10, pajust = 0.05)$plot
    
    save_ggplot(plot = wp_dotplot, 
                title = "WP_GSEA_Dotplot", 
                folder = file.path(output_dir, "01_Individual_plots"), 
                width = 10, height = 6, 
                units = "in")
    #ggsave(filename = file.path(output_dir, "01_Individual_plots", "10_WP_GSEA_Dotplot.png"), 
    #       wp_dotplot, dpi = 300, height = 6, width = 10)
    
    # volcanoplot
    print("3.WP GSEA-VolcanoPlot")
    volcano_wp <- volcanoGsea_v1(data = gse.wp, pvalue.cutoff = 0.05)
    save_ggplot(plot = volcano_wp, 
                title = "WP_GSEA_Volcanoplot", 
                folder = file.path(output_dir, "01_Individual_plots"), 
                width = 10, height = 10, 
                units = "in")
    #ggsave(filename = file.path(output_dir, "01_Individual_plots", "12_WP_GSEA_Volcanoplot.png"),
    #       volcano_wp, dpi = 300, height = 10, width = 10)
    
    # Combined plot
    print("4.WP Combined Plot")
    wp_combined <-  volcano_wp / (wp_dotplot + barplot_wp)  +
      plot_annotation(tag_levels = 'A', 
                      tag_prefix = "Fig. ",
                      title = "Gene Set Enrichment Analysis of WP Terms",
                      subtitle = "VolcanoPlot (Fig. A) \nDotplot of gene ratios of enriched terms (Fig. B) \nBarplot of NES of enriched terms (Fig. C)"
    ) & theme(plot.tag = element_text(size = 16),
              plot.subtitle = element_text(size = 18, hjust = 0.5),
              plot.title = element_text(
                family = "lobster", 
                size = 20,
                hjust = 0.5,
                face = "bold", 
                color = "#2a475e"))

    save_ggplot(plot = wp_combined, 
                title = "WikiPathway_GSEA_Combined_Plot", 
                folder = output_dir, 
                width = 20, height = 13, 
                units = "in")
    #ggsave(file.path(output_dir, "04_WikiPathway_GSEA_Combined_Plot.png"), 
    #       wp_combined, dpi = 300, width = 20, height = 13)
    
  } # If key for wikipathway check for visualizations
  
  cat("\n")
  print("Saving the Excel file")
  saveWorkbook(wb, file = file.path(output_dir, paste0("Enrichment_GSEA_Report_",title,".xlsx")))
  print_centered_note("ENDF OF THE ENRICHMENT ANALYSIS SCRIPT ")
} # Llave del final de la función
