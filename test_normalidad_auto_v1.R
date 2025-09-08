########## 19/05/2025 ##########
# ===================================== #
# ==== Normality distribution test ==== #
# ===================================== #
test_normalidad_auto <- function(df, 
                                 alpha = 0.05, 
                                 plot = TRUE, 
                                 guardar_plot = FALSE, 
                                 carpeta = "plots_normalidad",
                                 individual_folder = "plots_normalidad") {
  source("~/Documentos/09_scripts_R/Automate_Saving_ggplots.R")
  
  if (!is.data.frame(df)) stop("El input debe ser un data.frame.")
  
  if (guardar_plot && !dir.exists(carpeta)) {
    dir.create(carpeta)
  }
  
  if (!requireNamespace("nortest", quietly = TRUE)) {
    install.packages("nortest")
  }
  
  qq_plots_list <- list()
  
  resultados <- lapply(names(df), function(var) {
    x <- df[[var]]
    if (!is.numeric(x)) return(NULL)
    
    x <- na.omit(x)
    n <- length(x)
    
    if (n < 3) return(NULL)
    
    if (n <= 5000) {
      test <- shapiro.test(x)
      metodo <- "Shapiro-Wilk"
    } else {
      test <- nortest::ad.test(x)
      metodo <- "Anderson-Darling"
    }
    
    # Crear ggplot QQ plot
    qq_plot <- ggplot(data.frame(x = x), aes(sample = x)) +
      stat_qq() +
      stat_qq_line(colour = "red") +
      labs(
        title = paste0(var, " (", metodo, ")\np = ", signif(test$p.value, 3)),
        x = "Theoretic Quantile", y = "Sample Quantile"
      ) +
      theme_minimal()
    
    save_ggplot(plot = qq_plot, title = paste0("QQ_Plot_", var, "_", metodo), 
                folder = individual_folder, width = 2000, height = 2000)
    
    # Retornar ambos resultados y plot
    list(
      summary = data.frame(
        Variable = var,
        N = n,
        Method = metodo,
        `p-valor` = test$p.value,
        Interpretation = ifelse(test$p.value > alpha, 
                                "Normally Distributed", 
                                "Not Normally Distributed")
      ),
      plot = qq_plot
    )
  })
  
  # Separar resultados y plots
  resultados_df <- do.call(rbind, lapply(resultados, function(x) x$summary))
  qq_plots_list <- lapply(resultados, function(x) x$plot)
  names(qq_plots_list) <- names(df)[!sapply(qq_plots_list, is.null)]
  
  rownames(resultados_df) <- NULL
  
  # Si se quiere mostrar todos los QQ plots juntos
  if (plot && length(qq_plots_list) > 0) {
    ncol_fig <- min(3, length(qq_plots_list))
    qq_comb <- patchwork::wrap_plots(qq_plots_list, ncol = ncol_fig)
    save_ggplot(qq_comb, title = "QQ_Plots", folder = carpeta, width = 5000, height = 3000)
  }
  return(resultados_df)
}
