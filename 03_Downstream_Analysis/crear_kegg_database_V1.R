# All KEGG ID pathways and it's description
procesar_kegg_pathways <- function(path_txt,
                                   output_csv = NULL,
                                   organismo = "hsa") {
  
  library(KEGGREST)
  
  # Instrucciones para descargar el archivo .txt con la información de KEGG para una especie
  # determinada:
  # Insertar en el buscador el siguiente link y cambiar el nombre de la especie por el código 
  # de tres letras que le corresponda:
  # https://rest.kegg.jp/list/pathway/hsa
  
  # Leer archivo con pathways
  kegg_names <- read.delim(path_txt, header = FALSE)
  colnames(kegg_names) <- c("ID", "Pathway")
  
  # Obtener enlaces entre genes y pathways para el organismo
  mapping <- keggLink("pathway", organismo)
  
  # Crear dataframe con ENTREZID y SYMBOL
  if(organismo == "hsa"){
    library(org.Hs.eg.db)
    org.db <- org.Hs.eg.db
  }
  if(organismo == "mmu"){
    library(org.Mm.eg.db)
    org.db <- org.Mm.eg.db
  }else{stop("Please download the package from the db for your specie")}
  
  df <- data.frame(
    Pathway = gsub("path:", "", mapping),
    ENTREZID = gsub(paste0(organismo, ":"), "", names(mapping)),
    SYMBOL = AnnotationDbi::mapIds(
      org.db,
      keys = gsub(paste0(organismo, ":"), "", names(mapping)),
      column = "SYMBOL",
      keytype = "ENTREZID",
      multiVals = "first"  # O "list" si quieres mantener todos los símbolos
    ),
    stringsAsFactors = FALSE
  )
  
  # Agrupar por pathway y concatenar genes
  df <- df %>%
    group_by(Pathway) %>%
    summarise(
      ENTREZID = paste(ENTREZID, collapse = ", "),
      SYMBOL = paste(SYMBOL, collapse = ", "),
      .groups = "drop"
    )
  
  # Unir con los nombres de pathways
  kegg_db <- kegg_names %>%
    left_join(df, by = c("ID" = "Pathway"))
  
  # Guardar si se indica un archivo de salida
  if (!is.null(output_csv)) {
    write.csv(kegg_db, file = output_csv, row.names = FALSE)
  }
  
  return(kegg_db)
}