save_ggplot <- function(plot, title, folder, width = 2000, height = 3000, dpi = 300, 
                        limitsize = F, units = "px") {
  # Crear la carpeta si no existe
  if (!dir.exists(folder)) {
    dir.create(folder, recursive = TRUE)
  }
  
  # Obtener la lista de archivos existentes en la carpeta
  files <- list.files(folder, pattern = "^\\d+_", full.names = TRUE)
  
  # Extraer los números ya usados
  existing_numbers <- as.numeric(gsub("_(.*)", "", basename(files)))
  
  # Definir el siguiente número consecutivo
  next_number <- ifelse(length(existing_numbers) == 0, 1, max(existing_numbers, na.rm = TRUE) + 1)
  
  # Crear el nombre del archivo con numeración
  filename <- sprintf("%02d_%s.png", next_number, title)
  
  # Ruta completa del archivo
  file_path <- file.path(folder, filename)
  
  # Guardar el gráfico con los parámetros indicados
  ggsave(filename = file_path, 
         plot = plot, 
         dpi = dpi, 
         units = units, 
         width = width, 
         height = height, 
         bg = "white",
         limitsize = limitsize)
  
  # Mensaje de confirmación
  message("Gráfico guardado en: ", file_path)
}
