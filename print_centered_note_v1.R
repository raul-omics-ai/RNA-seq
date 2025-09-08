print_centered_note <- function(note, width = 80, separator = paste0(rep("-", 79), collapse = "")) {
  # Calcular espacios en blanco para centrar el texto
  white_spaces <- paste0(rep(" ", ((width - nchar(note)) / 2) + 0.1), collapse = "")
  
  # Imprimir mensaje formateado
  cat("\n")
  print(separator)
  print(paste0(white_spaces, note, white_spaces))
  print(separator)
}
