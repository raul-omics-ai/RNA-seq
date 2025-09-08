# 19/03/2025

# Función para calcular la moda
calcular_moda <- function(x) {
  ux <- unique(na.omit(x))
  ux[which.max(tabulate(match(x, ux)))]
}
