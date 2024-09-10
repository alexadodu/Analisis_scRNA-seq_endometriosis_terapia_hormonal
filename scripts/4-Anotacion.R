### Author: Madalina Alexandra Dodu

# Script para anotación de clústeres en un objeto Seurat. Este script realiza la 
# anotación de clústeres en un conjunto de datos Seurat combinado. Primero, 
# une las capas del conjunto de datos para la anotación de clústeres y guarda el 
# objeto resultante. Luego, realiza un análisis de marcadores de clústeres utilizando 
# la prueba de Wilcoxon, identifica los principales 10 marcadores y guarda los 
# resultados en un archivo Excel.

## *** PAQUETES ***
install.packages("openxlsx")
install.packages('Seurat')
install.packages("devtools")
install.packages("tidyverse")

library(devtools)
install_github("immunogenomics/presto")

library("openxlsx")
library('Seurat')
library(tidyverse)
library("presto")

## *** ANOTACIÓN DE CLÚSTERES ***

# Combinar diferentes capas del conjunto de datos Seurat en un solo objeto
joined <- JoinLayers(combined_dataset) # Presto da error si no se realiza esta unión

# Guardar el objeto unido como un archivo RDS para su posterior uso
saveRDS(joined, file = "data/processed/joined_clustered.rds")

# Análisis de marcadores de clústeres utilizando la prueba de Wilcoxon. Los grupos se especifican mediante el parámetro `group_by`, que en este caso es "RNA_snn_res.0.4".
markers <- wilcoxauc(joined, group_by = "RNA_snn_res.0.4")

# Identificación de los principales 10 marcadores de cada clúster
# Se seleccionan los 10 principales marcadores que están presentes en al menos el 70% de las observaciones dentro de cada grupo.
top10_markers <- top_markers(markers, n = 10, pct_in_min = 70)

# Guardar los resultados de los principales marcadores en un archivo Excel para su revisión y análisis
write.xlsx(top10_markers, file = "results/top10_markers.xlsx")

# Eliminación de clústeres de mala calidad, en este caso, clúster 15
Idents(object = joined) <- "RNA_snn_res.0.4"
joined <- subset(joined, idents = 15, invert = TRUE)

# Importación de la tabla con la anotación
annot_df <- read.xlsx("results/top10_ann.xlsx", sheet = 2)

# ANOTACIÓN
levels(Idents(joined)) <- as.numeric(levels(Idents(joined)))

if(all(levels(Idents(joined)) %in% annot_df $cluster)) {
  # Crear un vector para los nuevos niveles de acuerdo con el orden en annot_df
  new_levels <- annot_df $anotation_abv[match(levels(Idents(joined)), annot_df$cluster)]
  levels(Idents(joined)) <- new_levels
} else {
  stop("Los niveles de los clústeres no coinciden con el DataFrame.")
}

# Adición de la anotación a una nueva columna denominada `cell_type`
joined$cell_type <- Idents(joined)

# Guarda la lista de objetos Seurat en un archivo RDS para su uso posterior.
saveRDS(joined, file = "data/processed/joined_clustered.rds")