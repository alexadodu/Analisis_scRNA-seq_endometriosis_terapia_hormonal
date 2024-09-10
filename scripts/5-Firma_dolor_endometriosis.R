### Author: Madalina Alexandra Dodu

# Script para generar una firma génica asociada al dolor en dos conjuntos de 
# datos de scRNA-seq (A y B) y caracterización de tipos celulares con firma 
# del dolor y genes asociados a este síntoma. El objetivo es identificar los 
# genes relacionados con el dolor presentes en diferentes tipos celulares 
# dentro de estos datasets, utilizando la puntuación de módulos UCell. Además, 
# se realiza un análisis de expresión diferencial entre los tipos celulares clave 
# y se extraen los genes más significativos de cada uno. El script también procesa 
# los datos para escalar las expresiones y eliminar clústeres, lo que optimiza la 
# detección de los genes de dolor en tipos celulares específicos. Los resultados 
# obtenidos incluyen los 5 genes más representativos para cada tipo celular de 
# ambos datasets, que pueden estar involucrados en la señalización de dolor en el 
# contexto de la endometriosis.

## *** PAQUETES ***
install.packages("openxlsx")
install.packages('Seurat')
install.packages("tidyverse")

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("UCell")

library("openxlsx")
library('Seurat')
library(tidyverse)
library(UCell)
library(RColorBrewer)

## *** GENERACIÓN DE FIRMA DEL DOLOR ***

# Cargar genes asociados al dolor desde un archivo Excel, dividiendo aquellos genes separados por "and"
pain_genes <- read.xlsx("data/raw/pain_genes.xlsx", sheet = 1)
pain_genes <- pain_genes$Gene
pain_genes <- strsplit(pain_genes, "and")  # Separación de genes en una lista
pain_genes <- unlist(pain_genes)  # Descomponer en un vector
pain_genes <- trimws(pain_genes)  # Remover espacios en blanco innecesarios

# Se añaden genes los adicionales de Wistrom et. al 2022 
aditional_genes <- c("SRP14", "BMF", "GDAP1", "MLLT10", "BSN", "FSHB")
pain_genes <- c(pain_genes, aditional_genes)

# Se asegura que no haya duplicados y que todos los genes estén en mayúsculas
pain_genes <- sort(pain_genes)
pain_genes <- unique(pain_genes)
pain_genes <- toupper(pain_genes)

# Encontrar los genes relacionados con el dolor que están presentes en los datasets A y B
A_pain_genes_in_seurat <- intersect(pain_genes, rownames(A_joined))  # Comparar genes con las lecturas de A
B_pain_genes_in_seurat <- intersect(pain_genes, rownames(B_joined))  # Comparar genes con las lecturas de B

## *** PUNTUACIÓN DEL DOLOR CON UCELL ***

# Aplicar UCell para calcular las puntuaciones del módulo de dolor basado en los genes seleccionados
A_scores <- AddModuleScore_UCell(A_joined, 
                                 features=list(A_pain_genes_in_seurat),  # Genes relacionados al dolor
                                 name = "PainSignature_UCell")

# Aplicar el mismo análisis para el dataset B
B_scores <- AddModuleScore_UCell(B_joined, 
                                 features=list(B_pain_genes_in_seurat), 
                                 name = "PainSignature_UCell")

# Crear una nueva columna que combine las condiciones y el estado de tratamiento, para análisis posteriores
A_scores@meta.data$condition_treatment <- paste(A_scores@meta.data$condition, 
                                                A_scores@meta.data$TreatmentState, 
                                                sep = "_")
B_scores@meta.data$condition_treatment <- paste(B_scores@meta.data$condition, 
                                                B_scores@meta.data$TreatmentState, 
                                                sep = "_")

## VISUALIZACIÓN DE FIRMA DEL DOLOR/TIPO CELULAR
A1 <- VlnPlot(A_scores, features = "signature_1PainSignature_UCell",
              pt.size=0, cols = cell_colors) +
  ggtitle("Expresión de la firma de genes del dolor por tipo celular (Dataset A)") +
  theme_article() +
  RotatedAxis() +
  labs(x = "Tipo celular")

B1 <- VlnPlot(B_scores, features = "signature_1PainSignature_UCell", 
              pt.size=0, cols = cell_colors, raster=FALSE) +
  ggtitle("Expresión de la firma de genes del dolor por tipo celular (Dataset B)") +
  theme_article() +
  RotatedAxis() +
  labs(x = "Tipo celular")

## *** CARACTERIZACIÓN DE TIPOS CELULARES CON FIRMA DE DOLOR – DATASET A ***

# Definir los tipos celulares para el análisis de expresión diferencial de genes 
# relacionados con el dolor en dataset A
A_painGenes_typeCells_expression <- c("Cél. estromales endom. (sec.)",  # Subtipos celulares clave
                                      "Monocitos", 
                                      "Macrófagos",
                                      "Fibroblastos estrom.")

# Inicializar una lista para almacenar los resultados de la expresión diferencial por tipo celular
A_DifferentialResults_pain <- list()

# Realizar la escalación de datos tras la selección de las características más variables, 
# después de haber eliminado anteriormente los clústeres de baja calidad
A_scores <- FindVariableFeatures(A_scores, selection.method = "vst", nfeatures = 2000)
A_scores <- ScaleData(A_scores, features = rownames(A_scores))

# Guardar el objeto con las puntuaciones actualizadas
saveRDS(A_scores, file = "data/processed/A_pain_scores_UCell.rds")

## *** CARACTERIZACIÓN DE GENES DE LA VÍA DEL DOLOR – DATASET A ***

# Bucle que realiza el análisis de expresión diferencial por cada tipo celular en el dataset A
for (cell_type in A_painGenes_typeCells_expression) {
  A_DifferentialResults_pain[[cell_type]] <- FindMarkers(A_scores, ident.1 = cell_type, 
                                                         features = A_pain_genes_in_seurat)
}

# Extraer los 5 genes más significativos de cada tipo celular
top_genes_pain_A <- unique(unlist(lapply(A_DifferentialResults_pain, 
                                         function(x) head(rownames(x[order(x$p_val_adj),]), 5))))

# Filtrar las células de interés para los análisis posteriores y escalar los datos de los genes seleccionados
A_scores_Pain_Cells <- subset(A_scores, idents = A_painGenes_typeCells_expression)
A_scores_Pain_Cells <- ScaleData(A_scores_Pain_Cells, features = top_genes_pain_A)


## *** CARACTERIZACIÓN DE TIPOS CELULARES CON FIRMA DE DOLOR – DATASET B ***

# Repetir los pasos anteriores para el dataset B, ajustando a los tipos celulares específicos de B
Idents(B_scores) <- "cell_type"
B_painGenes_typeCells_expression <- c("Cél. estromales endom. (sec.)",  
                                      "Cél. estrom. decid.", 
                                      "Cél. estromales endom. (menstr.)",  
                                      "Cél. estromales endom. (prolif.)",  
                                      "Monocitos", 
                                      "Fibroblastos estrom. C7+",
                                      "Fibroblastos estrom.")

# Inicializar lista para almacenar resultados del dataset B
B_DifferentialResults_pain <- list()

# Escalar los datos del dataset B tras la selección de las características más variables
B_scores <- FindVariableFeatures(B_scores, selection.method = "vst", nfeatures = 2000)
B_scores <- ScaleData(B_scores)

# Guardar el objeto actualizado para dataset B
saveRDS(B_scores, file = "data/processed/B_pain_scores_UCell.rds")

## *** CARACTERIZACIÓN DE GENES DE LA VÍA DEL DOLOR – DATASET B ***

# Bucle para realizar el análisis de expresión diferencial en cada tipo celular de dataset B
for (cell_type in B_painGenes_typeCells_expression) {
  B_DifferentialResults_pain[[cell_type]] <- FindMarkers(B_scores, ident.1 = cell_type, 
                                                         features = B_pain_genes_in_seurat)
}

# Extraer los 5 genes más significativos de cada tipo celular
top_genes_pain_B <- unique(unlist(lapply(B_DifferentialResults_pain, 
                                         function(x) head(rownames(x[order(x$p_val_adj),]), 5))))

# Filtrar y escalar los datos de los genes seleccionados en células relevantes para el dataset B
B_scores_Pain_Cells <- subset(B_scores, idents = B_painGenes_typeCells_expression)
B_scores_Pain_Cells <- ScaleData(B_scores_Pain_Cells, features = top_genes_pain_B)

# HEATMAPS para ver expresión de genes 
heatmapA_pain.Cells <- DoHeatmap(A_scores_Pain_Cells, features = top_genes_pain_A, 
                                 group.bar = TRUE, size = 3, group.colors = cell_colors) 
heatmapA_pain.Cells <- heatmapA_pain.Cells +
  guides(color = FALSE) +
  labs(fill = "Expresión",
       y = "Genes marcadores relacionados con el dolor",
       title = "Expresión en las célulares con firma de genes del dolor (Dataset A)") +
  scale_fill_viridis(option="magma") +
  theme(text = element_text(face = "bold"))

heatmapB_pain.Cells <- DoHeatmap(B_scores_Pain_Cells, features = top_genes_pain_B, 
                                 group.bar = TRUE, size = 3, group.colors = cell_colors) 
heatmapB_pain.Cells <- heatmapB_pain.Cells +
  guides(color = FALSE) +
  labs(fill = "Expresión",
       y = "Genes marcadores relacionados con el dolor",
       title = "Expresión en células con firma de genes del dolor (Dataset B)") +
  scale_fill_viridis(option="magma") +
  theme(text = element_text(face = "bold"))

# DOTPLOTS para visualización de expresión por tratamiento y condición
A_dotplot <- DotPlot(A_scores_Pain_Cells, features = top_genes_pain_A,
                     group.by = "cell_type", split.by = "condition_treatment",
                     cols = "Spectral") +
  labs(
    y = "Tipo celular_Condición_Tratamiento",
    x = "Genes asociados al dolor"
  ) + theme_article() +
  theme(
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 10, face = "bold"),
  ) + RotatedAxis() +
  guides(size = guide_legend(title = "% expresado (céls.)"),
         colour = guide_colourbar("Expresión media"))
print(A_dotplot)

B_dotplot <- DotPlot(B_scores_Pain_Cells, features = top_genes_pain_B,
                     group.by = "cell_type", split.by = "condition_treatment",
                     cols = "Spectral") +
  labs(
    y = "Tipo celular_Condición_Tratamiento",
    x = "Genes asociados al dolor"
  ) + theme_article() +
  theme(
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 10, face = "bold"),
  ) + RotatedAxis() +
  guides(size = guide_legend(title = "% expresado (céls.)"),
         colour = guide_colourbar("Expresión media"))
print(B_dotplot)