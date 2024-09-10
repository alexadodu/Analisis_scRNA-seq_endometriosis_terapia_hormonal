### Author: Madalina Alexandra Dodu

# Script para normalización y detección del ciclo celular, combinación de datos y clustering 
# en datos de scRNA-seq. Este anexo describe el proceso de normalización, anotación del 
# ciclo celular, combinación de datos de múltiples muestras, y el posterior análisis de 
# agrupación de datos de scRNA-seq.

## *** PAQUETES ***
install.packages('Seurat')
install.packages("tidyverse")

library('Seurat')
library(tidyverse)

## *** NORMALIZACIÓN Y DETERMINACIÓN DEL CICLO CELULAR ***
# Cargar la lista de objetos Seurat previamente filtrados y sin dobletes
seurat_list <- readRDS(file = "data/processed/seurat_list_filtered_noDoublets.rds")

# Cargar las listas de genes asociados a las fases del ciclo celular (S y G2M) (vienen con Seurat)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Bucle para normalizar los datos y determinar la fase del ciclo celular de cada célula
for (i in 1:length(seurat_list)) {
  # Normalización de los datos de expresión génica
  seurat_list[[i]] <- NormalizeData(seurat_list[[i]], verbose = TRUE)
  
  # Anotación del ciclo celular basado en las fases S y G2M
  seurat_list[[i]] <- CellCycleScoring(seurat_list[[i]], g2m.features = g2m.genes, s.features = s.genes)
}

## *** COMBINACIÓN DE LOS DATOS DE LAS MUESTRAS ***
# Combinar todos los objetos Seurat de la lista en un único objeto para análisis conjunto
combined_dataset <- merge(seurat_list[[1]], y = seurat_list[2:length(seurat_list)], 
                          project = "Endometriosis_Dataset")

# Limpiar el metadata eliminando columnas innecesarias creadas durante la eliminación de dobletes
columns_to_remove <- c(grep("^pANN", colnames(combined_dataset@meta.data), value = TRUE), 
                      grep("^DF", colnames(combined_dataset@meta.data), value = TRUE))

# Eliminar las columnas de la metadata
combined_dataset@meta.data <- combined_dataset@meta.data[, !(colnames(combined_dataset@meta.data) 
                                                        %in% columns_to_remove)]

# Guardar el objeto combinado y limpio
saveRDS(combined_dataset, file = "data/processed/seurat_merged_postQC_cc.rds")

## *** CLUSTERING ***
# Identificación de las características genéticas más variables para su uso en clustering
combined_dataset <- FindVariableFeatures(combined_dataset, selection.method = "vst", nfeatures = 2000)

# Escalar los datos para que los valores sean comparables entre sí
combined_dataset <- ScaleData(combined_dataset, features = rownames(combined_dataset)) # Dataset A
combined_dataset <- ScaleData(combined_dataset) # Dataset B

# Realizar un Análisis de Componentes Principales (PCA) para reducir la dimensionalidad del dataset
combined_dataset <- RunPCA(combined_dataset, features = VariableFeatures(combined_dataset))

# Guardar el objeto Seurat después de PCA para posibles análisis posteriores
saveRDS(combined_dataset, file = "data/processed/seurat_merged_postPCA.rds")

# Volver a cargar el objeto Seurat post-PCA si es necesario
combined_dataset <- readRDS("data/processed/seurat_merged_postPCA.rds")

# Realizar el clustering en base a los resultados del PCA
ElbowPlot(combined_dataset) # Se utiliza el codo para determinar el número óptimo de PCs para el clustering

# Encontrar vecinos cercanos y realizar clustering de las células en base a 10 dimensiones principales
combined_dataset <- FindNeighbors(combined_dataset, dims = 1:10, verbose = FALSE)
combined_dataset <- FindClusters(combined_dataset, 
                                 resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))

# Ejecutar UMAP para visualizar los clústeres en un espacio de menor dimensión
combined_dataset <- RunUMAP(combined_dataset, dims = 1:10, reduction = "pca", verbose = FALSE)

# Guardar el objeto Seurat final con los clústeres identificados
saveRDS(combined_dataset, file = "data/processed/seurat_merged_clustered.rds")