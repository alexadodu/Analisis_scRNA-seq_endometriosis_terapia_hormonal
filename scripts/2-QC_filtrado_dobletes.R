### Author: Madalina Alexandra Dodu

# Script para filtrado de calidad y eliminación de dobletes en datos scRNA-seq. 
# Este script describe el proceso de filtrado y pre-procesamiento de datos scRNA-seq 
# mediante Seurat, con el objetivo de asegurar la calidad de los datos y eliminar 
# células potencialmente problemáticas como los dobletes. Primero, se realiza un 
# filtrado inicial basado en la cantidad de genes y la proporción de ARN mitocondrial 
# para mantener únicamente células de alta calidad. Luego, se implementa un análisis 
# más detallado para identificar y remover dobletes, utilizando la herramienta 
# DoubletFinder tras realizar una reducción de dimensionalidad mediante PCA y UMAP. 
# Este procedimiento garantiza la limpieza de los datos antes de su análisis posterior.

## *** PAQUETES ***
install.packages('Seurat')
install.packages("tidyverse")
install.packages("devtools")

library(devtools)
install_github('chris-mcginnis-ucsf/DoubletFinder')

library(DoubletFinder)
library('Seurat')
library(tidyverse)


## *** CONTROL DE CALIDAD ***

## FILTRADO DE CALIDAD
# --- DATASET A ---
# Se filtra cada objeto Seurat en la lista, eliminando células con menos de 500 genes detectados,
# menos de 1000 UMI, más de 100000 UMI, o con un porcentaje de genes mitocondriales mayor al 25%.
seurat_list <- lapply(seurat_list, function(x) {
  subset(x, subset = nGene >= 500 & nUMI >= 1000 & nUMI <= 100000 & percent.mt < 25)
})

# --- DATASET B ---
# Se filtra cada objeto Seurat en la lista, eliminando células con menos de 200 genes detectados y un porcentaje 
# mayor al 20%.
seurat_list <- lapply(seurat_list, function(x) {
  subset(x, subset = nGene >= 200 & percent.mt < 20)
})

## FILTRADO DE DOBLETES 
# Se inicia un bucle para iterar a través de cada objeto Seurat en la lista con el fin de identificar y eliminar 
# dobletes.
for (i in 1:length(seurat_list)) {
  # Imprime el índice de la muestra actual para seguimiento.
  print(paste0("Sample ", i))
  
  # Pre-procesamiento del objeto Seurat siguiendo el flujo de trabajo estándar de Seurat.
  sample <- NormalizeData(seurat_list[[i]])  # Normaliza los datos de expresión génica.
  sample <- FindVariableFeatures(sample)     # Identifica los genes con mayor variabilidad.
  sample <- ScaleData(sample)                # Escala los datos para que los valores sean comparables.
  sample <- RunPCA(sample)                   # Realiza un Análisis de Componentes Principales (PCA).
  
  # Construcción del gráfico de vecinos y agrupamiento.
  sample <- FindNeighbors(object = sample, dims = 1:20)  # Encuentra vecinos similares basado en las dimensiones 
especificadas.
  sample <- FindClusters(object = sample)                # Agrupa las células en clústeres basados en similitudes.
  sample <- RunUMAP(sample, dims = 1:20)                 # Ejecuta UMAP para visualizar los clústeres en un espacio de menor dimensión.

  # Identificación de pK óptimo para DoubletFinder (sin verdad de terreno)
  sweep.list <- paramSweep(sample, PCs = 1:20, sct = FALSE)   # Realiza una búsqueda de parámetros para DoubletFinder.
  sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)       # Resumen de las estadísticas de la búsqueda de parámetros.
  bcmvn <- find.pK(sweep.stats)                               # Identifica el valor óptimo de pK para maximizar la métrica BC.
  
  # Se selecciona el pK óptimo basado en la distribución del coeficiente de bimodalidad.
  pK <- bcmvn %>%
    filter(BCmetric == max(BCmetric)) %>%
    select(pK)
  pK <- as.numeric(as.character(pK[[1]]))
  
  ## Estimación de la proporción de dobletes homotípicos
  annotations <- sample@meta.data$seurat_clusters  # Extrae las anotaciones de clústeres.
  homotypic.prop <- modelHomotypic(annotations)    # Modela la proporción de dobletes homotípicos.
  nExp.poi <- round(0.075 * nrow(sample@meta.data)) ## Supone una tasa de formación de dobletes del 7.5% - ajustar según el dataset.
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop)) # Ajusta la estimación de dobletes esperados.

  # Corre DoubletFinder para identificar dobletes
  sample <- doubletFinder(seu = sample, 
                          PCs = 1:20, 
                          pN = 0.25,
                          pK = pK,
                          nExp = nExp.poi.adj, reuse.pANN = FALSE, sct = FALSE)
  
  # Actualiza el metadata con los resultados de DoubletFinder
  metadata <- sample@meta.data
  colnames(metadata)[44] <- "doublet_finder"  # Renombra la columna para identificar dobletes.
  sample@meta.data <- metadata 
  
  # Filtra las células etiquetadas como singletes.
  last_col_name <- colnames(metadata)[ncol(metadata)]
  is_singlet <- FetchData(sample, vars = last_col_name) == "Singlet"
  
  # Subconjunto de células singletes y guardar el resultado en la lista original.
  singlets <- subset(sample, cells = which(is_singlet))
  seurat_list[[i]] <- singlets
  remove(singlets)
}

# Guarda la lista de objetos Seurat filtrados y sin dobletes
saveRDS(seurat_list, file = "data/processed/seurat_list_filtered_noDoublets.rds")

# Limpia el entorno eliminando las variables que ya no son necesarias
rm(bcmvn, is_singlet, metadata, sample, sweep.stats, sweep.list, annotations, 
  homotypic.prop, i, last_col_name, nExp.poi, nExp.poi.adj, pK, current_sample)