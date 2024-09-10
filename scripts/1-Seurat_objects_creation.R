### Author: Madalina Alexandra Dodu

# Script para preparación de los datos de scRNA-seq con Seurat.  Este script detalla 
# el proceso de creación y manipulación de objetos Seurat a partir de matrices de conteo 
# de scRNA-seq para el análisis de células individuales. El código carga los datos en 
# formato H5, crea objetos Seurat, calcula porcentajes de genes mitocondriales y ribosomales, 
# y añade metadatos adicionales para cada muestra. Finalmente, los objetos Seurat procesados 
# se guardan para su posterior análisis, facilitando la visualización y el estudio de 
# la expresión génica específica.

## *** PAQUETES ***
install.packages("openxlsx")
install.packages('Seurat')
install.packages("hdf5r")
install.packages("tidyverse")

library(hdf5r)
library("openxlsx")
library('Seurat')
library(tidyverse)


## *** PREPARACIÓN DE LOS OBJETOS SEURAT ***

# Especifica el directorio donde están tus archivos HDF5.
files <- list.files(path = "data/raw/nombre_directorio", pattern = "*.h5", full.names = TRUE)

# Extrae los nombres de los archivos sin ruta ni extensión, para utilizarlos como nombres de muestra.
sample_names <- sapply(basename(files), function(x) sub("\\.h5$", "", x))

# Carga el archivo de metadatos desde un archivo Excel.
metadata <- read.xlsx("metadata/metadata.xlsx", sheet = 1)
ids <- metadata$id  # Extrae los IDs de las muestras desde la metadata.

# Crea una lista vacía para almacenar los objetos Seurat que se generarán.
seurat_list <- list()

# Bucle para cargar y procesar cada archivo H5 individualmente.
for (i in seq_along(files)) { 
  file <- files[i]
  sample_name <- sample_names[i]
  id <- ids[i]
  
  # Lee los datos del archivo HDF5 y crea un objeto Seurat a partir de ellos.
  seurat_data <- Read10X_h5(file)
  seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   min.features = 200, 
                                   project = sample_name)
  
  # Renombra las celdas del objeto Seurat para asegurar que cada celda tenga un nombre único.
  cell_ids <- colnames(seurat_obj)
  new_cell_ids <- paste(id, cell_ids, sep = "_")
  seurat_obj <- RenameCells(seurat_obj, new.names = new_cell_ids) 
  
  # Añade el objeto Seurat a la lista de objetos.
  seurat_list <- append(seurat_list, list(seurat_obj))
}

# Limpieza de memoria: Elimina las variables que ya no se necesitan para liberar espacio.
rm(seurat_data, seurat_obj, file, files, i, sample_name, sample_names, cell_ids, id, ids, new_cell_ids)

# Bucle para calcular porcentajes de genes mitocondriales y ribosomales, y añadir metadata a cada objeto Seurat.
for (i in seq_along(seurat_list)) {
  
  # Extrae el objeto Seurat actual de la lista.
  current_sample <- seurat_list[[i]]
  
  # Calcula el porcentaje de genes mitocondriales en el objeto Seurat.
  current_sample$percent.mt <- PercentageFeatureSet(object = current_sample, pattern = "^MT-")
  
  # Calcula el porcentaje de genes ribosomales en el objeto Seurat.
  current_sample$percent.Ribosomal <- PercentageFeatureSet(object = current_sample, pattern = "^RP[LS]")
  
  # Obtiene los nombres de las células del objeto Seurat.
  cell_ids <- colnames(current_sample)
  
  # Crea etiquetas de muestra basadas en los nombres de las células.
  sample_labels <- sapply(strsplit(cell_ids, "_"), `[`, 1)

  # Crea un dataframe de metadatos con los IDs de las células.
  metadata_seurat <- data.frame(row.names = cell_ids, id = sample_labels)

  # Guarda los nombres originales de las filas para usarlos más adelante.
  original_row_names <- row.names(metadata_seurat)
  
  # Combina los metadatos originales con los nuevos metadatos basados en los IDs.
  merged_metadata <- merge(metadata_seurat, metadata, by = "id", all.x = TRUE)
  
  # Restaura los nombres de las filas originales en los metadatos combinados.
  row.names(merged_metadata) <- original_row_names
  
  # Añade los metadatos combinados al objeto Seurat.
  current_sample <- AddMetaData(current_sample, metadata = merged_metadata)
  
  # Renombra las columnas de 'nCount_RNA' a 'nUMI' y 'nFeature_RNA' a 'nGene' para mayor claridad.
  colnames(current_sample@meta.data)[colnames(current_sample@meta.data) == "nCount_RNA"] <- "nUMI"
  colnames(current_sample@meta.data)[colnames(current_sample@meta.data) == "nFeature_RNA"] <- "nGene"
  
  # Guarda el objeto Seurat modificado de nuevo en la lista.
  seurat_list[[i]] <- current_sample
}

# Visualiza los metadatos del primer objeto Seurat para verificación.
View(seurat_list[[1]]@meta.data)

# Limpieza de memoria: Elimina variables temporales para liberar espacio.
rm(i, metadata, merged_metadata, metadata_seurat, current_sample, cell_ids, original_row_names, sample_labels)

# Guarda la lista de objetos Seurat en un archivo RDS para su uso posterior.
saveRDS(seurat_list, file = "data/processed/seurat_list_raw.rds")