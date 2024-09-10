### Author: Madalina Alexandra Dodu

## *** PAQUETES ***
install.packages("openxlsx")
install.packages('Seurat')
install.packages("tidyverse")

install.packages("devtools")
devtools::install_github("thomasp85/patchwork")

library("openxlsx")
library('Seurat')
library(tidyverse)
library(RColorBrewer)
library("patchwork")

## ELBOW PLOTS
elbowA <- ElbowPlot(A_joined_raw) + theme_article() +
    ggtitle("Dataset A") 
    # annotate("rect", xmin = 6, xmax = 12.5, ymin = 3 , ymax = 5.5, alpha = 0, color="pink") # añade recuadro

elbowB <- ElbowPlot(B_joined_raw) + theme_article() +
  ggtitle("Dataset B") 
  # annotate("rect", xmin = 5.5, xmax = 12.5, ymin = 2.2 , ymax = 4.2, alpha = 0, color="pink") # añade recuadro

elbow_plots <- (elbowA | elbowB) + 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 14, face = "bold"))

# Guardar los plots en alta resolución
ggsave(filename = "results/elbow_plots.png", 
       plot = elbow_plots, 
       width = 40,          # Ancho en centímetros
       height = 16,         # Alto en centímetros
       units = "cm",        # Unidades en centímetros
       dpi = 300) 


## VISUALIZACIÓN UMAP
A_combined_dataset <- readRDS(file = "data/processed/A_seurat_merged_clustered.rds")
Idents(object = A_joined) <- "RNA_snn_res.0.4"

B_combined_dataset <- readRDS(file = "data/processed/B_seurat_merged_clustered.rds")
Idents(object = B_joined) <- "RNA_snn_res.0.4"

# Ordena los clústeres para que aparezcan luego ordenados en la leyenda
cluster_numeric <- as.numeric(as.character(Idents(B_joined)))
B_joined$RNA_snn_res.0.4 <- factor(as.character(cluster_numeric),
                                             levels = sort(unique(cluster_numeric)))

A_umap <- DimPlot(A_joined,
                  reduction = "umap",
                  label = TRUE,
                  label.size = 5,
                  raster=FALSE) + 
  labs(
    x = "UMAP_1",              # Nombre del eje X
    y = "UMAP_2",              # Nombre del eje Y
  ) +
  ggtitle("Clústeres en el dataset A") +
  theme_article() +
  theme(
    legend.text = element_text(size = 11)  # Tamaño del texto de la leyenda
    ) 
print(A_umap)

B_umap <- DimPlot(B_joined,
                  reduction = "umap",
                  label = TRUE,
                  label.size = 5,
                  raster=FALSE) + 
  labs(
    x = "UMAP_1",              # Nombre del eje X
    y = "UMAP_2",              # Nombre del eje Y
    ) +
  ggtitle("Clústeres en el dataset B") +
  theme_article() +
  theme(
    legend.text = element_text(size = 11),  # Tamaño del texto de la leyenda
    ) 
print(B_umap)


## UMAPS FASE DEL CICLO CELULAR
A_plot_phases <- DimPlot(A_joined,
                         reduction = "umap", group.by = "Phase") + 
  labs(
    x = "UMAP_1",              # Nombre del eje X
    y = "UMAP_2",               # Nombre del eje Y
  ) +
  ggtitle("Fases del ciclo celular (dataset A)") +
  theme_article() + 
  theme(
    legend.text = element_text(size = 11),  # Tamaño del texto de la leyenda
  )
print(A_plot_phases)

B_plot_phases <- DimPlot(B_joined,
                         reduction = "umap", group.by = "Phase",
                         raster=FALSE) +
  labs(
    x = "UMAP_1",              # Nombre del eje X
    y = "UMAP_2",              # Nombre del eje Y
  ) +
  ggtitle("Fases del ciclo celular (Dataset B)") +
  theme_article() + 
  theme(
    legend.text = element_text(size = 11),  # Tamaño del texto de la leyenda
    )
print(B_plot_phases)


## UMAPS CON ANOTACIÓN
# Variables con colores para tener los tipos celulares/condiciones siempre del mismo color
cell_colors <- c(
  # endoteliales
  "Cél. endoteliales" = "#FF3333",
  
  #estromales
  "Cél. estromales endom. (sec.)" = "#F871FF",
  "Cél. estromales endom." = "#B83779FF",
  "Cél. estromales endom. (prolif.)" = "#EE2B6D",
  "Cél. glandulares" = "#FAB0B0",
  "Cél. estrom. decid." = "#822581FF",
  "Cél. estromales endom. (menstr.)" = "#FF65AC",
  "Cél. estromales" = "#B83779FF",
  
  #inmunes
  "Cél. T" = "#1AD1D3FF",
  "Cél. γ/δ & NK" = "#00CCFF",
  "Monocitos" = "#B684FF",
  "Macrófagos" = "#8B93FF",
  "Cél. B" = "#2B748EFF",
  "Cél. dendríticas" = "#4770E8FF",
  "Cél. plasmáticas" = "#42A0FF",
  "Cél. Inmunes" = "#00B5ED",
  "Cél. inmunes" = "#00B5ED",
  
  #epiteliales
  "Cél. no ciliadas" = "#FAF400",
  "Cél. epiteliales" = "#FFD521",
  "Cél. no ciliadas SPDEF+" = "#FEA832FF",
  "Cél. ciliadas" = "#FB8222FF",
  "Cél. epitel. prolif." = "#FEAB30",
  
  #musculares
  "Cél. musc. lisas" = "#DF6767",
  
  #fibroblastos
  "Fibroblastos estrom." = "#85C441",
  "Fibroblastos" = "#1B8E00",
  "Fibroblastos estrom. C7+" = "#33CC33"
)

treatment_colors <- c(
  "Control_N" = "#26828EFF",
  "Control_Y" = "#B4DE2CFF",
  "Endometriosis_N" = "#721F81FF",
  "Endometriosis_Y" = "#D6456CFF"
)

A_joined <- readRDS("data/processed/A_joined_clustered.rds")
B_joined <- readRDS("data/processed/B_joined_clustered.rds")

Idents(object = A_joined) <- "cell_type"
Idents(object = B_joined) <- "cell_type"

A_umap_ann <- DimPlot(A_joined, reduction = "umap", cols = cell_colors) + 
  labs(
    x = "UMAP_1",              # Nombre del eje X
    y = "UMAP_2"               # Nombre del eje Y
  ) +
  ggtitle("Dataset A") +
  theme_article() +
  theme(
    legend.text = element_text(size = 11)  # Tamaño del texto de la leyenda
  )
A_umap_ann <- LabelClusters(A_umap_ann, id = "ident", fontface = "bold")
print(A_umap_ann)

B_umap_ann <- DimPlot(B_joined, reduction = "umap", cols = cell_colors, raster=FALSE) + 
  labs(
    x = "UMAP_1",              # Nombre del eje X
    y = "UMAP_2"               # Nombre del eje Y
  ) +
  ggtitle("Dataset B") +
  theme_article() +
  theme(
    legend.text = element_text(size = 11)  # Tamaño del texto de la leyenda
  )
B_umap_ann <- LabelClusters(B_umap_ann, id = "ident",  fontface = "bold")
print(B_umap_ann)