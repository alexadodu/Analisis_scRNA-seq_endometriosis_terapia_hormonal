# Análisis de datos de single cell RNA-seq para evaluar la eficacia de la terapia hormonal en el tratamiento del dolor en la endometriosis

-----------------
La endometriosis es una enfermedad crónica caracterizada por dolor pélvico, que afecta la calidad de vida de muchas mujeres. El tratamiento hormonal ha sido ampliamente utilizado para manejar los síntomas, aunque la comprensión de su efecto a nivel molecular sigue siendo limitada. Este trabajo analiza la eficacia de la terapia hormonal en el tratamiento del dolor en pacientes con endometriosis mediante el uso de datos de secuenciación de RNA de células individuales. Para evaluar los cambios en la expresión génica, se procesaron dos conjuntos de datos con Seurat, sumando un total de 60 muestras. Se identificaron firmas de dolor en células estromales, inmunes y fibroblastos. Los resultados muestran una reducción en la expresión de varios genes asociados al dolor, como *SPARC* y *MDK*, en pacientes bajo tratamiento hormonal. Las limitaciones computacionales de 16GB de RAM impidieron la integración de ambos conjuntos de datos y el preprocesamiento de lecturas crudas, lo que podría haber afectado los resultados. Se recomienda realizar estudios adicionales con enfoques integrados para obtener conclusiones más robustas.

**Palabras clave:** Endometriosis, scRNA-seq, terapia hormonal, dolor, expresión génica, SPARC, MDK, Seurat.

*Este trabajo ha sido realizado como Trabajo de Fin de Máster en Bioinformática en la Universidad Internacional de Valencia (VIU).*

--------------

## Contenido
- En la carpeta `scripts` se hallan todos los códigos empleados para el procesamiento de datos de scRNA-seq de pacientes con endometriosis para visualizar el impacto del tratamiento hormonal.
- El documento del TFM, que sirve como guía para visualizar los resultados obtenidos se hallan en `documento_tfm`.
- Todas las figuras empleadas se encuentran en `imagenes`.

-----
## Datos empleados
1. Fonseca, M. A. S., Haro, M., Wright, K. N., Lin, X., Abbasi, F., Sun, J., Hernandez, L., Orr, N. L., Hong, J., Choi-Kuaea, Y., Maluf, H. M., Balzer, B. L., Fishburn, A., Hickey, R., Cass, I., Goodridge, H. S., Truong, M., Wang, Y., Pisarska, M. D., … Lawrenson, K. (2023). Single-cell transcriptomic analysis of endometriosis. Nature Genetics, 55(2), 255. https://doi.org/10.1038/S41588-022-01254-1
2. Tan, Y., Flynn, W. F., Sivajothi, S., Luo, D., Bozal, S. B., Davé, M., Luciano, A. A., Robson, P., Luciano, D. E., & Courtois, E. T. (2022). Single cell analysis of endometriosis reveals a coordinated transcriptional program driving immunotolerance and angiogenesis across eutopic and ectopic tissues. Nature Cell Biology, 24(8), 1306. https://doi.org/10.1038/S41556-022-00961-5

----
----

# Analysis of Single-Cell RNA-seq Data to Evaluate the Efficacy of Hormonal Therapy in the Treatment of Endometriosis-Related Pain

-----------------
Endometriosis is a chronic disease characterized by pelvic pain, severely impacting the quality of life of many women. Hormonal therapy has been widely used to manage symptoms, although its molecular effects remain poorly understood. This work analyzes the efficacy of hormonal therapy in treating pain in endometriosis patients using single-cell RNA sequencing (scRNA-seq) data. Two datasets, comprising a total of 60 samples, were processed using Seurat to evaluate changes in gene expression. Pain-related signatures were identified in stromal, immune, and fibroblast cells. Results show a reduction in the expression of several pain-associated genes, such as *SPARC* and *MDK*, in patients under hormonal treatment. Computational limitations, including 16GB of RAM, prevented the integration of both datasets and the preprocessing of raw reads, which may have impacted the results. Additional studies with integrated approaches are recommended to obtain more robust conclusions.

**Keywords:** Endometriosis, scRNA-seq, hormonal therapy, pain, gene expression, SPARC, MDK, Seurat.

*This work was carried out as part of a Master's Thesis in Bioinformatics at the International University of Valencia (VIU).*

--------------

## Content
- The `scripts` folder contains all the code used to process scRNA-seq data from endometriosis patients to visualize the impact of hormonal treatment.
- The thesis document, which serves as a guide to understanding the obtained results, can be found in the `documento_tfm` folder.
- All figures used are located in the `imagenes` folder.

-----
## Data from
1. Fonseca, M. A. S., Haro, M., Wright, K. N., Lin, X., Abbasi, F., Sun, J., Hernandez, L., Orr, N. L., Hong, J., Choi-Kuaea, Y., Maluf, H. M., Balzer, B. L., Fishburn, A., Hickey, R., Cass, I., Goodridge, H. S., Truong, M., Wang, Y., Pisarska, M. D., … Lawrenson, K. (2023). Single-cell transcriptomic analysis of endometriosis. Nature Genetics, 55(2), 255. https://doi.org/10.1038/S41588-022-01254-1
2. Tan, Y., Flynn, W. F., Sivajothi, S., Luo, D., Bozal, S. B., Davé, M., Luciano, A. A., Robson, P., Luciano, D. E., & Courtois, E. T. (2022). Single cell analysis of endometriosis reveals a coordinated transcriptional program driving immunotolerance and angiogenesis across eutopic and ectopic tissues. Nature Cell Biology, 24(8), 1306. https://doi.org/10.1038/S41556-022-00961-5