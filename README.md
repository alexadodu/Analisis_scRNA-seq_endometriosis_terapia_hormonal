# Análisis de datos de single cell RNA-seq para evaluar la eficacia de la terapia hormonal en el tratamiento del dolor en la endometriosis

-----------------
## Resumen
La endometriosis es un síndrome inflamatorio dependiente de estrógenos que afecta al 10% de las mujeres a nivel mundial. Se caracteriza por la presencia de tejido similar al endometrial fuera de su ubicación normal, denominado endometrio ectópico o lesiones endometriales, que puede causar inflamación y el consecuente dolor crónico y cíclico (dismenorrea) en las pacientes. El tratamiento principal de este síntoma a largo plazo es la terapia hormonal. Si bien los mecanismos que generan el dolor no están completamente comprendidos, se ha identificado una serie de genes involucrados en los procesos de inflamación y nocicepción. En este trabajo, se analizó el impacto del tratamiento hormonal en la expresión de genes asociados al dolor en pacientes con endometriosis mediante el análisis de datos de secuenciación de ARN de célula única o single cell (scRNA-seq). El objetivo principal del estudio fue analizar el impacto de la terapia hormonal en la expresión de estos genes en diferentes tipos celulares del endometrio eutópico y ectópico. Para ello, se procesaron y analizaron con Seurat 60 muestras de dos conjuntos de datos públicos de scRNA-seq, uno proveniente de pacientes con tratamiento hormonal y otro compuesto mayoritariamente por muestras sin tratamiento. 

Los resultados muestran una mayor expresión de genes asociados al dolor en células estromales, monocitos, macrófagos y fibroblastos estromales, tipos celulares que están implicados en el mecanismo del dolor. El tratamiento hormonal tuvo un impacto significativo en la reducción de la expresión de *SPARC*, *MDK*, *PTN* y *PRRX1* en células inmunes y estromales. Además, la reducción en la expresión del receptor de estrógenos *ESR1* en las células estromales sugiere una disminución de la inflamación mediada por estrógenos, lo que está relacionado con una menor percepción del dolor. También se observó una reducción en la expresión de la mayoría de los genes asociados al dolor en las células estromales durante la fase menstrual, lo que apoya la eficacia del tratamiento hormonal en la disminución de la dismenorrea. 

En conclusión, el análisis revela que el tratamiento hormonal modula la expresión de genes relacionados con la percepción del dolor en células estromales, monocitos y macrófagos, lo que podría explicar en parte la mejora de los síntomas en pacientes bajo esta terapia.


**Palabras clave:** endometriosis, scRNA-seq, terapia hormonal, dolor, dismenorrea, expresión génica, *SPARC*, *MDK*, *PNT*, *PRRX1*, *ESR1*, Seurat.

*Este trabajo ha sido realizado como Trabajo de Fin de Máster en Bioinformática en la Universidad Internacional de Valencia (VIU).*

## Contenido
- En la carpeta `scripts` se hallan todos los códigos empleados para el procesamiento de datos de scRNA-seq de pacientes con endometriosis para visualizar el impacto del tratamiento hormonal.
- El documento del TFM, que sirve como guía para visualizar los resultados obtenidos se hallan en `documento_tfm`.
- Todas las figuras empleadas se encuentran en `imagenes`. Las referencias correspondientes que se han empleado en la elbaoración de las figuras se hallan en el trabajo final.

## Datos empleados
### scRNA-seq
- Fonseca, M. A. S., Haro, M., Wright, K. N., Lin, X., Abbasi, F., Sun, J., Hernandez, L., Orr, N. L., Hong, J., Choi-Kuaea, Y., Maluf, H. M., Balzer, B. L., Fishburn, A., Hickey, R., Cass, I., Goodridge, H. S., Truong, M., Wang, Y., Pisarska, M. D., … Lawrenson, K. (2023). Single-cell transcriptomic analysis of endometriosis. *Nature Genetics*, 55(2), 255. https://doi.org/10.1038/S41588-022-01254-1
    GEO: [GSE179640](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179640)
- Tan, Y., Flynn, W. F., Sivajothi, S., Luo, D., Bozal, S. B., Davé, M., Luciano, A. A., Robson, P., Luciano, D. E., & Courtois, E. T. (2022). Single cell analysis of endometriosis reveals a coordinated transcriptional program driving immunotolerance and angiogenesis across eutopic and ectopic tissues. *Nature Cell Biology*, 24(8), 1306. https://doi.org/10.1038/S41556-022-00961-5
    GEO: [GSE213216](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213216)

### Genes asociados al dolor
- Rahmioglu, N., Mortlock, S., Ghiasi, M., Møller, P. L., Stefansdottir, L., Galarneau, G., Turman, C., Danning, R., Law, M. H., Sapkota, Y., Christofidou, P., Skarp, S., Giri, A., Banasik, K., Krassowski, M., Lepamets, M., Marciniak, B., Nõukas, M., Perro, D., … Zondervan, K. T. (2023). The genetic basis of endometriosis and comorbidity with other pain and inflammatory conditions. *Nature Genetics*, 55(3), 423–436. https://doi.org/10.1038/S41588-023-01323-Z
- Wistrom, E., Chase, R., Smith, P. R., & Campbell, Z. T. (2022). A compendium of validated pain genes. *Wires Mechanisms of Disease*, 14(6), e1570. https://doi.org/10.1002/WSBM.1570

----
----

# Analysis of Single-Cell RNA-seq Data to Evaluate the Efficacy of Hormonal Therapy in the Treatment of Endometriosis-Related Pain

-----------------
## Summary
Endometriosis is an estrogen-dependent inflammatory syndrome that affects 10% of women worldwide. It is characterized by the presence of tissue similar to the endometrium outside its normal location, known as ectopic endometrium or endometrial lesions, which can cause inflammation and result in chronic and cyclical pain (dysmenorrhea) in patients. The primary long-term treatment for this symptom is hormonal therapy. Although the mechanisms that generate pain are not fully understood, a series of genes involved in inflammation and nociception processes have been identified. In this study, the impact of hormonal treatment on the expression of pain-associated genes in endometriosis patients was analyzed using single-cell RNA sequencing (scRNA-seq) data. The main objective of the study was to examine the impact of hormonal therapy on the expression of these genes across different cell types in both eutopic and ectopic endometrial tissues. To achieve this, 60 samples from two public scRNA-seq datasets were processed and analyzed with Seurat, one derived from patients under hormonal treatment and the other mainly composed of samples without treatment.

The results show higher expression of pain-associated genes in stromal cells, monocytes, macrophages, and stromal fibroblasts, which are cell types implicated in pain mechanisms. Hormonal treatment significantly reduced the expression of *SPARC*, *MDK*, *PTN*, and *PRRX1* in immune and stromal cells. Additionally, the reduction in the expression of the estrogen receptor *ESR1* in stromal cells suggests a decrease in estrogen-mediated inflammation, which is associated with a lower perception of pain. A reduction in the expression of most pain-associated genes in stromal cells during the menstrual phase was also observed, supporting the effectiveness of hormonal treatment in reducing dysmenorrhea.

In conclusion, the analysis reveals that hormonal treatment modulates the expression of genes related to pain perception in stromal cells, monocytes, and macrophages, which may partly explain the improvement in symptoms for patients undergoing this therapy.

**Keywords:** endometriosis, scRNA-seq, hormonal therapy, pain, dysmenorrhea, gene expression, *SPARC*, *MDK*, *PNT*, *PRRX1*, *ESR1*, Seurat.

*This work was conducted as a Master’s Thesis in Bioinformatics at the International University of Valencia (VIU).*

## Content
- The `scripts` folder contains all the codes used for processing scRNA-seq data from endometriosis patients to visualize the impact of hormonal treatment.
- The Master’s Thesis document, which serves as a guide to visualize the obtained results, is located in `documento_tfm`.
- All the figures used are located in `imagenes`. The corresponding references employed in the creation of the figures are found in the final thesis.

## Data used
### scRNA-seq
- Fonseca, M. A. S., Haro, M., Wright, K. N., Lin, X., Abbasi, F., Sun, J., Hernandez, L., Orr, N. L., Hong, J., Choi-Kuaea, Y., Maluf, H. M., Balzer, B. L., Fishburn, A., Hickey, R., Cass, I., Goodridge, H. S., Truong, M., Wang, Y., Pisarska, M. D., … Lawrenson, K. (2023). Single-cell transcriptomic analysis of endometriosis. *Nature Genetics*, 55(2), 255. https://doi.org/10.1038/S41588-022-01254-1  
    GEO: [GSE179640](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179640)
- Tan, Y., Flynn, W. F., Sivajothi, S., Luo, D., Bozal, S. B., Davé, M., Luciano, A. A., Robson, P., Luciano, D. E., & Courtois, E. T. (2022). Single cell analysis of endometriosis reveals a coordinated transcriptional program driving immunotolerance and angiogenesis across eutopic and ectopic tissues. *Nature Cell Biology*, 24(8), 1306. https://doi.org/10.1038/S41556-022-00961-5  
    GEO: [GSE213216](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213216)

### Pain-associated genes
- Rahmioglu, N., Mortlock, S., Ghiasi, M., Møller, P. L., Stefansdottir, L., Galarneau, G., Turman, C., Danning, R., Law, M. H., Sapkota, Y., Christofidou, P., Skarp, S., Giri, A., Banasik, K., Krassowski, M., Lepamets, M., Marciniak, B., Nõukas, M., Perro, D., … Zondervan, K. T. (2023). The genetic basis of endometriosis and comorbidity with other pain and inflammatory conditions. *Nature Genetics*, 55(3), 423–436. https://doi.org/10.1038/S41588-023-01323-Z
- Wistrom, E., Chase, R., Smith, P. R., & Campbell, Z. T. (2022). A compendium of validated pain genes. *Wires Mechanisms of Disease*, 14(6), e1570. https://doi.org/10.1002/WSBM.1570