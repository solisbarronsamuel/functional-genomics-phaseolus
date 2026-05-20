# Cambios en el bacterioma y fungoma de la rizosfera asociados a bacterias fijadoras de nitrógeno en *Phaseolus vulgaris*

Samuel Solís Barrón

## 1 Objetivo
Comparar el bacterioma y fungoma de la rizosfera asociados a bacterias fijadoras de nitrógeno en *Phaseolus vulgaris*.

## 2 Metodología

Este proyecto es una réplcia de resultados de Influence of organic plant breeding on the rhizosphere microbiome of common bean (*Phaseolus vulgaris L.*). (Park et al., 2023)

### 2.1 Obtención de los datos

Los datos para el análisis se obtuvieron de NCBI BioProject: RJNA988238 (Bacteria) y  PPRJNA989655 (Hongo).

### 2.2 Procesamiento de variantes de secuencia de amplicón (ASV)

El procesamiento del filtrado de calidad de los ASV (Q ≥ 33), procesamiento de lecturas paired-end, eliminación de las secuencias de primers y quimeras ocasionadas por el PCR utilizará DADA2 (Sahil et al., 2025). Las lecturas tendrán una filtración de 260 pb en dirección forward y 240 pb en dirección reverse, utilizando los parámetros maxN = 0, maxEE = c(2,2), truncQ = 2 y rm.phix = TRUE (Mataranyika et al., 2024). Los ASV conservados serán ≥ 270 pb (Sahil et al., 2025). La asignación taxonómica se realizará mediante alineamiento contra la base de datos SILVA v131.1 y UNITE_v2025, los niveles taxonómicos van desde el filo hasta el género conservando ≥ 99% de identidad (Overgaard et al., 2022). La rarefacción se realizó para estandarizar la comparabilidad para los cálculos de diversidad (Mataranyika et al., 2024).

### 2.3 Análisis estadístico

Los análisis de diversidad alfa y beta se realizarán en R utilizando el paquete phyloseq. Las diferencias entre grupos en la diversidad alfa utilizará las pruebas de Kruskal–Wallis o de rangos de Wilcoxon, mientras que los valores de p serán ajustados por comparaciones múltiples mediante el método Benjamini–Hochberg FDR (Yurgel et al., 2026). Asimismo, las diferencias entre grupos se analizarán mediante PERMANOVA (función adonis) con un máximo de 999 permutaciones (Sahil et al., 2025). Se elaborará un mapa de calor (heatmap) a partir de los coeficientes de correlación calculados utilizando el paquete pheatmap, empleando la correlación de Spearman, con el propósito de estimar la relación del bacterioma respecto a la muestra. La gráfica de la composición taxonómica se realizará con el paquete ggplot2 (Ghotbi et al., 2025).

## 3 Referencias

1. Ghotbi, M., Ghotbi, M., Kuzyakov, Y., & Horwath, W. R. (2025). Management and rhizosphere microbial associations modulate genetic-driven nitrogen fate. Agriculture, Ecosystems & Environment, 378, 109308. https://doi.org/10.1016/j.agee.2024.109308
2. Mataranyika, P. N., Bez, C., Venturi, V., Chimwamurombe, P. M., & Uzabakiriho, J. D. (2024). Rhizospheric, seed, and root endophytic-associated bacteria of drought-tolerant legumes grown in arid soils of Namibia. Heliyon. https://doi.org/10.1016/j.heliyon.2024.e36718
3. Overgaard, C. K., Tao, K., Zhang, S., Christensen, B. T., Blahovska, Z., Radutoiu, S., Kelly, S., & Dueholm, M. K. D. (2022). Application of ecosystem-specific reference databases for increased taxonomic resolution in soil microbial profiling. Frontiers in Microbiology, 13, 942396. https://doi.org/10.3389/fmicb.2022.942396
4. Park, H. E., Nebert, L., King, R. M., Busby, P., & Myers, J. R. (2023). Influence of organic plant breeding on the rhizosphere microbiome of common bean (Phaseolus vulgaris L.). Frontiers in Plant Science, 14. https://doi.org/10.3389/fpls.2023.1251919
5. Sahil, R., Pal, V., Kharat, A. S., & Jain, M. (2025). A multi-omics meta-analysis of rhizosphere microbiome reveals growth-promoting marker bacteria at different stages of legume development. Plant, Cell & Environment. https://doi.org/10.1111/pce.15429
6. Yurgel, S. N., Miklas, P. N., & Porter, L. D. (2026). Nitrogen fixation, crop production and bacterial communities of common bean cultivars: A 77-year breeding perspective. Plant and Soil, 520, 1-17. https://doi.org/10.1007/s11104-025-08257-x