# ─────────────────────────────────────────────────────────────────────────────

# phyloseq

# ─────────────────────────────────────────────────────────────────────────────

# 1. CARGA DE LIBRERIAS

library(phyloseq)      # estructuras de datos para microbiomas.
library(vegan)         # métodos de ecología de comunidades multivariada.
library(ggplot2)       # visualización de gráficos.
library(RColorBrewer)  # paletas de colores para gráficos.
library(scales)        # funciones auxiliares para formatos de ejes.
library(knitr)         # genera tablas formateadas en HTML/PDF dentro de reportes.

# Se establece una semilla aleatoria para que sea reproducibles entre ejecuciones.
set.seed(123)

# ─────────────────────────────────────────────────────────────────────────────

# 2.  CONSTRUCCIÓN DE OBJETOS PHYLOSEQ


# phyloseq agrupa tres (o más) tablas de datos en un único objeto:
#   • otu_table   – matriz de abundancias (muestras × ASVs cuando
#                   taxa_are_rows = FALSE)
#   • tax_table   – cadenas de taxonomía por ASV
#   • sample_data – metadatos por muestra (similar a data.frame)

# 'metadata' debe estar cargado previamente.

# ──────────────────────────────────────────────────────────────────────────────
# LOAD OBJECTS FROM DADA2
# ──────────────────────────────────────────────────────────────────────────────

# Load the chimera-free seqtabs (rows are plain SRR accessions at this point)
seqtab.nochim16S <- readRDS("03_Results/rds/16S/seqtab.nochim16S.RDS")
seqtab.nochimITS <- readRDS("03_Results/rds/ITS/seqtab.nochimITS.RDS")

# Load taxonomy (DECIPHER output)
taxid16S <- readRDS("03_Results/rds/16S/taxid16S-decipher.RDS")
taxidITS  <- readRDS("03_Results/rds/ITS/taxidITS_decipher.RDS")

# Reload metadata (only needed to build Host_Run labels and sample_data)
metadata16S <- read.csv("01_RawData/csv/16S/SraRunTable16S.csv", stringsAsFactors = FALSE)
metadataITS <- read.csv("01_RawData/csv/ITS/SraRunTableITS.csv", stringsAsFactors = FALSE)
metadata16S$amplicon <- "16S"
metadataITS$amplicon <- "ITS"
metadata <- rbind(metadata16S, metadataITS)

# ──────────────────────────────────────────────────────────────────────────────
# RENAME SEQTAB ROWS: SRR → Host_Run
# seqtab row names come from FASTQ filenames (plain SRR IDs, e.g. "SRR25070549").
# sample_data row names must match, so both are set to Host_Run format.
# ──────────────────────────────────────────────────────────────────────────────

meta16S_sub <- metadata[metadata$amplicon == "16S", ]
idx16S      <- match(rownames(seqtab.nochim16S), meta16S_sub$Run)
rownames(seqtab.nochim16S) <- paste(meta16S_sub$Host[idx16S],
                                    meta16S_sub$Run[idx16S],
                                    sep = "_")

metaITS_sub <- metadata[metadata$amplicon == "ITS", ]
idxITS      <- match(rownames(seqtab.nochimITS), metaITS_sub$Run)
rownames(seqtab.nochimITS) <- paste(metaITS_sub$Host[idxITS],
                                    metaITS_sub$Run[idxITS],
                                    sep = "_")

# ──────────────────────────────────────────────────────────────────────────────
# BUILD SAMPLE DATA
# Row names must be identical to the seqtab row names set above.
# ──────────────────────────────────────────────────────────────────────────────

samdf16S <- metadata[metadata$amplicon == "16S", ]
samdfITS  <- metadata[metadata$amplicon == "ITS", ]

rownames(samdf16S) <- paste(samdf16S$Host, samdf16S$Run, sep = "_")
rownames(samdfITS) <- paste(samdfITS$Host, samdfITS$Run, sep = "_")

# ──────────────────────────────────────────────────────────────────────────────
# CONSTRUCT PHYLOSEQ OBJECTS
# ──────────────────────────────────────────────────────────────────────────────

ps16S <- phyloseq(
  otu_table(seqtab.nochim16S, taxa_are_rows = FALSE),
  tax_table(taxid16S),
  sample_data(samdf16S)
)
taxa_names(ps16S) <- paste0("ASV", seq_len(ntaxa(ps16S)))
cat("=== Bacterioma (16S) ===\n")
ps16S

ps_ITS <- phyloseq(
  otu_table(seqtab.nochimITS, taxa_are_rows = FALSE),
  tax_table(taxidITS),
  sample_data(samdfITS)
)
taxa_names(ps_ITS) <- paste0("ASV", seq_len(ntaxa(ps_ITS)))
cat("\n=== Fungoma (ITS) ===\n")
ps_ITS

# ─────────────────────────────────────────────────────────────────────────────

# 3.  RESUMEN DEL OBJETO

# ntaxa(): devuelve el número de taxa (ASVs/OTUs) en el objeto phyloseq.
cat("Número de ASVs:", ntaxa(ps16S), "\n")

# nsamples(): devuelve el número de muestras.
cat("Número de muestras:", nsamples(ps16S), "\n")

# rank_names(): devuelve los nombres de columna de los rangos taxonómicos
# almacenados en tax_table. collapse = " > " los formatea como jerarquía legible.
cat("Rangos taxonómicos:", paste(rank_names(ps16S), collapse = " > "), "\n")

# sample_variables(): devuelve los nombres de columna del sample_data.
cat("Variables de metadatos:", paste(sample_variables(ps16S), collapse = ", "), "\n")

# kable(): formatea un data.frame como tabla Markdown/HTML con título.
# data.frame(tax_table(ps16S)) convierte el objeto taxonomyTable a un
# data.frame estándar. [1:6, ] selecciona las primeras 6 filas (primeros 6 ASVs).
kable(data.frame(tax_table(ps16S))[1:6, ],
      caption = "Primeros 6 ASVs y su clasificación taxonómica")


# ─────────────────────────────────────────────────────────────────────────────

# 4.  PREPROCESAMIENTO Y FILTRADO

# Distribución de lecturas por muestra

# sample_sums() calcula el total de lecturas por muestra cuando taxa_are_rows = FALSE).

ss16S <- sample_sums(ps16S)
summary(ss16S)   # mín, Q1, mediana, Q3, máx

# Determinar cuántos tipos de amplicón distintos hay.

n_st   <- length(unique(sample_data(ps16S)$amplicon))

# brewer.pal() devuelve un vector de n códigos HEX de la paleta "Set1".
# max(..., 3) asegura al menos 3 colores (brewer.pal requiere n ≥ 3).

colores <- brewer.pal(max(n_st, 3), "Set1")

# Se construye un data.frame que mapea nombres de muestra a sus conteos totales.

df_reads <- data.frame(
  Sample     = names(ss16S),
  TotalReads = ss16S
)

# Gráfico de barras de lecturas totales por muestra, ordenadas de forma ascendente.
# reorder(Sample, TotalReads) ordena las muestras en el eje x por conteo.
# coord_flip() intercambia los ejes para que los nombres largos aparezcan
# horizontalmente en el eje y.
# scale_y_continuous(labels = comma) formatea las marcas con separador de miles.
# scale_fill_gradientn() aplica un gradiente continuo de color (similar a Viridis).

ggplot(df_reads, aes(x = reorder(Sample, TotalReads),
                     y = TotalReads,
                     fill = TotalReads)) +
  geom_bar(stat = "identity") +          
  coord_flip() +
  scale_y_continuous(labels = comma) +
  scale_fill_gradientn(
    colors = c("#440154", "#31688e", "#35b779", "#fde725")  
  ) +
  labs(title = "Total de Lecturas por Muestra (16S)",
       x = "Muestra",
       y = "Número de lecturas") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")


# Filtrado por prevalencia y abundancia

cat("ASVs originales:", ntaxa(ps16S), "\n")

# Elimina ASVs con menos de 10 lecturas totales
# taxa_sums() devuelve las lecturas totales por ASV en todas las muestras.
# prune_taxa() conserva solo los taxa que cumplen la condición lógica.
# El umbral de 10 es comúnmente usado para excluir ruido de secuenciación.

ps16S_f <- prune_taxa(taxa_sums(ps16S) > 10, ps16S)
cat("ASVs tras umbral de abundancia (> 10 lecturas):", ntaxa(ps16S_f), "\n")

# Eliminar ASVs presentes en menos del 5 % de las muestras.
# Umbral de prevalencia: un ASV debe aparecer en al menos el 5 % de las muestras.

prev_thr <- 0.05 * nsamples(ps16S_f)

# Se extrae la matriz como matriz R estándar (pierde la clase phyloseq).

asv <- as(otu_table(ps16S_f), "matrix")

# Si las muestras son filas (taxa_are_rows = FALSE), se transpone para que
# las filas sean taxa. Esto garantiza que rowSums cuente por taxa.

if (!taxa_are_rows(ps16S_f)) asv <- t(asv)   # ahora siempre taxa × muestras

# Se cuenta en cuántas muestras cada ASV tiene al menos una lectura.

taxa_prev <- rowSums(asv > 0)

# Se conservan solo los ASVs cuya prevalencia supera el umbral.

ps16S_f <- prune_taxa(taxa_prev >= prev_thr, ps16S_f)

cat("ASVs tras filtrado de prevalencia (≥ 5 % muestras):", ntaxa(ps16S_f), "\n")
cat("ASVs removidos en total:", ntaxa(ps16S) - ntaxa(ps16S_f), "\n")


# Rarefacción
# La rarefacción submuestrea cada muestra a la misma profundidad de
# secuenciación (min_depth) para que las estimaciones de diversidad no estén
# confundidas por diferencias en el tamaño de la biblioteca.

# min() de sample_sums da el menor número de lecturas entre todas las muestras.

min_depth <- min(sample_sums(ps16S_f))
cat("Profundidad mínima:", min_depth, "\n")

# rarefy_even_depth() submuestrea aleatoriamente las lecturas de cada muestra
# sin reemplazo (replace = FALSE) hasta exactamente sample.size lecturas.
# rngseed: semilla aleatoria para reproducibilidad.
# trimOTUs = TRUE: elimina los ASVs que llegan a cero tras la rarefacción.
# verbose = FALSE: suprime el mensaje informativo sobre muestras eliminadas.

ps16S_rare <- rarefy_even_depth(
  ps16S_f,
  sample.size = min_depth,
  rngseed     = 123,
  replace     = FALSE,
  trimOTUs    = TRUE,
  verbose     = FALSE
)

cat("Muestras tras rarefacción:", nsamples(ps16S_rare), "\n")
cat("ASVs tras rarefacción:",     ntaxa(ps16S_rare),    "\n")


# Abundancias relativas
# Se transforman los conteos crudos a proporciones para que las muestras con
# diferente tamaño de biblioteca sean comparables en análisis composicionales
# (p. ej., Bray-Curtis, gráficos de barras, PCoA).
# transform_sample_counts() aplica una función arbitraria f(x) a cada vector
# de muestra x. Aquí f(x) = x / sum(x) convierte conteos en proporciones.

ps16S_rel <- transform_sample_counts(ps16S_f, function(x) x / sum(x))

# Verificación: todas las sumas de muestra deben ser iguales a 1.0.

head(round(sample_sums(ps16S_rel), 4))


# ─────────────────────────────────────────────────────────────────────────────

# 5.  ANÁLISIS DE DIVERSIDAD

# Diversidad Alfa
# La diversidad alfa cuantifica la diversidad DENTRO de una sola muestra.
# Índices utilizados:
#   Observed   – conteo de ASVs presentes
#   Chao1      – estimador no paramétrico de riqueza total
#   Shannon    – balancea riqueza y equitatividad (H' = -Σ pᵢ ln pᵢ)
#   Simpson    – probabilidad de que dos lecturas elegidas al azar
#                pertenezcan a ASVs diferentes
#   InvSimpson – 1/Simpson; valor mayor = más diverso

# estimate_richness() calcula todos los índices solicitados sobre el objeto
# rarefaccionado para garantizar una comparación justa entre muestras.
alpha16S <- estimate_richness(
  ps16S_rare,
  measures = c("Observed", "Chao1", "Shannon", "Simpson", "InvSimpson")
)

# Se guardan los nombres de muestra como columna para fusionado/graficado posterior.
alpha16S$Sample <- rownames(alpha16S)

# Se muestran las primeras 8 filas como tabla formateada.
kable(head(alpha16S[, c("Sample","Observed","Chao1","Shannon","Simpson")], 8),
      digits  = 3,
      caption = "Índices de diversidad alfa — bacterioma 16S")

# plot_richness() es la función integrada de phyloseq para graficar diversidad alfa.
# x = "Run" coloca las muestras en el eje x agrupadas por ID de corrida.
# measures selecciona qué índices mostrar como facetas.
# geom_point() superpone puntos individuales sobre el gráfico de tiras por defecto.
plot_richness(ps16S_rare, x = "Run",
              measures = c("Shannon", "InvSimpson")) +
  geom_point(size = 3, color = "#e41a1c") +
  labs(title = "Diversidad Alfa — Bacterioma (16S)",
       x = NULL, y = "Índice de Diversidad") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# ==== 5b.  Diversidad Beta ====
# La diversidad beta mide las diferencias en COMPOSICIÓN microbiana ENTRE muestras.

# --- Disimilitud de Bray-Curtis ---
# phyloseq::distance() calcula la disimilitud por pares.
# method = "bray": Bray-Curtis (0 = idénticos, 1 = completamente diferentes).
# Se usa el objeto de abundancias relativas para que el tamaño de la biblioteca
# no sesgue las distancias.
dist_bray16S <- phyloseq::distance(ps16S_rel, method = "bray")
cat("Dimensiones Bray-Curtis:", attr(dist_bray16S, "Size"), "x",
    attr(dist_bray16S, "Size"), "\n")

# --- PCoA (Análisis de Coordenadas Principales) ---
# ordinate() envuelve varios métodos de ordenación. PCoA (= MDS) encuentra ejes
# que maximizan la varianza explicada en la matriz de distancias.
ord_bray16S <- ordinate(ps16S_rel, method = "PCoA", distance = dist_bray16S)

# Se extraen los autovalores para calcular el % de varianza explicada por cada eje.
eig16S <- ord_bray16S$values$Eigenvalues

# Cada autovalor se divide entre el total para obtener proporciones; × 100.
var16S <- round(eig16S / sum(eig16S) * 100, 1)
cat("Eje 1:", var16S[1], "% | Eje 2:", var16S[2], "%\n")

# plot_ordination() crea un ggplot con las puntuaciones de muestra en los ejes.
# paste0() inserta el % de varianza en las etiquetas de eje como contexto.
plot_ordination(ps16S_rel, ord_bray16S) +
  geom_point(size = 4, alpha = 0.9, color = "#377eb8") +
  labs(title = "PCoA — Bray-Curtis (16S)",
       x = paste0("PCoA1 [", var16S[1], "%]"),
       y = paste0("PCoA2 [", var16S[2], "%]")) +
  theme_bw(base_size = 12)

# --- NMDS (Escalamiento Multidimensional No Métrico) ---
# NMDS es una ordenación basada en rangos: preserva el ORDEN de las distancias
# en lugar de sus valores exactos, lo que lo hace más robusto ante
# no-linealidades. Un valor de estrés < 0.20 indica representación 2D adecuada.
ord_nmds16S <- ordinate(ps16S_rel, method = "NMDS", distance = "bray")
cat("Stress NMDS:", round(ord_nmds16S$stress, 4), "\n")

# annotate() añade una capa de texto en una posición fija del gráfico.
# Inf / -Inf posiciona el texto en la esquina derecha inferior.
# hjust / vjust ajustan la alineación respecto a ese punto de anclaje.
plot_ordination(ps16S_rel, ord_nmds16S) +
  geom_point(size = 4, alpha = 0.9, color = "#e41a1c") +
  annotate("text", x = Inf, y = -Inf,
           label = paste("Stress =", round(ord_nmds16S$stress, 3)),
           hjust = 1.1, vjust = -0.5, size = 3.5) +
  labs(title = "NMDS — Bray-Curtis (16S)") +
  theme_bw(base_size = 12)


# ------------------------------------------------------------------------------
# 6.  VISUALIZACIONES
# ------------------------------------------------------------------------------

# ==== 6a.  Gráfico de barras taxonómico ====

# tax_glom() aglomera todos los ASVs que comparten la misma asignación a nivel
# de Filo, sumando sus conteos. NArm = TRUE elimina ASVs sin asignación de Filo.
ps16S_phylum <- tax_glom(ps16S_rel, taxrank = "Phylum", NArm = TRUE)

# Se ordenan los filos por abundancia total descendente y se toman los 10 mejores.
top10_ph <- names(sort(taxa_sums(ps16S_phylum), decreasing = TRUE))[1:10]

# prune_taxa() conserva solo los 10 filos más abundantes.
ps16S_top10 <- prune_taxa(top10_ph, ps16S_phylum)

# psmelt() convierte un objeto phyloseq a un data.frame 'largo' compatible
# con ggplot2 (una fila por combinación muestra × taxón).
df_bar16S <- psmelt(ps16S_top10)

# Gráfico de barras apiladas: cada barra = una muestra, colores = Filo.
# position = "stack" apila los segmentos; percent_format() formatea el eje y en %.
# scale_fill_brewer() usa una paleta cualitativa con suficientes colores distintos.
ggplot(df_bar16S, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack", width = 0.9) +
  scale_fill_brewer(palette = "Paired") +
  scale_y_continuous(labels = percent_format()) +
  labs(title = "Composición Taxonómica — Phylum (16S)",
       x = NULL, y = "Abundancia Relativa", fill = "Phylum") +
  theme_bw(base_size = 9) +
  theme(axis.text.x     = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
        legend.position = "bottom",
        legend.key.size = unit(0.4, "cm"))


# ==== 6b.  Mapa de calor taxonómico ====

# Se seleccionan los 30 ASVs con mayor abundancia relativa total.
top30_16S   <- names(sort(taxa_sums(ps16S_rel), decreasing = TRUE))[1:30]
ps16S_top30 <- prune_taxa(top30_16S, ps16S_rel)

# plot_heatmap() organiza muestras y taxa mediante ordenación NMDS para que
# los perfiles similares se agrupen visualmente.
# taxa.label = "Genus": etiqueta las filas con nombres de género en vez de IDs.
# low/high: extremos del gradiente de color; na.value: color para taxonomía ausente.
plot_heatmap(ps16S_top30, method = "NMDS", distance = "bray",
             taxa.label = "Genus",
             low = "#f7fbff", high = "#08306b", na.value = "white") +
  labs(title = "Heatmap — Top 30 ASVs (16S)") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# ------------------------------------------------------------------------------
# 7.  ANÁLISIS ESTADÍSTICOS
# ------------------------------------------------------------------------------

# ==== 7a.  PERMANOVA ====
# PERMANOVA (vegan::adonis2) prueba si la composición microbiana difiere
# significativamente entre grupos definidos por una variable de metadatos.
# Particiona una matriz de distancias en componentes dentro y entre grupos,
# y usa permutaciones para evaluar la significancia.

# Se extrae el metadata en un data.frame estándar (requerido por adonis2).
meta16S <- data.frame(sample_data(ps16S_rel))

# Se convierte el objeto dist a una matriz simétrica.
mat_bray16S <- as.matrix(dist_bray16S)

# Se alinea el orden de filas/columnas de la matriz con el orden de filas del
# metadata. Un orden no coincidente produciría pruebas de permutación incorrectas.
mat_bray16S <- mat_bray16S[rownames(meta16S), rownames(meta16S)]

# set.seed() garantiza que la prueba de permutación sea reproducible.
# La fórmula 'mat_bray16S ~ Run' prueba si 'Run' explica la disimilitud.
# AJUSTAR ~ Run a la variable de agrupación real de tus metadatos
# (p. ej., ~ grupo_tratamiento, ~ sitio, etc.).
# permutations = 999 significa 999 permutaciones aleatorias + el valor
# observado = 1000 en total.
set.seed(123)
perm16S <- vegan::adonis2(
  mat_bray16S ~ Run,
  data = meta16S, permutations = 999
)
print(perm16S)


# ==== 7b.  Kruskal-Wallis sobre índice Shannon ====
# Prueba no paramétrica de una vía: verifica si la diversidad Shannon difiere
# entre niveles de una variable de agrupación sin asumir normalidad.

# Se adjunta la variable 'Run' de sample_data a la tabla de diversidad alfa.
# Esto vincula cada valor alfa con el grupo de metadatos correspondiente.
alpha16S$Run <- sample_data(ps16S_rare)$Run

# kruskal.test() requiere una fórmula y el data.frame que contiene ambas variables.
kw16S <- kruskal.test(Shannon ~ Run, data = alpha16S)
cat("=== Kruskal-Wallis: Shannon (16S) ===\n")
print(kw16S)

# Si p < 0.05 puedes realizar comparaciones por pares con Wilcoxon
# (pairwise.wilcox.test con p.adjust.method = "BH") para identificar qué
# grupos difieren — consulta el tutorial phyloseq2 para un ejemplo completo.


# ==== 7c.  Abundancia diferencial (DESeq2) ====
# DESeq2 modela los conteos crudos con una distribución binomial negativa y
# prueba si cada ASV es significativamente más o menos abundante entre dos grupos.
# Este bloque está envuelto en requireNamespace() para que el script se ejecute
# incluso si DESeq2 no está instalado.

if (requireNamespace("DESeq2", quietly = TRUE)) {
  library(DESeq2)

  # USO: descomenta y adapta las líneas siguientes a tu diseño experimental.
  # Sustituye 'grupo' por tu columna de agrupación real, y "GrupoA"/"GrupoB"
  # por los niveles que quieres comparar.

  # ps_sub16S <- subset_samples(ps16S, grupo %in% c("GrupoA", "GrupoB"))
  # ps_sub16S <- prune_taxa(taxa_sums(ps_sub16S) > 5, ps_sub16S)
  #
  # phyloseq_to_deseq2() convierte el objeto phyloseq en un DESeqDataSet.
  # La fórmula de diseño ~ grupo indica a DESeq2 qué variable probar.
  # dds16S    <- phyloseq_to_deseq2(ps_sub16S, ~ grupo)
  #
  # estimateSizeFactors con type = "poscounts" maneja mejor los datos de
  # microbioma con muchos ceros que el enfoque de media geométrica por defecto.
  # dds16S    <- DESeq2::estimateSizeFactors(dds16S, type = "poscounts")
  #
  # DESeq() ajusta el modelo binomial negativo y realiza pruebas de Wald.
  # fitType = "mean" es un estimador de dispersión simple adecuado para datos
  # dispersos.
  # dds16S    <- DESeq2::DESeq(dds16S, fitType = "mean", quiet = TRUE)
  #
  # results() extrae fold changes y p-valores; cooksCutoff = FALSE conserva
  # todos los resultados sin filtrar valores atípicos (recomendado para
  # datos de microbioma).
  # res16S    <- DESeq2::results(dds16S, cooksCutoff = FALSE)
  #
  # Se filtran los ASVs significativos: se eliminan filas donde padj es NA
  # (no probado) y se conservan solo los con p-valor ajustado < 0.05 (FDR ≤ 5 %).
  # res_sig16S <- res16S[which(!is.na(res16S$padj) & res16S$padj < 0.05), ]
  # cat("ASVs diferencialmente abundantes (FDR < 5%):", nrow(res_sig16S), "\n")
}


# ------------------------------------------------------------------------------
# 8.  EXPORTACIÓN DE RESULTADOS
# ------------------------------------------------------------------------------

# saveRDS() guarda los objetos phyloseq filtrados y transformados en disco.
# Pueden recargarse en una nueva sesión de R con readRDS().
saveRDS(ps16S_f,    "03_Results/rds/16S/ps16S_filtrado.RDS")
saveRDS(ps16S_rel,  "03_Results/rds/16S/ps16S_rel.RDS")
saveRDS(ps16S_rare, "03_Results/rds/16S/ps16S_rare.RDS")

# write.csv() exporta cada tabla de datos como CSV de texto plano.
# as.data.frame() es necesario porque otu_table y tax_table son objetos S4,
# no data.frames estándar; write.csv() no puede manejarlos directamente.
write.csv(as.data.frame(otu_table(ps16S_f)),   "03_Results/csv/16S/tabla_asv.csv")
write.csv(as.data.frame(tax_table(ps16S_f)),   "03_Results/csv/16S/taxonomia.csv")
write.csv(as.data.frame(sample_data(ps16S_f)), "03_Results/csv/16S/metadatos.csv")

# alpha16S es un data.frame estándar, por lo que write.csv() funciona directamente.
write.csv(alpha16S, "03_Results/csv/16S/diversidad_alfa.csv")