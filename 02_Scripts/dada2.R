# ─────────────────────────────────────────────────────────────────────────────

# PIPELINE DE ASV — 16S & 18S/ITS (DADA2)

# ─────────────────────────────────────────────────────────────────────────────

# 1. INSTALACIÓN DE PAQUETES

# requireNamespace() comprueba si "BiocManager" ya está instalado sin cargarlo.
# El operador ! niega el resultado, por lo que la condición es verdadera cuando
# el paquete NO está disponible.
# quietly = TRUE suprime mensajes de advertencia durante la comprobación.

# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

# BiocManager::install(version = ...) actualiza o instala el núcleo de
# Bioconductor a la versión indicada (3.23). Esto garantiza compatibilidad
# entre todos los paquetes de Bioconductor que se usarán a continuación.

# BiocManager::install(version = "3.23")

# BiocManager::install() con un vector de nombres instala varios paquetes
# de Bioconductor de una sola vez:
#   dada2      — pipeline principal de denoising de amplicones.
#   ShortRead  — lectura e inspección de archivos FASTQ.
#   Biostrings — manipulación de cadenas de ADN .
#   phyloseq   — análisis y visualizar datos de metagenómica.
#   DESeq2     — pruebas de abundancia diferencial.
#   DECIPHER   — asignación taxonómica con IdTaxa.
# force = TRUE fuerza la reinstalación aunque ya estén instalados.

# BiocManager::install(
#  c("dada2", "ShortRead", "Biostrings", "phyloseq", "DESeq2", "DECIPHER"),
#  force = TRUE)

# install.packages() instala paquetes del repositorio CRAN:
#   ggplot2      — visualización de datos.
#   vegan        — métricas de diversidad ecológica.
#   scales       — transformaciones de ejes en gráficos.
#   RColorBrewer — paletas de colores categóricas.
#   ggpubr       — figuras listas para publicación basadas en ggplot2.
#   knitr        — generación de reportes reproducibles.

# install.packages(c("ggplot2", "vegan", "scales",
#                   "RColorBrewer", "ggpubr", "knitr"))

# ─────────────────────────────────────────────────────────────────────────────

# 2. CARGA DE LIBRERIAS

# library() carga cada paquete instalado en la sesión activa de R.

library(dada2)        # funciones del pipeline DADA2.
library(ShortRead)    # lectura/escritura de FASTQ.
library(DESeq2)       # abundancia diferencial.
library(Biostrings)   # operaciones sobre cadenas de ADN.
library(phyloseq)     # estructuras de datos para microbiomas.
library(vegan)        # ecología de comunidades.
library(ggplot2)      # visualización de gráficos.
library(RColorBrewer) # paletas de colores para gráficos.
library(ggpubr)       # utilidades de publicación para ggplot2.
library(scales)       # funciones auxiliares para formatos de ejes.
library(knitr)        # genera tablas formateadas en HTML/PDF dentro de reportes.
library(DECIPHER)     # clasificador IdTaxa.


# ─────────────────────────────────────────────────────────────────────────────

# 3. METADATOS DE MUESTRAS

# read.csv() lee un archivo de texto delimitado por comas y lo convierte en
# un data.frame de R.
# stringsAsFactors = FALSE evita que R convierta automáticamente las columnas.
# de texto en factores, lo cual facilita la manipulación posterior como cadenas.

# La ruta apunta a la tabla de metadatos descargada del
# NCBI SRA para el experimento de 16S (bacterioma) e ITS (fungoma).

metadata16S <- read.csv("01_RawData/csv/16S/SraRunTable16S.csv",
                        stringsAsFactors = FALSE)

metadataITS <- read.csv("01_RawData/csv/ITS/SraRunTableITS.csv",
                        stringsAsFactors = FALSE)

# El operador $ accede a una columna.
# Aquí se añade la columna "amplicon" a cada tabla con el valor correspondiente,
# para poder distinguir el origen de cada muestra después de unir ambas tablas.

metadata16S$amplicon <- "16S"
metadataITS$amplicon <- "ITS"

# rbind() ("row bind") apila dos data.frames uno encima del otro por filas,
# combinando ambas tablas de metadatos en una sola. Las columnas deben coincidir.
# Resultado: un único data.frame con todas las muestras de ambos amplicones.

metadata <- rbind(metadata16S, metadataITS)

# cat() imprime texto en la consola concatenando sus argumentos.
# sum(metadata$amplicon == "16S") cuenta cuántas filas tienen el valor "16S"
# en la columna amplicon (TRUE cuenta como 1, FALSE como 0).
# \n es el carácter de nueva línea.

cat("Muestras 16S:", sum(metadata$amplicon == "16S"), "\n")
cat("Muestras ITS:", sum(metadata$amplicon == "ITS"),  "\n")

# head() muestra las primeras 6 filas del data.frame en la consola.
# Sirve como verificación rápida de que la unión fue correcta.

head(metadata)

# write.csv() exporta el data.frame combinado como archivo CSV en el directorio

write.csv(metadata, "03_Results/csv/metadata/metadata.csv")

# Indexación lógica: metadata$Run devuelve el vector de accesiones SRR,
# y el corchete [] lo filtra conservando solo las que pertenecen al amplicón 16S.
# Este vector se usará para localizar los archivos FASTQ correspondientes.

srr16S <- metadata$Run[metadata$amplicon == "16S"]
srrITS <- metadata$Run[metadata$amplicon == "ITS"]

# ─────────────────────────────────────────────────────────────────────────────

# 4. RUTAS A LOS ARCHIVOS FASTQ


# Se define la ruta raíz donde están los FASTQ crudos del 16S como cadena de texto.

path16S <- "01_RawData/fastq/16S/"

# list.files() lista todos los archivos en ese directorio: verificación de que
# la ruta existe y contiene los archivos esperados antes de continuar.

list.files(path16S)

# Lo mismo para ITS/18S.

pathITS <- "01_RawData/fastq/ITS/"

list.files(pathITS)

# list.files() con pattern = "_1.fastq" selecciona solo los archivos de lecturas
# forward (R1). 
# full.names = TRUE devuelve la ruta completa, no solo el nombre.
# sort() ordena alfabéticamente para garantizar que los pares R1/R2 queden
# alineados en el mismo orden (muestra i de fnFs corresponde a muestra i de fnRs).

fnFs16S <- sort(list.files(path16S, pattern = "_1.fastq", full.names = TRUE))
fnRs16S <- sort(list.files(path16S, pattern = "_2.fastq", full.names = TRUE))

fnFsITS <- sort(list.files(pathITS, pattern = "_1.fastq", full.names = TRUE))
fnRsITS <- sort(list.files(pathITS, pattern = "_2.fastq", full.names = TRUE))

# Vector de caracteres con los IDs SRR de las muestras control
# que deben excluirse del análisis para evitar sesgo.

remove_ids <- c(
  "SRR25070675",   # control 16S
  "SRR25070676",   # control 16S
  "SRR25100830",   # control ITS
  "SRR25100831"    # control ITS
)

# paste(..., collapse = "|") une los IDs con el separador "|" (operador OR en
# expresiones regulares), creando un patrón como "ID1|ID2|ID3|ID4".
# Así, grepl() podrá encontrar cualquiera de los IDs en una sola búsqueda.

patron <- paste(remove_ids, collapse = "|")

# grepl(patron) devuelve TRUE para cada ruta que contiene algún ID de
# control. El operador ! invierte el resultado, y la indexación [] conserva
# solo las rutas que NO corresponden a controles.

fnFs16S <- fnFs16S[!grepl(patron, fnFs16S)]
fnRs16S <- fnRs16S[!grepl(patron, fnRs16S)]
fnFsITS <- fnFsITS[!grepl(patron, fnFsITS)]
fnRsITS <- fnRsITS[!grepl(patron, fnRsITS)]

# basename() extrae solo el nombre del archivo (sin la ruta).
# strsplit(..., "_") divide el nombre en partes usando "_" como separador.
# `[`(1) extrae el primer elemento (el acceso SRR), que es el nombre de muestra.
# sapply() aplica esto a cada elemento del vector y devuelve un vector de nombres.

samples.names16S <- sapply(strsplit(basename(fnFs16S), "_"), `[`, 1)
samples.namesITS <- sapply(strsplit(basename(fnFsITS), "_"), `[`, 1)

# length() devuelve cuántos elementos tiene el vector; se reporta para confirmar
# que la exclusión de controles fue correcta.

cat("Muestras 16S:", length(samples.names16S), "\n")
cat("Muestras ITS:", length(samples.namesITS),  "\n")

# Imprime los vectores de nombres en la consola para inspección visual.

samples.names16S
samples.namesITS

# saveRDS() serializa un objeto de R en un archivo binario .RDS en el directorio.
# Esto permite recuperarlo exactamente igual en sesiones futuras sin recalcular.

saveRDS(samples.names16S, file = "03_Results/rds/16S/samples.names16S.RDS")
saveRDS(samples.namesITS, file = "03_Results/rds/ITS/samples.namesITS.RDS")

# readRDS() carga el objeto guardado de vuelta a la sesión.
# Útil para retomar el análisis desde este punto sin volver a ejecutar lo anterior.

samples.names16S <- readRDS("03_Results/rds/16S/samples.names16S.RDS")
samples.namesITS <- readRDS("03_Results/rds/ITS/samples.namesITS.RDS")


# ─────────────────────────────────────────────────────────────────────────────

# 5. ELIMINACIÓN DE PRIMERS CON CUTADAPT — 16S V4-V5 (515F / 926R)

# Script utilizado de la página web: https://benjjneb.github.io/dada2/tutorial.html
# Los primers fueron seleccionado de acuerdo a lo publicado en este paper: 
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0212355

# Definir secuencias de primers:

# Secuencia del primer forward 515F para la región V4-V5 del gen 16S rRNA.
# Las letras IUPAC ambiguas (Y, M) representan mezclas de bases en la posición.

FWD.16S <- "GTGYCAGCMGCCGCGGTAA"
REV.16S <- "CCGYCAATTYMTTTRAGTTT"

# Función auxiliar: genera las cuatro orientaciones posibles de un primer
# (Forward, Complemento, Reverso, Reverso-Complemento) como cadenas de texto.
# Esto es necesario porque los primers pueden aparecer en cualquier orientación
# dentro de las lecturas según la química de la librería.

allOrients <- function(primer) {
  require(Biostrings)          # asegura que Biostrings esté disponible
  dna <- DNAString(primer)     # convierte la cadena de texto en objeto DNAString
                               # para poder aplicar operaciones biológicas
  orients <- c(
    Forward    = dna,
    Complement = Biostrings::complement(dna),       # complemento base a base (5'→3')
    Reverse    = Biostrings::reverse(dna),           # secuencia leída al revés
    RevComp    = Biostrings::reverseComplement(dna) # reverso-complemento (artefacto más común)
  )
  return(sapply(orients, toString))  # convierte cada DNAString a texto plano
}

# Se calculan todas las orientaciones de cada primer para usarlas en la búsqueda.
FWD.orients.16S <- allOrients(FWD.16S)
REV.orients.16S <- allOrients(REV.16S)
FWD.orients.16S   # imprime para verificar las cuatro orientaciones

# Filtrar y contar primers antes de cortar

# file.path() construye rutas de archivo de forma portable.
# Se crea la ruta de salida para lecturas filtradas en el subdirectorio "filtN".

fnFs16S.filtN <- file.path(path16S, "filtN", basename(fnFs16S))
fnRs16S.filtN <- file.path(path16S, "filtN", basename(fnRs16S))

# filterAndTrim() filtra y recorta lecturas paired-end.
# Aquí solo se usa para eliminar lecturas con bases N (maxN = 0), lo cual es
# un requisito previo al aprendizaje del modelo de error de DADA2.
# multithread = TRUE acelera el proceso usando múltiples núcleos (usar FALSE en Windows).

filterAndTrim(fnFs16S, fnFs16S.filtN,
              fnRs16S, fnRs16S.filtN,
              maxN = 0, multithread = TRUE)

# Función que cuenta cuántas lecturas de un archivo FASTQ contienen un primer.
# vcountPattern() cuenta coincidencias del patrón en cada secuencia.
# sread(readFastq(fn)) lee las secuencias del archivo FASTQ.
# fixed = FALSE activa la interpretación de códigos IUPAC ambiguos en el patrón.
# sum(nhits > 0) cuenta cuántas lecturas tienen al menos una coincidencia.

primerHits <- function(primer, fn) {
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

# rbind() construye una tabla donde cada fila es una combinación dirección/orientación.
# sapply() aplica primerHits a cada orientación del primer y devuelve un vector.
# Si los primers están presentes (antes de cortar), se esperan cuentas > 0.

rbind(
  FWD.ForwardReads = sapply(FWD.orients.16S, primerHits, fn = fnFs16S.filtN[[1]]),
  FWD.ReverseReads = sapply(FWD.orients.16S, primerHits, fn = fnRs16S.filtN[[1]]),
  REV.ForwardReads = sapply(REV.orients.16S, primerHits, fn = fnFs16S.filtN[[1]]),
  REV.ReverseReads = sapply(REV.orients.16S, primerHits, fn = fnRs16S.filtN[[1]])
)

# Ejecutar Cutadapt

# Ruta al ejecutable de Cutadapt instalado en el sistema.
# Cambiar si Cutadapt está instalado en otra ubicación.

cutadapt <- "/usr/bin/cutadapt"

# system2() ejecuta un comando del sistema operativo desde R.
# args = "--version" imprime la versión instalada como verificación.

system2(cutadapt, args = "--version")

# Se construye la ruta al directorio de salida para lecturas recortadas del 16S.

path.cut.16S <- file.path(path16S, "cutadapt")

# dir.exists() verifica si el directorio ya existe; si no, dir.create() lo crea.
# Esto evita errores si el directorio ya fue creado en una ejecución anterior.

if (!dir.exists(path.cut.16S)) dir.create(path.cut.16S)

# Se construyen las rutas de los archivos de salida recortados (mismo nombre que
# la entrada pero dentro del directorio "cutadapt").

fnFs16S.cut <- file.path(path.cut.16S, basename(fnFs16S))
fnRs16S.cut <- file.path(path.cut.16S, basename(fnRs16S))

# dada2:::rc() calcula el reverso-complemento de la secuencia del primer.
# El triple ::: accede a una función interna (no exportada) del paquete dada2.
# Los RC de los primers son necesarios para recortar artefactos de read-through
# (cuando la lectura es más larga que el amplicón y atraviesa el otro primer).

FWD.RC.16S <- dada2:::rc(FWD.16S)
REV.RC.16S <- dada2:::rc(REV.16S)

# paste() construye cadenas con los argumentos de Cutadapt separados por espacio:
# -g: recortar el primer en el extremo 5' de R1 (forward)
# -a: recortar el RC del reverse en el extremo 3' de R1 (artefacto read-through)
# -G y -A: análogos para R2 (reverse)

R1.flags.16S <- paste("-g", FWD.16S, "-a", REV.RC.16S)
R2.flags.16S <- paste("-G", REV.16S, "-A", FWD.RC.16S)

# Bucle for sobre todas las muestras: se ejecuta Cutadapt una vez por muestra.
# seq_along(fnFs16S) genera la secuencia 1, 2, ..., n donde n es el nº de muestras.
# En cada iteración i, se procesan el par R1[i] / R2[i] correspondiente.

for (i in seq_along(fnFs16S)) {
  system2(cutadapt,
          args = c(R1.flags.16S, R2.flags.16S,
                   "-n", 2,               # permitir hasta 2 ocurrencias de adaptadores por lectura
                   "-o", fnFs16S.cut[i],  # archivo de salida R1 recortado
                   "-p", fnRs16S.cut[i],  # archivo de salida R2 recortado
                   fnFs16S.filtN[i],      # archivo de entrada R1 filtrado de Ns
                   fnRs16S.filtN[i]))     # archivo de entrada R2 filtrado de Ns
}

# Verificación post-corte: todas las cuentas deben ser 0 o cercanas a 0.
# Si quedan cuentas altas, puede indicar un problema con los flags de Cutadapt.

rbind(
  FWD.ForwardReads = sapply(FWD.orients.16S, primerHits, fn = fnFs16S.cut[[1]]),
  FWD.ReverseReads = sapply(FWD.orients.16S, primerHits, fn = fnRs16S.cut[[1]]),
  REV.ForwardReads = sapply(REV.orients.16S, primerHits, fn = fnFs16S.cut[[1]]),
  REV.ReverseReads = sapply(REV.orients.16S, primerHits, fn = fnRs16S.cut[[1]])
)

# Se actualizan los vectores de rutas para apuntar a los archivos recortados,
# que serán la entrada para el resto del pipeline DADA2.

cutFs16S <- sort(list.files(path.cut.16S, pattern = "_1.fastq", full.names = TRUE))
cutRs16S <- sort(list.files(path.cut.16S, pattern = "_2.fastq", full.names = TRUE))

# Función auxiliar: extrae el nombre de muestra del nombre de archivo recortado.

get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]

# unname() elimina los nombres del vector resultante (los nombres serían las rutas).

sample.names16S <- unname(sapply(cutFs16S, get.sample.name))

head(sample.names16S)   # verificar los primeros nombres

# Guardar los vectores de rutas recortadas para poder reanudar la sesión.

saveRDS(cutFs16S, file = "03_Results/rds/16S/cutFs16S.RDS")
saveRDS(cutRs16S, file = "03_Results/rds/16S/cutRs16S.RDS")

# Punto de reanudación: cargar desde disco sin necesidad de reejecutar Cutadapt.

cutFs16S <- readRDS("03_Results/rds/16S/cutFs16S.RDS")
cutRs16S <- readRDS("03_Results/rds/16S/cutRs16S.RDS")


# ─────────────────────────────────────────────────────────────────────────────

# 6. ELIMINACIÓN DE PRIMERS CON CUTADAPT — 18S V4


# Definir primers

FWD.18S <- "CCAGCASCYGCGGTAATTCC"  # primer forward V4-18S
REV.18S <- "ACTTTCGTTCTTGATYRA"    # primer reverse V4-18S

# Misma función allOrients definida anteriormente para el 16S.

FWD.orients.18S <- allOrients(FWD.18S)
REV.orients.18S <- allOrients(REV.18S)
FWD.orients.18S   # inspeccionar

# Filtrar Ns y verificar

# Rutas de salida para las lecturas ITS filtradas de Ns.

fnFs18S.filtN <- file.path(pathITS, "filtN", basename(fnFsITS))
fnRs18S.filtN <- file.path(pathITS, "filtN", basename(fnRsITS))

# Eliminar lecturas con bases N del conjunto ITS.
# multithread = FALSE es la opción segura en Windows.

filterAndTrim(fnFsITS, fnFs18S.filtN,
              fnRsITS, fnRs18S.filtN,
              maxN = 0, multithread = FALSE)

# Verificar presencia de primers antes de cortar (se esperan cuentas > 0).

rbind(
  FWD.ForwardReads = sapply(FWD.orients.18S, primerHits, fn = fnFs18S.filtN[[1]]),
  FWD.ReverseReads = sapply(FWD.orients.18S, primerHits, fn = fnRs18S.filtN[[1]]),
  REV.ForwardReads = sapply(REV.orients.18S, primerHits, fn = fnFs18S.filtN[[1]]),
  REV.ReverseReads = sapply(REV.orients.18S, primerHits, fn = fnRs18S.filtN[[1]])
)

# Ejecutar Cutadapt y verificar

# Directorio de salida para lecturas ITS recortadas.

path.cut.18S <- file.path(pathITS, "cutadapt")
if (!dir.exists(path.cut.18S)) dir.create(path.cut.18S)

# Rutas de salida de los archivos recortados para el 18S/ITS.

fnFs18S.cut <- file.path(path.cut.18S, basename(fnFsITS))
fnRs18S.cut <- file.path(path.cut.18S, basename(fnRsITS))

# Calcular reversos-complementos de los primers del 18S.

FWD.RC.18S <- dada2:::rc(FWD.18S)
REV.RC.18S <- dada2:::rc(REV.18S)

# Construir los flags de Cutadapt para el 18S.

R1.flags.18S <- paste("-g", FWD.18S, "-a", REV.RC.18S)
R2.flags.18S <- paste("-G", REV.18S, "-A", FWD.RC.18S)

# Ejecutar Cutadapt para cada muestra ITS.

for (i in seq_along(fnFsITS)) {
  system2(cutadapt,
          args = c(R1.flags.18S, R2.flags.18S,
                   "-n", 2,
                   "-o", fnFs18S.cut[i],
                   "-p", fnRs18S.cut[i],
                   fnFs18S.filtN[i],
                   fnRs18S.filtN[i]))
}

# Confirmar que los primers fueron eliminados (todas las cuentas deben ser ~0).

rbind(
  FWD.ForwardReads = sapply(FWD.orients.18S, primerHits, fn = fnFs18S.cut[[1]]),
  FWD.ReverseReads = sapply(FWD.orients.18S, primerHits, fn = fnRs18S.cut[[1]]),
  REV.ForwardReads = sapply(REV.orients.18S, primerHits, fn = fnFs18S.cut[[1]]),
  REV.ReverseReads = sapply(REV.orients.18S, primerHits, fn = fnRs18S.cut[[1]])
)

# Actualizar vectores de rutas para apuntar a los archivos recortados del 18S.

cutFs18S <- sort(list.files(path.cut.18S, pattern = "_1.fastq", full.names = TRUE))
cutRs18S <- sort(list.files(path.cut.18S, pattern = "_2.fastq", full.names = TRUE))

# Extraer nombres de muestra de los archivos 18S recortados.

sample.names18S <- unname(sapply(cutFs18S, get.sample.name))
head(sample.names18S)

saveRDS(cutFs18S, file = "03_Results/rds/ITS/cutFsITS.RDS")
saveRDS(cutRs18S, file = "03_Results/rds/ITS/cutRsITS.RDS")

# Punto de reanudación.
cutFs18S <- readRDS("03_Results/rds/ITS/cutFsITS.RDS")
cutRs18S <- readRDS("03_Results/rds/ITS/cutRsITS.RDS")

# ─────────────────────────────────────────────────────────────────────────────

# 7. PERFILES DE CALIDAD PHRED

# file.info()$size devuelve el tamaño en bytes de cada archivo.
# Se filtran archivos < 50 bytes porque están vacíos o prácticamente vacíos,
# lo que causaría errores al intentar graficar su perfil de calidad.

# Ayuda del chat porque dos plots no se generaron y marcaba error.

cutFs16S.nonempty <- cutFs16S[file.info(cutFs16S)$size > 20]
cutRs16S.nonempty <- cutRs16S[file.info(cutRs16S)$size > 20]


# Para diagnóstico adicional: contar el número de lecturas en los primeros 20 archivos.
# length(readFastq(f)) devuelve cuántas lecturas tiene el archivo f.

sapply(cutFs16S[1:20], function(f) length(readFastq(f)))

# pdf() abre un dispositivo gráfico y redirige todas las gráficas al archivo PDF.
# width y height controlan el tamaño en pulgadas.
# plotQualityProfile() genera el perfil de calidad PHRED de los primeros 20 archivos.
# Estos perfiles se usan para decidir los parámetros de truncLen en el paso siguiente.
# dev.off() cierra el dispositivo gráfico y termina de escribir el PDF.

# 16S — lecturas forward
pdf("03_Results/phred/16S/quality-forward-16S.pdf", width = 20, height = 12)
plotQualityProfile(cutFs16S[1:20])
dev.off()

# 16S — lecturas reverse
pdf("03_Results/phred/16S/quality-reverse-16S.pdf", width = 20, height = 12)
plotQualityProfile(cutRs16S.nonempty[1:20])
dev.off()

# Sigue el error: Error in density.default(qscore): 'x' contains missing values
# lectura dentro del archivo FASTQ tiene longitud 0 después del recorte de Cutadapt

# ITS — lecturas forward
pdf("03_Results/phred/ITS/quality-forward-ITS.pdf", width = 20, height = 12)
plotQualityProfile(cutFs18S[1:20])
dev.off()

# ITS — lecturas reverse
pdf("03_Results/phred/ITS/quality-reverse-ITS.pdf", width = 20, height = 12)
plotQualityProfile(cutRs18S[1:20])
dev.off()

# ─────────────────────────────────────────────────────────────────────────────

# 8. FILTRADO DE CALIDAD

# Se construyen las rutas de salida para las lecturas filtradas por calidad,
# usando un subdirectorio "filtered" dentro de la carpeta de lecturas recortadas.

filtFs16S <- file.path(path.cut.16S, "filtered", basename(cutFs16S))
filtRs16S <- file.path(path.cut.16S, "filtered", basename(cutRs16S))
filtFsITS <- file.path(path.cut.18S, "filtered", basename(cutFs18S))
filtRsITS <- file.path(path.cut.18S, "filtered", basename(cutRs18S))

# names() asigna nombres a los elementos del vector de rutas.
# DADA2 requiere que los vectores de archivos filtrados tengan los nombres de
# muestra para asociar correctamente los resultados a cada muestra.

names(filtFs16S) <- samples.names16S
names(filtRs16S) <- samples.names16S
names(filtFsITS) <- samples.namesITS
names(filtRsITS) <- samples.namesITS

# filterAndTrim() aplica el filtrado de calidad completo a las lecturas 16S.
# El resultado asv16S es una tabla con el conteo de lecturas antes y después
# del filtrado para cada muestra (útil para el tracking posterior).

truncLen = c(240, 160) # trunca R1 a 240 pb y R2 a 160 pb
                       #  (ajustar según los perfiles PHRED del paso anterior)
   maxN = 0            # elimina cualquier lectura que contenga una N
   maxEE = c(2, 2)     # descarta lecturas con más de 2 errores esperados
                       #  (errores esperados = suma de probabilidades de error por base)
   truncQ = 2          # trunca en la primera base con calidad ≤ 2
   rm.phix = TRUE      # elimina lecturas que mapean contra PhiX (control interno)
   compress = TRUE     # escribe los archivos de salida en formato gzip
   multithread = TRUE  # usa múltiples núcleos (FALSE en Windows)

asv16S <- filterAndTrim(
  cutFs16S, filtFs16S,
  cutRs16S, filtRs16S,
  truncLen    = c(240, 160),
  maxN        = 0,
  maxEE       = c(2, 2),
  truncQ      = 2,
  rm.phix     = TRUE,
  compress    = TRUE,
  multithread = TRUE
)

# Filtrado de calidad para lecturas ITS/18S con los mismos parámetros.
# CORRECCIÓN: el original pasaba fnFsITS/fnRsITS (archivos crudos) como entrada,
# en lugar de los archivos recortados cutFs18S/cutRs18S. Corregido aquí.
asvITS <- filterAndTrim(
  cutFs18S, filtFsITS,
  cutRs18S, filtRsITS,
  truncLen    = c(240, 160),
  maxN        = 0,
  maxEE       = c(2, 2),
  truncQ      = 2,
  rm.phix     = TRUE,
  compress    = TRUE,
  multithread = TRUE
)

# Guardar las tablas de resumen del filtrado.
saveRDS(asv16S, file = "03_Results/rds/16S/asv16S.RDS")
saveRDS(asvITS, file = "03_Results/rds/ITS/asvITS.RDS")

# Punto de reanudación.
asv16S <- readRDS("03_Results/rds/16S/asv16S.RDS")
asvITS <- readRDS("03_Results/rds/ITS/asvITS.RDS")


# ─────────────────────────────────────────────────────────────────────────────
# 9. APRENDIZAJE DEL MODELO DE ERRORES
# ─────────────────────────────────────────────────────────────────────────────
# Las líneas activas cargan modelos de error pre-calculados desde disco.
# Descomentar los pares learnErrors() + saveRDS() para calcularlos desde cero.

# learnErrors() estima las tasas de error de secuenciación de Illumina a partir
# de las lecturas filtradas. DADA2 modela la probabilidad de que una base
# observada sea el resultado de un error a partir de la base "real".
# Este modelo es fundamental para distinguir errores de variación biológica real.

# Modelo de error forward del 16S (cargado desde RDS pre-calculado)
errF16S <- readRDS("03_Results/rds/16S/errF16S.RDS")
# errF16S <- learnErrors(filtFs16S, multithread = TRUE)
# saveRDS(errF16S, file = "03_Results/rds/16S/errF16S.RDS")

# Modelo de error reverse del 16S
errR16S <- readRDS("03_Results/rds/16S/errR16S.RDS")
# errR16S <- learnErrors(filtRs16S, multithread = TRUE)
# saveRDS(errR16S, file = "03_Results/rds/16S/errR16S.RDS")

# Modelo de error forward del ITS
errFITS <- readRDS("03_Results/rds/ITS/errFITS.RDS")
# errFITS <- learnErrors(filtFsITS, multithread = TRUE)
# saveRDS(errFITS, file = "03_Results/rds/ITS/errFITS.RDS")

# Modelo de error reverse del ITS
errRITS <- readRDS("03_Results/rds/ITS/errRITS.RDS")
# errRITS <- learnErrors(filtRsITS, multithread = TRUE)
# saveRDS(errRITS, file = "03_Results/rds/ITS/errRITS.RDS")

# plotErrors() visualiza el modelo de error aprendido.
# Los puntos son las tasas observadas; la línea es el ajuste del modelo.
# Un buen modelo muestra puntos que siguen aproximadamente la línea.
# nominalQ = TRUE superpone la tasa de error esperada según el score PHRED nominal.

# png() abre un dispositivo de imagen PNG (más liviano que PDF para figuras únicas).
png("03_Results/figures/16S/errores-forward-16S.png")
plotErrors(errF16S, nominalQ = TRUE)
dev.off()

png("03_Results/figures/16S/errores-reverse-16S.png")
plotErrors(errR16S, nominalQ = TRUE)
dev.off()

png("03_Results/figures/ITS/errores-forward-ITS.png")
plotErrors(errFITS, nominalQ = TRUE)
dev.off()

png("03_Results/figures/ITS/errores-reverse-ITS.png")
plotErrors(errRITS, nominalQ = TRUE)
dev.off()


# ─────────────────────────────────────────────────────────────────────────────
# 10. INFERENCIA DE ASVs Y FUSIÓN DE PARES
# ─────────────────────────────────────────────────────────────────────────────

# dada() aplica el algoritmo de denoising de DADA2 a las lecturas filtradas,
# usando el modelo de error aprendido para distinguir variantes reales de errores.
# Se ejecuta por separado para lecturas forward y reverse.
# El resultado es un objeto dada con las variantes de secuencia de amplicón (ASVs).
dadaFs16S <- dada(filtFs16S, err = errF16S, multithread = TRUE)
dadaRs16S <- dada(filtRs16S, err = errR16S, multithread = TRUE)
dadaFsITS <- dada(filtFsITS, err = errFITS, multithread = TRUE)
dadaRsITS <- dada(filtRsITS, err = errRITS, multithread = TRUE)

# mergePairs() fusiona las lecturas forward y reverse denoiseadas en secuencias
# de longitud completa. Requiere que las lecturas se superpongan suficientemente.
# Se carga desde RDS; descomentar para recalcular.
# mergers16S <- mergePairs(dadaFs16S, filtFs16S, dadaRs16S, filtRs16S, verbose = TRUE)
# saveRDS(mergers16S, file = "03_Results/rds/16S/mergers16S.RDS")
mergers16S <- readRDS("03_Results/rds/16S/mergers16S.RDS")

# makeSequenceTable() construye la tabla de secuencias (ASV table):
# filas = muestras, columnas = secuencias de ASV, valores = número de lecturas.
# Esta es la tabla central del análisis, equivalente a una OTU table.
seqtab16S <- makeSequenceTable(mergers16S)

# Mismo procedimiento para ITS.
# mergersITS <- mergePairs(dadaFsITS, filtFsITS, dadaRsITS, filtRsITS, verbose = TRUE)
# saveRDS(mergersITS, file = "03_Results/rds/ITS/mergersITS.RDS")
mergersITS <- readRDS("03_Results/rds/ITS/mergersITS.RDS")
seqtabITS <- makeSequenceTable(mergersITS)

# dim(seqtab16S) mostraría [n_muestras × n_ASVs] — útil para inspección.
# saveRDS(seqtab16S, file = "03_Results/rds/16S/seqtab16S.RDS")
seqtab16S <- readRDS("03_Results/rds/16S/seqtab16S.RDS")

# dim(seqtabITS)
# saveRDS(seqtabITS, file = "03_Results/rds/ITS/seqtabITS.RDS")
seqtabITS <- readRDS("03_Results/rds/ITS/seqtabITS.RDS")


# ─────────────────────────────────────────────────────────────────────────────
# 11. ELIMINACIÓN DE QUIMERAS
# ─────────────────────────────────────────────────────────────────────────────

# removeBimeraDenovo() elimina secuencias quiméricas (bimeras) de la tabla de ASVs.
# Las quimeras son artefactos de PCR formados por la fusión de dos secuencias reales.
# method = "consensus" requiere que la quimera sea reconocida como tal en la
# mayoría de las muestras antes de eliminarla.
# verbose = TRUE imprime el número de bimeras detectadas y eliminadas.
seqtab.nochim16S <- removeBimeraDenovo(seqtab16S, method = "consensus",
                                        multithread = TRUE, verbose = TRUE)
# Calcular la fracción de lecturas conservadas después de eliminar quimeras.
# sum(seqtab.nochim16S) / sum(seqtab16S): si este valor es bajo (< 0.7),
# puede indicar un problema con los primers o el protocolo de librería.
cat("Fraccion de lecturas 16S conservadas:",
    round(sum(seqtab.nochim16S) / sum(seqtab16S), 4), "\n")

seqtab.nochimITS <- removeBimeraDenovo(seqtabITS, method = "consensus",
                                        multithread = TRUE, verbose = TRUE)
cat("Fraccion de lecturas ITS conservadas:",
    round(sum(seqtab.nochimITS) / sum(seqtabITS), 4), "\n")


# ─────────────────────────────────────────────────────────────────────────────
# 12. TABLA DE SEGUIMIENTO DE LECTURAS
# ─────────────────────────────────────────────────────────────────────────────

# Función auxiliar: getUniques() extrae las secuencias únicas de un objeto DADA2;
# sum() suma sus abundancias para obtener el total de lecturas en ese objeto.
getN <- function(x) sum(getUniques(x))

# cbind() ("column bind") une columnas para construir la tabla de seguimiento 16S.
# Cada columna representa el conteo de lecturas en una etapa del pipeline:
#   asv16S             — lecturas crudas (input) y después del filtrado de calidad
#   sapply(dadaFs16S)  — después del denoising forward
#   sapply(dadaRs16S)  — después del denoising reverse
#   sapply(mergers16S) — después de la fusión de pares
#   rowSums(...)       — después de la eliminación de quimeras
track16S <- cbind(
  asv16S,
  sapply(dadaFs16S,  getN),
  sapply(dadaRs16S,  getN),
  sapply(mergers16S, getN),
  rowSums(seqtab.nochim16S)
)
# colnames() y rownames() asignan etiquetas a columnas y filas de la matriz.
colnames(track16S) <- c("input", "filtered", "denoisedF",
                         "denoisedR", "merged", "nonchim")
rownames(track16S) <- samples.names16S

head(track16S)
write.csv(track16S, "03_Results/csv/16S/track16S.csv")

# Lo mismo para ITS.
trackITS <- cbind(
  asvITS,
  sapply(dadaFsITS,  getN),
  sapply(dadaRsITS,  getN),
  sapply(mergersITS, getN),
  rowSums(seqtab.nochimITS)
)
colnames(trackITS) <- c("input", "filtered", "denoisedF",
                         "denoisedR", "merged", "nonchim")
rownames(trackITS) <- samples.namesITS

head(trackITS)
write.csv(trackITS, "03_Results/csv/ITS/trackITS.csv")


# ─────────────────────────────────────────────────────────────────────────────
# 13. ASIGNACIÓN TAXONÓMICA — Bayesiano ingenuo de DADA2 (assignTaxonomy / addSpecies)
# ─────────────────────────────────────────────────────────────────────────────

# ── 16S — base de datos SILVA ────────────────────────────────────────────────

# assignTaxonomy() clasifica las secuencias de ASV contra una base de datos de
# referencia usando el clasificador bayesiano ingenuo de Wang et al. (2007).
# Devuelve una matriz con la clasificación hasta género para cada ASV.
# El archivo .fa.gz es el conjunto de entrenamiento de SILVA v138.1.
taxa16S <- assignTaxonomy(seqtab.nochim16S,
                           "01_RawData/silva/silva_nr99_v138.1_train_set.fa.gz",
                           multithread = TRUE)

# addSpecies() extiende la clasificación hasta nivel de especie mediante
# asignación exacta de la secuencia completa contra la base de datos de especies.
taxa16S <- addSpecies(taxa16S,
                       "01_RawData/silva/silva_species_assignment_v138.1.fa.gz")

saveRDS(taxa16S, file = "03_Results/rds/16S/taxa16S.RDS")

# Crear una copia sin nombres de fila (secuencias de ADN) para imprimir limpio.
taxa16S.print <- taxa16S
rownames(taxa16S.print) <- NULL  # NULL elimina los nombres de fila
head(taxa16S.print)

# ── ITS — base de datos UNITE ────────────────────────────────────────────────
# Descargar desde https://unite.ut.ee/repository.php
# multithread = FALSE porque la base de UNITE es generalmente más pequeña.
taxa_ITS <- assignTaxonomy(seqtab.nochimITS,
                             "01_RawData/unite_general_release.fasta.gz",
                             multithread = FALSE)
saveRDS(taxa_ITS, file = "03_Results/rds/ITS/taxaITS.RDS")


# ─────────────────────────────────────────────────────────────────────────────
# 14. ASIGNACIÓN TAXONÓMICA CON DECIPHER (IdTaxa)
# ─────────────────────────────────────────────────────────────────────────────

# ── 14a. 16S — conjunto de entrenamiento SILVA ────────────────────────────────

# getSequences() extrae las secuencias de los ASVs de la tabla seqtab.
# DNAStringSet() convierte el vector de secuencias en un objeto de Biostrings,
# que es el formato requerido por IdTaxa().
dna16S <- DNAStringSet(getSequences(seqtab.nochim16S))

# load() carga un archivo .RData que contiene el objeto "trainingSet" de SILVA.
# Este objeto es el clasificador pre-entrenado de DECIPHER para 16S.
load("01_RawData/silva/SILVA_SSU_r138_2019.RData")

# IdTaxa() clasifica las secuencias usando el enfoque probabilístico de DECIPHER,
# que generalmente es más preciso que el clasificador bayesiano a nivel de género.
# strand = "top" analiza solo la cadena tal como está (sin probar el reverso-complemento).
# processors = NULL usa todos los núcleos disponibles.
ids16S <- IdTaxa(dna16S, trainingSet,
                 strand = "top", processors = NULL, verbose = FALSE)

# Vector de rangos taxonómicos esperados en el resultado de IdTaxa.
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")

# IdTaxa devuelve una lista; este sapply() la convierte en una matriz:
# match(ranks, x$rank) alinea los rangos del resultado con el orden esperado.
# x$taxon[m] extrae el nombre del taxón en cada rango en el orden correcto.
# Las asignaciones "unclassified_*" se reemplazan por NA para consistencia.
# t() transpone la matriz para que quede como filas = ASVs, columnas = rangos.
taxid16S <- t(sapply(ids16S, function(x) {
  m    <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))

# Nombrar columnas (rangos) y filas (secuencias de ASV) de la matriz resultante.
colnames(taxid16S) <- ranks
rownames(taxid16S) <- getSequences(seqtab.nochim16S)

saveRDS(taxid16S, file = "03_Results/rds/16S/taxid16S-decipher.RDS")
taxid16S <- readRDS("03_Results/rds/16S/taxid16S-decipher.RDS")

# ── 14b. ITS — conjunto de entrenamiento UNITE ────────────────────────────────

# Mismo procedimiento que para el 16S, usando el clasificador UNITE para hongos.
dnaITS <- DNAStringSet(getSequences(seqtab.nochimITS))

# Carga el objeto trainingSet de UNITE v2025 (reemplaza el trainingSet de SILVA).
load("01_RawData/silva/UNITE_v2025.RData")

idsITS <- IdTaxa(dnaITS, trainingSet,
                 strand = "top", processors = NULL, verbose = FALSE)

# Misma conversión de lista a matriz que para el 16S.
taxidITS <- t(sapply(idsITS, function(x) {
  m    <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))

colnames(taxidITS) <- ranks
rownames(taxidITS) <- getSequences(seqtab.nochimITS)

saveRDS(taxidITS, file = "03_Results/rds/ITS/taxidITS_decipher.RDS")