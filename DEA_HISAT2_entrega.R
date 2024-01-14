library(DESeq2)
library(pheatmap)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)
library(VennDiagram)
library(ComplexHeatmap)
library(org.Hs.eg.db)
library(circlize)
library(grid)
library(tidyverse)
library(msigdbr)
library(clusterProfiler)
library(fgsea)
library(magrittr)
library(knitr)

################################################################################
############ Importación de los datos y preparación de la matriz ###############
################################################################################
setwd("/home/usr/MT/data/Quant_HISAT2")  # Cambia el directorio de trabajo si es necesario

# Cogemos una lista de todos los directorios
filesH = list.files(path="/Users/martatoledocastillo/Creative Cloud Files/Desktop/QuantHISAT/", pattern="txt")

# Unimos dos archivos y generamos una matriz
file1H = read.table(filesH[1], col.names=c("transcript",files[1]))
file2H = read.table(filesH[2], col.names=c("transcript",files[2]))
count_data_h = merge (file1H, file2H, by=c("transcript"))

# Con un loop unimos todos los archivos en la matriz
for(i in 3:length(filesH))
{
  file = read.table(filesH[i], col.names = c("transcript", sub("\\.txt", "", filesH[i])))
  count_data_h <- merge(count_data_h, file, by=c("transcript"))
}

# Se hacen las modificaciones necesarias para que sea legible
colnames(count_data_h) <- gsub("\\.txt", "", colnames(count_data_h))
count_data_h <- tail(count_data_h, -5)
row.names(count_data_h) <- count_data_h$transcript
count_data_h$transcript <- NULL

# Cargamos la información de las muestras
fenotipos = read.csv2("home/usr/MT/data/Quant_HISAT2/labels.csv", header = TRUE)
rownames(fenotipos) <- fenotipos$sample
fenotipos$sample <- NULL
fenotipos <- fenotipos %>% arrange(fenotipos$phenotype)
orden<- rownames(fenotipos)
# Convierte los nombres de columna en índices de columna
count_data_h <- count_data_h[, orden, drop=FALSE]


# Comprobación: las columnas de la información de las muestras coincide con las filas de la matriz de cuantificación
if ( all(colnames(count_data_h) %in% rownames(fenotipos))) {
  cat("Todas las muestras están presentes tanto en los datos de cuantificación
      como de subtipo")
} else {
  cat("Falta alguna muestra en alguno de los dos data frames")
}

if(all(colnames(count_data_h) == rownames(fenotipos))) {
  cat("Los datos de cuantificación y de fenotipo siguen el mismo orden")
} else {
  cat("Las muestras no siguen el mismo orden en los datos de cuantificación que
      en los datos de subtipos")
}


################################################################################
################################################################################
######################## EXPLORACIÓN DATOS SUBTIPOS ############################
################################################################################
################################################################################

################# Creamos el objeto DESeq para estudiar subtipos ###############
dds_h <- DESeqDataSetFromMatrix(countData = count_data_h, colData = fenotipos, design = ~ phenotype)
# Indicamos que no hay un valor de referencia
dds_h <- DESeq(dds_h, test="LRT", reduced = ~ 1)

########################### Exploración 1 vs 1 #################################

res.TN.LumB_h <- results(dds_h, contrast=c("phenotype", "TN", "LumB"), alpha = 0.05)
res.TN.LumB_h <- res.TN.LumB_h[order(res.TN.LumB_h$padj, abs(res.TN.LumB_h$log2FoldChange)),]
summary(res.TN.LumB_h)
sigs.TN.LumB_h <- na.omit(res.TN.LumB_h)
sigs.TN.LumB_h <- sigs.TN.LumB_h[which(sigs.TN.LumB_h$padj <= 0.05 &
                                    abs(sigs.TN.LumB_h$log2FoldChange) > 1) &
                              sigs.TN.LumB_h$baseMean > 50, ]
sigs.df.TN.LumB_h <- as.data.frame(sigs.TN.LumB_h)

sigs.df.TN.LumB_h$symbol <- mapIds(org.Hs.eg.db, keys = rownames(sigs.df.TN.LumB_h),
                                 keytype = "ENSEMBL",
                                 column = "SYMBOL")

res.TN.HER2_h <- results(dds_h, contrast=c("phenotype", "TN", "HER2"), alpha = 0.05)
res.TN.HER2_h <- res.TN.HER2_h[order(res.TN.HER2_h$padj, abs(res.TN.HER2_h$log2FoldChange)),]
summary(res.TN.HER2_h)
sigs.TN.HER2_h <- na.omit(res.TN.HER2_h)
sigs.TN.HER2_h <- sigs.TN.HER2_h[which(sigs.TN.HER2_h$padj <= 0.05 &
                                    abs(sigs.TN.HER2_h$log2FoldChange) > 1) &
                              sigs.TN.HER2_h$baseMean > 50, ]
sigs.df.TN.HER2_h <- as.data.frame(sigs.TN.HER2_h)

sigs.df.TN.HER2_h$symbol <- mapIds(org.Hs.eg.db, keys = rownames(sigs.df.TN.HER2_h),
                                 keytype = "ENSEMBL",
                                 column = "SYMBOL")

res.LumB.HER2_h <- results(dds_h, contrast=c("phenotype", "HER2", "LumB"), alpha = 0.05)
res.LumB.HER2_h <- res.LumB.HER2_h[order(res.LumB.HER2_h$pvalue, abs(res.LumB.HER2_h$log2FoldChange)),]
summary(res.LumB.HER2_h)
sigs.LumB.HER2_h <- na.omit(res.LumB.HER2_h)
sigs.LumB.HER2_h <- sigs.LumB.HER2_h[which(sigs.LumB.HER2_h$padj <= 0.05 &
                                    abs(sigs.LumB.HER2_h$log2FoldChange) > 1) &
                              sigs.LumB.HER2_h$baseMean > 50, ]
sigs.df.LumB.HER2_h <- as.data.frame(sigs.LumB.HER2_h)

sigs.df.LumB.HER2_h$symbol <- mapIds(org.Hs.eg.db, keys = rownames(sigs.df.LumB.HER2_h),
                                 keytype = "ENSEMBL",
                                 column = "SYMBOL")

res_h <- results(dds_h,alpha = 0.05)


################################################################################
################################################################################
######################## EXPLORACIÓN DATOS RESPUESTA ###########################
################################################################################
################################################################################

################### Exploración resultados: 3 subtipos #########################

######################## Generamos objeto DESeq2 ###############################
dds_res_h <- DESeqDataSetFromMatrix(countData = count_data_h, colData = fenotipos, design = ~ response)
dds_res_h <- DESeq(dds_res_h)

############################ Obtención resultado ###############################
res_res_h <- results(dds_res_h, contrast = c("response", "no", "si"), alpha = 0.05)
sigs_res_h <- na.omit(res_res_h)
sigs_res_h<- sigs_res[which(sigs_res_h$padj <= 0.05 & abs(sigs_res_h$log2FoldChange) > 1) & sigs_res_h$baseMean > 50, ]
sigs_resh.df <- as.data.frame(sigs_res_h)
sigs_resh.df$symbol <- mapIds(org.Hs.eg.db, keys = rownames(sigs_resh.df), keytype = "ENSEMBL",
                         column = "SYMBOL")
df.top_res <- sigs_resh.df[order(sigs_resh.df$log2FoldChange, decreasing = TRUE),]
