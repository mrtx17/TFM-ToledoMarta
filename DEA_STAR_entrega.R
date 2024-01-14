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
setwd("/home/usr/MT/data/Quant_STAR")  # Cambia el directorio de trabajo

# Cogemos una lista de todos los directorios
files <- list.files(path="/Users/martatoledocastillo/Creative Cloud Files/Desktop/Quant/", pattern="txt")

# Unimos dos archivos y generamos una matriz
file1 <- read.table(files[1], col.names=c("transcript",files[1]))
file2 <- read.table(files[2], col.names=c("transcript",files[2]))
count_data <- merge (file1, file2, by=c("transcript"))

# Con un loop unimos todos los archivos en la matriz
for(i in 3:length(files))
{
  file <- read.table(files[i], col.names = c("transcript", sub("\\.txt", "", files[i])))
  count_data <- merge(count_data, file, by=c("transcript"))
}

# Se hacen las modificaciones necesarias para que sea legible
colnames(count_data) <- gsub("\\.txt", "", colnames(count_data))
count_data <- tail(count_data, -5)
row.names(count_data) <- count_data$transcript
count_data$transcript <- NULL

# Cargamos la información de las muestras
fenotipos = read.csv2("home/usr/MT/data/Quant_STAR/labels.csv", header = TRUE)
rownames(fenotipos) <- fenotipos$sample
fenotipos$sample <- NULL
fenotipos <- fenotipos %>% arrange(fenotipos$response)
orden<- rownames(fenotipos)
# Convierte los nombres de columna en índices de columna
count_data <- count_data[, orden, drop=FALSE]


# Comprobación: las columnas de la información de las muestras coincide con las filas de la matriz de cuantificación
if ( all(colnames(count_data) %in% rownames(fenotipos))) {
  cat("Todas las muestras están presentes tanto en los datos de cuantificación
      como de subtipo")
} else {
  cat("Falta alguna muestra en alguno de los dos data frames")
}

if(all(colnames(count_data) == rownames(fenotipos))) {
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
dds_obj <- DESeqDataSetFromMatrix(countData = count_data, colData = fenotipos,
                                  design = ~ phenotype)
ddsLRT <- DESeq(dds_obj, test="LRT", reduced = ~ 1)

########################### Heatmap de subtipos ################################
# Se obtienen los resultados
res <- results(ddsLRT)
# Se depuran los datos suprimiendo aquellos que sean NA y aquellos que no concuerden
# con ciertos criterios de p ajustados, Fold Change y la media de bases
sigs <- na.omit(res)
sigs<- sigs[which(sigs$padj <= 0.05 & abs(sigs$log2FoldChange) > 1) & sigs$baseMean > 50, ]
sigs.df <- as.data.frame(sigs)
# Obtenemos el símbolo de cada tránscrito para mejor interpretación
sigs.df$symbol <- mapIds(org.Hs.eg.db, keys = rownames(sigs.df), keytype = "ENSEMBL",
                         column = "SYMBOL")
# Ordenamos los datos por Fold Change de forma que será útil en siguientes pasos
df.top <- sigs.df[order(sigs.df$log2FoldChange, decreasing = TRUE),]

# Normalizamos los datos y generamos una matriz
ntd <- normTransform(ddsLRT)
mat <- assay(ntd)[rownames(df.top), rownames(fenotipos)]
symbol <- ifelse(is.na(df.top$symbol), rownames(df.top), df.top$symbol)
rownames(mat) <- symbol
mat_n <- na.omit(rownames(mat))
baseMean <- rowMeans(mat)
# Creamos una matriz con z-score, estandarizando la expresión génica
mat.scaled <- t(apply(mat, 1, scale))
colnames(mat.scaled) <- fenotipos$phenotype

# Creamos un vector que mantiene las 20 primeras y últimas filas (con fold mayor
# y manor)
rows_keep <- c(seq(1:20), seq((nrow(mat.scaled)-19), nrow(mat.scaled)))

# Generación de un primer Heatmap agrupando por subtipo
png(file= "home/usr/MT/data/Quant_STAR/RNASEQ-TFM/file_heatmap_nocluster.png", res=120, width=9000, height = 6000)
ha <- HeatmapAnnotation(summary = anno_summary(gp = gpar(fill =2),
                                               height = unit(10, "cm")))
h1 <- Heatmap(mat.scaled[rows_keep, ], cluster_rows = F,
              column_labels = colnames(mat.scaled),
              name="Z-score",
              cluster_columns = F,
              column_names_gp = grid::gpar(fontsize = 40),
              row_names_gp = grid::gpar(fontsize = 30),
heatmap_legend_param = list(labels_gp = gpar(fontsize = 50), title_gp = gpar(fontsize = 50)))
draw(h1)
dev.off()

# Generación de un segundo Heatmap agrupando por expresión
pdf(file= "home/usr/MT/data/Quant_STAR/RNASEQ-TFM/file_heatmap_cluster.pdf", width = 10)
ha <- HeatmapAnnotation(summary = anno_summary(gp = gpar(fill =2),
                                               height = unit(10, "cm")))
h2 <- Heatmap(mat.scaled[rows_keep, ], cluster_rows = F,
              column_labels = colnames(mat.scaled),
              name="Z-score",
              cluster_columns = T,
              column_names_gp = grid::gpar(fontsize = 4),
              row_names_gp = grid::gpar(fontsize = 3),
heatmap_legend_param = list(labels_gp = gpar(fontsize = 5), title_gp = gpar(fontsize = 5)))
h2
draw(h2)
dev.off()

################################################################################
############################ EXPLORACIÓN RESULTADOS  ###########################
################################################################################
res_LRT_tb <- res %>%
  data.frame() %>%
  rownames_to_column(var="transcript") %>%
  as_tibble()
sigsLRT_gene <- res_LRT_tb %>% filter(padj < 0.05) %>%
  filter(abs(log2FoldChange) >= 1)

clustering_sig_genes <- sigsLRT_gene %>% arrange(log2FoldChange)

ntd <- normTransform(ddsLRT)
mat <- assay(ntd)
clustering_rlog <- mat[clustering_sig_genes$transcript, ]

#Agrupación de los genes por expresión y ordenador por log2FoldChange para ver su expresión en cada subtipo
clusters <- degPatterns(clustering_rlog, metadata = fenotipos,
                        time ="phenotype", col=NULL)
# Generación listado de datos de expresión
expression_data <- clusters$normalized[, c(1,3,4)]
clustering_sig_genes$symbol <- mapIds(org.Hs.eg.db,
                                   keys = clustering_sig_genes$transcript,
                                   keytype = "ENSEMBL",
                                   column = "SYMBOL")
expression_data$symbol <- mapIds(org.Hs.eg.db,
                                      keys = expression_data$genes,
                                      keytype = "ENSEMBL",
                                      column = "SYMBOL")
expression_data <- expression_data[order(expression_data$value,
                                         decreasing = TRUE), ]
expression_data <- na.omit(expression_data)


################################################################################
####################### EXPLORACIÓN RESULTADOS 1 VS 1 ##########################
################################################################################

# TN vs LumB
res.TN.LumB <- results(ddsLRT, contrast=c("phenotype", "TN", "LumB"), alpha = 0.05)
res.TN.LumB <- res.TN.LumB[order(res.TN.LumB$padj,
                                 abs(res.TN.LumB$log2FoldChange)),]
summary(res.TN.LumB)
sigs.TN.LumB <- na.omit(res.TN.LumB)
sigs.TN.LumB<- sigs.TN.LumB[which(sigs.TN.LumB$padj <= 0.05 &
                            abs(sigs.TN.LumB$log2FoldChange) > 1), ]
sigs.df.TN.LumB <- as.data.frame(sigs.TN.LumB)

sigs.df.TN.LumB$symbol <- mapIds(org.Hs.eg.db, keys = rownames(sigs.df.TN.LumB),
                                 keytype = "ENSEMBL",
                                 column = "SYMBOL")


df.top.TN.LumB <- sigs.df.TN.LumB[order(sigs.df.TN.LumB$log2FoldChange,
                                        decreasing = TRUE),]
df.top.TN.LumB <- na.omit(df.top.TN.LumB)
write.csv(df.top.TN.LumB, "home/usr/MT/data/Quant_STAR/RNASEQ-TFM/df.top.TN.LumB.csv", row.names=FALSE)

# TN vs HER2
res.TN.HER2 <- results(ddsLRT, contrast=c("phenotype", "TN", "HER2"), alpha = 0.05)
res.TN.HER2 <- res.TN.HER2[order(res.TN.HER2$padj,
                                 abs(res.TN.HER2$log2FoldChange)),]
summary(res.TN.HER2)
sigs.TN.HER2 <- na.omit(res.TN.HER2)
sigs.TN.HER2<- sigs.TN.HER2[which(sigs.TN.HER2$padj <= 0.05 &
                                    abs(sigs.TN.HER2$log2FoldChange) > 1), ]
sigs.df.TN.HER2 <- as.data.frame(sigs.TN.HER2)

sigs.df.TN.HER2$symbol <- mapIds(org.Hs.eg.db, keys = rownames(sigs.df.TN.HER2),
                                 keytype = "ENSEMBL",
                                 column = "SYMBOL")


df.top.TN.HER2 <- sigs.df.TN.HER2[order(sigs.df.TN.HER2$log2FoldChange,
                                        decreasing = TRUE),]
df.top.TN.HER2 <- na.omit(df.top.TN.HER2)
write.csv(df.top.TN.LumB, "home/usr/MT/data/Quant_STAR/RNASEQ-TFM/df.top.TN.LumB.csv", row.names=FALSE)


#HER2 vs LumB
res.LumB.HER2 <- results(ddsLRT, contrast=c("phenotype", "HER2", "LumB"), alpha = 0.05)
res.LumB.HER2 <- res.LumB.HER2[order(res.LumB.HER2$pvalue,
                                     abs(res.LumB.HER2$log2FoldChange)),]
summary(res.LumB.HER2)

sigs.LumB.HER2 <- na.omit(res.LumB.HER2)
sigs.LumB.HER2<- sigs.LumB.HER2[which(sigs.LumB.HER2$padj <= 0.05 &
                                    abs(sigs.LumB.HER2$log2FoldChange) > 1), ]
sigs.df.LumB.HER2 <- as.data.frame(sigs.LumB.HER2)

sigs.df.LumB.HER2$symbol <- mapIds(org.Hs.eg.db,
                                   keys = rownames(sigs.df.LumB.HER2),
                                   keytype = "ENSEMBL",
                                   column = "SYMBOL")


df.top.LumB.HER2 <- sigs.df.LumB.HER2[order(sigs.df.LumB.HER2$log2FoldChange,
                                            decreasing = TRUE),]
df.top.LumB.HER2 <- na.omit(df.top.LumB.HER2)
write.csv(df.top.TN.LumB, "home/usr/MT/data/Quant_STAR/RNASEQ-TFM/df.top.LumB.HER2", row.names=FALSE)




################################################################################
################################################################################
#################### DEA: Respondedoras vs. No Respondedoras ###################
################################################################################
################################################################################

################### ExploraciÓn resultados: 3 subtipos #########################

######################## Generamos objeto DESeq2 ###############################
dds_res <- DESeqDataSetFromMatrix(countData = count_data, colData = fenotipos, design = ~ response)
dds_res <- DESeq(dds_res)

############################ Obtención resultado ###############################
res_res <- results(dds_res, contrast = c("response", "no", "si"), alpha = 0.05)
sigs_res <- na.omit(res_res)
sigs_res<- sigs_res[which(sigs_res$padj <= 0.05 & abs(sigs_res$log2FoldChange) > 1) & sigs_res$baseMean > 50, ]
sigs_res.df <- as.data.frame(sigs_res)
sigs_res.df$symbol <- mapIds(org.Hs.eg.db, keys = rownames(sigs_res.df), keytype = "ENSEMBL",
                         column = "SYMBOL")
df.top_res <- sigs_res.df[order(sigs_res.df$log2FoldChange, decreasing = TRUE),]

############################### Generamos Heatmaps #############################
ntd_res <- normTransform(dds_res)
mat_res <- assay(ntd_res)[rownames(df.top_res), rownames(fenotipos)]
symbol_res <- ifelse(is.na(df.top_res$symbol), rownames(df.top_res), df.top_res$symbol)
rownames(mat_res) <- symbol_res
mat_n_res <- na.omit(rownames(mat_res))
baseMean <- rowMeans(mat_res)
mat.scaled_r <- t(apply(mat_res, 1, scale))
colnames(mat.scaled_r) <- fenotipos$response
# Creamos un vector que mantiene las 20 primeras y últimas filas (con fold mayor
# y manor)
rows_keep_r <- c(seq(1:20), seq((nrow(mat.scaled_r)-19), nrow(mat.scaled_r)))
# Heatmap de la expresión ordenado por respuesta
pdf(file= "home/usr/MT/data/Quant_STAR/RNASEQ-TFM/response_heatmap_nocluster.pdf", width = 10, height = 6)
hr1 <- Heatmap(mat.scaled_r[rows_keep_r, ], cluster_rows = F,
              column_labels = colnames(mat.scaled_r),
              name="Z-score",
              cluster_columns = F,
              column_names_gp = grid::gpar(fontsize = 6),
              row_names_gp = grid::gpar(fontsize = 6))
draw(hr1)
dev.off()
# Heatmap de la expresión ordenado por expresión
pdf(file= "home/usr/MT/data/Quant_STAR/RNASEQ-TFM/RNASEQ-TFM/response_heatmap_yescluster.pdf", width = 10, height = 6)
hr2 <- Heatmap(mat.scaled_r[rows_keep_r, ], cluster_rows = F,
              column_labels = colnames(mat.scaled_r),
              name="Z-score",
              cluster_columns = T,
              column_names_gp = grid::gpar(fontsize = 6),
              row_names_gp = grid::gpar(fontsize = 6))
draw(hr2)
dev.off()


################################################################################
###########################       HER2        ##################################
################################################################################
# Tomamos la matriz solo de HER2
fenotipos <- fenotipos %>% arrange(fenotipos$phenotype)
orden<- rownames(fenotipos)
count_data <- count_data[, orden, drop = FALSE]
fenotipo_her <- fenotipos %>% filter(fenotipos$phenotype == 'HER2')
fenotipo_her
orden<- rownames(fenotipo_her)
count_data_HER2 <- count_data[, orden, drop = TRUE]
count_data_HER2
# Creamos el objeto DESeq2
dds_her <- DESeqDataSetFromMatrix(countData = count_data_HER2,
                                  colData = fenotipo_her, design = ~ response)
dds_her <- DESeq(dds_her)
# Exploramos resultados
res_her <- results(dds_her)
summary(res_her)
sigs_her <- na.omit(res_her)
sigs_her<- sigs_her[which(sigs_her$padj <= 0.05 & abs(sigs_her$log2FoldChange) > 1), ]
sigs_her <- as.data.frame(sigs_her)
sigs_her$symbol <- mapIds(org.Hs.eg.db, keys = rownames(sigs_her),
                         keytype = "ENSEMBL",
                         column = "SYMBOL")
df.top_her <- sigs_her[order(sigs_her$log2FoldChange, decreasing = TRUE),]



################################################################################
###############################       TN        ################################
################################################################################
# Tomamos la matriz solo de TN
fenotipos <- fenotipos %>% arrange(fenotipos$phenotype)
orden<- rownames(fenotipos)
count_data <- count_data[, orden, drop = FALSE]
fenotipo_tn <- fenotipos %>% filter(fenotipos$phenotype == 'TN')
fenotipo_tn
orden<- rownames(fenotipo_tn)
count_data_TN <- count_data[, orden, drop = TRUE]
count_data_TN
# Creamos el objeto DESeq2
dds_tn <- DESeqDataSetFromMatrix(countData = count_data_TN,
                                  colData = fenotipo_tn, design = ~ response)
dds_tn <- DESeq(dds_tn)
# Exploramos resultados
res_tn <- results(dds_tn)
summary(res_tn)
sigs_tn <- na.omit(res_tn)
sigs_tn<- sigs_tn[which(sigs_tn$padj <= 0.05 & abs(sigs_tn$log2FoldChange) > 1), ]
sigs_tn <- as.data.frame(sigs_tn)
sigs_tn$symbol <- mapIds(org.Hs.eg.db, keys = rownames(sigs_tn),
                         keytype = "ENSEMBL",
                         column = "SYMBOL")
df.top_tn <- sigs_tn[order(sigs_tn$log2FoldChange, decreasing = TRUE),]


################################################################################
#############################       LumB       #################################
################################################################################
# Tomamos la matriz solo de LumB
fenotipos <- fenotipos %>% arrange(fenotipos$phenotype)
orden<- rownames(fenotipos)
count_data <- count_data[, orden, drop = FALSE]
fenotipo_LumB <- fenotipos %>% filter(fenotipos$phenotype == 'LumB')
fenotipo_LumB
orden_LumB <- rownames(fenotipo_LumB)
count_data_LumB <- count_data[, orden_LumB, drop = TRUE]
count_data_LumB
# Creamos el objeto DESeq2
dds_LumB <- DESeqDataSetFromMatrix(countData = count_data_LumB,
                                   colData = fenotipo_LumB,
                                   design = ~ response)
dds_LumB <- DESeq(dds_LumB)
# Exploramos resultados
res_LumB <- results(dds_LumB)
summary(res_LumB)
sigs_LumB <- na.omit(res_LumB)
sigs_LumB<- sigs_LumB[which(sigs_LumB$padj <= 0.05 & abs(sigs_LumB$log2FoldChange) > 1), ]
sigs_LumB <- as.data.frame(sigs_LumB)
sigs_LumB$symbol <- mapIds(org.Hs.eg.db, keys = rownames(sigs_LumB),
                          keytype = "ENSEMBL",
                          column = "SYMBOL")
df.top_LumB <- sigs_LumB[order(sigs_LumB$log2FoldChange, decreasing = TRUE),]
