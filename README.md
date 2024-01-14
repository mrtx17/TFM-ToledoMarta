En esta carpeta se encuentran los códigos utilizados durante el análisis transcriptómico de exosomas de biopsia líquida, realizado para el Trabajo de Fin de Máster de la UOC y la UB, y los datos clínicos de las pacientes. Tanto las muestras como los datos clínicos pertenecen a Atrys Health.

También se encuentran los datos obtenidos durante el análisis de expresión diferencial y el análisis de enriquecimiento de conjunto de genes, tanto entre los subtipos y las diferentes respuestas a la quimioterapia neoadyuvante.

A continuación se describe todo el contenido:


- 1Preprocessing.sh: script en bash a través del cual se hace el pre-procesamiento.
- 2AligHISAT2.sh: script en bash a través del cual se alinea con HISAT2.
- 2AligSTAR.sh: script en bash a través del cual se alinea con STAR.
- 3Quant.sh: script en bash que permite la cuantificación de las secuencias alineadas.
- DEA_HISAT2_entrega.R: script en R a través del cual se genera el análisis de expresión diferencial de los datos obtenidos por HISAT2.
- DEA_STAR_entrega.R: script en R a través del cual se genera el análisis de expresión diferencial de los datos obtenidos por STAR.
- DEG-LumBvsHER2.STAR.csv: tabla que recoge los genes expresados diferencialmente entre el subtipo luminal B y HER2+.
- DEG-TNvsHER2.STAR.csv: tabla que recoge los genes expresados diferencialmente entre el subtipo triple negativo y HER2+.
- DEG-TNvsLumB.STAR.csv: tabla que recoge los genes expresados diferencialmente entre el subtipo triple negativo y luminal B.
- Datos clínicos RNAseq.csv: tabla que recoge los datos clínicos y demográficos de las pacientes.
- GSEA-Novs.Si-STAR.csv: tabla que recoge los conjuntos de genes enriquecidos diferencialmente en respondedoras frente a no respondedoras.
- GSEA-TNvsHER2-STAR.csv: tabla que recoge los conjuntos de genes enriquecidos diferencialmente en subtipo TN frente a HER2+.
- GSEA-TNvsLumB-STAR.csv: tabla que recoge los conjuntos de genes enriquecidos diferencialmente en subtipo TN frente a luminal B.
- GSEA-TNvsREST-STAR.csv: tabla que recoge los conjuntos de genes enriquecidos diferencialmente en subtipo TN frente al resto.
- STAR-listado_genes_Z-score.csv: tabla que recoge los genes expresados diferencialmente entre los tres fenotipos, con el Z-score.
- df.top.response_STAR.csv: tabla que recoge los genes expresados diferencilamente entre respondedoras y no respondedoras en los tres subtipos.
- df.top_LumB.response_STAR.csv: tabla que recoge los genes expresados diferencilamente entre respondedoras y no respondedoras en el subtipo luminal B.
- df.top_her.response_STAR.csv: tabla que recoge los genes expresados diferencilamente entre respondedoras y no respondedoras en el subtipo HER2+.
