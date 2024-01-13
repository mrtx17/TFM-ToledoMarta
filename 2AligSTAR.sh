#!/bin/bash

#Directorio donde están los datos pre-procesados
home_dir="/home/usr/MT/data/processed"
#Directorio donde están los índices de STAR
STAR_dir="/home/usr/resources/STAR"
#Directorio final de los datos alineados por STAR
alig_STAR="/home/usr/MT/data/allignedSTAR"

cd ${home_dir}
echo "El directorio actual es $(pwd)"

for file in $(find *RNA_fwd.fq.gz -type f);
do
	SECONDS=0
	FICHERO="${file}"
	NOMBRE="${FICHERO%%.RNA*}"
	echo "La alineación de ${NOMBRE} con STAR comienza"

	STAR --runThreadN 24 --genomeDir ${STAR_dir} --readFilesCommand zcat \
	--readFilesIn ${NOMBRE}.RNA_fwd.fq.gz ${NOMBRE}.RNA_rev.fq.gz \
	--outFileNamePrefix ${alig_STAR}/${NOMBRE} \
	--outSAMtype BAM SortedByCoordinate 

	duration=$SECONDS
	echo "$(($duration/60)) minutos y $((duration % 60)) segundos"

	done
