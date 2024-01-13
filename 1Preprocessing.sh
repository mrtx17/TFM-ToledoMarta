#!/bin/bash

# Definición de directorios
home="/home/usr/MT/data/" # Directorio de los datos en bruto
bbdukres="/usr/local/bin/bbmap/resources" # Directorio datos referencia adaptadores
bbduk="/usr/local/bin/bbmap/" # Directorio para llegar a bbduk
sortmerna="/home/usr/apps/sortmernadb" # Directorio para referencia de rRNA
fin_dir="/home/usr/MT/data/processed" # directorio final

# Cambio al directorio principal
cd ${home}
echo "El directorio actual es $(pwd)"

# Iteración sobre archivos *_merge_R1.fastq.gz
for file in $(find *_merge_R1.fastq.gz -type f);
do
	SECONDS=0
	FICHERO="${file}"
	NOMBRE="${FICHERO%%_merge_R1*}"
	echo "${NOMBRE} se procesará ahora"

	## Análisis de calidad de las secuencias
	echo "El análisis de calidad comienza"
	fastqc -t 6 --noextract \
	${NOMBRE}_merge_R1.fastq ${NOMBRE}_merge_R2.fastq

	## Eliminación de adaptadores con BBDuk
	echo "El procesamiento por BBDuk comienza"
	${bbduk}/bbduk.sh in1=${NOMBRE}_merge_R1.fastq \
	in2=${NOMBRE}_merge_R2.fastq \
	ref=${bbdukres}/adapters.fa \
	out1=clean1.${NOMBRE}.fastq \
	out2=clean2.${NOMBRE}.fastq \
	tbo tpe mink=20 ktrim=r

	## Análisis de calidad post-BBDuk
	echo "El análisis de calidad comienza"
	fastqc -t 6 --noextract \
	clean1.${NOMBRE}.fastq clean2.${NOMBRE}.fastq

	## Eliminación de secuencias con calidad <5 y longitud < 20 con CutAdapt
	echo "El procesamiento por CutAdapt comienza"
	cutadapt -q 5,5 --pair-filter=any -m 20 \
	--output=out.${NOMBRE}.1_Q5.fastq \
	--paired-output=out.${NOMBRE}.2_Q5.fastq \
	clean1.${NOMBRE}.fastq \
	clean2.${NOMBRE}.fastq

	## Análisis de calidad post-CutAdapt
	echo "El análisis de calidad comienza"
	fastqc -t 6 --noextract \
	out.${NOMBRE}.1_Q5.fastq out.${NOMBRE}.2_Q5.fastq

	## Pre-procesamiento: eliminación de rRNA con SortMeRNA
	echo "El procesamiento de RNA por sortmeRNA comienza"
	sortmerna --workdir ${sortmerna} --ref ${sortmerna}/smr_v4.3_fast_db.fasta \
	--reads out.${NOMBRE}.1_Q5.fastq.gz --reads out.${NOMBRE}.2_Q5.fastq.gz \
	--fastx -out2 --aligned ${home}/${NOMBRE}.rRNA --other ${fin_dir}/${NOMBRE}.RNA --threads 24

	## Eliminación de archivos temporales
	rm -rf ${sortmerna}/kvdb/

	## Análisis de las secuencias post-sortmeRNA
	echo "El análisis de calidad comienza"
	cd ${fin_dir}
	fastqc -t 6 --noextract \
	${NOMBRE}.RNA_fwd.fq.gz ${NOMBRE}.RNA_rev.fq.gz

	## Generación de informes con MultiQC
	multiqc .

	## Calcular y mostrar el tiempo de ejecución
	duration=$SECONDS
	echo "$(($duration/60)) minutos y $((duration % 60)) segundos"

done
