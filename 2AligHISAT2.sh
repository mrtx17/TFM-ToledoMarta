#!/bin/bash

#Directorio donde están los datos pre-procesados
home_dir="/home/usr/MT/data/processed"
#Directorio donde están los índices de HISAT2
HISAT_dir="/home/usr/resources/hisat2/genome_snp_tran"
#Directorio final de los datos alineados por HISAT2
alig_HISAT="/home/usr/MT/data/alignedHISAT"

cd ${home_dir}
echo "El directorio actual es $(pwd)"

for file in $(find *RNA_fwd.fq.gz -type f);
do
	SECONDS=0
	FICHERO="${file}"
	NOMBRE="${FICHERO%%.RNA*}"
	echo "La alineación de ${NOMBRE} con HISAT2 comienza"

	hisat2 -x ${HISAT_dir} --threads 24 \
	-1 ${NOMBRE}_RNA_fwd.fq.gz -2 ${NOMBRE}_RNA_rev.fq.gz --rna-stradness FR \
	-S ${alig_HISAT}/${NOMBRE}_Haligned.sam --summary-file ${alig_HISAT}/${NOMBRE}_summary.txt

	samtools sort ${alig_HISAT}/${NOMBRE}_Haligned.sam -o ${alig_HISAT}/${NOMBRE}_Haligned.bam
	samtools sort -n ${alig_HISAT}/${NOMBRE}_Haligned.bam -o${alig_HISAT}/${NOMBRE}_Haligned.sorted.bam

	rm ${alig_HISAT}/${NOMBRE}_Haligned.sam
	rm ${alig_HISAT}/${NOMBRE}_Haligned.bam

	duration=$SECONDS
	echo "$(($duration/60)) minutos y $((duration % 60)) segundos"

done
