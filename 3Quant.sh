#!/bin/bash

#Directorio donde están los datos pre-procesados
home_dir_STAR="/home/usr/MT/data/processed"
#Directorio donde están los datos pre-procesados
home_dir_HISAT2="/home/usr/MT/data/alignedHISAT"
#Directorio donde están los índices de STAR
gtf_file="/home/usr/resources/Homo_sapiens.GRCh38.110.chr.gtf"
#
fin_dir="/home/usr/MT/data/Quant"


cd ${home_dir_HISAT2}
echo "El directorio actual es $(pwd)"

for file in $(find *Haligned.sorted.bam -type f);
do
	SECONDS=0
	FICHERO="${file}"
	NOMBRE="${FICHERO%%_Haligned*}"
	echo "La cuantificación de ${NOMBRE} con HTSeq comienza"

	htseq-count ${NOMBRE}_Haligned.sorted.bam ${gtf_file} \
	--format bam --order name --mode union \
	--nonunique all -n 24 --idattr gene_id \
	--supplementary-alignments ignore --secondary-alignments score > ${fin_dir}/${NOMBRE}_quantHT.txt

	duration=$SECONDS
	echo "$(($duration/60)) minutos y $((duration % 60)) segundos"
done

cd ${home_dir_STAR}
echo "El directorio actual es $(pwd)"

for file in $(find *SortedByCoordinate.bam -type f);
do
		SECONDS=0
		FICHERO="${file}"
		NOMBRE="${FICHERO%%SortedByCoordinate*}"
		echo "La cuantificación de ${NOMBRE} con HTSeq comienza"

		htseq-count ${NOMBRE}_Haligned.sorted.bam ${gtf_file} \
		--format bam --order name --mode union \
		--nonunique all -n 24 --idattr gene_id \
		--supplementary-alignments ignore --secondary-alignments score > ${fin_dir}/${NOMBRE}_quantSTAR.txt

		duration=$SECONDS
		echo "$(($duration/60)) minutos y $((duration % 60)) segundos"
done
