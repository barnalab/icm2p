#!/bin/bash
#SBATCH -n 1
#SBATCH -t 12:0:0
#SBATCH --mem-per-cpu=16G
#SBATCH -o logs/%x.out
#SBATCH -e logs/%x.err

#Dependencies
#---
#cutdapt
#bowtie2
#bbmap
#samtools
#bamtools

source activate m2env

READ1=$1
READ2=$2
OUTPATH=$3
OUTPREFIX=$4

MINPHRED=30
INDEX=index/reference

mkdir -p $OUTPATH

#Trim illumina adapter sequences
cutadapt -j 1 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -A GATCGTCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT -o $OUTPATH/$OUTPREFIX"_R1_trimmed.fastq" -p $OUTPATH/$OUTPREFIX"_R2_trimmed.fastq" $READ1 $READ2
rm $READ1
rm $READ2

#Merge paired reads to one
bbmerge.sh in1=$OUTPATH/$OUTPREFIX"_R1_trimmed.fastq" in2=$OUTPATH/$OUTPREFIX"_R2_trimmed.fastq" out=$OUTPATH/$OUTPREFIX".fastq" outu=$OUTPATH/$OUTPREFIX".unmerged"
rm $OUTPATH/$OUTPREFIX"_R1_trimmed.fastq"
rm $OUTPATH/$OUTPREFIX"_R2_trimmed.fastq"

#Bowtie2
bowtie2 -p 1 --end-to-end --very-sensitive --mp 3,1 --rdg 5,1 --rfg 5,1 --dpad 30 --no-unal -x $INDEX -U $OUTPATH/$OUTPREFIX".fastq" | samtools view -bS - > $OUTPATH/$OUTPREFIX".bam"
rm $OUTPATH/$OUTPREFIX".fastq"
rm $OUTPATH/$OUTPREFIX".unmerged"

#Split to each reference
bamtools split -in $OUTPATH/$OUTPREFIX".bam" -reference -refPrefix ""
rm $OUTPATH/$OUTPREFIX".bam"

#For each reference sequence
for BAM in `ls $OUTPATH/$OUTPREFIX"".*.bam`
do
	REF=$(basename $BAM)
	REF=${REF%.bam}
	REF=${REF#$OUTPREFIX"."}

	samtools view $BAM > $OUTPATH/$OUTPREFIX"."$REF".sam"
	rm $BAM

	#ShapeMapper2 mutation parser
	./shapemapper_mutation_parser -i $OUTPATH/$OUTPREFIX"."$REF".sam" -o $OUTPATH/$OUTPREFIX"."$REF".mutParse" -m 0 --max_internal_match 0 --min_qual $MINPHRED --max_paired_fragment_length 10000

	rm $OUTPATH/$OUTPREFIX"."$REF".sam"
	mkdir -p $OUTPATH/$REF
   
	#Make correlated mutation matrix 
	python m2matrix.py --ref $INDEX".fa" --refname $REF --mutParse $OUTPATH/$OUTPREFIX"."$REF".mutParse" --outprefix $OUTPATH/$REF/$OUTPREFIX"."$REF
	rm $OUTPATH/$OUTPREFIX"."$REF".mutParse"
done
