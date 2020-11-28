#!/bin/bash

OUTPATH=$1
OUTPREFIX=$2

for INPREFIX in "${@:3}"
do
    READ1=$INPREFIX"_R1.fastq.gz"
    READ2=$INPREFIX"_R2.fastq.gz"

    mkdir -p $OUTPATH

    zcat $READ1 | split -l 4000000 --filter='cat > $FILE.fastq' - $OUTPATH/$OUTPREFIX"_R1."
    zcat $READ2 | split -l 4000000 --filter='cat > $FILE.fastq' - $OUTPATH/$OUTPREFIX"_R2."
done

for READ in `ls $OUTPATH/*R1.*.fastq`
do
    SPLIT=$(basename $READ)
    SPLIT=${SPLIT%.fastq}
    SPLIT=${SPLIT#$OUTPREFIX"_R1."}

    READ1=$OUTPATH/$OUTPREFIX"_R1."$SPLIT".fastq"
    READ2=$OUTPATH/$OUTPREFIX"_R2."$SPLIT".fastq"

	mkdir -p logs
	p0=$(sbatch --parsable --job-name=$OUTPREFIX"_"$SPLIT ./p0.sh $READ1 $READ2 $OUTPATH $OUTPREFIX"_"$SPLIT);
done
