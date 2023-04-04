#!/bin/bash

# set reference genome file path
ref_genome="NC_001477.fasta"

# gunzip all fastq.gz files in folder
for gz_file in *.fastq.gz
do
    gunzip "$gz_file"
done

# loop through all fastq files in folder
for fastq_file in *.fastq
do
    # generate output SAM file name based on input fastq file name
    sam_file="${fastq_file%.fastq}_minimap2.sam"

    # run minimap2 to align fastq file to reference genome and output SAM file
    minimap2 -a "$ref_genome" "$fastq_file" > "$sam_file"

    # generate output BAM file name based on input SAM file name
    bam_file="${sam_file%.sam}_sorted.bam"

    # sort SAM file and output as BAM file
    samtools sort "$sam_file" -o "$bam_file"

    # generate index file for BAM file
    samtools index "$bam_file"

    # generate output consensus fasta file name based on input BAM file name
    consensus_file="${bam_file%.bam}.consensus.fasta"

    # generate consensus fasta file from BAM file
    samtools consensus -f fasta "$bam_file" -o "$consensus_file"
done
