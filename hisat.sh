#!/bin/bash
#ref_genome="/"
ref_gtf=/mnt/raid/illumina/AlexR/transcriptomes/splcing.detection.comparsion/data/refgenome/Mus_musculus.GRCm38.98.gtf
#extract_

echo "Mapping started!"
ls *.fq | parallel 'hisat2 -p 8 --dta -x /mnt/raid/illumina/AlexR/transcriptomes/splcing.detection.comparsion/data/refgenome/hs2-index/index {} -S {.}.sam'
echo "SAM-BAM conversion started"
ls *.sam | parallel 'samtools sort -@ 16 -o {.}.bam {}'
echo "Assembling transcripts"
ls *.bam | parallel 'stringtie -p 16 -G /mnt/raid/illumina/AlexR/transcriptomes/splcing.detection.comparsion/data/refgenome/Mus_musculus.GRCm38.98.gtf -o {.}.gtf -l {.} {}'
echo "Merging transcripts"
mkdir merged
stringtie --merge -p 8 -G /mnt/raid/illumina/AlexR/transcriptomes/splcing.detection.comparsion/data/refgenome/Mus_musculus.GRCm38.98.gtf -o merged/merged.gtf *.gtf
echo "Executing gffcompare!"
cd merged
gffcompare -r /mnt/raid/illumina/AlexR/transcriptomes/splcing.detection.comparsion/data/refgenome/Mus_musculus.GRCm38.98.gtf -o merged merged.gtf
cd ..
echo "Creating ballgown counts tables..."
mkdir ballgown
ls *.bam | parallel 'stringtie -e -B -p 16 -G merged/merged.gtf -o ballgown/{.}/{.}.gtf {.}.bam'
echo "Done!"
