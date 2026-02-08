#!/bin/bash

#######Each species was processed using the same workflow/pipeline.
###Genome assemble
hifiasm -o Species -t 40 -l 2 -n 4 \
./Species.hifi.fastq.gz \
2>hifi.assemble.log

###Genome annotation
##BRAKER
singularity exec -B ${PWD}:${PWD} ${BRAKER_SIF} braker.pl \
  --genome=Species_hifiasm_masked.fasta \
  --bam=Species_sort.bam \ 
  --workingdir=./  \
  --GENEMARK_PATH=gmes_linux_64_4 \
  --threads 40 \
  --gm_max_intergenic 10000 \
  --skipOptimize \
  --AUGUSTUS_CONFIG_PATH=augustus/config \
  --AUGUSTUS_BIN_PATH=augustus/bin \
  --AUGUSTUS_SCRIPTS_PATH=augustus/scripts \
  --busco_lineage eukaryota_odb10 &> braker.log

##Gemoma
##Homologous species were selected differently for each species to assist in annotation.

java -Xms300G -Xmx350G -jar GeMoMa-1.9.jar CLI GeMoMaPipeline \
        threads=32 \
        GeMoMa.Score=ReAlign \
        AnnotationFinalizer.r=NO o=true p=false \
        GAF.tf=false \
        t=Species_genome.fasta \
        s=own i=Human a=human-chm13.gff g=human-chm13.fa

##Stringtie

hisat2-build Species_genome.fasta Species_genome

hisat2 -p 20 --dta -x Species_genome -1 ../Species_1_clean.fq.gz -2 ../Species_2_clean.fq.gz -S ./Species.sam

samtools view  -@ 20 -b Species.sam >Species.bam

samtools sort  -@ 20 Species.bam -o Species_sort.bam

samtools index Species_sort.bam

samtools flagstat Species_sort.bam >stat.txt

stringtie -p 20 -o stringtie.gtf Species_sort.bam

mkdir transdecoder
cd transdecoder
gtf_genome_to_cdna_fasta.pl ../stringtie.gtf Species_hifiasm_masked.fasta > transcripts.fasta
gtf_to_alignment_gff3.pl ../stringtie.gtf > transcripts.gff3
TransDecoder.LongOrfs -t transcripts.fasta
TransDecoder.Predict -t transcripts.fasta
cdna_alignment_orf_to_genome_orf.pl transcripts.fasta.transdecoder.gff3 transcripts.gff3 transcripts.fasta > transcripts.fasta.transdecoder.genome.gff3

###EVM
EVidenceModeler \
--sample_id Species \
--genome Species_hifiasm_masked.fasta \
--gene_predictions predict.gff3 \
--protein_alignments protein.gff3 \
--segmentSize 100000 \
--overlapSize 10000 \
--weights weight.txt