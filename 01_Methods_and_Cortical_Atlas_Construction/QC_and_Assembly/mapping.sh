#!/bin/bash

source mobivision-v3.0/source.sh

###Each species was processed using the same workflow/pipeline.
mobivision mkindex -n Spe_genome \
-f genome.fasta \
-g genomic.gtf

mobivision quantify -i Spe_genome \
-f readsfile/ \
-t 16 \
-o ./outs