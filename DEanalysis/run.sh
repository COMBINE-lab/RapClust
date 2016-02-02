#!/bin/sh

echo "Enter clust file"

read clustfile

adir=$(dirname ${clustfile})

echo "Enter species"

read species

python GetContig2Clust.py $clustfile ${adir}/rapclust_clusters.flat

cfile=rapclust_clusters.flat

Rscript GetDEGenes.R truth $adir contig2cuffGene.tsv $species
Rscript GetDEGenes.R sailfish $adir $cfile $species
Rscript GetDEGenes.R corset $adir corset-counts.txt $species

python DEcomp.py --method sailfish --adir $adir
python DEcomp.py --method corset --adir $adir
