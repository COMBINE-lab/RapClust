#!/bin/sh

echo "Enter clust file"

read clustfile

#/mnt/scratch3/avi/clustering/logFoldChange/sailfish/full/graph.clust.filt

adir=$(dirname ${clustfile})

python GetContig2Clust.py $clustfile ${adir}/rapclust_clusters.flat

cfile=rapclust_clusters.flat

#Rscript GetDEGenes.R truth $adir contig2cuffGene.tsv
Rscript GetDEGenes.R sailfish $adir $cfile
#Rscript GetDEGenes.R corset $adir corset-counts.txt

python DEcomp.py --method sailfish --adir $adir
python DEcomp.py --method corset --adir $adir
