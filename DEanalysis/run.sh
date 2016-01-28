#!/bin/sh

echo "Enter clust file"

read clustfile

#/mnt/scratch3/avi/clustering/logFoldChange/sailfish/full/graph.clust.filt

python GetContig2Clust.py $clustfile /home/laraib/clust/DE_analysis/contig2clust.tsv

Rscript GetDEGenes.R truth
Rscript GetDEGenes.R sailfish
Rscript GetDEGenes.R corset

python DEcomp.py --method sailfish
python DEcomp.py --method corset
