#!/bin/bash

echo "Type species name"

read species

Rscript GetDEGenes.R truth /mnt/scratch3/avi/clustering/data/$species/trinity/sailfish/sailfish_quant/ /mnt/scratch3/avi/clustering/data/corsetData/${species^}-Trinity/contigs2genes.disambiguous.txt ./ $species
Rscript GetDEGenes.R sailfish /mnt/scratch3/avi/clustering/data/$species/trinity/sailfish/sailfish_quant/ /mnt/scratch3/avi/clustering/data/$species/trinity/sailfish/results/mag.flat.clust ./ $species
Rscript GetDEGenes.R corset /mnt/scratch3/avi/clustering/data/$species/trinity/sailfish/sailfish_quant/ /mnt/scratch3/avi/clustering/data/corsetData/${species^}-Trinity/corset-counts.txt ./ $species

python DEcomp.py --method sailfish --clustfile /mnt/scratch3/avi/clustering/data/$species/trinity/sailfish/results/mag.flat.clust --contig2gene /mnt/scratch3/avi/clustering/data/corsetData/${species^}-Trinity/contigs2genes.disambiguous.txt --adir ./ --odir ./
python DEcomp.py --method corset --clustfile /mnt/scratch3/avi/clustering/data/corsetData/${species^}-Trinity/corset-clusters.txt --contig2gene /mnt/scratch3/avi/clustering/data/corsetData/${species^}-Trinity/contigs2genes.disambiguous.txt --adir ./ --odir ./
