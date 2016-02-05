#!/bin/bash

echo "Type species name"

read species

#Rscript GetDEGenes.R truth /mnt/scratch3/avi/clustering/data/$species/trinity/sailfish/sailfish_quant/ /mnt/scratch3/avi/clustering/data/corsetData/${species^}-Trinity/contigs2genes.disambiguous.txt ./$species $species
#Rscript GetDEGenes.R sailfish /mnt/scratch3/avi/clustering/data/$species/trinity/sailfish/sailfish_quant/ /mnt/scratch3/avi/clustering/data/$species/trinity/sailfish/results/mag.flat.clust ./$species $species
#Rscript GetDEGenes.R corset /mnt/scratch3/avi/clustering/data/$species/trinity/sailfish/sailfish_quant/ /mnt/scratch3/avi/clustering/data/corsetData/${species^}-Trinity/corset-counts.txt ./$species $species
Rscript GetDEGenes.R cdhit /mnt/scratch3/avi/clustering/data/$species/trinity/sailfish/sailfish_quant/ /mnt/scratch3/avi/clustering/data/$species/trinity/cd_hit/0.95.flat ./$species $species

#python DEcomp.py --method sailfish --clustfile /mnt/scratch3/avi/clustering/data/$species/trinity/sailfish/results/mag.flat.clust --contig2gene /mnt/scratch3/avi/clustering/data/corsetData/${species^}-Trinity/contigs2genes.disambiguous.txt --adir ./$species --odir ./$species
#python DEcomp.py --method corset --clustfile /mnt/scratch3/avi/clustering/data/corsetData/${species^}-Trinity/corset-clusters.txt --contig2gene /mnt/scratch3/avi/clustering/data/corsetData/${species^}-Trinity/contigs2genes.disambiguous.txt --adir ./$species --odir ./$species
python DEcomp.py --method cdhit --clustfile /mnt/scratch3/avi/clustering/data/$species/trinity/cd_hit/0.95.flat --contig2gene /mnt/scratch3/avi/clustering/data/corsetData/${species^}-Trinity/contigs2genes.disambiguous.txt --adir ./$species --odir ./$species

#python ProcessCurves.py --dset $species --label $species/sailfish_de_recovery.scored-label RapClust --label $species/corset_de_recovery.scored-label corset --label $species/cdhit_de_recovery.scored-label cdhit
