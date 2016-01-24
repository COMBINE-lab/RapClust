#!/bin/sh

Rscript GetDEGenes.R truth
Rscript GetDEGenes.R sailfish
#Rscript GetDEGenes.R corset

python DEcomp.py sailfish
python DEcomp.py corset
