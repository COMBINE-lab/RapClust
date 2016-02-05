# RapClust :  Accurate, Fast and Lightweight Clustering of de novo Transcriptomes using Fragment Equivalence Classes
============

RapClust is a tool for clustering contigs from *de novo* transcriptome assemblies.  RapClust is designed to be run downstream of the [Sailfish](https://github.com/kingsfordgroup/sailfish) or [Salmon](https://github.com/COMBINE-lab/salmon) tools for rapid transcript-level quantification.  Specifically, RapClust relies on the *fragment equivalence classes* computed by these tools in order to determine how seqeunce is shared across the transcriptome, and how reads map to potentially-related contigs across different conditions.  RapClust is heavily inspired by the approach of [Corset](https://github.com/Oshlack/Corset), and one of the main goals of the tool is to make the same type of high-quality *de novo* transcriptome clusterings available in the new breed of ultra-fast expression analysis pipelines.  RapClust achieves its speed partly by replacing traditional alignment with the novel concept of [quasi-mapping](https://github.com/COMBINE-lab/RapMap), which yeilds sufficient information for highly-accurate quantification and clustering orders of magnitude faster than standard read alignment tools.  RapClust also achieves its speed and accuracy by relying on the concise (yet surprisingly rich) information contained in fragment-level equivalence classes, and the mapping ambiguity graph they induce.

## Using RapClust
=================

RapClust is written in Python, is easy to use, and adds only marginal runtime to what is already required for rapid *de novo* transcriptome quantification (typically a few seconds or minutes, even for very large experiments with many reads and samples).  RapClust is compatible with both Sailfish and Salmon; in the instructions below, we explain how to use it with Sailfish.  There are two main steps involved in running RapClust:

    1. Run Sailfish on each sample in your experiment, passing it the `--dumpEq` option.  This will tell Sailfish to dump a representation of the fragment equivalence classes that it computed during quasi-mapping of each sample.  Apart from this additional option, Sailfish should be run normally (i.e. passing in whatever other options are appropriate for your samples).
    
    2. Run RapClust, providing it with a configuration file that describes the experimental setup of your samples, and where the Sailfish quantification results have been written.
    
Let's illustrate this pipeline with a particular example, the following experimental data from the [Trapnell et al. paper](http://www.nature.com/nbt/journal/v31/n1/full/nbt.2450.html):

Accession | Condition | Replicate
----------|-----------|----------
SRR493366 | scramble  | 1
SRR493367	| scramble  | 2
SRR493368	| scramble  | 3
SRR493369	| HOXA1KD	  | 1
SRR493370	| HOXA1KD	  | 2
SRR493371 | HOXA1KD   | 3

We'll assume that the raw read files reside in the directory `reads`.  Assuming that you've already built the index on the transcriptome you wish to quantify, a typical run of Sailfish on this data would look something like.

```
> parallel -j 6 "samp={}; sailfish quant -i index -l IU -1 <(gunzip -c reads/{$samp}_1.fq.gz) -2 <(gunzip -c reads/{$samp}_2.fq.gz) -o {$samp}_quant --dumpEq -p 4" ::: SRR493366 SRR493367 SRR493368 SRR4933669 SRR493370 SRR493371
```

This will quantify each sample, and write the result to the directory `samplename_quant`.  Given this setup, we're now ready to run RapClust.  First, we have to make an appropriate config file.  We'll use the following:

```
conditions:
    - Control
    - HOXA1 Knockdown
samples:
    Control:
        - SRR493366_quant
        - SRR493367_quant
        - SRR493368_quant
    HOXA1 Knockdown:
        - SRR493369_quant
        - SRR493370_quant
        - SRR493371_quant
outdir: human_rapclust
```

you can place this in a file called `config.yaml`.  RapClust uses [YAML](http://yaml.org/) to specify its configuration files.  The setup here is hopefully self-explanatory.  There configuration file must contain the following three entries; `conditions`, `samples`, and `outdir`.  The `conditions` entry lists the conditions present in the sample. The `samples` entry is a nested dictionary of lists; there is a key corrseponding to each condition listed in the `conditions` entry, and the value associated with this key is a list of quantification directories of the samples for this condition.  Finally, the `outdir` entry specifies where the RapClust output and intermediate files should be stored.  Given the above, we can run RapClust as:

```
> ./RapClust.py --config config.yaml
```

This will process the samples, generate the mapping ambiguity graph, filter it according to the conditions, and cluster the resuling graph (RapClust uses [MCL](http://micans.org/mcl/) internally for clustering).  Once RapClust is finished, the `human_rapclust` directory should exist.  It will contain the following files:

`mag.clust, mag.filt.net,  mag.flat.clust,  mag.net,  stats.json`

The most important file for downstream processing is `mag.flat.clust`.  It contains the computed cluster information in a "transcript-to-gene" mapping formation that is compatible with downstream tools like [tximport](https://github.com/mikelove/tximport).  The otherfiles may be useful for exploration, but they are more intended for RapClusts's internal use (e.g. `mag.filt.net` contains the filtered mapping ambiguity graph that is used for clustering).


## Citations:
=============

Differential analysis of gene regulation at transcript resolution with RNA-seq by Cole Trapnell, David G Henderickson, Martin Savageau, Loyal Goff, John L Rinn and Lior Pachter, Nature Biotechnology 31, 46–53 (2013).

Stijn van Dongen. Graph Clustering by Flow Simulation. PhD thesis, University of Utrecht, 2000

Charlotte Soneson, Michael I Love, and Mark D Robinson. Differential analyses for rna-seq: transcript-level estimates improve gene-level inferences. F1000Research, 4, 2015.
