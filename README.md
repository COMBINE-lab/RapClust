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


## Citations:
=============
