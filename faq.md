---
layout: page
title: FAQ 
subtitle: Frequently Asked Questions about RapClust 
---

**Q**: What is RapClust?

  - RapClust is a tool for clustering contigs in *de novo* transcriptomes.  By clustering related contigs, 
    we hope to group together contigs deriving from the same transcript or gene.  Given a good clustering, 
    one can obtain more robust and accurate estimates of "group" expression, which can also improve downstream
    analysis (e.g. differential expression).


**Q**: How does RapClust work?

  - RapClust works by assessing the similarity of contigs based on the fraction of reads that multi-map between them.
    The multi-mapping fragments induce a *mapping ambiguity graph*, which can then be partitioned to yield clusters of 
    contigs that are highly self-similar.  A full manuscript describing the methodology behind RapClust will be available
    soon.

