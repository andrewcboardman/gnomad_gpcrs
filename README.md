# workflow for gnomad re-analysis
## Introduction
* gnomadIC (gnomAD Inference of Constraint) is a workflow for re-analysing results from the gnomad dataset using custom annotations and aggregation groups.
The constraint computation pipeline is written in [Hail 0.2](https://hail.is) and is based on the workflow from the gnomAD flagship manuscript [Karczewski et al., 2019](https://www.biorxiv.org/content/10.1101/531210v2)
* Analysis can be run using a set of chosen genes using custom_constraint_analysis.py. You can either keep standard settings (aggegrating missense variants using provided PolyPhen2 labels) or add your own annotations.
* Custom annotations must be linked to HGVSP ids for variants of your choice. These can then be used to modify the annotation column before aggregation.

## Standard workflow
* Models for the expected number of variants based on coverage, nucleotide context and methylation status are downloaded from the gnomAD gsutil bucket.
* Observed and possible variant sets are downloaded from the same source and stored. They can be extracted using the script `extract_variant_tables.py`
* The possible variant set is annotated with the expected frequency of each variant.
* Both the observed and expected variant counts are aggregated for each possible protein effect and the results are merged.
* 
