# workflow for gnomad re-analysis
* gnomADIC (gnomAD inference) is a workflow for re-analysing results from the gnomad dataset using custom aggregation groups.The constraint computation pipeline is written in [Hail 0.2](https://hail.is) and is based on the workflow from the gnomAD flagship manuscript [Karczewski et al., 2019](https://www.biorxiv.org/content/10.1101/531210v2)
* Using gsutil, download the prop-observed table from the URLs given below ; this is about 4GB. Google cloud accounts come with 1TB of free download from public repos like this one. Put it into`data/standard/prop_observed_standard`. If you want to use differently pre-processed datasets, these can be found in the same Google cloud bucket.
* To run the analysis, use the script `src/finalise_constraint_custom.py` with selected parameters. Aggregation should take about 25 minutes on a desktop.
* Current workflow allows reannotation at the level of gene-transcript-consequence. 
* Functions for reannotation and aggregation with custom labels have been added but are not fully functional yet and should not be used.

URLs for standard data download:
`gs://gnomad-public/papers/2019-flagship-lof/v1.0/standard/prop_observed_standard.txt.bgz`, `gs://gnomad-public/papers/2019-flagship-lof/v1.0/standard/prop_observed_standard_x.txt.bgz` and `gs://gnomad-public/papers/2019-flagship-lof/v1.0/standard/prop_observed_standard_y.txt.bgz`