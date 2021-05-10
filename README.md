# workflow for gnomad re-analysis
* this is a workflow for re-analysing results from the gnomad dataset using custom aggregation groups.The constraint computation pipeline is written in [Hail 0.2](https://hail.is) and is based on the workflow from the gnomAD flagship manuscript [Karczewski et al., 2019](https://www.biorxiv.org/content/10.1101/531210v2)
* Using gsutil, download the prop-observed table from gs://gnomad-public/papers/2019-flagship-lof/v1.0/standard/prop_observed_standard.txt.bgz 
    - this is about 4GB. google cloud accounts come with 1TB of data downloads
    - Can be replaced with any of the other results if you so choose. Put it into`data/standard/prop_observed_standard`
* then run the script `src/finalise_constraint_custom.py` with selected parameters
* the code in `constraint`, `constraint_utils`, and `gnomad` is provided for reference and will be removed in future versions