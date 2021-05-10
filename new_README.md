# workflow for gnomad re-analysis
* this is a workflow for re-analysing results from the gnomad dataset using a different aggregation 
* Using gsutil, download the prop-observed table from gs://gnomad-public/papers/2019-flagship-lof/v1.0/standard/prop_observed_standard.txt.bgz 
    - this is about 4GB. google cloud accounts come with 1TB of data downloads
    - Can be replaced with any of the other results if you so choose. Put it into`data/standard/prop_observed_standard`
* then run the script `constraint/finalise_constraint_custom.py` with selected parameters