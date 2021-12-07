tail -n +45 data/GPCR_synthetic_variants_output.txt | awk '{if(substr(4,1,11) != IMPACT=MODI) {print -bash}}' - > data/GPCR_synthetic_coding_variants.txt
