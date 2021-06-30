import hail as hl; hl.init()

types = {'context': hl.tstr,
    'ref': hl.tstr,
    'alt': hl.tstr,
    'methylation_level': hl.tint32,
    'exome_coverage': hl.tint32,
    'variant_count': hl.tint64,
    'downsampling_counts_global': hl.tarray(hl.tint64),
    'downsampling_counts_afr': hl.tarray(hl.tint64),
    'downsampling_counts_amr': hl.tarray(hl.tint64),
    'downsampling_counts_eas': hl.tarray(hl.tint64),
    'downsampling_counts_nfe': hl.tarray(hl.tint64),
    'downsampling_counts_sas': hl.tarray(hl.tint64),
    'mu_snp': hl.tfloat64,
    'possible_variants': hl.tint64,
    'transition': hl.tbool,
    'cpg': hl.tbool,
    'variant_type': hl.tstr,
    'variant_type_model': hl.tstr}

ht = hl.import_table('data/prop_observed_by_coverage_no_common_pass_filtered_bins.txt.bgz',types=types)
ht.write('data/prop_observed_by_coverage_no_common_pass_filtered_bins.ht')
    