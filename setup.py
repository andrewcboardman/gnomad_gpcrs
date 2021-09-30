import os

def setup_paths(model):
    root = './data'
    # Google storage paths 
    gs_paths = dict(
        exomes_path = 'gs://gcp-public-data--gnomad/papers/2019-flagship-lof/v1.0/model/exomes_processed.ht/',
        context_path = 'gs://gcp-public-data--gnomad/resources/context/grch37_context_vep_annotated.ht/',
        mutation_rate_path = 'gs://gcp-public-data--gnomad/papers/2019-flagship-lof/v1.0/model/mutation_rate_methylation_bins.ht',
        po_coverage_path = 'gs://gcp-public-data--gnomad/papers/2019-flagship-lof/v1.0/model/prop_observed_by_coverage_no_common_pass_filtered_bins.ht'
    )
    # Local input paths
    local_input_paths = dict(
        exomes_local_path = f'{root}/exomes.ht',
        context_local_path = f'{root}/context.ht',
        mutation_rate_local_path = f'{root}/mutation_rate_methylation_bins.ht',
        po_coverage_local_path = f'{root}/prop_observed_by_coverage_no_common_pass_filtered_bins.ht',
        coverage_models_local_path = f'{root}/coverage_models.pkl'
    )
    # Local output paths
    output_subdir = f'{root}/{model}'
    if not os.path.isdir(output_subdir):
        os.mkdir(output_subdir)
    local_output_paths = dict(
        possible_variants_ht_path = f'{output_subdir}/possible_transcript_pop.ht',
        po_output_path = f'{output_subdir}/prop_observed.ht',
        finalized_output_path = f'{output_subdir}/constraint.ht',
        summary_output_path = f'{output_subdir}/constraint_final.ht'
    )
    paths = {**gs_paths, **local_input_paths, **local_output_paths}

    return paths