#!/usr/bin/env python3
import argparse

from itertools import product
import hail as hl
from typing import Dict, List, Optional, Set, Tuple
from gnomad.utils.vep import *
import hail as hl
import argparse
from itertools import product
import hail as hl
from typing import Dict, List, Optional, Set, Tuple, Any
from gnomad_pipeline.utils import *

# set pops
POPS = ('global', 'afr', 'amr', 'eas', 'nfe', 'sas')
HIGH_COVERAGE_CUTOFF = 40

def preprocess_exomes_and_genomes():
    '''Add context, coverage and methylation data to genome and exome hail tables'''
    # This is not ready yet; 
    # context, coverage and methylation data are too large to download
    # URLs stored here for the moment; need to access these on cloud and filter for required genes
    exomes_path = 'gs://gcp-public-data--gnomad/release/2.1/ht/exomes/gnomad.exomes.r2.1.sites.ht/'
    coverage_path_exomes = 'gs://gcp-public-data--gnomad/release/2.1/coverage/exomes/gnomad.exomes.r2.1.coverage.ht'
    context_path_vep = 'gs://gcp-public-data--gnomad/resources/context/grch37_context_vep_annotated.ht/'
    methylation_path = 'gs://gnomad-public/resources/grch37/methylation_sites/methylation.ht'

    split_context_mt(vep_context_ht_path, {'exomes': coverage_ht_path('exomes'), 'genomes': coverage_ht_path('genomes')},
                         methylation_sites_ht_path(), context_ht_path, args.overwrite)
    pre_process_data(get_gnomad_public_data('genomes'), context_ht_path, processed_genomes_ht_path, args.overwrite)
    pre_process_data(get_gnomad_public_data('exomes'), context_ht_path, processed_exomes_ht_path, args.overwrite)


def prepare_exomes(exome_ht: hl.Table, groupings: List, impose_high_af_cutoff_upfront: bool = True) -> hl.Table:

    # Manipulate VEP annotations and explode by them
    exome_ht = add_most_severe_csq_to_tc_within_ht(exome_ht)
    exome_ht = exome_ht.transmute(transcript_consequences=exome_ht.vep.transcript_consequences)
    exome_ht = exome_ht.explode(exome_ht.transcript_consequences)
    
    # Annotate variants with grouping variables. 
    exome_ht, grouping = annotate_constraint_groupings(exome_ht,groupings) # This function needs to be adapted
    exome_ht = exome_ht.select(
        'context', 'ref', 'alt', 'methylation_level', 'freq', 'pass_filters', *groupings)

    # Filter by allele count
    # Likely to need to adapt this function as well
    af_cutoff = 0.001
    freq_index = exome_ht.freq_index_dict.collect()[0][dataset]

    def keep_criteria(ht):
        crit = (ht.freq[freq_index].AC > 0) & ht.pass_filters & (ht.coverage > 0)
        if impose_high_af_cutoff_upfront:
            crit &= (ht.freq[freq_index].AF <= af_cutoff)
        return crit

    exome_ht = exome_ht.filter(keep_criteria(exome_ht))
    return exome_ht


def get_proportion_observed(
        exome_ht: hl.Table, 
        possible_variants_ht: hl.Table,
        groupings: List,
        dataset: str = 'gnomad', 
        impose_high_af_cutoff_upfront: bool = True, 
        half_cutoff = False) -> hl.Table:
    '''Filter exomes and aggregate by grouping variables'''
    # Filter variant table and annotate with grouping variables
    exome_ht = prepare_exomes(exome_ht, groupings)  
    # Aggregate on grouping variables    
    # Should probably adapt this function
    ht = count_variants(exome_ht, additional_grouping=groupings, partition_hint=2000, force_grouping=True,
                        count_downsamplings=POPS, impose_high_af_cutoff_here=not impose_high_af_cutoff_upfront)
    ht = ht.join(possible_variants_ht, 'outer')
    agg_expr = {
        'variant_count': hl.agg.sum(ht.variant_count),
        'adjusted_mutation_rate': hl.agg.sum(ht.adjusted_mutation_rate),
        'possible_variants': hl.agg.sum(ht.possible_variants),
        'expected_variants': hl.agg.sum(ht.expected_variants),
        'mu': hl.agg.sum(ht.mu)
    }
    for pop in POPS:
        agg_expr[f'adjusted_mutation_rate_{pop}'] = hl.agg.array_sum(ht[f'adjusted_mutation_rate_{pop}'])
        agg_expr[f'expected_variants_{pop}'] = hl.agg.array_sum(ht[f'expected_variants_{pop}'])
        agg_expr[f'downsampling_counts_{pop}'] = hl.agg.array_sum(ht[f'downsampling_counts_{pop}'])
    ht = ht.group_by(*grouping).partition_hint(1000).aggregate(**agg_expr)
    return ht.annotate(obs_exp=ht.variant_count / ht.expected_variants)

def finalize_dataset(po_ht: hl.Table, keys: Tuple[str] = ('gene', 'transcript', 'canonical'),
                     n_partitions: int = 1000) -> hl.Table:
    '''aggregate variants to calculate constraint metrics and significance'''
    # Tidy this code a bit
    # This function aggregates over genes in all cases, as XG spans PAR and non-PAR X
    po_ht = po_ht.repartition(n_partitions).persist()

    # Getting classic LoF annotations (no LOFTEE)
    classic_lof_annotations = hl.literal({'stop_gained', 'splice_donor_variant', 'splice_acceptor_variant'})
    lof_ht_classic = po_ht.filter(classic_lof_annotations.contains(po_ht.annotation) &
                                  ((po_ht.modifier == 'HC') | (po_ht.modifier == 'LC')))
    lof_ht_classic = collapse_lof_ht(lof_ht_classic, keys, False)
    lof_ht_classic = lof_ht_classic.rename({x: f'{x}_classic' for x in list(lof_ht_classic.row_value)})

    # Getting all LoF annotations (LOFTEE HC + OS)
    lof_ht_classic_hc = po_ht.filter((po_ht.modifier == 'HC') | (po_ht.modifier == 'OS'))
    lof_ht_classic_hc = collapse_lof_ht(lof_ht_classic_hc, keys, False)
    lof_ht_classic_hc = lof_ht_classic_hc.rename({x: f'{x}_with_os' for x in list(lof_ht_classic_hc.row_value)})

    # Getting all LoF annotations (LOFTEE HC)
    lof_ht = po_ht.filter(po_ht.modifier == 'HC')
    lof_ht = collapse_lof_ht(lof_ht, keys, False)

    # Aggregate missense variants
    mis_ht = po_ht.filter(po_ht.annotation == 'missense_variant')
    agg_expr = {
        'obs_mis': hl.agg.sum(mis_ht.variant_count),
        'exp_mis': hl.agg.sum(mis_ht.expected_variants),
        'oe_mis': hl.agg.sum(mis_ht.variant_count) / hl.agg.sum(mis_ht.expected_variants),
        'mu_mis': hl.agg.sum(mis_ht.mu),
        'possible_mis': hl.agg.sum(mis_ht.possible_variants)
    }
    for pop in POPS:
        agg_expr[f'exp_mis_{pop}'] = hl.agg.array_sum(mis_ht[f'expected_variants_{pop}'])
        agg_expr[f'obs_mis_{pop}'] = hl.agg.array_sum(mis_ht[f'downsampling_counts_{pop}'])
    mis_ht = mis_ht.group_by(*keys).aggregate(**agg_expr)

    # Aggregate pphen and non-pphen missense variants
    mis_pphen_ht = po_ht.filter(po_ht.modifier == 'probably_damaging')
    mis_pphen_ht = mis_pphen_ht.group_by(*keys).aggregate(obs_mis_pphen=hl.agg.sum(mis_pphen_ht.variant_count),
                                                          exp_mis_pphen=hl.agg.sum(mis_pphen_ht.expected_variants),
                                                          oe_mis_pphen=hl.agg.sum(mis_pphen_ht.variant_count) / hl.agg.sum(mis_pphen_ht.expected_variants),
                                                          possible_mis_pphen=hl.agg.sum(mis_pphen_ht.possible_variants))

    mis_non_pphen_ht = po_ht.filter(po_ht.modifier != 'probably_damaging')
    mis_non_pphen_ht = mis_non_pphen_ht.group_by(*keys).aggregate(obs_mis_non_pphen=hl.agg.sum(mis_non_pphen_ht.variant_count),
                                                          exp_mis_non_pphen=hl.agg.sum(mis_non_pphen_ht.expected_variants),
                                                          oe_mis_non_pphen=hl.agg.sum(mis_non_pphen_ht.variant_count) / hl.agg.sum(mis_non_pphen_ht.expected_variants),
                                                          possible_mis_non_pphen=hl.agg.sum(mis_non_pphen_ht.possible_variants))

            
    # Aggregate synonymous variants                                
    syn_ht = po_ht.filter(po_ht.annotation == 'synonymous_variant').key_by(*keys)
    agg_expr = {
        'obs_syn': hl.agg.sum(syn_ht.variant_count),
        'exp_syn': hl.agg.sum(syn_ht.expected_variants),
        'oe_syn': hl.agg.sum(syn_ht.variant_count) / hl.agg.sum(syn_ht.expected_variants),
        'mu_syn': hl.agg.sum(syn_ht.mu),
        'possible_syn': hl.agg.sum(syn_ht.possible_variants)
    }
    for pop in POPS:
        agg_expr[f'exp_syn_{pop}'] = hl.agg.array_sum(syn_ht[f'expected_variants_{pop}'])
        agg_expr[f'obs_syn_{pop}'] = hl.agg.array_sum(syn_ht[f'downsampling_counts_{pop}'])
    syn_ht = syn_ht.group_by(*keys).aggregate(**agg_expr)

    # join constraint metrics
    ht = lof_ht_classic.annotate(
        **mis_ht[lof_ht_classic.key], 
        **mis_pphen_ht[lof_ht_classic.key],
        **mis_non_pphen_ht[lof_ht_classic.key],
        **syn_ht[lof_ht_classic.key], 
        **lof_ht[lof_ht_classic.key],
        **lof_ht_classic_hc[lof_ht_classic.key]
        )

    # calculate confidence intervals
    syn_cis = oe_confidence_interval(ht, ht.obs_syn, ht.exp_syn, prefix='oe_syn')
    mis_cis = oe_confidence_interval(ht, ht.obs_mis, ht.exp_mis, prefix='oe_mis')
    lof_cis = oe_confidence_interval(ht, ht.obs_lof, ht.exp_lof, prefix='oe_lof')
    mis_pphen_cis = oe_confidence_interval(ht, ht.obs_mis_pphen, ht.exp_mis_pphen, prefix='oe_mis_pphen')
    mis_non_pphen_cis = oe_confidence_interval(ht, ht.obs_mis_non_pphen, ht.exp_mis_non_pphen, prefix='oe_mis_non_pphen')

    # join confidence intervals
    ht = ht.annotate(
        **syn_cis[ht.key], 
        **mis_cis[ht.key], 
        **lof_cis[ht.key],
        **mis_pphen_cis[ht.key],
        **mis_non_pphen_cis[ht.key]
        )
    # Calculate significance
    # Z score calculation not feasible with partial dataset
    # Need to include flagging of issues in constraint calculations
    #ht = calculate_all_z_scores(ht)
    return ht

def run_tests(ht):
    """Tests loading of autosome po table"""
    # Incorporate more tests to check that aggregation by new variant annotations is legit
    ht = prepare_ht(ht)
    ht.show()

def main(args):
    # Set paths for data access based on command line parameters
    root = './data'
    
    context_ht_path = f'{root}/context/Homo_sapiens_assembly19.fasta.snps_only.vep_20181129.ht'
    processed_genomes_ht_path = f'{root}/model/genomes_processed.ht'
    processed_exomes_ht_path = f'{root}/model/exomes_processed.ht'
    mutation_rate_ht_path = f'{root}/model/mutation_rate_methylation_bins.ht'
    po_coverage_ht_path = f'{root}/model/prop_observed_by_coverage_no_common_pass_filtered_bins.ht'
    po_ht_path = f'{root}/{{subdir}}/prop_observed_{{subdir}}.ht'
    raw_constraint_ht_path = f'{root}/{{subdir}}/constraint_{{subdir}}.ht'
    final_constraint_ht_path = f'{root}/{{subdir}}/constraint_final_{{subdir}}.ht'
    possible_variants_ht_path = f'{root}/model/possible_data/possible_transcript_pop_{args.model}.ht'
    
    variants_table_path = f'{root}/gnomad_v2.1.1_gpcr_variants_unfiltered.tsv'

    po_output_path = po_ht_path.format(subdir=args.model)
    output_path = raw_constraint_ht_path.format(subdir=args.model)
    final_path = final_constraint_ht_path.format(subdir=args.model)

    # Sets method for aggregation, will need to be changed for custom analysis
    MODEL_KEYS = {
        'worst_csq': ['gene'],
        'tx_annotation': ['gene', 'expressed'],
        'standard': ['gene', 'transcript', 'canonical']
    }
    
    if args.test:
        ht = load_variant_table(variants_table_path)
        run_tests(ht)

    if args.get_proportion_observed:
        # Build a model for methylation-dependent mutation rate and apply it to get proportion of variants observed
        # Also need to incorporate genomes and v3 if possible

        print('Running aggregation of variants by grouping variables')
        # Tables of observed mutations in exomes
        full_exome_ht = hl.read_table(processed_exomes_ht_path)
        full_exome_ht = prepare_ht(full_exome_ht, args.trimers) 

        # filter into X, Y and autosomal regions
        exome_ht = full_exome_ht.filter(full_exome_ht.locus.in_autosome_or_par())
        exome_x_ht = hl.filter_intervals(full_exome_ht, [hl.parse_locus_interval('X')])
        exome_x_ht = exome_x_ht.filter(exome_x_ht.locus.in_x_nonpar())
        exome_y_ht = hl.filter_intervals(full_exome_ht, [hl.parse_locus_interval('Y')])
        exome_y_ht = exome_y_ht.filter(exome_y_ht.locus.in_y_nonpar())

        # Modelling results of estimated mutation rates for genes, coverage, methylation level and base context
        possible_variants_ht = hl.read_table(possible_variants_ht_path)

        # Set chosen groupings to aggregate on
        groupings = ['gene','annotation','modifier']

        # Apply model; aggregate by chosen groupings & get proportion observed; write to file
        po_exome_ht, po_exome_x_ht, po_exome_y_ht = \
            [get_proportion_observed(ht, possible_variants_ht, groupings) for ht in (exome_ht, exome_x_ht, exome_y_ht)]
        po_exome_ht.write(po_output_path, overwrite=args.overwrite)
        po_exome_x_ht.write(po_output_path.replace('.ht', '_x.ht'), overwrite=args.overwrite)
        po_exome_y_ht.write(po_output_path.replace('.ht', '_y.ht'), overwrite=args.overwrite)

       
    if args.aggregate:
        print('Running aggregation by gene')
        # read PO hail tables for autosomes, X and Y chromosomes and join them
        print(f'Reading hail table from {po_output_path}')
        # Autosomes
        ht = load_or_import(po_output_path, args.overwrite)
        # X chromosome
        ht_x = load_or_import(po_output_path.replace('.ht','_x.ht'), args.overwrite)
        # Y chromosome
        ht_y = load_or_import(po_output_path.replace('.ht','_y.ht'), args.overwrite)
        # Combine into one table
        ht = ht.union(ht_x).union(ht_y)
        # group by gene/transcript and calculate summary stats
        ht = finalize_dataset(ht, keys=MODEL_KEYS[args.model])
        # write hail table to output path
        ht.write(output_path, args.overwrite)
        hl.read_table(output_path).export(output_path.replace('.ht', '.txt.bgz'))

    if args.summarise:
        print('Finalising summary stats')
        # write summary stats to output path
        ht = load_or_import(output_path, args.overwrite)
        mut_types = ('lof', 'mis', 'syn','mis_pphen','mis_non_pphen')
        output_var_types = zip(('obs', 'exp', 'oe', 'oe', 'oe'),
                                ('', '', '', '_lower', '_upper'))
        output_vars = product(mut_types,output_var_types)
        ht.select(
            'gene','transcript','canonical',
            *[f'{t}_{m}{ci}' for m, (t, ci) in output_vars],
            #*[f'{m}_z' for m in mut_types[:3]], 
            'pLI', 'pRec', 'pNull', 
            #gene_issues=ht.constraint_flag
        ).select_globals().write(final_path, overwrite=args.overwrite)
        hl.read_table(final_path).export(final_path.replace('.ht', '.txt.bgz'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--test', help='Run tests without actually requesting data',action='store_true')
    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--trimers', help='Use trimers instead of heptamers', action='store_true')
    parser.add_argument('--dataset', help='Which dataset to use (one of gnomad, non_neuro, non_cancer, controls)', default='gnomad')
    parser.add_argument('--model', help='Which model to apply (one of "standard", "syn_canonical", or "worst_csq" for now)', default='standard')
    parser.add_argument('--skip_af_filter_upfront', help='Skip AF filter up front (to be applied later to ensure that it is not affecting population-specific constraint): not generally recommended', action='store_true')
    parser.add_argument('--get_proportion_observed',help='Apply model to calculate proportion observed')
    parser.add_argument('--aggregate', help='Aggregate p_obs table', action='store_true')
    parser.add_argument('--summarise', help='Report summary stats', action='store_true')
    args = parser.parse_args()
    main(args)