#!/usr/bin/python

import argparse
from itertools import product
import pickle
from typing import Dict, List, Optional, Set, Tuple, Any
import hail as hl
import pandas as pd
from gnomad_pipeline.utils import *

# Summary of pipeline steps
# Get gene list & parameters 
# Pull relevant data from 
#   (a) exomes/genomes variant call table 
#   (b) site context table 
#   (c) mutation rate table
# Pre-process exomes table (and possibly others)
# Pull all of coverage model table (prop_observed_by_coverage) and build a coverage model
# Calculate all possible variants in relevant genes with expected rate of observation
# Calculate proportion of variants observed by category
# Finalise and calculate summary stats for release

def get_gene_intervals():
    # Get Ensembl gene intervals from file
    df_intervals = pd.read_csv('data/Ensembl_Grch37_gpcr_genome_locations.csv')
    df_intervals = df_intervals[['HGNC symbol','Grch37 symbol','Grch37 chromosome','Grch37 start bp','Grch37 end bp']]
    df_intervals['locus_interval_txt'] = df_intervals['Grch37 chromosome'] + ':'  + \
        + df_intervals['Grch37 start bp'].map(str) + '-' + df_intervals['Grch37 end bp'].map(str)
    return list(df_intervals['locus_interval_txt'].map(hl.parse_locus_interval))


def get_exomes(exomes_path, gene_intervals):
    # Path to full exome table
    full_exome_ht = hl.read_table(exomes_path)
    # Select relevant columns to avoid getting too much data 
    full_exome_ht = full_exome_ht.select(
        full_exome_ht.freq,
        full_exome_ht.vep.transcript_consequences,
        full_exome_ht.context,
        full_exome_ht.methylation,
        full_exome_ht.coverage
    )
    # Filter by gene
    exome_ht = hl.filter_intervals(full_exome_ht, gene_intervals)
    # Prepare exomes
    exome_ht = prepare_ht(exome_ht)
    return exome_ht


def get_context(context_ht_path, gene_intervals):
    full_context_ht = hl.read_table(context_ht_path)
    full_context_ht = full_context_ht.select(
        full_context_ht.vep.transcript_consequences,
        full_context_ht.context,
        full_context_ht.methylation,
        full_context_ht.coverage
    )
    # Filter this table in similar way to exomes ht
    selected_context_ht = hl.filter_intervals(full_context_ht,gene_intervals)
    selected_context_ht = prepare_ht(selected_context_ht)
    return selected_context_ht


def prepare_exomes_and_context(
        exome_ht: hl.Table,
        context_ht: hl.Table,
        dataset: str = 'gnomad', 
        impose_high_af_cutoff_upfront: bool = True
    ) -> hl.Table:
    '''Prepare exome mutation data for modelling. Untested'''
    # For standard analysis, assume the worst predicted transcript consequence takes place
    exome_ht = add_most_severe_csq_to_tc_within_ht(exome_ht)
    exome_ht = exome_ht.transmute(transcript_consequences=exome_ht.vep.transcript_consequences)
    exome_ht = exome_ht.explode(exome_ht.transcript_consequences)
    
    # Annotate variants with grouping variables. 
    exome_ht, grouping = annotate_constraint_groupings(exome_ht)
    exome_ht = exome_ht.select('context', 'ref', 'alt', 'methylation_level', 'freq', 'pass_filters', *grouping)

    # annotate and filter context table in same way. Exome coverage or genome coverage?
    context_ht, grouping = annotate_constraint_groupings(context_ht)
    context_ht = context_ht.filter(hl.is_defined(context_ht.exome_coverage)).select(
        'context', 'ref', 'alt', 'methylation_level', *grouping)

    # Filter by allele count > 0, allele fraction < cutoff, coverage > 0
    af_cutoff = 0.001
    freq_index = exome_ht.freq_index_dict.collect()[0][dataset]
    def keep_criteria(ht):
        crit = (ht.freq[freq_index].AC > 0) & ht.pass_filters & (ht.coverage > 0)
        if impose_high_af_cutoff_upfront:
            crit &= (ht.freq[freq_index].AF <= af_cutoff)
        return crit
    exome_ht = exome_ht.filter(keep_criteria(exome_ht))

    return exome_ht, context_ht, grouping


def get_possible_variants(context_ht,mutation_ht, coverage_model, plateau_models, grouping, possible_file,
        half_cutoff = False, HIGH_COVERAGE_CUTOFF = 40, POPS = ('global', 'afr', 'amr', 'eas', 'nfe', 'sas')):
    '''Compute table of possible variants with needed properties. Untested'''
    ht = count_variants(context_ht, additional_grouping=grouping, partition_hint=2000, force_grouping=True)
    ht = annotate_with_mu(ht, mutation_ht)
    ht = ht.transmute(possible_variants=ht.variant_count)
    ht = annotate_variant_types(ht.annotate(mu_agg=ht.mu_snp * ht.possible_variants))
    model = hl.literal(plateau_models.total)[ht.cpg]
    cov_cutoff = (HIGH_COVERAGE_CUTOFF / half_cutoff) if half_cutoff else HIGH_COVERAGE_CUTOFF
    ann_expr = {
        'adjusted_mutation_rate': ht.mu_agg * model[1] + model[0],
        'coverage_correction': hl.case()
            .when(ht.coverage == 0, 0)
            .when(ht.coverage >= cov_cutoff, 1)
            .default(coverage_model[1] * hl.log10(ht.coverage) + coverage_model[0])
    }
    for pop in POPS:
        pop_model = hl.literal(plateau_models[pop])
        slopes = hl.map(lambda f: f[ht.cpg][1], pop_model)
        intercepts = hl.map(lambda f: f[ht.cpg][0], pop_model)
        ann_expr[f'adjusted_mutation_rate_{pop}'] = ht.mu_agg * slopes + intercepts
    ht = ht.annotate(**ann_expr)
    ann_expr = {
        'expected_variants': ht.adjusted_mutation_rate * ht.coverage_correction,
        'mu': ht.mu_agg * ht.coverage_correction
    }
    for pop in POPS:
        ann_expr[f'expected_variants_{pop}'] = ht[f'adjusted_mutation_rate_{pop}'] * ht.coverage_correction
    ht = ht.annotate(**ann_expr)
    ht.write(possible_file, True)
    return ht
    

def get_proportion_observed(
        exome_ht: hl.Table, 
        context_ht: hl.Table,
        mutation_ht: hl.Table,
        coverage_model, plateau_models,
        possible_variants_ht_path,
        grouping: List,
        impose_high_af_cutoff_upfront: bool = True) -> hl.Table:
    '''Aggregate by grouping variables'''

    # Get possible variants for this context table
    possible_variants_ht = get_possible_variants(context_ht,mutation_ht, coverage_model, plateau_models, grouping, possible_variants_ht_path)

    # Count observed variants by grouping
    ht = count_variants(exome_ht, additional_grouping=grouping, partition_hint=2000, force_grouping=True,
                        count_downsamplings=False, impose_high_af_cutoff_here=not impose_high_af_cutoff_upfront,
                        POPS = ('global', 'afr', 'amr', 'eas', 'nfe', 'sas'))

    # Match with expected frequency of variants by grouping
    ht = ht.join(possible_variants_ht, 'outer')

    # Aggregate all groupings across populations
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
    # ht = calculate_all_z_scores(ht)
    return ht

def run_tests(exomes_path, context_path, po_observed_by_coverage_path):
    """Tests loading of autosome po table"""
    print('Test access to exomes and filtering')
  

    # Also need to test: storing in files and reloading, pre-processing, model usage, etc.

def main(args):
    # Set paths for data access based on command line parameters
    POPS = ('global', 'afr', 'amr', 'eas', 'nfe', 'sas') # set pops

    # New paths to public access bucket
    exomes_path = 'gs://gcp-public-data--gnomad/papers/2019-flagship-lof/v1.0/model/exomes_processed.ht/'
    context_path_vep = 'gs://gcp-public-data--gnomad/resources/context/grch37_context_vep_annotated.ht/'
    mutation_rate_ht_path = 'gs://gcp-public-data--gnomad/papers/2019-flagship-lof/v1.0/model/mutation_rate_methylation_bins.ht'
    po_coverage_ht_path = 'gs://gcp-public-data--gnomad/papers/2019-flagship-lof/v1.0/model/prop_observed_by_coverage_no_common_pass_filtered_bins.ht'

    # Paths for local storage
    root = './data'
    exomes_local_path = f'{root}/exomes.ht'
    context_local_path = f'{root}/context.ht'
    mutation_rate_local_path = f'{root}/context.ht'
    po_coverage_local_path = f'{root}/prop_observed_by_coverage_no_common_pass_filtered_bins.ht'
    
    po_ht_path = f'{root}/{{subdir}}/prop_observed_{{subdir}}.ht'
    raw_constraint_ht_path = f'{root}/{{subdir}}/constraint_{{subdir}}.ht'
    final_constraint_ht_path = f'{root}/{{subdir}}/constraint_final_{{subdir}}.ht'
    possible_variants_ht_path = f'{root}/model/possible_data/possible_transcript_pop_{args.model}.ht'

    po_output_path = po_ht_path.format(subdir=args.model)
    output_path = raw_constraint_ht_path.format(subdir=args.model)
    final_path = final_constraint_ht_path.format(subdir=args.model)

    hl.init(quiet=True)
    if args.test:
        # Currently tests whether it can read data from exomes sites table
        #test_genes = ['ADRB1']
        test_intervals = get_gene_intervals()
        # print(test_intervals)
        # ht_exomes_test = get_exomes(exomes_path, test_intervals)
        # ht_exomes_test.show()
        # ht_context_test = get_context(context_path,test_intervals)
        # ht_context_test.show()
       
 
        coverage_ht = hl.read_table(po_coverage_local_path)
        coverage_model, plateau_models = build_models(coverage_ht, args.trimers, True) 
        with open('data/coverage_models.pkl','wb') as fid:
            pickle.dump((coverage_model,plateau_models),fid)
        
        print('Ran coverage model building')

        #run_tests(exomes_path,context_path_vep, po_coverage_local_path)

    if args.get_data:
        print('Getting data from Google Cloud...')
        # Get intervals for chosen genes from Ensembl
        gene_intervals = get_gene_intervals()
        # Get relevant exomes data 
        full_exome_ht = get_exomes(exomes_path, gene_intervals)

        # filter into X, Y and autosomal regions
        exome_ht = full_exome_ht.filter(full_exome_ht.locus.in_autosome_or_par())
        exome_x_ht = hl.filter_intervals(full_exome_ht, [hl.parse_locus_interval('X')])
        exome_x_ht = exome_x_ht.filter(exome_x_ht.locus.in_x_nonpar())
        exome_y_ht = hl.filter_intervals(full_exome_ht, [hl.parse_locus_interval('Y')])
        exome_y_ht = exome_y_ht.filter(exome_y_ht.locus.in_y_nonpar())
        exome_ht.write(exomes_local_path,overwrite=args.overwrite)
        exome_x_ht.write(exomes_local_path.replace('.ht', '_x.ht'),overwrite=args.overwrite)
        exome_y_ht.write(exomes_local_path.replace('.ht', '_y.ht'),overwrite=args.overwrite)

        # Do the same for context table
        full_context_ht = get_context(context_path_vep,gene_intervals)
        context_ht = full_context_ht.filter(full_context_ht.locus.in_autosome_or_par())
        context_x_ht = hl.filter_intervals(full_context_ht, [hl.parse_locus_interval('X')])
        context_x_ht = context_x_ht.filter(context_x_ht.locus.in_x_nonpar())
        context_y_ht = hl.filter_intervals(full_context_ht, [hl.parse_locus_interval('Y')])
        context_y_ht = context_y_ht.filter(context_y_ht.locus.in_y_nonpar())
        context_ht.write(context_local_path,overwrite=args.overwrite)
        context_x_ht.write(context_local_path.replace('.ht', '_x.ht'),overwrite=args.overwrite)
        context_y_ht.write(context_local_path.replace('.ht', '_y.ht'),overwrite=args.overwrite)

        # Get methylation-dependent mutation rate
        mutation_rate_ht = hl.read_table(mutation_rate_ht_path).select('mu_snp')
        mutation_rate_ht.write(mutation_rate_local_path,overwrite=args.overwrite)

        # Build model for observation rate based on coverage
        coverage_ht = hl.read_table(po_coverage_local_path)
        coverage_model, plateau_models = build_models(coverage_ht, args.trimers, True) 
        with open('data/coverage_models.pkl','wb') as fid:
            pickle.dump((coverage_model,plateau_models),fid)
        



    if args.aggregate:
        print('Running aggregation of variants by grouping variables')
        # Set chosen grouping variables to aggregate on
        grouping = ['gene','annotation','modifier']

        # aggregate by chosen groupings & get proportion observed; write to file
        input_hts = zip((exome_ht, exome_x_ht, exome_y_ht), (context_ht, context_x_ht, context_y_ht))
        output_hts = []
        for exome_ht_, context_ht_ in input_hts:
            output_hts.append(get_proportion_observed(exome_ht_,context_ht_,grouping))
        po_exome_ht, po_exome_x_ht, po_exome_y_ht = output_hts

        # Write proportion observed to file        
        po_exome_ht.write(po_output_path, overwrite=args.overwrite)
        po_exome_x_ht.write(po_output_path.replace('.ht', '_x.ht'), overwrite=args.overwrite)
        po_exome_y_ht.write(po_output_path.replace('.ht', '_y.ht'), overwrite=args.overwrite)

       
    if args.finalise:
        print('Joining tables and running aggregation by gene')
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
    parser.add_argument('--get_data',help='Extract data from google cloud',action='store_true')
    parser.add_argument('--aggregate', help='Aggregate by grouping variables', action='store_true')
    parser.add_argument('--finalise',help='Calculate estimates',action='store_true')
    parser.add_argument('--summarise', help='Report summary stats', action='store_true')
    args = parser.parse_args()
    main(args)