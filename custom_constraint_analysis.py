#!/usr/bin/python

import argparse
from itertools import product
import hail as hl
from typing import Dict, List, Optional, Set, Tuple, Union
#from gnomad.utils.vep import *
import argparse
from itertools import product
from typing import Dict, List, Optional, Set, Tuple, Any
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

def get_gene_intervals(genes):
    # Amend this to use Ensembl gene intervals from file
    return [hl.parse_locus_interval('1:START-50K')]

def get_exomes(exomes_path, genes, gene_intervals):
    # Path to full exome table
    full_exome_ht = hl.read_table(exomes_path)
    # limit to first 100 rows for now
    full_exome_ht.head(100)
    # Select relevant columns to avoid getting too much data 
    exome_ht = full_exome_ht.select(
        full_exome_ht.freq,
        full_exome_ht.vep.transcript_consequences
    )
    # Filter by gene
    selected_exomes_ht = hl.filter_intervals(exome_ht, gene_intervals)
    # Prepare exomes
    selected_exomes_ht = prepare_ht(selected_exomes_ht)
    return selected_exomes_ht


def get_context(context_ht_path, gene_intervals):
    full_context_ht = prepare_ht(hl.read_table(context_ht_path), args.trimers).head(5)
    # Filter this table in similar way to exomes ht
    selected_context_ht = hl.filter_intervals(full_context_ht,gene_intervals)
    return selected_context_ht


def get_mutation_rate(mutation_rate_ht_path,gene_intervals):
    mutation_rate_ht = hl.read_table(mutation_rate_ht_path).select('mu_snp').head(5)
    # don't know whether i need to filter this table, take the head to be clear
    return mutation_rate_ht


def prepare_for_modelling(
        exome_ht: hl.Table,
        dataset: str = 'gnomad', 
        impose_high_af_cutoff_upfront: bool = True
    ) -> hl.Table:
    '''Prepare exome mutation data for modelling'''
    # Manipulate VEP annotations and explode by them - why?
    exome_ht = add_most_severe_csq_to_tc_within_ht(exome_ht)
    exome_ht = exome_ht.transmute(transcript_consequences=exome_ht.vep.transcript_consequences)
    exome_ht = exome_ht.explode(exome_ht.transcript_consequences)
    
    # Annotate variants with grouping variables. 
    # This seems to need some things that you don't have (ref, alt, methylation level, pass_filters)
    exome_ht, groupings = annotate_constraint_groupings(exome_ht)
    exome_ht = exome_ht.select('context', 'ref', 'alt', 'methylation_level', 'freq', 'pass_filters', *groupings)

    context_ht = context_ht.filter(hl.is_defined(context_ht.exome_coverage)).select(
        'context', 'ref', 'alt', 'methylation_level', *grouping)

    # Filter by allele count > 0, allele fraction < cutoff, coverage > 0
    # what does this code do? look at what you get from annotation
    af_cutoff = 0.001
    freq_index = exome_ht.freq_index_dict.collect()[0][dataset]
    def keep_criteria(ht):
        crit = (ht.freq[freq_index].AC > 0) & ht.pass_filters & (ht.coverage > 0)
        if impose_high_af_cutoff_upfront:
            crit &= (ht.freq[freq_index].AF <= af_cutoff)
        return crit
    exome_ht = exome_ht.filter(keep_criteria(exome_ht))

    return exome_ht


def build_coverage_model(po_coverage_ht_path):
    '''Build linear regression model for coverage, not sure of function'''
    coverage_ht = hl.read_table(po_coverage_ht_path)
    coverage_x_ht = hl.read_table(po_coverage_ht_path.replace('.ht', '_x.ht'))
    coverage_y_ht = hl.read_table(po_coverage_ht_path.replace('.ht', '_y.ht'))

    coverage_model, plateau_models = build_models(coverage_ht, args.trimers, True) # Have not read this function
    _, plateau_x_models = build_models(coverage_x_ht, args.trimers, True)
    _, plateau_y_models = build_models(coverage_y_ht, args.trimers, True)

    return coverage_model, plateau_models


def get_possible_variants(context_ht,mutation_ht, coverage_model, plateau_models, grouping, possible_file,
        half_cutoff = False, HIGH_COVERAGE_CUTOFF = 40, POPS = ('global', 'afr', 'amr', 'eas', 'nfe', 'sas')):
    '''Compute table of possible variants with needed properties'''
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
        possible_variants_ht: hl.Table,
        grouping: List,
        impose_high_af_cutoff_upfront: bool = True) -> hl.Table:
    '''Aggregate by grouping variables'''

    # Count observed variants by grouping
    ht = count_variants(exome_ht, additional_grouping=grouping, partition_hint=2000, force_grouping=True,
                        count_downsamplings=False, impose_high_af_cutoff_here=not impose_high_af_cutoff_upfront,
                        POPS = ('global', 'afr', 'amr', 'eas', 'nfe', 'sas'))
    # Match with expected frequency of variants by grouping
    ht = ht.join(possible_variants_ht, 'outer')

    # Aggregate all groupings by 
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

def run_tests(exomes_path):
    """Tests loading of autosome po table"""
    print('Test access to exomes and filtering')
    test_genes = ['ADRB1']
    test_intervals = get_gene_intervals(test_genes)
    ht_exomes_test = get_exomes(exomes_path, test_genes, test_intervals)
    ht_exomes_test.show()

    # Also need to test: storing in files and reloading, pre-processing, model usage, etc.

def main(args):
    # Set paths for data access based on command line parameters
    POPS = ('global', 'afr', 'amr', 'eas', 'nfe', 'sas') # set pops

    # New paths to public access bucket
    exomes_path = 'gs://gcp-public-data--gnomad/release/2.1/ht/exomes/gnomad.exomes.r2.1.sites.ht/'
    coverage_path_exomes = 'gs://gcp-public-data--gnomad/release/2.1/coverage/exomes/gnomad.exomes.r2.1.coverage.ht'
    context_path_vep = 'gs://gcp-public-data--gnomad/resources/context/grch37_context_vep_annotated.ht/'
    methylation_path = 'gs://gnomad-public/resources/grch37/methylation_sites/methylation.ht'


    root = 'gs://gcp-public-data--gnomad/papers/2019-flagship-lof/v1.0/model/'
    
    processed_genomes_ht_path = f'{root}/model/genomes_processed.ht'
    processed_exomes_ht_path = f'{root}/model/exomes_processed.ht'

    coverage_autosomes = f'{root}/prop_observed_by_coverage_no_common_pass_filtered_bins.ht'
    coverage_x = f'{root}/prop_observed_by_coverage_no_common_pass_filtered_bins_x.ht/'
    coverage_y = f'{root}/prop_observed_by_coverage_no_common_pass_filtered_bins_y.ht/' 

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
        # Currently tests whether it can read data from exomes sites table
        run_tests(exomes_path)

    if args.get_data:
        print('Getting data from Google cloud...')
        # Get intervals for chosen genes from Ensembl
        gene_intervals = get_gene_intervals()
        # Get relevant exomes, genomes, context and mutation rate data
        exomes_ht = get_exomes(exomes_path, gene_intervals)
        # genomes_ht = get_genomes(processed_genomes_ht_path)
        context_ht = get_context(context_path_vep,gene_intervals)
        mutation_rate_ht = get_mutation_rate(mutation_rate_ht_path,gene_intervals)

        # Get coverage table and build a model for methylation-dependent mutation rate
        coverage_model, plateau_models = build_coverage_model(po_coverage_ht_path)

        print('Pre-processing variant data...')
        full_exome_ht = prepare_ht(full_exome_ht, args.trimers) # what does this do? does it need the other tables?

        # filter into X, Y and autosomal regions
        exome_ht = full_exome_ht.filter(full_exome_ht.locus.in_autosome_or_par())
        exome_x_ht = hl.filter_intervals(full_exome_ht, [hl.parse_locus_interval('X')])
        exome_x_ht = exome_x_ht.filter(exome_x_ht.locus.in_x_nonpar())
        exome_y_ht = hl.filter_intervals(full_exome_ht, [hl.parse_locus_interval('Y')])
        exome_y_ht = exome_y_ht.filter(exome_y_ht.locus.in_y_nonpar())

        exome_ht.write(exomes_ht_path)
        exome_x_ht.write(exomes_x_ht_path)
        exome_y_ht.write(exomes_y_ht_path)

    if args.aggregate:
        print('Running aggregation of variants by grouping variables')
        # Set chosen groupings to aggregate on
        groupings = ['gene','annotation','modifier']

        # Filter variant table and annotate with grouping variables
        exome_ht = prepare_exomes(exome_ht, groupings) # do this for X and Y?
        
        # Calculate possible variants
        possible_variants_ht = get_possible_variants(groupings) # whatever is needed to do this

        # aggregate by chosen groupings & get proportion observed; write to file
        po_exome_ht, po_exome_x_ht, po_exome_y_ht = \
            [get_proportion_observed(ht, possible_variants_ht, groupings) # adjust to make up for lack of possible_variants \ 
                for ht in (exome_ht, exome_x_ht, exome_y_ht)]

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
    parser.add_argument('--get_data',help='Extract data from google cloud')
    parser.add_argument('--aggregate', help='Aggregate by grouping variables', action='store_true')
    parser.add_argument('--finalise',help='Calculate estimates',action='store_true')
    parser.add_argument('--summarise', help='Report summary stats', action='store_true')
    args = parser.parse_args()
    main(args)