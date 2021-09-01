#!/usr/bin/python

import argparse
from itertools import product
import pickle
import random
from typing import Dict, List, Optional, Set, Tuple, Any
import hail as hl
import pandas as pd
from gnomad_pipeline.utils import *
from data_loader import * 
from constraint_analysis import *
from finalise_results import *

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


def run_tests(paths):
    """Run tests of functions for randomly selected gene"""
    print('Test access to exomes and filtering')       
    test_intervals = get_gene_intervals()
    test_intervals = random.sample(test_intervals,k=1)
    print(f'Testing genes: {test_intervals}')
    # test_data = get_data(test_intervals)
    # print(f'Test data loaded: {test_data.summarise()}')
    # test_results_provisional = get_proportion_observed(test_data)
    # print(f'Test results filtered: {test_results_provisional.summarise()}')
    # test_results_final = finalize_dataset(test_results_provisional)
    # print(f'Test results finalised: {test_results_final.summarise()}')
    # print(test_results_final)


def main(args):
    # Set chosen grouping variables to aggregate on
    grouping = ['gene','annotation','modifier']
    POPS = ('global', 'afr', 'amr', 'eas', 'nfe', 'sas') # set pops

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
        mutation_rate_local_path = f'{root}/context.ht', #???
        po_coverage_local_path = f'{root}/prop_observed_by_coverage_no_common_pass_filtered_bins.ht',  
        possible_variants_ht_path = f'{root}/model/possible_data/possible_transcript_pop_{args.model}.ht',
        po_output_path = f'{root}/{{subdir}}/prop_observed_{{subdir}}.ht'.format(subdir=args.model),
        finalized_output_path = f'{root}/{{subdir}}/constraint_{{subdir}}.ht'.format(subdir=args.model),
        summary_output_path = f'{root}/{{subdir}}/constraint_final_{{subdir}}.ht'.format(subdir=args.model)
    )
    paths = {**gs_paths, **local_input_paths}

    
    hl.init(quiet=True)

    if args.test:
        run_tests(paths)

    if args.get_data:
        # Load gene intervals
        gene_intervals = get_gene_intervals()
        print('Getting data from Google Cloud...')
        data = get_data(paths, gene_intervals)
    
    if args.aggregate:
        print('Running aggregation of variants by grouping variables')
        aggregate_by_groupings(paths, data)
       
    if args.finalise:
        print('Joining tables and running aggregation by gene')
        # read PO hail tables for autosomes, X and Y chromosomes and join them
        finalize_dataset(paths, data)

    if args.summarise:
        print('Finalising summary stats')
        # write summary stats to output path
        summarise(paths)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--test', help='Run tests without actually requesting data',action='store_true')
    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--trimers', help='Use trimers instead of heptamers', action='store_true')
    parser.add_argument('--dataset', help='Which dataset to use (one of gnomad, non_neuro, non_cancer, controls)', default='gnomad')
    parser.add_argument('--model', help='Which model to apply (one of "standard", "syn_canonical", or "worst_csq" for now) - warning not implemented', default='standard')
    parser.add_argument('--skip_af_filter_upfront', help='Skip AF filter up front (to be applied later to ensure that it is not affecting population-specific constraint): not generally recommended', action='store_true')
    parser.add_argument('--get_data',help='Extract data from google cloud',action='store_true')
    parser.add_argument('--aggregate', help='Aggregate by grouping variables', action='store_true')
    parser.add_argument('--finalise',help='Calculate estimates',action='store_true')
    parser.add_argument('--summarise', help='Report summary stats', action='store_true')
    args = parser.parse_args()
    main(args)