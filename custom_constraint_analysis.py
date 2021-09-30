#!/usr/bin/python

import argparse
import pickle
import random
import hail as hl
import pandas as pd
from itertools import product
from typing import Dict, List, Optional, Set, Tuple, Any
from gnomad_pipeline.utils import *
from setup import *
from data_loader import * 
from constraint_analysis import *
from finalise_results import *
from summarise_results import *

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

def run_tests(paths):
    """Run tests of functions for randomly selected gene"""
    print('Test access to exomes and filtering')       
    test_intervals = get_gene_intervals()
    test_intervals = test_intervals.sample(n=1,random_state=0)
    symbols = test_intervals['HGNC symbol'].unique()
    print(f'Testing genes: {symbols}')
    test_data = get_data(paths, list(test_intervals['interval']))
    #print(f'Test data loaded: {test_data.summarise()}')
    #test_data = aggregate_by_groupings(paths, test_data)
    #print(f'Test results filtered: {test_results_provisional.summarise()}')
    #test_data = finalize_dataset(paths, test_data)
    #test_summary = summarise_test(test_data)
    # print(f'Test results finalised: {test_results_final.summarise()}')
    # print(test_results_final)


def main(args):
    # Set chosen grouping variables to aggregate on
    POPS = ('global', 'afr', 'amr', 'eas', 'nfe', 'sas') # set pops

    # setup paths
    paths = setup_paths(args.model)
    # setup data 
    data = {}
    
    hl.init(quiet=True)

    if args.test:
        run_tests(paths)

    if args.get_data:
        # Load gene intervals
        gene_intervals = get_gene_intervals()
        print('Getting data from Google Cloud...')
        # Load observed mutation data, possible mutation data, and associated data
        data = get_data(paths, gene_intervals)
    
    if args.aggregate:
        if len(data) == 0:
            print('Loading data...')
            data = load_data_to_aggregate(paths)
        print('Running aggregation of variants by grouping variables...')
        aggregate_by_groupings(paths, data)
       
    if args.finalise:
        if len(data) == 0:
            print('Loading data...')
            data = load_data_to_finalise(paths)
        print('Joining tables and running aggregation by gene')
        # read PO hail tables for autosomes, X and Y chromosomes and join them
        # Aggregate by final markers and calculate type-level metrics
        data = finalize_dataset(paths, data)       

    if args.summarise:
        if len(data) == 0:
            data = load_data_to_summarise(paths)
        print('Calculating summary stats')
        # write summary stats to output path
        data = summarise(paths, data)


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