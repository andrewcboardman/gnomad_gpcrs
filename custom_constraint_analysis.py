#!/usr/bin/python

import argparse
import hail as hl
import gnomadIC

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


def run_tasks(tasks, paths, test = False):
    '''Runs all requested tasks in specified path'''
    data = {}
    
    if 'load_data' in tasks:
        # Load gene intervals
        gene_intervals = gnomadIC.get_gene_intervals(test)
        # If in test mode only load 1 gene
        print('Getting data from Google Cloud...')
        # Load observed mutation data, possible mutation data, and associated data
        data = gnomadIC.load_data(paths, gene_intervals)
        print('Data loaded successfully!')
    
    if 'aggregate' in tasks:
        if len(data) == 0:
            print('Loading data...')
            data = gnomadIC.load_data_to_aggregate(paths)
        print('Running aggregation of variants by grouping variables...')
        data = gnomadIC.aggregate(paths, data)
        print('Aggregated variants successfully!')
       
    if 'estimate' in tasks:
        if len(data) == 0:
            print('Loading data...')
            data = gnomadIC.load_data_to_estimate(paths)
        print('Joining tables and running aggregation by gene')
        data = gnomadIC.estimate(paths, data)


def main(args):
    '''Controls whether to setup in test mode or not, and generates a run ID if not in test mode'''
    # Initialise Hail, setting output 
    hl.init(log='hail_logs/test_load_data.log', quiet=True)

    # Setup paths
    if args.test:
        print('Running in test mode: you can relax, sit back and enjoy the ride')
        paths = gnomadIC.setup_paths('test')
    else:
        run_ID = f'{args.dataset}_{args.model}'
        print(f'Running without test mode active: THIS IS NOT A DRILL. \n Run ID: {run_ID}')        
        paths = gnomadIC.setup_paths(run_ID)
    
    # Run chosen tasks
    run_tasks(args.tasks, paths = paths, test=args.test)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--test', help='Run tests',action='store_true')
    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--dataset', help='Which dataset to use (one of gnomad, non_neuro, non_cancer, controls)', default='gnomad')
    parser.add_argument('--model', help='Which model to apply (one of "standard", "syn_canonical", or "worst_csq" for now) - warning not implemented', default='standard')
    parser.add_argument('--tasks', nargs='+', help='Which tasks to perform (select from load_data, aggregate, estimate)')
    args = parser.parse_args()
    main(args)