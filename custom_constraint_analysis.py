#!/usr/bin/python
import os
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

def setup_paths(run_ID):
    root = './data'
    # Google storage paths 
    gs_paths = dict(
        exomes_path = 'gs://gcp-public-data--gnomad/papers/2019-flagship-lof/v1.0/model/exomes_processed.ht/',
        context_path = 'gs://gcp-public-data--gnomad/resources/context/grch37_context_vep_annotated.ht/',
        mutation_rate_path = 'gs://gcp-public-data--gnomad/papers/2019-flagship-lof/v1.0/model/mutation_rate_methylation_bins.ht',
        po_coverage_path = 'gs://gcp-public-data--gnomad/papers/2019-flagship-lof/v1.0/model/prop_observed_by_coverage_no_common_pass_filtered_bins.ht'
    )
    # Local paths
    output_subdir = f'{root}/{run_ID}'
    if not os.path.isdir(output_subdir):
        os.mkdir(output_subdir)
    local_paths = dict(
        # models - shared between runs
        mutation_rate_local_path = f'{root}/mutation_rate_methylation_bins.ht',
        po_coverage_local_path = f'{root}/prop_observed_by_coverage_no_common_pass_filtered_bins.ht',
        coverage_models_local_path = f'{root}/coverage_models.pkl',
        # outputs - specific to run
        exomes_local_path = f'{output_subdir}/exomes.ht',
        context_local_path = f'{output_subdir}/context.ht',        
        possible_variants_ht_path = f'{output_subdir}/possible_transcript_pop.ht',
        po_output_path = f'{output_subdir}/prop_observed.ht',
        finalized_output_path = f'{output_subdir}/constraint.ht',
        summary_output_path = f'{output_subdir}/constraint_final.ht'
    )
    paths = {**gs_paths, **local_paths}

    return paths


def run_tasks(tasks, paths, model, test = False):
    '''Runs all requested tasks in specified path'''
    data = {}
    
    if 'load_data' in tasks:
        # Load gene intervals
        gene_intervals = gnomadIC.get_gene_intervals(test)
        # If in test mode only load 1 gene
        print('Getting data from Google Cloud...')
        # Load observed mutation data, possible mutation data, and annotations
        data = gnomadIC.get_data(paths, gene_intervals, model)
        print('Data loaded successfully!')
  
    if 'aggregate' in tasks:
        print('Running aggregation of variants by grouping variables...')
        data = gnomadIC.aggregate(paths, data, model)
        print('Aggregated variants successfully!')



def main(args):
    '''Controls whether to setup in test mode or not, and generates a run ID if not in test mode'''
    # Initialise Hail, setting output 
    hl.init(log='hail_logs/test_load_data.log', quiet=args.quiet)

    # Setup paths
    if args.test:
        run_ID = 'test'
        print('Running in test mode: Relax, sit back and enjoy the ride')
    else:
        run_ID = f'{args.dataset}_{args.model}'  
        print(f'Running without test mode active: THIS IS NOT A DRILL. \n Run ID: {run_ID}')
    paths = setup_paths(run_ID)
    
    # Run chosen tasks
    run_tasks(args.tasks, paths = paths, model = args.model, test=args.test)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--test', help='Run tests',action='store_true',default=False)
    parser.add_argument('-q','--quiet',help='Run in quiet mode',action='store_true',default=False)
    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--dataset', help='Which dataset to use (one of gnomad, non_neuro, non_cancer, controls)', default='gnomad')
    parser.add_argument('--model', nargs= '+', help='Which model to apply (one of "standard", "syn_canonical", or "worst_csq" for now) - warning not implemented', default='standard')
    parser.add_argument('--tasks', nargs='+', help='Which tasks to perform (select from load_data, aggregate)')
    args = parser.parse_args()
    main(args)