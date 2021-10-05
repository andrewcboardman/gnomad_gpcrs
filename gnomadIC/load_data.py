import argparse
import pickle
import os
import hail as hl
import pandas as pd
from typing import Dict, List, Optional, Set, Tuple, Any
import utils.utils as utils

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

def get_gene_intervals(test=False):
    
    # Get Ensembl gene intervals from file
    df_intervals = pd.read_csv('data/Ensembl_Grch37_gpcr_genome_locations.csv')
    df_intervals = df_intervals[['HGNC symbol','Grch37 symbol','Grch37 chromosome','Grch37 start bp','Grch37 end bp']]
    df_intervals['locus_interval_txt'] = df_intervals['Grch37 chromosome'] + ':'  + \
        + df_intervals['Grch37 start bp'].map(str) + '-' + df_intervals['Grch37 end bp'].map(str)
    df_intervals['interval'] = df_intervals['locus_interval_txt'].map(hl.parse_locus_interval)

    if test:
        df_intervals = df_intervals.sample(n=1,random_state=0)
        print(f"{str(df_intervals['HGNC symbol'].values[0])} chosen as test gene")
        return df_intervals
    else:
        return df_intervals



def get_data(paths, gene_intervals, overwrite=True, trimers=False):
    '''
    This is the new master function for loading all necessary data for constraint analysis on the given genes
    Paths are passed in from the main program. 
    The exomes and context data should always be downloaded as new gene intervals are passed in. 
    The mutation rate by methylation and proportion observed by coverage tables are stored locally.
    They should be downloaded if not present but the control flow to do this isn't yet implemented 
    '''

    data = {}

    # Get relevant exomes & context data
    full_exome_ht = get_exomes(paths['exomes_path'], gene_intervals)    
    full_context_ht = get_context(paths['context_path'], gene_intervals)

    # Preprocess by adding context
    annotated_exome_ht = pre_process_exomes(full_exome_ht, full_context_ht)

    # Prepare tables with constraint groupings
    prepped_exome_ht, prepped_context_ht, grouping = prepare_exomes_and_context(annotated_exome_ht,full_context_ht)

    # Split tables into autosomes, X and Y chromosomes and write to file
    exome_ht, exome_x_ht, exome_y_ht = filter_exomes(prepped_exome_ht, paths['exomes_local_path'], overwrite)
    context_ht,context_x_ht, context_y_ht = filter_context(prepped_context_ht, paths['context_local_path'], overwrite)

    # Add tables to data dictionary
    data.update(dict(zip(
        ('exome_ht', 'exome_x_ht', 'exome_y_ht', 'context_ht', 'context_x_ht', 'context_y_ht','grouping'),
        (exome_ht, exome_x_ht, exome_y_ht, context_ht, context_x_ht, context_y_ht,grouping)
    )))
    
    # Get table for mutation rate if it doesn't exist  
    if os.path.isdir(paths['mutation_rate_local_path']):
        data['mutation_rate_ht'] = hl.read_table(paths['mutation_rate_local_path'])
    else:
        mutation_rate_ht = hl.read_table(paths['mutation_rate_path']).select('mu_snp')
        mutation_rate_ht.write(paths['mutation_rate_local_path'],overwrite=overwrite)
        data['mutation_rate_ht'] = mutation_rate_ht
    
    # Get coverage models if they don't exist
    if os.path.isfile(paths['coverage_models_local_path']):
        # Load local coverage models
        with open('data/coverage_models.pkl','rb') as fid:
            coverage_model, plateau_models = pickle.load(fid)
    else:
        if os.path.isdir(paths['po_coverage_local_path']):
            # Load local coverage table
           coverage_ht = hl.read_table(paths['po_coverage_local_path'])
        else:
            # Download coverage table
            coverage_ht = hl.read_table(paths['po_coverage_path'])
            coverage_ht.write(paths['po_coverage_local_path'])
        # Build models from coverage table
        coverage_model, plateau_models = utils.build_models(coverage_ht, trimers=True) 
        with open('data/coverage_models.pkl','wb') as fid:
            pickle.dump((coverage_model,plateau_models),fid)  

    data.update(dict(zip(
        ('coverage_model', 'plateau_models'),
        (coverage_model, plateau_models)
    )))

    return data


def get_exomes(exomes_path, gene_intervals):
    # Path to full exome table
    full_exome_ht = hl.read_table(exomes_path)
    # Select relevant columns to avoid getting too much data 
    full_exome_ht = full_exome_ht.select('freq','vep','context','methylation','coverage','filters'
    )
    # Filter by gene
    exome_ht = hl.filter_intervals(full_exome_ht, gene_intervals)
    # Prepare exomes
    exome_ht = utils.prepare_ht(exome_ht, trimer=True)
    return exome_ht


def get_context(context_ht_path, gene_intervals):
    full_context_ht = hl.read_table(context_ht_path)
    full_context_ht = full_context_ht.select('vep','context','methylation','coverage')
    # Filter this table in similar way to exomes ht
    selected_context_ht = hl.filter_intervals(full_context_ht,gene_intervals)
    selected_context_ht = utils.prepare_ht(selected_context_ht, trimer=True)
    return selected_context_ht


def pre_process_exomes(exomes_ht: hl.Table, context_ht: hl.Table,
                      overwrite: bool = False) -> None:
    context_ht = context_ht.annotate(vep=context_ht.vep.drop('colocated_variants'))
    return exomes_ht.annotate(**context_ht[exomes_ht.key], pass_filters=hl.len(exomes_ht.filters) == 0)#.write(output_ht_path, overwrite)


def prepare_exomes_and_context(
        exome_ht: hl.Table,
        context_ht: hl.Table,
        dataset: str = 'gnomad', 
        impose_high_af_cutoff_upfront: bool = True
    ) -> hl.Table:
    '''Prepare exome mutation data for modelling '''
    # For standard analysis, assume the worst predicted transcript consequence takes place
    exome_ht = utils.add_most_severe_csq_to_tc_within_ht(exome_ht)
    exome_ht = exome_ht.transmute(transcript_consequences=exome_ht.vep.transcript_consequences)  
    exome_ht = exome_ht.explode(exome_ht.transcript_consequences)
    # Annotate variants with grouping variables
    exome_ht, grouping = utils.annotate_constraint_groupings(exome_ht)
    exome_ht = exome_ht.select('context', 'ref', 'alt', 'methylation_level', 'freq', 'pass_filters', *grouping)

    # annotate and filter context table in same way. Exome coverage or genome coverage?
    context_ht = utils.add_most_severe_csq_to_tc_within_ht(context_ht)
    context_ht = context_ht.transmute(transcript_consequences=context_ht.vep.transcript_consequences)  
    context_ht = context_ht.explode(context_ht.transcript_consequences)
    context_ht, grouping = utils.annotate_constraint_groupings(context_ht)
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


def filter_exomes(full_exome_ht, exomes_local_path, overwrite):
    # filter into X, Y and autosomal regions
    exome_ht = full_exome_ht.filter(full_exome_ht.locus.in_autosome_or_par())
    exome_x_ht = hl.filter_intervals(full_exome_ht, [hl.parse_locus_interval('X')])
    exome_x_ht = exome_x_ht.filter(exome_x_ht.locus.in_x_nonpar())
    exome_y_ht = hl.filter_intervals(full_exome_ht, [hl.parse_locus_interval('Y')])
    exome_y_ht = exome_y_ht.filter(exome_y_ht.locus.in_y_nonpar())
    exome_ht.write(exomes_local_path,overwrite=overwrite)
    exome_x_ht.write(exomes_local_path.replace('.ht', '_x.ht'),overwrite=overwrite)
    exome_y_ht.write(exomes_local_path.replace('.ht', '_y.ht'),overwrite=overwrite)
    return exome_ht, exome_x_ht, exome_y_ht


def filter_context(full_context_ht, context_local_path, overwrite):
    context_ht = full_context_ht.filter(full_context_ht.locus.in_autosome_or_par())
    context_x_ht = hl.filter_intervals(full_context_ht, [hl.parse_locus_interval('X')])
    context_x_ht = context_x_ht.filter(context_x_ht.locus.in_x_nonpar())
    context_y_ht = hl.filter_intervals(full_context_ht, [hl.parse_locus_interval('Y')])
    context_y_ht = context_y_ht.filter(context_y_ht.locus.in_y_nonpar())
    context_ht.write(context_local_path,overwrite=overwrite)
    context_x_ht.write(context_local_path.replace('.ht', '_x.ht'),overwrite=overwrite)
    context_y_ht.write(context_local_path.replace('.ht', '_y.ht'),overwrite=overwrite)
    return context_ht, context_x_ht, context_y_ht


def load_data_to_aggregate(paths):
    exomes_local_path = paths['exomes_local_path']
    context_local_path = paths['context_local_path']
    mutation_rate_local_path = paths['mutation_rate_local_path']

    data = dict(
        zip(
            ['exome_ht','exome_x_ht','exome_y_ht',
            'context_ht','context_x_ht','context_y_ht',
            'mutation_rate_ht'],
            [hl.read_table(path) for path in \
                (
                    exomes_local_path, 
                    exomes_local_path.replace('.ht', '_x.ht'), 
                    exomes_local_path.replace('.ht', '_y.ht'),
                    context_local_path, 
                    context_local_path.replace('.ht', '_x.ht'), 
                    context_local_path.replace('.ht', '_y.ht'),
                    mutation_rate_local_path
                )
            ]
        )
    )

    # Get coverage models
    with open('data/coverage_models.pkl','rb') as fid:
        coverage_model, plateau_models = pickle.load(fid)
    data.update(dict(zip(
        ('coverage_model', 'plateau_models'),
        (coverage_model, plateau_models)
    )))
    # Adjust this to allow custom grouping analysis
    data['grouping'] = ['annotation','modifier','transcript', 'gene','canonical', 'coverage']
    return data


def load_data_to_finalise(paths):
    data = dict(zip(
        ('po_ht','po_x_ht','po_y_ht'),
        (hl.read_table(x) for x in (
                paths['po_output_path'], 
                paths['po_output_path'].replace('.ht','_x.ht'),
                paths['po_output_path'].replace('.ht','_y.ht')
            )
        )
    ))
    return data


def print_data_loader_test_summary(test_data):
    return test_data['exome_ht'].summarize()


def run_data_loader_test(print_summary=True):
    hl.init(log='hail_logs/test_load_data.log')
    paths = setup_paths('test')
    test_intervals = get_gene_intervals(test=True)
    #print(list(test_intervals['interval']))
    test_data = get_data(paths, list(test_intervals['interval']))
    if print_summary:
        print_data_loader_test_summary(test_data)
    return paths, test_data


def main(args):
    if args.test:
        print('Running tests...')
        run_data_loader_test(print_summary=True)
    else:
        print('Please run this script from custom_constraint_analysis.py')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--test', help='Run tests',action='store_true')
    args = parser.parse_args()
    main(args)