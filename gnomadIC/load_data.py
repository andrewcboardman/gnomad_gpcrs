import argparse
import pickle
import os

from numpy.lib import utils
import hail as hl
import pandas as pd
from typing import Dict, List, Optional, Set, Tuple, Any
from .utils import utils

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
        return df_intervals['interval'].tolist()
    else:
        return df_intervals['interval'].tolist()


def get_table(path, columns, intervals):
   # Path to full exome table
    full_ht = hl.read_table(path)
    # Select relevant columns to avoid getting too much data 
    full_ht = full_ht.select(*columns)
    # Filter by gene
    ht = hl.filter_intervals(full_ht, intervals)
    # Prepare exomes
    ht = utils.prepare_ht(ht, trimer=True)
    return ht 


def annotate_exomes(exomes_ht: hl.Table, context_ht: hl.Table,
                      overwrite: bool = False) -> None:
    context_ht = context_ht.annotate(
        vep=context_ht.vep.drop('colocated_variants')
        )
    exomes_ht = exomes_ht.annotate(**context_ht[exomes_ht.key], pass_filters=hl.len(exomes_ht.filters) == 0)
    return exomes_ht


def filter_exomes(exome_ht, af_cutoff=0.001, dataset: str = 'gnomad', 
        impose_high_af_cutoff_upfront: bool = True):
    # Filter by allele count > 0, allele fraction < cutoff, coverage > 0

    freq_index = exome_ht.freq_index_dict.collect()[0][dataset]
    def keep_criteria(ht):
        crit = (ht.freq[freq_index].AC > 0) & ht.pass_filters & (ht.coverage > 0)
        if impose_high_af_cutoff_upfront:
            crit &= (ht.freq[freq_index].AF <= af_cutoff)
        return crit
    exome_ht = exome_ht.filter(keep_criteria(exome_ht))
    return exome_ht


def prepare_table_vep(ht, columns):
    ht = utils.add_most_severe_csq_to_tc_within_ht(ht)
    ht = ht.transmute(transcript_consequences=ht.vep.transcript_consequences)  
    ht = ht.explode(ht.transcript_consequences)
    ht, _ = utils.annotate_constraint_groupings(ht)
    return ht


def prepare_exomes_and_context(
        exome_ht: hl.Table,
        context_ht: hl.Table,
        
    ) -> hl.Table:
    '''Prepare exome mutation data for modelling '''
    # For standard analysis, assume the worst predicted transcript consequence takes place
    exome_ht = prepare_table_vep(exome_ht)
    context_ht = prepare_table_vep(context_ht)


    # Annotate variants with grouping variables
    exome_ht, grouping = utils.annotate_constraint_groupings(exome_ht)
    exome_ht = exome_ht.select('context', 'ref', 'alt', 'methylation_level', 'freq', 'pass_filters', *grouping)    
    context_ht, grouping = utils.annotate_constraint_groupings(context_ht)
    context_ht = context_ht.filter(hl.is_defined(context_ht.exome_coverage)).select(
        'context', 'ref', 'alt', 'methylation_level', *grouping)

    exome_ht = filter_exomes(exome_ht)

    return exome_ht, context_ht, grouping


def split_table(full_ht):
    # filter into X, Y and autosomal regions for separate aggregation
    autosome_ht = full_ht.filter(full_ht.locus.in_autosome_or_par())
    x_ht = hl.filter_intervals(full_ht, [hl.parse_locus_interval('X')])
    x_ht = x_ht.filter(x_ht.locus.in_x_nonpar())
    y_ht = hl.filter_intervals(full_ht, [hl.parse_locus_interval('Y')])
    y_ht = y_ht.filter(y_ht.locus.in_y_nonpar())
    return (autosome_ht, x_ht, y_ht)


def get_data(paths, gene_intervals, overwrite=True, trimers=False):
    '''
    This is the new master function for loading all necessary data for constraint analysis on the given genes
    Paths are passed in from the main program. 
    The exomes and context data should always be downloaded as new gene intervals are passed in. 
    The mutation rate by methylation and proportion observed by coverage tables are stored locally.
    They should be downloaded if not present but the control flow to do this isn't yet implemented 

    Changes that could be made 
    - refactor some of the functions in here to serve exomes & context 
    - save whole exomes and context tables to avoid splitting data into files that don't need to exist
    '''

    # Get relevant exomes & context data
    full_exome_ht = get_table(
        paths['exomes_path'], 
        (
            'vep',
            'context',
            'methylation',
            'coverage',
            'freq',
            'filters'
        ),
        gene_intervals
    )
    full_context_ht = get_table(
        paths['context_path'], 
        (
            'vep',
            'context',
            'methylation',
            'coverage'
        ),
        gene_intervals
    )

    # Add context to exomes and filter
    annotated_exome_ht = annotate_exomes(
        full_exome_ht, 
        full_context_ht
    )

    # Prepare tables including all necessary information
    prepped_exome_ht = prepare_table_vep(
        annotated_exome_ht,
        (
            'context', 
            'ref', 
            'alt', 
            'methylation_level', 
            'annotation',
            'modifier',
            'transcript',
            'gene',
            'canonical',
            'coverage',
            'freq', 
            'pass_filters'
        )
    )

    prepped_context_ht = prepare_table_vep(
        full_context_ht,
        (
            'context', 
            'ref', 
            'alt', 
            'methylation_level', 
            'annotation',
            'modifier',
            'transcript',
            'gene',
            'canonical', 
            'coverage'
        )
    )

    # filter exome variants by frequency
    filtered_exome_ht = filter_exomes(prepped_exome_ht)

    # Write to file
    filtered_exome_ht.write(paths['exomes_local_path'], overwrite=overwrite)
    prepped_context_ht.write(paths['context_local_path'], overwrite=overwrite)
     
    data = load_data_to_aggregate(paths)
    return data


def load_models(paths):
    # Get table for mutation rate if it doesn't exist  
    if os.path.isdir(paths['mutation_rate_local_path']):
        mutation_rate_ht = hl.read_table(paths['mutation_rate_local_path'])
    else:
        mutation_rate_ht = hl.read_table(paths['mutation_rate_path']).select('mu_snp')
        mutation_rate_ht.write(paths['mutation_rate_local_path'])
    
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

    models = {
        'mutation_rate_ht': mutation_rate_ht,
        'coverage_model': coverage_model,
        'plateau_model': plateau_models
    }
    return models

def load_data_to_aggregate(paths):
    # Load data
    exome_ht = hl.read_table(paths['exomes_local_path'])
    context_ht = hl.read_table(paths['context_local_path'])
    exome_auto_ht, exome_x_ht, exome_y_ht = split_table(exome_ht)
    context_auto_ht, context_x_ht, context_y_ht = split_table(context_ht)

    # filter into X, Y and autosomal regions
    data = dict(zip(
        [
            'exome_ht','exome_x_ht','exome_y_ht',
            'context_ht','context_x_ht','context_y_ht',
        ],
        [
            exome_auto_ht, exome_x_ht, exome_y_ht,
            context_auto_ht, context_x_ht, context_y_ht,
        ]
    ))

    # Adjust this to allow custom grouping analysis
    data['grouping'] = [
        'annotation',
        'modifier',
        'transcript', 
        'gene',
        'canonical',
        'coverage'
        ]
    
    models = load_models(paths)
    data.update(models)
    return data


def load_data_to_estimate(paths):
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
