import argparse
import pickle
import os

from numpy.lib import utils
import hail as hl
import pandas as pd
from typing import Dict, List, Optional, Set, Tuple, Any
from .utils import utils

def get_mutation_annotations(model):    
    # Get custom annotations for mutations from file
    df_annotations = pd.read_csv('data/mutation_annotations.csv')
    if model == 'standard':
        df_annotations = df_annotations['gene','transcript','pos','ref_aa','alt_aa','standard_csq']
        return hl.Table.from_pandas(df_annotations, key=['locus','alleles'])
    else:
        raise NotImplementedError('Model not implemented!')


def get_table(path, intervals, model, additional_fields = [], trimer=True):
   # Path to full exome table
    full_ht = hl.read_table(path)
    # Select relevant fields to avoid getting too much data
    fields = ['vep','context', 'methylation','coverage'] + additional_fields
    full_ht = full_ht.select(*fields)
    # Filter by gene
    ht = hl.filter_intervals(full_ht, intervals)
    # Prepare exomes
    ht = utils.prepare_ht(ht, trimer=trimer)
    # Extract relevant parts of VEP struct and set as groupings for annotation join
    ht = utils.add_most_severe_csq_to_tc_within_ht(ht)
    ht = ht.transmute(transcript_consequences=ht.vep.transcript_consequences) 
    ht = ht.explode(ht.transcript_consequences)
    ht, groupings = utils.annotate_constraint_groupings(ht, model)
    # Extract coverage as exome coverage
    ht = ht.transmute(coverage=ht.exome_coverage)

    # Select fields for modelling
    final_fields = ['ref', 'alt', 'context', 'methylation_level', 'coverage'] + \
        groupings + additional_fields
    ht = ht.select(*final_fields)
    return ht, groupings


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


def get_data(paths, gene_intervals, model, overwrite=True, trimer=True):
    '''
    This is the new master function for loading all necessary data for constraint analysis on the given genes
    Paths are passed in from the main program. 
    The exomes and context data should always be downloaded as new gene intervals are passed in. 
    The mutation rate by methylation and proportion observed by coverage tables are stored locally.
    They should be downloaded if not present but the control flow to do this isn't yet implemented 
    '''
    # Prepare context table by filtering on gene intervals and selecting correct VEP annotations
    context_ht, groupings = get_table(paths['context_path'], gene_intervals, model, trimer=trimer)

    # Get exomes data by filtering on gene intervals & selecting correct VEP annotations
    exome_ht, _ = get_table(paths['exomes_path'], gene_intervals, model, additional_fields= ['freq', 'filters'], trimer=trimer)

    # Do extra filtering of exomes
    exome_ht = exome_ht.annotate(pass_filters = hl.len(exome_ht.filters)==0)
    exome_ht = filter_exomes(exome_ht)  

    # Write to file
    exome_ht.write(paths['exomes_local_path'], overwrite=overwrite)
    context_ht.write(paths['context_local_path'], overwrite=overwrite)

    data = {
        'exome_ht': exome_ht,
        'context_ht':context_ht,
        'groupings':groupings
    }

    return data