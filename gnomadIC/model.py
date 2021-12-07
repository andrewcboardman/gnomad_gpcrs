
from itertools import product
import argparse
import hail as hl
from typing import List
from .utils import utils
from .data import *
import os

HIGH_COVERAGE_CUTOFF = 40
POPS = ('global', 'afr', 'amr', 'eas', 'nfe', 'sas')
  

def split_table(path, full_ht):
    # filter into X, Y and autosomal regions for separate aggregation
    auto_path = path.replace('.ht','_auto.ht')
    if not os.path.isdir(auto_path):
        auto_ht = full_ht.filter(full_ht.locus.in_autosome_or_par())
        auto_ht.write(auto_path)
    auto_ht = hl.read_table(auto_path)

    x_path = path.replace('.ht','_x.ht')
    if not os.path.isdir(x_path):
        x_ht = hl.filter_intervals(full_ht, [hl.parse_locus_interval('X')])
        x_ht = x_ht.filter(x_ht.locus.in_x_nonpar())
        x_ht.write(x_path)
    x_ht = hl.read_table(x_path)

    y_path = path.replace('.ht','_y.ht')
    if not os.path.isdir(y_path):
        y_ht = hl.filter_intervals(full_ht, [hl.parse_locus_interval('Y')])
        y_ht = y_ht.filter(y_ht.locus.in_y_nonpar())
        y_ht.write(y_path)
    y_ht = hl.read_table(y_path)

    return {'auto':auto_ht, 'x':x_ht, 'y': y_ht}


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
        with open(paths['coverage_models_local_path'],'rb') as fid:
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
        with open(paths['coverage_models_local_path'],'wb') as fid:
            pickle.dump((coverage_model,plateau_models),fid)  

    models = {
        'mutation_rate_ht': mutation_rate_ht,
        'coverage_model': coverage_model,
        'plateau_models': plateau_models
    }
    return models

def preprocess(paths, data, grouping, model):
    # Add extra annotations based on VEP 
    # annotations_ht = hl.get_annotations(model)
    # data['context_ht'] = data['context_ht'].annotate(**annotations_ht[data['context_ht'].hgvsp])
    # data['exome_ht'] = data['exome_ht'].annotate(**annotations_ht[data['exome_ht'].hgvsp])
    # Split data; load models; modify grouping
    data.update({
        'exomes': split_table(paths['exomes_local_path'],data['exome_ht']),
        'context': split_table(paths['context_local_path'],data['context_ht']),
        'models': load_models(paths),
        'grouping': grouping
    })
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


def get_expected_variants(ht, models, grouping, possible_path, pops=False):
    '''Compute table of possible variants with needed properties'''
    # Apply model to calculated expected variants
    print('Calculating expected variants')
    ht = ht.annotate(variant_count=hl.literal(1))
    ht = utils.annotate_expected_mutations(ht, models['mutation_rate_ht'], models['plateau_models'], models['coverage_model'], pops = pops)

    # Count possible variants by context, ref, alt & grouping - need to expand list of groupings to keep this from destroying information
    agg_expr = {
        'expected_variants': hl.agg.sum(ht.expected_variants),
        'possible_variants': hl.agg.sum(ht.possible_variants),
        'adjusted_mutation_rate': hl.agg.sum(ht.adjusted_mutation_rate),      
        'raw_mutation_rate': hl.agg.sum(ht.mu)
    }
    ht = ht.group_by(*grouping).aggregate(**agg_expr)
    ht.write(possible_path, overwrite=True)
    return ht


def get_proportion_observed(
        exome_ht: hl.Table, 
        expected_variants_ht: hl.Table,
        grouping: List,
        proportion_variants_observed_ht_path,
        impose_high_af_cutoff_upfront: bool = True,
        pops = False, overwrite=True) -> hl.Table:
    '''Aggregate by grouping variables'''

    # Count observed variants by grouping - expand grouping to include hgvsp to prevent information loss
    agg_expr = {
        'observed_variants': hl.agg.count()
    }
    observed_variants_ht = exome_ht.group_by(*grouping).aggregate(**agg_expr)
    obs_path = proportion_variants_observed_ht_path.replace('.ht','_raw.ht')
    observed_variants_ht.write(obs_path,overwrite=overwrite)

    # Merge observed variants with expected variants
    observed_variants_ht = hl.read_table(obs_path)
    # Why does this join lead to variable duplication? Think it's because of missing values
    observed_variants_ht = observed_variants_ht.join(expected_variants_ht, 'outer')
    observed_variants_ht.write(proportion_variants_observed_ht_path,overwrite=overwrite)
    return observed_variants_ht


def model(paths, data, model):
    '''
    This is the new master function for performing constraint analysis
    Possible variants for populations currently switched off
    '''
    # Get data if not given 
    if not data:
        print('Loading data...')
        data = {
            'exome_ht': hl.read_table(paths['exomes_local_path']),
            'context_ht': hl.read_table(paths['context_local_path'])
        }
    # Pre-process data
    grouping = [
        'annotation',
        'modifier',
        'transcript', 
        'gene',
        'canonical'
        ]
    data = preprocess(paths, data, grouping, model)


    # Loop over autosomes, x y: aggregate by chosen groupings & get proportion observed
    tables = ('auto','x','y')
    data['prop_observed'] = {}

    for table in tables:
        expected_variants_ht = get_expected_variants(
            data['context'][table], 
            data['models'],
            data['grouping'],
            paths['possible_variants_ht_path'].replace('.ht',f'_{table}.ht'),
            pops=False
        )

        prop_observed_ht = get_proportion_observed(
            data['exomes'][table],
            expected_variants_ht,
            data['grouping'],
            paths['po_output_path'].replace('.ht',f'_{table}.ht'), 
            overwrite=True)
        data['prop_observed'][table] = prop_observed_ht

    # Take union of answers and write to file
    data['prop_observed_ht'] = (
                data['prop_observed']['auto']
                .union(data['prop_observed']['x'])
                .union(data['prop_observed']['y'])
                )
    data['prop_observed_ht'].write(paths['po_output_path'], overwrite=True)

    
    return data