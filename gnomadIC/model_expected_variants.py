
from itertools import product
import argparse
import hail as hl
from typing import List
from .utils import utils
from .load_data import *

HIGH_COVERAGE_CUTOFF = 40
POPS = ('global', 'afr', 'amr', 'eas', 'nfe', 'sas')
  

def get_expected_variants(context_ht, mutation_rate_ht, coverage_model, plateau_models, grouping, possible_file):
    '''Compute table of possible variants with needed properties'''
    # Count possible variants by groupings - need to expand list of groupings to keep this from destroying information
    ht = utils.count_variants(context_ht, additional_grouping=grouping, partition_hint=2000, force_grouping=True)

    # Apply model to calculated expected variants
    ht = utils.annotate_expected_mutations(ht, mutation_rate_ht, plateau_models, coverage_model)
   
    # write to file
    ht.write(possible_file, True)

    return ht


def get_proportion_observed(
        exome_ht: hl.Table, 
        expected_variants_ht: hl.Table,
        grouping: List,
        proportion_variants_observed_ht_path,
        impose_high_af_cutoff_upfront: bool = True,
        pops = False, overwrite=True) -> hl.Table:
    '''Aggregate by grouping variables'''

    # Count observed variants by context, ref, alt + grouping
    observed_variants_ht = utils.count_variants(exome_ht, additional_grouping=grouping, partition_hint=2000, force_grouping=True,
                        count_downsamplings=[], impose_high_af_cutoff_here=not impose_high_af_cutoff_upfront)

    # Merge observed variants with expected variants - why outer join?
    observed_variants_ht = observed_variants_ht.join(expected_variants_ht, 'outer')

    # Now aggregate again - Now just on groupings and not on context, ref, alt
    agg_expr = {
        'observed_variants': hl.agg.sum(observed_variants_ht.variant_count),
        'expected_variants': hl.agg.sum(observed_variants_ht.expected_variants),
        'possible_variants': hl.agg.sum(observed_variants_ht.possible_variants),
        'adjusted_mutation_rate': hl.agg.sum(observed_variants_ht.adjusted_mutation_rate),      
        'raw_mutation_rate': hl.agg.sum(observed_variants_ht.mu)
    }
    if pops:
        for pop in POPS:
            agg_expr[f'adjusted_mutation_rate_{pop}'] = hl.agg.array_sum(ht[f'adjusted_mutation_rate_{pop}'])
            agg_expr[f'expected_variants_{pop}'] = hl.agg.array_sum(ht[f'expected_variants_{pop}'])
            agg_expr[f'downsampling_counts_{pop}'] = hl.agg.array_sum(ht[f'downsampling_counts_{pop}'])
        
    observed_variants_ht = observed_variants_ht.group_by(*grouping).partition_hint(1000).aggregate(**agg_expr)

    # Take ratio of observed/expected mutations
    observed_variants_ht = observed_variants_ht.annotate(obs_exp=observed_variants_ht.variant_count / observed_variants_ht.expected_variants)
    observed_variants_ht.write(proportion_variants_observed_ht_path,overwrite)
    return observed_variants_ht


def summarise_prop_observed(po_ht):
    keys = ('gene', 'transcript', 'canonical')
    
    classic_lof_annotations = hl.literal({'stop_gained', 'splice_donor_variant', 'splice_acceptor_variant'})
    filtered_tables = dict(
        lof_ht_classic = pos.ht.filter(classic_lof_annotations.contains(po_ht.annotation) & ((po_ht.modifier == 'HC') | (po_ht.modifier == 'LC'))),
        lof_ht_hc = po_ht.filter(po_ht.modifier == 'HC'),
        lof_ht_classic_hc = po_ht.filter((po_ht.modifier == 'HC') | (po_ht.modifier == 'OS')),
        mis_ht = po_ht.filter(po_ht.annotation == 'missense_variant'),
        mis_pphen_ht = po_ht.filter(po_ht.modifier == 'probably_damaging'),
        mis_non_pphen_ht = po_ht.filter(po_ht.modifier != 'probably_damaging'),
        syn_ht = po_ht.filter(po_ht.annotation == 'synonymous_variant').key_by(*keys)
    )

    for i , (table_name, table) in enumerate(filtered_tables):
        if table_name.startswith('lof'):
            table_agg = utils.collapse_lof_ht(table)
        else:
            agg_expr = {
                'obs': hl.agg.sum(table.variant_count),
                'exp': hl.agg.sum(table.expected_variants),
                'oe': hl.agg.sum(table.variant_count) / hl.agg.sum(table.expected_variants),
                'mutation_rate': hl.agg.sum(table.mu),
                'possible': hl.agg.sum(table.possible_variants)
            }
            table_agg = table.groupby(*keys).agg(**agg_expr)
        # add confidence intervals
        table_agg_cis = utils.oe_confidence_interval(table_agg)
        table_agg.annotate(**table_agg_cis[table_agg.key])
        # label names
        table_agg = utils.label_metrics(table_agg, table_name)
        # join tables
        if i == 0:
            finalised_ht = table_agg
        else:
            finalised_ht = finalised_ht.annotate(**table_agg[finalised_ht.key])

    mut_types = (
        'lof', 'mis', 'syn','mis_pphen','mis_non_pphen'
        )
    output_var_types = zip(('obs', 'exp', 'oe', 'oe', 'oe'),
                            ('', '', '', '_lower', '_upper'))
    output_vars = product(mut_types,output_var_types)
    summary = (finalised_ht
        .select(
            *[f'{t}_{m}{ci}' for m, (t, ci) in output_vars],
            #gene_issues=data['finalised_ht'].constraint_flag 
            # This is not currently included but needs to be
        )
        .select_globals()
    )
    return summary


def aggregate(paths, data):
    '''
    This is the new master function for performing constraint analysis
    Possible variants for populations currently switched off
    '''
    # Loop over autosomes, x y: aggregate by chosen groupings & get proportion observed
    tables = ('auto','x','y')

    for table in tables:
        expected_variants_ht = get_expected_variants(
            data[f'context_{table}_ht'], 
            data['mutation_rate_ht'],
            data['coverage_model'], 
            data['plateau_models'],
            data['grouping'],
            paths['possible_variants_ht_path'].replace('.ht',f'_{table}.ht'),
            pops=False
        )

        prop_observed_ht = get_proportion_observed(
            data[f'exome_{table}_ht'],
            expected_variants_ht,
            data['grouping'],
            paths['po_output_path'].replace('.ht',f'_{table}.ht'))
        data[f'prop_observed_{table}'] = prop_observed_ht

    # Take union of answers and write to file
    data['prop_observed_ht'] = (
                data['prop_observed_auto_ht']
                .union(data['prop_observed_x_ht'])
                .union(data['prop_observed_y_ht'])
                )
    data['prop_observed_ht'].write(paths['po_output_path'], overwrite=True)

    data['summary_ht'] = summarise_prop_observed(data['prop_observed_ht'])
    

    return data