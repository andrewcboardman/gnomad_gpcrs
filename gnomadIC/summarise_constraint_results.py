import argparse
import hail as hl
from itertools import product
from typing import Dict, List, Optional, Set, Tuple, Any
import utils.utils as utils


def finalize_dataset(paths, data, pops = False):
    '''aggregate variants to calculate constraint metrics and significance'''
    # Z score calculation not feasible with partial dataset
    # Need to include flagging of issues in constraint calculations
    keys = ('gene', 'transcript', 'canonical')
    n_partitions = 1000


    # Take union of proportion observed tables
    po_ht = data['po_ht'].union(data['po_x_ht']).union(data['po_y_ht'])
    
    # This function aggregates over genes in all cases, as XG spans PAR and non-PAR X
    #po_ht = po_ht.repartition(n_partitions).persist()

    # Getting classic LoF annotations (no LOFTEE)
    classic_lof_annotations = hl.literal({'stop_gained', 'splice_donor_variant', 'splice_acceptor_variant'})
    lof_ht_classic = po_ht.filter(classic_lof_annotations.contains(po_ht.annotation) &
                                  ((po_ht.modifier == 'HC') | (po_ht.modifier == 'LC')))
    lof_ht_classic = utils.collapse_lof_ht(lof_ht_classic, keys, False)
    lof_ht_classic = lof_ht_classic.rename({x: f'{x}_classic' for x in list(lof_ht_classic.row_value)})

    # Getting all LoF annotations (LOFTEE HC + OS)
    lof_ht_classic_hc = po_ht.filter((po_ht.modifier == 'HC') | (po_ht.modifier == 'OS'))
    lof_ht_classic_hc = utils.collapse_lof_ht(lof_ht_classic_hc, keys, False)
    lof_ht_classic_hc = lof_ht_classic_hc.rename({x: f'{x}_with_os' for x in list(lof_ht_classic_hc.row_value)})

    # Getting all LoF annotations (LOFTEE HC)
    lof_ht = po_ht.filter(po_ht.modifier == 'HC')
    lof_ht = utils.collapse_lof_ht(lof_ht, keys, False)

    # Aggregate missense variants
    mis_ht = po_ht.filter(po_ht.annotation == 'missense_variant')
    agg_expr = {
        'obs_mis': hl.agg.sum(mis_ht.variant_count),
        'exp_mis': hl.agg.sum(mis_ht.expected_variants),
        'oe_mis': hl.agg.sum(mis_ht.variant_count) / hl.agg.sum(mis_ht.expected_variants),
        'mu_mis': hl.agg.sum(mis_ht.mu),
        'possible_mis': hl.agg.sum(mis_ht.possible_variants)
    }

    if pops:
        for pop in utils.POPS:
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
    if pops:
        for pop in utils.POPS:
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
    syn_cis = utils.oe_confidence_interval(ht, ht.obs_syn, ht.exp_syn, prefix='oe_syn')
    mis_cis = utils.oe_confidence_interval(ht, ht.obs_mis, ht.exp_mis, prefix='oe_mis')
    lof_cis = utils.oe_confidence_interval(ht, ht.obs_lof, ht.exp_lof, prefix='oe_lof')
    mis_pphen_cis = utils.oe_confidence_interval(ht, ht.obs_mis_pphen, ht.exp_mis_pphen, prefix='oe_mis_pphen')
    mis_non_pphen_cis = utils.oe_confidence_interval(ht, ht.obs_mis_non_pphen, ht.exp_mis_non_pphen, prefix='oe_mis_non_pphen')

    # join confidence intervals
    ht = ht.annotate(
        **syn_cis[ht.key], 
        **mis_cis[ht.key], 
        **lof_cis[ht.key],
        **mis_pphen_cis[ht.key],
        **mis_non_pphen_cis[ht.key]
        )

    data['finalised_ht'] = ht
    ht.write(paths['finalised_output_path'])

    # Calculate significance
    # ht = calculate_all_z_scores(ht)
    return data
    

def summarise(paths, data):
    mut_types = ('lof', 'mis', 'syn','mis_pphen','mis_non_pphen')
    output_var_types = zip(('obs', 'exp', 'oe', 'oe', 'oe'),
                            ('', '', '', '_lower', '_upper'))
    output_vars = product(mut_types,output_var_types)
    data['summary'] = (data['finalised_ht']
        .select(
            'gene','transcript','canonical',
            *[f'{t}_{m}{ci}' for m, (t, ci) in output_vars],
            gene_issues=data['finalised_ht'].constraint_flag
        )
        .select_globals()
    )
    data['summary'].write(
        paths['summary_output_path'], 
        overwrite=args.overwrite
    )
    data['summary'].export(paths['summary_output_path'].replace('.ht', '.txt.bgz'))


def run_output_test(print_summary=True):
    hl.init(log='hail_logs/test_summarise_constraint_results.log')
    paths = utils.setup_paths('test')
    test_data = utils.load_data_to_finalise(paths)
    test_data = finalize_dataset(paths, test_data)
    print('Test completed!')


def main(args):
    if args.test:
        print('Running tests...')
        run_output_test()
    else:
        print('Please run this script from custom_constraint_analysis.py')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--test', help='Run tests',action='store_true')
    args = parser.parse_args()
    main(args)