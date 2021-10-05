import argparse
import hail as hl
from typing import List
from .gnomad_pipeline.utils import *
from .load_data import *

HIGH_COVERAGE_CUTOFF = 40
POPS = ('global', 'afr', 'amr', 'eas', 'nfe', 'sas')


def aggregate_by_groupings(paths, data):
    '''This is the new master function for performing constraint analysis'''
    # aggregate by chosen groupings & get proportion observed; write to file
    name = 'proportion_observed_output'
    possible_variants_path = paths['possible_variants_ht_path']
    prop_variants_observed_path = paths['po_output_path']
    inputs = zip(
        (data['exome_ht'], data['exome_x_ht'], data['exome_y_ht']), 
        (data['context_ht'], data['context_x_ht'], data['context_y_ht']),
        (name, name + '_x', name + '_y'),
        (possible_variants_path, possible_variants_path.replace('.ht','_x.ht'), possible_variants_path.replace('.ht','_y.ht')),
        (prop_variants_observed_path, prop_variants_observed_path.replace('.ht','_x.ht'), prop_variants_observed_path.replace('.ht','_y.ht'))
    )
    for exome_ht_, context_ht_, name_, possible_variants_path_, prop_variants_observed_path_ in inputs:
        po_exome_ht_ = get_proportion_observed(
            exome_ht_, 
            context_ht_, 
            data['mutation_rate_ht'],
            data['coverage_model'], 
            data['plateau_models'],
            possible_variants_path_,
            prop_variants_observed_path_,
            data['grouping'])
        data[name_] = po_exome_ht_
    return data
    

def get_proportion_observed(
        exome_ht: hl.Table, 
        context_ht: hl.Table,
        mutation_ht: hl.Table,
        coverage_model, plateau_models,
        possible_variants_ht_path,
        proportion_variants_observed_ht_path,
        grouping: List,
        impose_high_af_cutoff_upfront: bool = True,
        pops = False, overwrite=True) -> hl.Table:
    '''Aggregate by grouping variables'''

    # Get possible variants for this context table

    possible_variants_ht = get_possible_variants(context_ht,mutation_ht, coverage_model, plateau_models, grouping, possible_variants_ht_path, pops = pops)

    # Save possible variants??

    # Count observed variants by grouping
    ht = count_variants(exome_ht, additional_grouping=grouping, partition_hint=2000, force_grouping=True,
                        count_downsamplings=[], impose_high_af_cutoff_here=not impose_high_af_cutoff_upfront)

    # Match with expected frequency of variants by grouping
    # This should not be renaming fields!!!
    print('Observed variants \n: ', ht.describe())
    print('Expected variants \n: ', possible_variants_ht.describe())

    ht = ht.join(possible_variants_ht, 'outer')

    # Aggregate all groupings 
    agg_expr = {
        'variant_count': hl.agg.sum(ht.variant_count),
        'adjusted_mutation_rate': hl.agg.sum(ht.adjusted_mutation_rate),
        'possible_variants': hl.agg.sum(ht.possible_variants),
        'expected_variants': hl.agg.sum(ht.expected_variants),
        'mu': hl.agg.sum(ht.mu)
    }
    if pops:
        for pop in POPS:
            agg_expr[f'adjusted_mutation_rate_{pop}'] = hl.agg.array_sum(ht[f'adjusted_mutation_rate_{pop}'])
            agg_expr[f'expected_variants_{pop}'] = hl.agg.array_sum(ht[f'expected_variants_{pop}'])
            agg_expr[f'downsampling_counts_{pop}'] = hl.agg.array_sum(ht[f'downsampling_counts_{pop}'])
        
    ht = ht.group_by(*grouping).partition_hint(1000).aggregate(**agg_expr)
    ht = ht.annotate(obs_exp=ht.variant_count / ht.expected_variants)
    ht.write(proportion_variants_observed_ht_path,overwrite)
    return ht


def get_possible_variants(context_ht, mutation_ht, coverage_model, plateau_models, grouping, possible_file,
        half_cutoff = False, pops = False):
    '''Compute table of possible variants with needed properties'''
    ht = count_variants(context_ht, additional_grouping=grouping, partition_hint=2000, force_grouping=True)
    ht = annotate_with_mu(ht, mutation_ht)
    ht = ht.transmute(possible_variants=ht.variant_count)
    ht = annotate_variant_types(ht.annotate(mu_agg=ht.mu_snp * ht.possible_variants))
    model = hl.literal(plateau_models.total)[ht.cpg]
    cov_cutoff = (HIGH_COVERAGE_CUTOFF / half_cutoff) if half_cutoff else HIGH_COVERAGE_CUTOFF
    ann_expr = {
        'adjusted_mutation_rate': ht.mu_agg * model[1] + model[0],
        'coverage_correction': hl.case()
            .when(ht.coverage == 0, 0)
            .when(ht.coverage >= cov_cutoff, 1)
            .default(coverage_model[1] * hl.log10(ht.coverage) + coverage_model[0])
    }
    if pops:
        for pop in POPS:
            pop_model = hl.literal(plateau_models[pop])
            slopes = hl.map(lambda f: f[ht.cpg][1], pop_model)
            intercepts = hl.map(lambda f: f[ht.cpg][0], pop_model)
            ann_expr[f'adjusted_mutation_rate_{pop}'] = ht.mu_agg * slopes + intercepts
    ht = ht.annotate(**ann_expr)
    ann_expr = {
        'expected_variants': ht.adjusted_mutation_rate * ht.coverage_correction,
        'mu': ht.mu_agg * ht.coverage_correction
    }
    if pops:
        for pop in POPS:
            ann_expr[f'expected_variants_{pop}'] = ht[f'adjusted_mutation_rate_{pop}'] * ht.coverage_correction
    ht = ht.annotate(**ann_expr)
    ht.write(possible_file, True)
    return ht
    



def run_analysis_test(print_summary=True):
    hl.init(log='hail_logs/test_model_expected_variants.log')
    paths = setup_paths('test')
    test_data = load_data_to_aggregate(paths)
    test_data = aggregate_by_groupings(paths, test_data)
    #return paths, test_data


def main(args):
    if args.test:
        print('Running tests...')
        run_analysis_test(print_summary=True)
    else:
        print('Please run this script from custom_constraint_analysis.py')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--test', help='Run tests',action='store_true')
    args = parser.parse_args()
    main(args)