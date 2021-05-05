#!/usr/bin/env python3


#from constraint_utils import generic
import pickle
import copy
import uuid
import hail as hl


def get_all_pop_lengths(ht, prefix: str = 'observed_', pops: List[str] = POPS, skip_assertion: bool = False):
    ds_lengths = ht.aggregate([hl.agg.min(hl.len(ht[f'{prefix}{pop}'])) for pop in pops])
    # temp_ht = ht.take(1)[0]
    # ds_lengths = [len(temp_ht[f'{prefix}{pop}']) for pop in pops]
    pop_lengths = list(zip(ds_lengths, pops))
    print('Found: ', pop_lengths)
    if not skip_assertion:
        assert ht.all(hl.all(lambda f: f, [hl.len(ht[f'{prefix}{pop}']) == length for length, pop in pop_lengths]))
    return pop_lengths


def collapse_lof_ht(lof_ht: hl.Table, keys: Tuple[str], calculate_pop_pLI: bool = False) -> hl.Table:
    agg_expr = {
        'obs_lof': hl.agg.sum(lof_ht.variant_count),
        'mu_lof': hl.agg.sum(lof_ht.mu),
        'possible_lof': hl.agg.sum(lof_ht.possible_variants),
        'exp_lof': hl.agg.sum(lof_ht.expected_variants)
    }
    for pop in POPS:
        agg_expr[f'exp_lof_{pop}'] = hl.agg.array_sum(lof_ht[f'expected_variants_{pop}'])
        agg_expr[f'obs_lof_{pop}'] = hl.agg.array_sum(lof_ht[f'downsampling_counts_{pop}'])
    lof_ht = lof_ht.group_by(*keys).aggregate(**agg_expr).persist()
    lof_ht = lof_ht.filter(lof_ht.exp_lof > 0)
    if calculate_pop_pLI:
        pop_lengths = get_all_pop_lengths(lof_ht, 'obs_lof_')
        print(pop_lengths)
        for pop_length, pop in pop_lengths:
            print(f'Calculating pLI for {pop}...')
            plis = []
            for i in range(8, pop_length):
                print(i)
                ht = lof_ht.filter(lof_ht[f'exp_lof_{pop}'][i] > 0)
                pli_ht = pLI(ht, ht[f'obs_lof_{pop}'][i], ht[f'exp_lof_{pop}'][i])
                plis.append(pli_ht[lof_ht.key])
            lof_ht = lof_ht.annotate(**{
                f'pLI_{pop}': [pli.pLI for pli in plis],
                f'pRec_{pop}': [pli.pRec for pli in plis],
                f'pNull_{pop}': [pli.pNull for pli in plis],
            })
    return lof_ht.annotate(
        **pLI(lof_ht, lof_ht.obs_lof, lof_ht.exp_lof)[lof_ht.key],
        oe_lof=lof_ht.obs_lof / lof_ht.exp_lof).key_by(*keys)

def oe_confidence_interval(ht: hl.Table, obs: hl.expr.Int32Expression, exp: hl.expr.Float32Expression,
                           prefix: str = 'oe', alpha: float = 0.05, select_only_ci_metrics: bool = True) -> hl.Table:
    ht = ht.annotate(_obs=obs, _exp=exp)
    oe_ht = ht.annotate(_range=hl.range(0, 2000).map(lambda x: hl.float64(x) / 1000))
    oe_ht = oe_ht.annotate(_range_dpois=oe_ht._range.map(lambda x: hl.dpois(oe_ht._obs, oe_ht._exp * x)))

    oe_ht = oe_ht.transmute(_cumulative_dpois=hl.cumulative_sum(oe_ht._range_dpois))
    max_cumulative_dpois = oe_ht._cumulative_dpois[-1]
    oe_ht = oe_ht.transmute(_norm_dpois=oe_ht._cumulative_dpois.map(lambda x: x / max_cumulative_dpois))
    oe_ht = oe_ht.transmute(
        _lower_idx=hl.argmax(oe_ht._norm_dpois.map(lambda x: hl.or_missing(x < alpha, x))),
        _upper_idx=hl.argmin(oe_ht._norm_dpois.map(lambda x: hl.or_missing(x > 1 - alpha, x)))
    )
    oe_ht = oe_ht.transmute(**{
        f'{prefix}_lower': hl.cond(oe_ht._obs > 0, oe_ht._range[oe_ht._lower_idx], 0),
        f'{prefix}_upper': oe_ht._range[oe_ht._upper_idx]
    })
    if select_only_ci_metrics:
        return oe_ht.select(f'{prefix}_lower', f'{prefix}_upper')
    else:
        return oe_ht.drop('_exp')


def pLI(ht: hl.Table, obs: hl.expr.Int32Expression, exp: hl.expr.Float32Expression) -> hl.Table:
    last_pi = {'Null': 0, 'Rec': 0, 'LI': 0}
    pi = {'Null': 1 / 3, 'Rec': 1 / 3, 'LI': 1 / 3}
    expected_values = {'Null': 1, 'Rec': 0.463, 'LI': 0.089}
    ht = ht.annotate(_obs=obs, _exp=exp)

    while abs(pi['LI'] - last_pi['LI']) > 0.001:
        last_pi = copy.deepcopy(pi)
        ht = ht.annotate(
            **{k: v * hl.dpois(ht._obs, ht._exp * expected_values[k]) for k, v in pi.items()})
        ht = ht.annotate(row_sum=hl.sum([ht[k] for k in pi]))
        ht = ht.annotate(**{k: ht[k] / ht.row_sum for k, v in pi.items()})
        pi = ht.aggregate({k: hl.agg.mean(ht[k]) for k in pi.keys()})

    ht = ht.annotate(
        **{k: v * hl.dpois(ht._obs, ht._exp * expected_values[k]) for k, v in pi.items()})
    ht = ht.annotate(row_sum=hl.sum([ht[k] for k in pi]))
    return ht.select(**{f'p{k}': ht[k] / ht.row_sum for k, v in pi.items()})


def calculate_z(input_ht: hl.Table, obs: hl.expr.NumericExpression, exp: hl.expr.NumericExpression, output: str = 'z_raw') -> hl.Table:
    ht = input_ht.select(_obs=obs, _exp=exp)
    ht = ht.annotate(_chisq=(ht._obs - ht._exp) ** 2 / ht._exp)
    return ht.select(**{output: hl.sqrt(ht._chisq) * hl.cond(ht._obs > ht._exp, -1, 1)})


def calculate_all_z_scores(ht: hl.Table) -> hl.Table:
    ht = ht.annotate(**calculate_z(ht, ht.obs_syn, ht.exp_syn, 'syn_z_raw')[ht.key])
    ht = ht.annotate(**calculate_z(ht, ht.obs_mis, ht.exp_mis, 'mis_z_raw')[ht.key])
    ht = ht.annotate(**calculate_z(ht, ht.obs_lof, ht.exp_lof, 'lof_z_raw')[ht.key])
    reasons = hl.empty_set(hl.tstr)
    reasons = hl.cond(hl.or_else(ht.obs_syn, 0) + hl.or_else(ht.obs_mis, 0) + hl.or_else(ht.obs_lof, 0) == 0, reasons.add('no_variants'), reasons)
    reasons = hl.cond(ht.exp_syn > 0, reasons, reasons.add('no_exp_syn'), missing_false=True)
    reasons = hl.cond(ht.exp_mis > 0, reasons, reasons.add('no_exp_mis'), missing_false=True)
    reasons = hl.cond(ht.exp_lof > 0, reasons, reasons.add('no_exp_lof'), missing_false=True)
    reasons = hl.cond(hl.abs(ht.syn_z_raw) > 5, reasons.add('syn_outlier'), reasons, missing_false=True)
    reasons = hl.cond(ht.mis_z_raw < -5, reasons.add('mis_too_many'), reasons, missing_false=True)
    reasons = hl.cond(ht.lof_z_raw < -5, reasons.add('lof_too_many'), reasons, missing_false=True)
    ht = ht.annotate(constraint_flag=reasons)
    sds = ht.aggregate(hl.struct(
        syn_sd=hl.agg.filter(
            ~ht.constraint_flag.contains('no_variants') &
            ~ht.constraint_flag.contains('syn_outlier') &
            ~ht.constraint_flag.contains('no_exp_syn') &
            hl.is_defined(ht.syn_z_raw),
            hl.agg.stats(ht.syn_z_raw)).stdev,
        mis_sd=hl.agg.filter(
            ~ht.constraint_flag.contains('no_variants') &
            ~ht.constraint_flag.contains('mis_outlier') &
            ~ht.constraint_flag.contains('no_exp_mis') &
            hl.is_defined(ht.mis_z_raw) & (ht.mis_z_raw < 0),
            hl.agg.explode(lambda x: hl.agg.stats(x), [ht.mis_z_raw, -ht.mis_z_raw])
        ).stdev,
        lof_sd=hl.agg.filter(
            ~ht.constraint_flag.contains('no_variants') &
            ~ht.constraint_flag.contains('lof_outlier') &
            ~ht.constraint_flag.contains('no_exp_lof') &
            hl.is_defined(ht.lof_z_raw) & (ht.lof_z_raw < 0),
            hl.agg.explode(lambda x: hl.agg.stats(x), [ht.lof_z_raw, -ht.lof_z_raw])
        ).stdev
    ))
    print(sds)
    ht = ht.annotate_globals(**sds)
    return ht.transmute(syn_z=ht.syn_z_raw / sds.syn_sd,
                        mis_z=ht.mis_z_raw / sds.mis_sd,
                        lof_z=ht.lof_z_raw / sds.lof_sd)

    
def reannotate_gene_regions(ht, annotation_flags):
    '''Function to add custom gene annotations'''
    return ht

def finalize_dataset(po_ht: hl.Table, keys: Tuple[str] = ('gene', 'transcript', 'canonical'),
                     n_partitions: int = 1000) -> hl.Table:
    # This function aggregates over genes in all cases, as XG spans PAR and non-PAR X
    po_ht = po_ht.repartition(n_partitions).persist()

    # Getting classic LoF annotations (no LOFTEE)
    classic_lof_annotations = hl.literal({'stop_gained', 'splice_donor_variant', 'splice_acceptor_variant'})
    lof_ht_classic = po_ht.filter(classic_lof_annotations.contains(po_ht.annotation) &
                                  ((po_ht.modifier == 'HC') | (po_ht.modifier == 'LC')))
    lof_ht_classic = collapse_lof_ht(lof_ht_classic, keys, False)
    lof_ht_classic = lof_ht_classic.rename({x: f'{x}_classic' for x in list(lof_ht_classic.row_value)})

    # Getting all LoF annotations (LOFTEE HC + OS)
    lof_ht_classic_hc = po_ht.filter((po_ht.modifier == 'HC') | (po_ht.modifier == 'OS'))
    lof_ht_classic_hc = collapse_lof_ht(lof_ht_classic_hc, keys, False)
    lof_ht_classic_hc = lof_ht_classic_hc.rename({x: f'{x}_with_os' for x in list(lof_ht_classic_hc.row_value)})

    # Getting all LoF annotations (LOFTEE HC)
    lof_ht = po_ht.filter(po_ht.modifier == 'HC')
    lof_ht = collapse_lof_ht(lof_ht, keys, False)

    mis_ht = po_ht.filter(po_ht.annotation == 'missense_variant')
    agg_expr = {
        'obs_mis': hl.agg.sum(mis_ht.variant_count),
        'exp_mis': hl.agg.sum(mis_ht.expected_variants),
        'oe_mis': hl.agg.sum(mis_ht.variant_count) / hl.agg.sum(mis_ht.expected_variants),
        'mu_mis': hl.agg.sum(mis_ht.mu),
        'possible_mis': hl.agg.sum(mis_ht.possible_variants)
    }
    for pop in POPS:
        agg_expr[f'exp_mis_{pop}'] = hl.agg.array_sum(mis_ht[f'expected_variants_{pop}'])
        agg_expr[f'obs_mis_{pop}'] = hl.agg.array_sum(mis_ht[f'downsampling_counts_{pop}'])
    mis_ht = mis_ht.group_by(*keys).aggregate(**agg_expr)

    pphen_mis_ht = po_ht.filter(po_ht.modifier == 'probably_damaging')
    pphen_mis_ht = pphen_mis_ht.group_by(*keys).aggregate(obs_mis_pphen=hl.agg.sum(pphen_mis_ht.variant_count),
                                                          exp_mis_pphen=hl.agg.sum(pphen_mis_ht.expected_variants),
                                                          oe_mis_pphen=hl.agg.sum(pphen_mis_ht.variant_count) / hl.agg.sum(pphen_mis_ht.expected_variants),
                                                          possible_mis_pphen=hl.agg.sum(pphen_mis_ht.possible_variants))
    syn_ht = po_ht.filter(po_ht.annotation == 'synonymous_variant').key_by(*keys)
    agg_expr = {
        'obs_syn': hl.agg.sum(syn_ht.variant_count),
        'exp_syn': hl.agg.sum(syn_ht.expected_variants),
        'oe_syn': hl.agg.sum(syn_ht.variant_count) / hl.agg.sum(syn_ht.expected_variants),
        'mu_syn': hl.agg.sum(syn_ht.mu),
        'possible_syn': hl.agg.sum(syn_ht.possible_variants)
    }
    for pop in POPS:
        agg_expr[f'exp_syn_{pop}'] = hl.agg.array_sum(syn_ht[f'expected_variants_{pop}'])
        agg_expr[f'obs_syn_{pop}'] = hl.agg.array_sum(syn_ht[f'downsampling_counts_{pop}'])
    syn_ht = syn_ht.group_by(*keys).aggregate(**agg_expr)

    ht = lof_ht_classic.annotate(**mis_ht[lof_ht_classic.key], **pphen_mis_ht[lof_ht_classic.key],
                                 **syn_ht[lof_ht_classic.key], **lof_ht[lof_ht_classic.key],
                                 **lof_ht_classic_hc[lof_ht_classic.key])
    syn_cis = oe_confidence_interval(ht, ht.obs_syn, ht.exp_syn, prefix='oe_syn')
    mis_cis = oe_confidence_interval(ht, ht.obs_mis, ht.exp_mis, prefix='oe_mis')
    lof_cis = oe_confidence_interval(ht, ht.obs_lof, ht.exp_lof, prefix='oe_lof')
    ht = ht.annotate(**syn_cis[ht.key], **mis_cis[ht.key], **lof_cis[ht.key])
    return calculate_all_z_scores(ht)

def run_tests():
    """Tests paths given have correct schema and functions work correctly"""
    print('Tests need to be implemented')

def main():

    root = 'gs://gnomad-public/papers/2019-flagship-lof/v1.0'
    po_ht_path = f'{root}/{{subdir}}/prop_observed_{{subdir}}.ht'
    raw_constraint_ht_path = f'{root}/{{subdir}}/constraint_{{subdir}}.ht'
    final_constraint_ht_path = f'{root}/{{subdir}}/constraint_final_{{subdir}}.ht'
    if args.dataset != 'gnomad':
        po_coverage_ht_path = po_coverage_ht_path.replace(root, root + f'/{args.dataset}')
    POPS = ('global', 'afr', 'amr', 'eas', 'nfe', 'sas')
    MODEL_KEYS = {
        'worst_csq': ['gene'],
        'tx_annotation': ['gene', 'expressed'],
        'standard': ['gene', 'transcript', 'canonical']
    }
    # Set paths based on input dataset and model
    po_output_path = po_ht_path.format(subdir=args.model)
    output_path = raw_constraint_ht_path.format(subdir=args.model)
    final_path = final_constraint_ht_path.format(subdir=args.model)
    if args.test:
        run_tests()
    
    else:
        if args.aggregate:
            print('Running aggregation')
            # read PO hail tables for autosomes, X and Y chromosomes and join them
            print(f'Reading hail table from {po_output_path}')
            ht = hl.read_table(po_output_path).union(
                hl.read_table(po_output_path.replace('.ht', '_x.ht'))
            ).union(
                hl.read_table(po_output_path.replace('.ht', '_y.ht'))
            )
            # Reannotate gene regions
            ht = reannotate_gene_regions(ht, {})
            # group by gene/transcript and calculate summary stats
            if args.model != 'syn_canonical':
                ht = finalize_dataset(ht, keys=MODEL_KEYS[args.model])
            # write hail table to output path
            ht.write(output_path, args.overwrite)
            hl.read_table(output_path).export(output_path.replace('.ht', '.txt.bgz'))
        if args.summarise:
            print('Finalising summary stats')
            # write summary stats to output path
            ht = hl.read_table(output_path)
            var_types = ('lof', 'mis', 'syn')
            ht.select(
                *[f'{t}_{v}{ci}' for v in var_types
                    for t, ci in zip(('obs', 'exp', 'oe', 'mu', 'oe', 'oe'),
                                    ('', '', '', '', '_lower', '_upper'))],
                *[f'{v}_z' for v in var_types], 'pLI', 'pRec', 'pNull', gene_issues=ht.constraint_flag
            ).select_globals().write(final_path, overwrite=args.overwrite)
            hl.read_table(final_path).export(final_path.replace('.ht', '.txt.bgz'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--test', help='Run tests without actually requesting data',action='store_true')
    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--dataset', help='Which dataset to use (one of gnomad, non_neuro, non_cancer, controls)', default='gnomad')
    parser.add_argument('--model', help='Which model to apply (one of "standard", "syn_canonical", or "worst_csq" for now)', default='standard')
    parser.add_argument('--aggregate', help='Get p_obs table and aggregate', action='store_true')
    parser.add_argument('--summarise', help='Report summary stats', action='store_true')
    args = parser.parse_args()
    main(args)