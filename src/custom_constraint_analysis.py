#!/usr/bin/env python3
import argparse
import os
import pickle
import copy
import uuid
from itertools import product
import hail as hl
from typing import Dict, List, Optional, Set, Tuple

# set pops
POPS = ('global', 'afr', 'amr', 'eas', 'nfe', 'sas')

def get_all_pop_lengths(ht, prefix: str = 'observed_', pops: List[str] = POPS, skip_assertion: bool = False):
    '''Get lengths of arrays for each pop'''
    ds_lengths = ht.aggregate([hl.agg.min(hl.len(ht[f'{prefix}{pop}'])) for pop in pops])
    temp_ht = ht.take(1)[0]
    ds_lengths = [len(temp_ht[f'{prefix}{pop}']) for pop in pops]
    pop_lengths = list(zip(ds_lengths, pops))
    print('Found: ', pop_lengths)
    if not skip_assertion:
        assert ht.all(hl.all(lambda f: f, [hl.len(ht[f'{prefix}{pop}']) == length for length, pop in pop_lengths]))
    return pop_lengths


def collapse_lof_ht(lof_ht: hl.Table, keys: Tuple[str], calculate_pop_pLI: bool = False) -> hl.Table:
    '''Aggregate lof variants in genes for each population'''
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
    '''Calculate CI for observed/expected ratio'''
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
    '''Calculate p(lof intolerant) - metric for constraint'''
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
    '''get significance of a constraint finding'''
    ht = input_ht.select(_obs=obs, _exp=exp)
    ht = ht.annotate(_chisq=(ht._obs - ht._exp) ** 2 / ht._exp)
    return ht.select(**{output: hl.sqrt(ht._chisq) * hl.cond(ht._obs > ht._exp, -1, 1)})


def calculate_all_z_scores(ht: hl.Table) -> hl.Table:
    '''calculate significance for all constraint findings in table'''

    # Raw z scores from chi squared distribution
    ht = ht.annotate(**calculate_z(ht, ht.obs_syn, ht.exp_syn, 'syn_z_raw')[ht.key])
    ht = ht.annotate(**calculate_z(ht, ht.obs_mis, ht.exp_mis, 'mis_z_raw')[ht.key])
    ht = ht.annotate(**calculate_z(ht, ht.obs_lof, ht.exp_lof, 'lof_z_raw')[ht.key])
    ht = ht.annotate(**calculate_z(ht, ht.obs_mis_pphen, ht.exp_mis_pphen, 'mis_pphen_z_raw')[ht.key])
    ht = ht.annotate(**calculate_z(ht, ht.obs_mis_non_pphen, ht.exp_mis_non_pphen, 'mis_non_pphen_z_raw')[ht.key])
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
        ).stdev,
        # mis_pphen_sd=hl.agg.filter(
        #     ~ht.constraint_flag.contains('no_variants') &
        #     ~ht.constraint_flag.contains('mis_outlier') &
        #     ~ht.constraint_flag.contains('no_exp_mis') &
        #     hl.is_defined(ht.mis_pphen_z_raw) & (ht.mis_pphen_z_raw < 0),
        #     hl.agg.explode(lambda x: hl.agg.stats(x), [ht.mis_pphen_z_raw, -ht.mis_pphen_z_raw])
        # ).stdev,
        # mis_non_pphen_sd=hl.agg.filter(
        #     ~ht.constraint_flag.contains('no_variants') &
        #     ~ht.constraint_flag.contains('mis_outlier') &
        #     ~ht.constraint_flag.contains('no_exp_mis') &
        #     hl.is_defined(ht.mis_non_pphen_z_raw) & (ht.mis_non_pphen_z_raw < 0),
        #     hl.agg.explode(lambda x: hl.agg.stats(x), [ht.mis_non_pphen_z_raw, -ht.mis_non_pphen_z_raw])
        # ).stdev
    ))
    print(sds)
    ht = ht.annotate_globals(**sds)
    # ht_z = ht.transmute(
    #     syn_z=ht.syn_z_raw / sds.syn_sd,
    #     mis_z=ht.mis_z_raw / sds.mis_sd,
    #     lof_z=ht.lof_z_raw / sds.lof_sd,
    #     #mis_pphen_z=ht.mis_pphen_z_raw / sds.mis_pphen_sd,
    #     #mis_non_pphen_z=ht.mis_non_pphen_z_raw / sds.mis_non_pphen_sd
    # )
    return ht_z


def finalize_dataset(po_ht: hl.Table, keys: Tuple[str] = ('gene', 'transcript', 'canonical'),
                     n_partitions: int = 1000) -> hl.Table:
    '''aggregate variants to calculate constraint metrics and significance'''
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

    # Aggregate missense variants
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
    for pop in POPS:
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
    syn_cis = oe_confidence_interval(ht, ht.obs_syn, ht.exp_syn, prefix='oe_syn')
    mis_cis = oe_confidence_interval(ht, ht.obs_mis, ht.exp_mis, prefix='oe_mis')
    lof_cis = oe_confidence_interval(ht, ht.obs_lof, ht.exp_lof, prefix='oe_lof')
    mis_pphen_cis = oe_confidence_interval(ht, ht.obs_mis_pphen, ht.exp_mis_pphen, prefix='oe_mis_pphen')
    mis_non_pphen_cis = oe_confidence_interval(ht, ht.obs_mis_non_pphen, ht.exp_mis_non_pphen, prefix='oe_mis_non_pphen')

    # join confidence intervals
    ht = ht.annotate(
        **syn_cis[ht.key], 
        **mis_cis[ht.key], 
        **lof_cis[ht.key],
        **mis_pphen_cis[ht.key],
        **mis_non_pphen_cis[ht.key]
        )
    # Calculate significance
    #ht = calculate_all_z_scores(ht)
    return ht

def run_tests(ht):
    """Tests loading of autosome po table"""
    # Incorporate more tests to check that aggregation by new variant annotations is legit
    ht.show()

def load_or_import(path, overwrite):
    # Specify input format to avoid coercion to string
    types = {
        'adjusted_mutation_rate_global': hl.expr.types.tarray(hl.tfloat64),
        'expected_variants_global': hl.expr.types.tarray(hl.tfloat64),
        'downsampling_counts_global': hl.expr.types.tarray(hl.tint32),
        'adjusted_mutation_rate_afr': hl.expr.types.tarray(hl.tfloat64),
        'expected_variants_afr': hl.expr.types.tarray(hl.tfloat64),
        'downsampling_counts_afr': hl.expr.types.tarray(hl.tint32),
        'adjusted_mutation_rate_amr': hl.expr.types.tarray(hl.tfloat64),
        'expected_variants_amr': hl.expr.types.tarray(hl.tfloat64),
        'downsampling_counts_amr': hl.expr.types.tarray(hl.tint32),
        'adjusted_mutation_rate_eas': hl.expr.types.tarray(hl.tfloat64),
        'expected_variants_eas': hl.expr.types.tarray(hl.tfloat64),
        'downsampling_counts_eas': hl.expr.types.tarray(hl.tint32),
        'adjusted_mutation_rate_nfe': hl.expr.types.tarray(hl.tfloat64),
        'expected_variants_nfe': hl.expr.types.tarray(hl.tfloat64),
        'downsampling_counts_nfe': hl.expr.types.tarray(hl.tint32),
        'adjusted_mutation_rate_sas': hl.expr.types.tarray(hl.tfloat64),
        'expected_variants_sas': hl.expr.types.tarray(hl.tfloat64),
        'downsampling_counts_sas': hl.expr.types.tarray(hl.tint32)
        }

    if os.path.isdir(path) and not overwrite:
            ht = hl.read_table(path)
    else:
        ht = hl.import_table(path.replace('.ht','.txt.bgz'),impute=True,types=types)
        ht.write(path,overwrite)
    return ht

def main(args):
    # Set paths for data access based on command line parameters
    root = './data'
    
    context_ht_path = f'{root}/context/Homo_sapiens_assembly19.fasta.snps_only.vep_20181129.ht'
    processed_genomes_ht_path = f'{root}/model/genomes_processed.ht'
    processed_exomes_ht_path = f'{root}/model/exomes_processed.ht'
    mutation_rate_ht_path = f'{root}/model/mutation_rate_methylation_bins.ht'
    po_coverage_ht_path = f'{root}/model/prop_observed_by_coverage_no_common_pass_filtered_bins.ht'
    po_ht_path = f'{root}/{{subdir}}/prop_observed_{{subdir}}.ht'
    raw_constraint_ht_path = f'{root}/{{subdir}}/constraint_{{subdir}}.ht'
    final_constraint_ht_path = f'{root}/{{subdir}}/constraint_final_{{subdir}}.ht'

    po_output_path = po_ht_path.format(subdir=args.model)
    output_path = raw_constraint_ht_path.format(subdir=args.model)
    final_path = final_constraint_ht_path.format(subdir=args.model)

    # Sets method for aggregation, will need to be changed for custom analysis
    MODEL_KEYS = {
        'worst_csq': ['gene'],
        'tx_annotation': ['gene', 'expressed'],
        'standard': ['gene', 'transcript', 'canonical']
    }
    
    if args.test:
        ht = load_or_import(po_output_path, args.overwrite)
        run_tests(ht)

    if args.get_proportion_observed:

        full_context_ht = prepare_ht(hl.read_table(context_ht_path), args.trimers)
       #full_genome_ht = prepare_ht(hl.read_table(processed_genomes_ht_path), args.trimers)
        full_exome_ht = prepare_ht(hl.read_table(processed_exomes_ht_path), args.trimers)

        context_ht = full_context_ht.filter(full_context_ht.locus.in_autosome_or_par())
        #genome_ht = full_genome_ht.filter(full_genome_ht.locus.in_autosome_or_par())
        exome_ht = full_exome_ht.filter(full_exome_ht.locus.in_autosome_or_par())

        context_x_ht = hl.filter_intervals(full_context_ht, [hl.parse_locus_interval('X')])
        context_x_ht = context_x_ht.filter(context_x_ht.locus.in_x_nonpar())
        context_y_ht = hl.filter_intervals(full_context_ht, [hl.parse_locus_interval('Y')])
        context_y_ht = context_y_ht.filter(context_y_ht.locus.in_y_nonpar())

        exome_x_ht = hl.filter_intervals(full_exome_ht, [hl.parse_locus_interval('X')])
        exome_x_ht = exome_x_ht.filter(exome_x_ht.locus.in_x_nonpar())
        exome_y_ht = hl.filter_intervals(full_exome_ht, [hl.parse_locus_interval('Y')])
        exome_y_ht = exome_y_ht.filter(exome_y_ht.locus.in_y_nonpar())

        # mutation rate
        mutation_ht = hl.read_table(mutation_rate_ht_path).select('mu_snp')

        # Get proportion observed tables
        coverage_ht = hl.read_table(po_coverage_ht_path)
        coverage_x_ht = hl.read_table(po_coverage_ht_path.replace('.ht', '_x.ht'))
        coverage_y_ht = hl.read_table(po_coverage_ht_path.replace('.ht', '_y.ht'))

        # Build models for mutation frequency based on context
        coverage_model, plateau_models = build_models(coverage_ht, args.trimers, True)
        _, plateau_x_models = build_models(coverage_x_ht, args.trimers, True)
        _, plateau_y_models = build_models(coverage_y_ht, args.trimers, True)

        get_proportion_observed(exome_ht, context_ht, mutation_ht, plateau_models,
                                coverage_model, recompute_possible=True,
                                custom_model=args.model, dataset=args.dataset,
                                impose_high_af_cutoff_upfront=not args.skip_af_filter_upfront
                                ).write(po_output_path, overwrite=args.overwrite)
        hl.read_table(po_output_path).export(po_output_path.replace('.ht', '.txt.bgz'))

        get_proportion_observed(exome_x_ht, context_x_ht, mutation_ht, plateau_x_models,
                                coverage_model, recompute_possible=True,
                                custom_model=args.model, dataset=args.dataset,
                                impose_high_af_cutoff_upfront=not args.skip_af_filter_upfront
                                ).write(po_output_path.replace('.ht', '_x.ht'), overwrite=args.overwrite)
        hl.read_table(po_output_path.replace('.ht', '_x.ht')).export(po_output_path.replace('.ht', '_x.txt.bgz'))

        get_proportion_observed(exome_y_ht, context_y_ht, mutation_ht, plateau_y_models,
                                coverage_model, recompute_possible=True,
                                custom_model=args.model, dataset=args.dataset,
                                impose_high_af_cutoff_upfront=not args.skip_af_filter_upfront
                                ).write(po_output_path.replace('.ht', '_y.ht'), overwrite=args.overwrite)
        hl.read_table(po_output_path.replace('.ht', '_y.ht')).export(po_output_path.replace('.ht', '_y.txt.bgz'))

    if args.aggregate:
        print('Running aggregation')
        # read PO hail tables for autosomes, X and Y chromosomes and join them
        print(f'Reading hail table from {po_output_path}')
        # Autosomes
        ht_autosomes = load_or_import(po_output_path, args.overwrite)
        # X chromosome
        ht_x = load_or_import(po_output_path.replace('.ht','_x.ht'), args.overwrite)
        # Y chromosome
        ht_y = load_or_import(po_output_path.replace('.ht','_y.ht'), args.overwrite)
        # Combine into one table
        ht = ht_autosomes.union(ht_x).union(ht_y)
        # group by gene/transcript and calculate summary stats
        if args.model != 'syn_canonical':
            ht = finalize_dataset(ht, keys=MODEL_KEYS[args.model])
        # write hail table to output path
        ht.write(output_path, args.overwrite)
        hl.read_table(output_path).export(output_path.replace('.ht', '.txt.bgz'))

    if args.summarise:
        print('Finalising summary stats')
        # write summary stats to output path
        ht = load_or_import(output_path, args.overwrite)
        mut_types = ('lof', 'mis', 'syn','mis_pphen','mis_non_pphen')
        output_var_types = zip(('obs', 'exp', 'oe', 'oe', 'oe'),
                                ('', '', '', '_lower', '_upper'))
        output_vars = product(mut_types,output_var_types)
        ht.select(
            'gene','transcript','canonical',
            *[f'{t}_{m}{ci}' for m, (t, ci) in output_vars],
            #*[f'{m}_z' for m in mut_types[:3]], 
            'pLI', 'pRec', 'pNull', 
            #gene_issues=ht.constraint_flag
        ).select_globals().write(final_path, overwrite=args.overwrite)
        hl.read_table(final_path).export(final_path.replace('.ht', '.txt.bgz'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--test', help='Run tests without actually requesting data',action='store_true')
    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--dataset', help='Which dataset to use (one of gnomad, non_neuro, non_cancer, controls)', default='gnomad')
    parser.add_argument('--model', help='Which model to apply (one of "standard", "syn_canonical", or "worst_csq" for now)', default='standard')
    parser.add_argument('--skip_af_filter_upfront', help='Skip AF filter up front (to be applied later to ensure that it is not affecting population-specific constraint): not generally recommended', action='store_true')
    parser.add_argument('--apply_model',help='Apply model to calculate proportion observed')
    parser.add_argument('--aggregate', help='Aggregate p_obs table', action='store_true')
    parser.add_argument('--summarise', help='Report summary stats', action='store_true')
    args = parser.parse_args()
    main(args)