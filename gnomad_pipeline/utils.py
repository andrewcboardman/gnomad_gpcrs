from gnomad.utils.vep import *
import hail as hl
import argparse
import os
import pickle
import copy
import uuid
from itertools import product
import hail as hl
from typing import Dict, List, Optional, Set, Tuple, Any

def annotate_with_mu(ht: hl.Table, mutation_ht: hl.Table, output_loc: str = 'mu_snp',
                     keys: Tuple[str] = ('context', 'ref', 'alt', 'methylation_level')) -> hl.Table:
    mu = hl.literal(mutation_ht.aggregate(hl.dict(hl.agg.collect(
        (hl.struct(**{k: mutation_ht[k] for k in keys}), mutation_ht.mu_snp)))))
    mu = mu.get(hl.struct(**{k: ht[k] for k in keys}))
    return ht.annotate(**{output_loc: hl.case().when(hl.is_defined(mu), mu).or_error('Missing mu')})

def build_models(coverage_ht: hl.Table, trimers: bool = False, weighted: bool = False, half_cutoff = False,
                 ) -> Tuple[Tuple[float, float], Dict[str, Tuple[float, float]]]:
    keys = ['context', 'ref', 'alt', 'methylation_level', 'mu_snp']

    cov_cutoff = (HIGH_COVERAGE_CUTOFF / half_cutoff) if half_cutoff else HIGH_COVERAGE_CUTOFF
    all_high_coverage_ht = coverage_ht.filter(coverage_ht.exome_coverage >= cov_cutoff)
    agg_expr = {
        'observed_variants': hl.agg.sum(all_high_coverage_ht.variant_count),
        'possible_variants': hl.agg.sum(all_high_coverage_ht.possible_variants)
    }
    for pop in POPS:
        agg_expr[f'observed_{pop}'] = hl.agg.array_sum(all_high_coverage_ht[f'downsampling_counts_{pop}'])
    high_coverage_ht = all_high_coverage_ht.group_by(*keys).aggregate(**agg_expr)

    high_coverage_ht = annotate_variant_types(high_coverage_ht, not trimers)
    plateau_models = build_plateau_models_pop(high_coverage_ht, weighted=weighted)

    high_coverage_scale_factor = all_high_coverage_ht.aggregate(
        hl.agg.sum(all_high_coverage_ht.variant_count) /
        hl.agg.sum(all_high_coverage_ht.possible_variants * all_high_coverage_ht.mu_snp))

    all_low_coverage_ht = coverage_ht.filter((coverage_ht.exome_coverage < cov_cutoff) &
                                             (coverage_ht.exome_coverage > 0))

    low_coverage_ht = all_low_coverage_ht.group_by(log_coverage=hl.log10(all_low_coverage_ht.exome_coverage)).aggregate(
        low_coverage_obs_exp=hl.agg.sum(all_low_coverage_ht.variant_count) /
                             (high_coverage_scale_factor * hl.agg.sum(all_low_coverage_ht.possible_variants * all_low_coverage_ht.mu_snp)))
    coverage_model = build_coverage_model(low_coverage_ht)
    # TODO: consider weighting here as well

    return coverage_model, plateau_models

    
def split_context_mt(raw_context_ht_path: str, coverage_ht_paths: Dict[str, str], methylation_ht_path: str,
                     context_ht_path: str, overwrite: bool = False) -> None:
    '''Pre-process gnomad release data with context, coverage and methylation data'''
    raw_context_ht = hl.split_multi_hts(hl.read_table(raw_context_ht_path))
    raw_context_ht.write(f'{raw_context_ht_path}.temp_split.ht', overwrite)
    raw_context_ht = hl.read_table(f'{raw_context_ht_path}.temp_split.ht')

    coverage_hts = {loc: hl.read_table(coverage_ht_path) for loc, coverage_ht_path in coverage_ht_paths.items()}
    coverage_hts = {loc: coverage_ht.drop('#chrom', 'pos') if '#chrom' in list(coverage_ht.row) else coverage_ht
                    for loc, coverage_ht in coverage_hts.items()}
    methylation_ht = hl.read_table(methylation_ht_path)
    gerp_ht = hl.read_table(gerp_annotations_path)

    raw_context_ht = raw_context_ht.annotate(
        methylation=methylation_ht[raw_context_ht.locus],
        coverage=hl.struct(**{loc: coverage_ht[raw_context_ht.locus] for loc, coverage_ht in coverage_hts.items()}),
        gerp=gerp_ht[raw_context_ht.locus].S)

    raw_context_ht.write(context_ht_path, overwrite)


def pre_process_data(ht: hl.Table, split_context_ht_path: str,
                     output_ht_path: str, overwrite: bool = False) -> None:
    context_ht = hl.read_table(split_context_ht_path).drop('a_index', 'was_split')
    context_ht = context_ht.annotate(vep=context_ht.vep.drop('colocated_variants'))
    ht.annotate(**context_ht[ht.key], pass_filters=hl.len(ht.filters) == 0).write(output_ht_path, overwrite)

def load_variant_table(file):
    ht = hl.import_table(file,impute=True)
    return ht

# Preparation of exomes table 

def reverse_complement_bases(bases: hl.expr.StringExpression) -> hl.expr.StringExpression:
    return hl.delimit(hl.range(bases.length() - 1, -1, -1).map(lambda i: flip_base(bases[i])), '')
    # return bases[::-1].map(lambda x: flip_base(x))


def flip_base(base: hl.expr.StringExpression) -> hl.expr.StringExpression:
    return (hl.switch(base)
            .when('A', 'T')
            .when('T', 'A')
            .when('G', 'C')
            .when('C', 'G')
            .default(base))

def collapse_strand(ht: Union[hl.Table, hl.MatrixTable]) -> Union[hl.Table, hl.MatrixTable]:
    collapse_expr = {
        'ref': hl.cond(((ht.ref == 'G') | (ht.ref == 'T')),
                       reverse_complement_bases(ht.ref), ht.ref),
        'alt': hl.cond(((ht.ref == 'G') | (ht.ref == 'T')),
                       reverse_complement_bases(ht.alt), ht.alt),
        'context': hl.cond(((ht.ref == 'G') | (ht.ref == 'T')),
                           reverse_complement_bases(ht.context), ht.context),
        'was_flipped': (ht.ref == 'G') | (ht.ref == 'T')
    }
    return ht.annotate(**collapse_expr) if isinstance(ht, hl.Table) else ht.annotate_rows(**collapse_expr)


def prepare_ht(ht, trimer: bool = False, annotate_coverage: bool = True):
    if trimer:
        ht = trimer_from_heptamer(ht)
    str_len = 3 if trimer else 7

    if isinstance(ht, hl.Table): 
        ht = ht.annotate(ref=ht.alleles[0], alt=ht.alleles[1])
        ht = ht.filter((hl.len(ht.ref) == 1) & (hl.len(ht.alt) == 1) & ht.context.matches(f'[ATCG]{{{str_len}}}'))
        ht = annotate_variant_types(collapse_strand(ht), not trimer)
    else: # handle case where ht is a matrix table
        ht = ht.annotate_rows(ref=ht.alleles[0], alt=ht.alleles[1])
        ht = ht.filter_rows((hl.len(ht.ref) == 1) & (hl.len(ht.alt) == 1) & ht.context.matches(f'[ATCG]{{{str_len}}}'))
        ht = annotate_variant_types(collapse_strand(ht), not trimer)
    annotation = {
        'methylation_level': hl.case().when(
            ht.cpg & (ht.methylation.MEAN > 0.6), 2
        ).when(
            ht.cpg & (ht.methylation.MEAN > 0.2), 1
        ).default(0)
    }
    if annotate_coverage:
        annotation['exome_coverage'] = ht.coverage.exomes.median
    return ht.annotate(**annotation) if isinstance(ht, hl.Table) else ht.annotate_rows(**annotation)

def trimer_from_heptamer(t: Union[hl.MatrixTable, hl.Table]) -> Union[hl.MatrixTable, hl.Table]:
    trimer_expr = hl.cond(hl.len(t.context) == 7, t.context[2:5], t.context)
    return t.annotate_rows(context=trimer_expr) if isinstance(t, hl.MatrixTable) else t.annotate(context=trimer_expr)

def annotate_variant_types(t: Union[hl.MatrixTable, hl.Table],
                           heptamers: bool = False) -> Union[hl.MatrixTable, hl.Table]:
    """
    Adds cpg, transition, and variant_type, variant_type_model columns
    """
    mid_index = 3 if heptamers else 1
    transition_expr = (((t.ref == "A") & (t.alt == "G")) | ((t.ref == "G") & (t.alt == "A")) |
                       ((t.ref == "T") & (t.alt == "C")) | ((t.ref == "C") & (t.alt == "T")))
    cpg_expr = (((t.ref == "G") & (t.alt == "A") & (t.context[mid_index - 1:mid_index] == 'C')) |
                ((t.ref == "C") & (t.alt == "T") & (t.context[mid_index + 1:mid_index + 2] == 'G')))
    if isinstance(t, hl.MatrixTable):
        t = t.annotate_rows(transition=transition_expr, cpg=cpg_expr)
    else:
        t = t.annotate(transition=transition_expr, cpg=cpg_expr)
    variant_type_expr = (hl.case()
                         .when(t.cpg, 'CpG')
                         .when(t.transition, 'non-CpG transition')
                         .default('transversion'))
    variant_type_model_expr = hl.cond(t.cpg, t.context, "non-CpG")
    if isinstance(t, hl.MatrixTable):
        return t.annotate_rows(variant_type=variant_type_expr, variant_type_model=variant_type_model_expr)
    else:
        return t.annotate(variant_type=variant_type_expr, variant_type_model=variant_type_model_expr)

# Annotation of exomes table using VEP information

def add_most_severe_csq_to_tc_within_ht(t):
    annotation = t.vep.annotate(transcript_consequences=t.vep.transcript_consequences.map(
        add_most_severe_consequence_to_consequence))
    return t.annotate_rows(vep=annotation) if isinstance(t, hl.MatrixTable) else t.annotate(vep=annotation)

def annotate_constraint_groupings(ht: Union[hl.Table, hl.MatrixTable],
                                  custom_model: str = None) -> Tuple[Union[hl.Table, hl.MatrixTable], List[str]]:
    """
    HT must be exploded against whatever axis

    Need to add `'coverage': ht.exome_coverage` here (which will get corrected out later)
    """
    groupings = {
        'annotation': ht.transcript_consequences.most_severe_consequence,
        'modifier': hl.case()
            .when(hl.is_defined(ht.transcript_consequences.lof),
                    ht.transcript_consequences.lof)
            .when(hl.is_defined(ht.transcript_consequences.polyphen_prediction),
                    ht.transcript_consequences.polyphen_prediction)
            .default('None'),
        'transcript': ht.transcript_consequences.transcript_id,
        'gene': ht.transcript_consequences.gene_symbol,
        'canonical': hl.or_else(ht.transcript_consequences.canonical == 1, False),
        'coverage': ht.exome_coverage
        }
    
    ht = ht.annotate(**groupings) if isinstance(ht, hl.Table) else ht.annotate_rows(**groupings)
    return ht, list(groupings.keys())

# Aggregation of variant counts

def count_variants(ht: hl.Table,
                   count_singletons: bool = False, count_downsamplings: Optional[List[str]] = (),
                   additional_grouping: Optional[List[str]] = (), partition_hint: int = 100,
                   omit_methylation: bool = False, return_type_only: bool = False,
                   force_grouping: bool = False, singleton_expression: hl.expr.BooleanExpression = None,
                   impose_high_af_cutoff_here: bool = False) -> Union[hl.Table, Any]:
    """
    Count variants by context, ref, alt, methylation_level
    """

    grouping = hl.struct(context=ht.context, ref=ht.ref, alt=ht.alt)
    if not omit_methylation:
        grouping = grouping.annotate(methylation_level=ht.methylation_level)
    for group in additional_grouping:
        grouping = grouping.annotate(**{group: ht[group]})

    if count_singletons:
        # singleton = hl.any(lambda f: (f.meta.size() == 1) & (f.meta.get('group') == 'adj') & (f.AC[1] == 1), ht.freq)
        if singleton_expression is None:
            singleton_expression = ht.freq[0].AC == 1

    if count_downsamplings or force_grouping:
        # Slower, but more flexible (allows for downsampling agg's)
        output = {'variant_count': hl.agg.count_where(ht.freq[0].AF <= 0.001) if impose_high_af_cutoff_here else hl.agg.count()}
        for pop in count_downsamplings:
            output[f'downsampling_counts_{pop}'] = downsampling_counts_expr(ht, pop, impose_high_af_cutoff=impose_high_af_cutoff_here)
        if count_singletons:
            output['singleton_count'] = hl.agg.count_where(singleton_expression)
            for pop in count_downsamplings:
                output[f'singleton_downsampling_counts_{pop}'] = downsampling_counts_expr(ht, pop, singleton=True)
        return ht.group_by(**grouping).partition_hint(partition_hint).aggregate(**output)
    else:
        agg = {'variant_count': hl.agg.counter(grouping)}
        if count_singletons:
            agg['singleton_count'] = hl.agg.counter(hl.agg.filter(singleton_expression, grouping))

        if return_type_only:
            return agg['variant_count'].dtype
        else:
            return ht.aggregate(hl.struct(**agg))

def downsampling_counts_expr(ht: Union[hl.Table, hl.MatrixTable], pop: str = 'global', variant_quality: str = 'adj',
                             singleton: bool = False, impose_high_af_cutoff: bool = False) -> hl.expr.ArrayExpression:
    indices = hl.zip_with_index(ht.freq_meta).filter(
        lambda f: (f[1].size() == 3) & (f[1].get('group') == variant_quality) &
                  (f[1].get('pop') == pop) & f[1].contains('downsampling')
    )
    sorted_indices = hl.sorted(indices, key=lambda f: hl.int(f[1]['downsampling'])).map(lambda x: x[0])
    # TODO: this likely needs to be fixed for aggregations that return missing (need to be 0'd out)

    def get_criteria(i):
        if singleton:
            return hl.int(ht.freq[i].AC == 1)
        elif impose_high_af_cutoff:
            return hl.int((ht.freq[i].AC > 0) & (ht.freq[i].AF <= 0.001))
        else:
            return hl.int(ht.freq[i].AC > 0)
    return hl.agg.array_sum(hl.map(get_criteria, sorted_indices))

# Further aggregation by gene level to finalise dataset

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

# Calculation of summary stats

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


def annotate_issues(ht: hl.Table) -> hl.Table:
    '''Annotate issues with constraint calculations'''
    reasons = hl.empty_set(hl.tstr)
    reasons = hl.cond(hl.or_else(ht.obs_syn, 0) + hl.or_else(ht.obs_mis, 0) + hl.or_else(ht.obs_lof, 0) == 0, reasons.add('no_variants'), reasons)
    reasons = hl.cond(ht.exp_syn > 0, reasons, reasons.add('no_exp_syn'), missing_false=True)
    reasons = hl.cond(ht.exp_mis > 0, reasons, reasons.add('no_exp_mis'), missing_false=True)
    reasons = hl.cond(ht.exp_lof > 0, reasons, reasons.add('no_exp_lof'), missing_false=True)
    ht = ht.annotate(constraint_flag=reasons)
    return ht


# functions for loading compressed tables from text form to re-run 

def load_or_import_po(path, overwrite):
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
