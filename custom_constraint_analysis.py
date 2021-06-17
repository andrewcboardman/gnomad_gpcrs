#!/usr/bin/python

import argparse
from itertools import product
import hail as hl
from typing import Dict, List, Optional, Set, Tuple
from gnomad.utils.vep import *
import hail as hl
import argparse
from itertools import product
import hail as hl
from typing import Dict, List, Optional, Set, Tuple, Any
from gnomad_pipeline.utils import *

POPS = ('global', 'afr', 'amr', 'eas', 'nfe', 'sas') # set pops
HIGH_COVERAGE_CUTOFF = 40
exomes_path = 'gs://gcp-public-data--gnomad/release/2.1/ht/exomes/gnomad.exomes.r2.1.sites.ht/'
coverage_path_exomes = 'gs://gcp-public-data--gnomad/release/2.1/coverage/exomes/gnomad.exomes.r2.1.coverage.ht'
context_path_vep = 'gs://gcp-public-data--gnomad/resources/context/grch37_context_vep_annotated.ht/'
methylation_path = 'gs://gnomad-public/resources/grch37/methylation_sites/methylation.ht'

# utility functions

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


# steps in pipeline

def prepare_exomes(exome_ht: hl.Table, groupings: List, impose_high_af_cutoff_upfront: bool = True) -> hl.Table:
    '''Prepare exome mutation data for modelling'''

    # Manipulate VEP annotations and explode by them
    exome_ht = add_most_severe_csq_to_tc_within_ht(exome_ht)
    exome_ht = exome_ht.transmute(transcript_consequences=exome_ht.vep.transcript_consequences)
    exome_ht = exome_ht.explode(exome_ht.transcript_consequences)
    
    # Annotate variants with grouping variables. 
    exome_ht, grouping = annotate_constraint_groupings(exome_ht,groupings) # This function needs to be adapted
    exome_ht = exome_ht.select(
        'context', 'ref', 'alt', 'methylation_level', 'freq', 'pass_filters', *groupings)

    # Filter by allele count
    # Likely to need to adapt this function as well
    af_cutoff = 0.001
    freq_index = exome_ht.freq_index_dict.collect()[0][dataset]

    def keep_criteria(ht):
        crit = (ht.freq[freq_index].AC > 0) & ht.pass_filters & (ht.coverage > 0)
        if impose_high_af_cutoff_upfront:
            crit &= (ht.freq[freq_index].AF <= af_cutoff)
        return crit

    exome_ht = exome_ht.filter(keep_criteria(exome_ht))
    return exome_ht

def build_coverage_model(po_coverage_ht_path):
    '''Build linear regression model for variants observed'''
    coverage_ht = hl.read_table(po_coverage_ht_path)
    coverage_x_ht = hl.read_table(po_coverage_ht_path.replace('.ht', '_x.ht'))
    coverage_y_ht = hl.read_table(po_coverage_ht_path.replace('.ht', '_y.ht'))

    coverage_model, plateau_models = build_models(coverage_ht, args.trimers, True)
    _, plateau_x_models = build_models(coverage_x_ht, args.trimers, True)
    _, plateau_y_models = build_models(coverage_y_ht, args.trimers, True)

    return coverage_model, plateau_models

def get_proportion_observed(
        exome_ht: hl.Table, 
        possible_variants_ht: hl.Table,
        groupings: List,
        dataset: str = 'gnomad', 
        impose_high_af_cutoff_upfront: bool = True, 
        half_cutoff = False) -> hl.Table:
    '''Aggregate by grouping variables'''
    # Filter variant table and annotate with grouping variables
    exome_ht = prepare_exomes(exome_ht, groupings)  
    

    # Aggregate on grouping variables    
    # Should probably adapt this function
    ht = count_variants(exome_ht, additional_grouping=groupings, partition_hint=2000, force_grouping=True,
                        count_downsamplings=POPS, impose_high_af_cutoff_here=not impose_high_af_cutoff_upfront)
    ht = ht.join(possible_variants_ht, 'outer')
    agg_expr = {
        'variant_count': hl.agg.sum(ht.variant_count),
        'adjusted_mutation_rate': hl.agg.sum(ht.adjusted_mutation_rate),
        'possible_variants': hl.agg.sum(ht.possible_variants),
        'expected_variants': hl.agg.sum(ht.expected_variants),
        'mu': hl.agg.sum(ht.mu)
    }
    for pop in POPS:
        agg_expr[f'adjusted_mutation_rate_{pop}'] = hl.agg.array_sum(ht[f'adjusted_mutation_rate_{pop}'])
        agg_expr[f'expected_variants_{pop}'] = hl.agg.array_sum(ht[f'expected_variants_{pop}'])
        agg_expr[f'downsampling_counts_{pop}'] = hl.agg.array_sum(ht[f'downsampling_counts_{pop}'])
    ht = ht.group_by(*grouping).partition_hint(1000).aggregate(**agg_expr)
    return ht.annotate(obs_exp=ht.variant_count / ht.expected_variants)

def finalize_dataset(po_ht: hl.Table, keys: Tuple[str] = ('gene', 'transcript', 'canonical'),
                     n_partitions: int = 1000) -> hl.Table:
    '''aggregate variants to calculate constraint metrics and significance'''
    # Tidy this code a bit
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
    # Z score calculation not feasible with partial dataset
    # Need to include flagging of issues in constraint calculations
    # ht = calculate_all_z_scores(ht)
    return ht

def run_tests(ht):
    """Tests loading of autosome po table"""
    # Need to test access to gnomad data in GS buckets and make sure this doesn't cost too much
    # Then test that input keys match 
    ht = prepare_ht(ht)
    ht.show()

def main(args):
    # Set paths for data access based on command line parameters
    root = 'gs://gcp-public-data--gnomad/papers/2019-flagship-lof/v1.0/model/'
    
    processed_genomes_ht_path = f'{root}/model/genomes_processed.ht'
    processed_exomes_ht_path = f'{root}/model/exomes_processed.ht'

    coverage_autosomes = f'{root}/prop_observed_by_coverage_no_common_pass_filtered_bins.ht'
    coverage_x = f'{root}/prop_observed_by_coverage_no_common_pass_filtered_bins_x.ht/'
    coverage_y = f'{root}/prop_observed_by_coverage_no_common_pass_filtered_bins_y.ht/' 

    mutation_rate_ht_path = f'{root}/model/mutation_rate_methylation_bins.ht'
    po_coverage_ht_path = f'{root}/model/prop_observed_by_coverage_no_common_pass_filtered_bins.ht'
    po_ht_path = f'{root}/{{subdir}}/prop_observed_{{subdir}}.ht'
    raw_constraint_ht_path = f'{root}/{{subdir}}/constraint_{{subdir}}.ht'
    final_constraint_ht_path = f'{root}/{{subdir}}/constraint_final_{{subdir}}.ht'
    possible_variants_ht_path = f'{root}/model/possible_data/possible_transcript_pop_{args.model}.ht'
    
    variants_table_path = f'{root}/gnomad_v2.1.1_gpcr_variants_unfiltered.tsv'

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
        # Currently tests whether it can read data from exomes sites table
        ht_exomes_path = 'gs://gcp-public-data--gnomad/release/2.1/ht/exomes/gnomad.exomes.r2.1.sites.ht/'
        ht_exomes = hl.read_table(ht_exomes_path)
        ht_exomes.describe()
        # ht = load_variant_table(variants_table_path)
        # run_tests(ht)

    if args.get_proportion_observed:
        # Build a model for methylation-dependent mutation rate and apply it to get proportion of variants observed
        # Also need to incorporate genomes and v3 if possible

        # get model for coverage
        coverage_model, plateau_models = build_coverage_model(po_coverage_ht_path)
        print('Running aggregation of variants by grouping variables')
        # Tables of observed mutations in exomes
        full_exome_ht = hl.read_table(processed_exomes_ht_path)
        full_exome_ht = prepare_ht(full_exome_ht, args.trimers) 

        # filter into X, Y and autosomal regions
        exome_ht = full_exome_ht.filter(full_exome_ht.locus.in_autosome_or_par())
        exome_x_ht = hl.filter_intervals(full_exome_ht, [hl.parse_locus_interval('X')])
        exome_x_ht = exome_x_ht.filter(exome_x_ht.locus.in_x_nonpar())
        exome_y_ht = hl.filter_intervals(full_exome_ht, [hl.parse_locus_interval('Y')])
        exome_y_ht = exome_y_ht.filter(exome_y_ht.locus.in_y_nonpar())

        # Modelling results of estimated mutation rates for genes, coverage, methylation level and base context
        possible_variants_ht = hl.read_table(possible_variants_ht_path)

        # Set chosen groupings to aggregate on
        groupings = ['gene','annotation','modifier']

        # Apply model; aggregate by chosen groupings & get proportion observed; write to file
        po_exome_ht, po_exome_x_ht, po_exome_y_ht = \
            [get_proportion_observed(ht, possible_variants_ht, groupings) for ht in (exome_ht, exome_x_ht, exome_y_ht)]
        po_exome_ht.write(po_output_path, overwrite=args.overwrite)
        po_exome_x_ht.write(po_output_path.replace('.ht', '_x.ht'), overwrite=args.overwrite)
        po_exome_y_ht.write(po_output_path.replace('.ht', '_y.ht'), overwrite=args.overwrite)

       
    if args.aggregate:
        print('Running aggregation by gene')
        # read PO hail tables for autosomes, X and Y chromosomes and join them
        print(f'Reading hail table from {po_output_path}')
        # Autosomes
        ht = load_or_import(po_output_path, args.overwrite)
        # X chromosome
        ht_x = load_or_import(po_output_path.replace('.ht','_x.ht'), args.overwrite)
        # Y chromosome
        ht_y = load_or_import(po_output_path.replace('.ht','_y.ht'), args.overwrite)
        # Combine into one table
        ht = ht.union(ht_x).union(ht_y)
        # group by gene/transcript and calculate summary stats
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

# Summary of pipeline steps
# Get gene list & parameters and pull relevant data from (a) exomes/genomes variant call table (b) site context table (c) mutation rate table
# Filter exome and other tables for use and annotate 
# Pull all of coverage model table (prop_observed_by_coverage) and build a coverage model
# Calculate all possible variants in relevant genes with expected rate of observation
# Calculate proportion of variants observed by category
# Finalise and calculate summary stats for release

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--test', help='Run tests without actually requesting data',action='store_true')
    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--trimers', help='Use trimers instead of heptamers', action='store_true')
    parser.add_argument('--dataset', help='Which dataset to use (one of gnomad, non_neuro, non_cancer, controls)', default='gnomad')
    parser.add_argument('--model', help='Which model to apply (one of "standard", "syn_canonical", or "worst_csq" for now)', default='standard')
    parser.add_argument('--skip_af_filter_upfront', help='Skip AF filter up front (to be applied later to ensure that it is not affecting population-specific constraint): not generally recommended', action='store_true')
    parser.add_argument('--get_proportion_observed',help='Apply model to calculate proportion observed')
    parser.add_argument('--aggregate', help='Aggregate p_obs table', action='store_true')
    parser.add_argument('--summarise', help='Report summary stats', action='store_true')
    args = parser.parse_args()
    main(args)