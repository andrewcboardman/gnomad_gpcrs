def add_most_severe_csq_to_tc_within_ht(t):
    annotation = t.vep.annotate(transcript_consequences=t.vep.transcript_consequences.map(
        add_most_severe_consequence_to_consequence))
    return t.annotate_rows(vep=annotation) if isinstance(t, hl.MatrixTable) else t.annotate(vep=annotation)


def take_one_annotation_from_tc_within_ht(t):
    annotation = t.vep.annotate(transcript_consequences=t.vep.transcript_consequences[0])
    return t.annotate_rows(vep=annotation) if isinstance(t, hl.MatrixTable) else t.annotate(vep=annotation)

def load_tx_expression_data(context=True):
    tx_ht_file = 'gs://gnomad-public/papers/2019-tx-annotation/pre_computed/all.possible.snvs.tx_annotated.021819.ht' if context \
        else 'gs://gnomad-public/papers/2019-tx-annotation/gnomad.exomes.r2.1.1.sites.tx_annotated.021319.ht'
    tx_ht = hl.read_matrix_table(tx_ht_file).rows()

    def process_expression_data(csq_expression):
        exprs_to_drop = ['ensg', 'csq', 'symbol', 'lof', 'lof_flag', 'mean_proportion']
        expression_data = csq_expression.drop(*exprs_to_drop)
        all_tissues = list(expression_data.values())
        expression_data_list = list(zip(list(expression_data), all_tissues))
        brain_tissues = [x[1] for x in expression_data_list if 'Brain' in x[0]]
        return csq_expression.select('ensg', 'csq', 'symbol', 'lof',
                                     mean_expression=hl.mean(hl.filter(lambda e: ~hl.is_nan(e), all_tissues), filter_missing=True),
                                     max_expression=hl.max(hl.filter(lambda e: ~hl.is_nan(e), all_tissues), filter_missing=True),
                                     mean_brain_expression=hl.mean(hl.filter(lambda k: ~hl.is_nan(k), brain_tissues), filter_missing=True),
                                     Brain_Cortex=csq_expression.Brain_Cortex
                                     )

    return tx_ht.annotate(tx_annotation=tx_ht.tx_annotation.map(process_expression_data))

def annotate_constraint_groupings(ht: Union[hl.Table, hl.MatrixTable],
                                  custom_model: str = None) -> Tuple[Union[hl.Table, hl.MatrixTable], List[str]]:
    """
    HT must be exploded against whatever axis

    Need to add `'coverage': ht.exome_coverage` here (which will get corrected out later)
    """
    if custom_model == 'worst_csq':
        groupings = {
            'annotation': ht.worst_csq_by_gene.most_severe_consequence,
            'modifier': hl.case()
                .when(hl.is_defined(ht.worst_csq_by_gene.lof),
                      ht.worst_csq_by_gene.lof)
                .when(hl.is_defined(ht.worst_csq_by_gene.polyphen_prediction),
                      ht.worst_csq_by_gene.polyphen_prediction)
                .default('None'),
            'gene': ht.worst_csq_by_gene.gene_symbol,
            'coverage': ht.exome_coverage
        }
    elif custom_model == 'tx_annotation':
        groupings = {
            'annotation': ht.tx_annotation.csq,
            'modifier': ht.tx_annotation.lof,
            'gene': ht.tx_annotation.symbol,
            'expressed': hl.case(missing_false=True).when(
                ht.tx_annotation.mean_expression >= 0.9, 'high').when(
                ht.tx_annotation.mean_expression > 0.1, 'medium').when(
                hl.is_defined(ht.tx_annotation.mean_expression), 'low').default('missing'),
            'coverage': ht.exome_coverage
        }
    else:
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
        if custom_model == 'splice_region':
            groupings['distance_splice'] = ht.transcript_consequences
    ht = ht.annotate(**groupings) if isinstance(ht, hl.Table) else ht.annotate_rows(**groupings)
    return ht, list(groupings.keys())


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

def prepare_ht(ht, trimer: bool = False, annotate_coverage: bool = True):
    if trimer:
        ht = trimer_from_heptamer(ht)
    str_len = 3 if trimer else 7

    if isinstance(ht, hl.Table):
        ht = ht.annotate(ref=ht.alleles[0], alt=ht.alleles[1])
        ht = ht.filter((hl.len(ht.ref) == 1) & (hl.len(ht.alt) == 1) & ht.context.matches(f'[ATCG]{{{str_len}}}'))
        ht = annotate_variant_types(collapse_strand(ht), not trimer)
    else:
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


def build_plateau_models_pop(ht: hl.Table, weighted: bool = False) -> Dict[str, Tuple[float, float]]:
    """
    Calibrates high coverage model (returns intercept and slope)
    """
    pop_lengths = get_all_pop_lengths(ht)
    agg_expr = {
        pop: [hl.agg.group_by(ht.cpg,
                             hl.agg.linreg(ht[f'observed_{pop}'][i] / ht.possible_variants, [1, ht.mu_snp],
                                           weight=ht.possible_variants if weighted else None)
                             ).map_values(lambda x: x.beta) for i in range(length)]
        for length, pop in pop_lengths
    }
    agg_expr['total'] = hl.agg.group_by(ht.cpg,
                                        hl.agg.linreg(ht.observed_variants / ht.possible_variants, [1, ht.mu_snp],
                                                      weight=ht.possible_variants if weighted else None)
                                        ).map_values(lambda x: x.beta)
    return ht.aggregate(hl.struct(**agg_expr))

# Model building
def build_coverage_model(coverage_ht: hl.Table) -> (float, float):
    """
    Calibrates coverage model (returns intercept and slope)
    """
    return tuple(coverage_ht.aggregate(hl.agg.linreg(coverage_ht.low_coverage_obs_exp, [1, coverage_ht.log_coverage])).beta)

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