def aggregate_by_groupings(paths_dict, data_dict):
    '''This is the new master function for performing constraint analysis'''
    # aggregate by chosen groupings & get proportion observed; write to file
    input_hts = zip(
        (exome_ht, exome_x_ht, exome_y_ht), 
        (context_ht, context_x_ht, context_y_ht),
        (
            po_output_path, 
            po_output_path.replace('.ht', '_x.ht'),
            po_output_path.replace('.ht', '_y.ht')
        )
    )
    for exome_ht_, context_ht_, po_output_path_ in input_hts:
        po_exome_ht_ = get_proportion_observed(exome_ht_,context_ht_,grouping)
        po_exome_ht_.write(po_output_path_, overwrite=args.overwrite)




def get_possible_variants(context_ht,mutation_ht, coverage_model, plateau_models, grouping, possible_file,
        half_cutoff = False, HIGH_COVERAGE_CUTOFF = 40, POPS = ('global', 'afr', 'amr', 'eas', 'nfe', 'sas')):
    '''Compute table of possible variants with needed properties. Untested'''
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
    for pop in POPS:
        ann_expr[f'expected_variants_{pop}'] = ht[f'adjusted_mutation_rate_{pop}'] * ht.coverage_correction
    ht = ht.annotate(**ann_expr)
    ht.write(possible_file, True)
    return ht
    

def get_proportion_observed(
        exome_ht: hl.Table, 
        context_ht: hl.Table,
        mutation_ht: hl.Table,
        coverage_model, plateau_models,
        possible_variants_ht_path,
        grouping: List,
        impose_high_af_cutoff_upfront: bool = True) -> hl.Table:
    '''Aggregate by grouping variables'''

    # Get possible variants for this context table
    possible_variants_ht = get_possible_variants(context_ht,mutation_ht, coverage_model, plateau_models, grouping, possible_variants_ht_path)

    # Save possible variants??
    #possible_variants_ht_path.write(possible_variants_ht_path)

    # Count observed variants by grouping
    ht = count_variants(exome_ht, additional_grouping=grouping, partition_hint=2000, force_grouping=True,
                        count_downsamplings=False, impose_high_af_cutoff_here=not impose_high_af_cutoff_upfront,
                        POPS = ('global', 'afr', 'amr', 'eas', 'nfe', 'sas'))

    # Match with expected frequency of variants by grouping
    ht = ht.join(possible_variants_ht, 'outer')

    # Aggregate all groupings across populations
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