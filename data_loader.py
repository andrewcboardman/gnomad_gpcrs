def get_data(paths, data, gene_intervals):
    '''This is the new master function for loading all necessary data for constraint analysis on the given genes'''
    # Get paths
    exomes_path = paths['exomes_path']
    exomes_local_path = paths['exomes_local_path']
    context_path = paths['context_path']
    context_local_path = paths['context_local_path']
    mutation_rate_path = paths['mutation_rate_path']
    mutation_rate_local_path = paths['mutation_rate_local_path']
    po_coverage_local_path = paths['po_coverage_local_path']

    # Get relevant exomes data, filter and split  
    full_exome_ht = get_exomes(exomes_path, gene_intervals)
    exome_ht, exome_x_ht, exome_y_ht = filter_exomes(full_exome_ht, exomes_local_path, args.overwrite)
    data.update(dict(zip(
        ('exome_ht', 'exome_x_ht', 'exome_y_ht'),
        (exome_ht, exome_x_ht, exome_y_ht)
    )))

    # Do the same for nucleotide context
    full_context_ht = get_context(context_path,gene_intervals)
    context_ht,context_x_ht, context_y_ht = filter_context(full_context_ht, exomes_local_path, args.overwrite)
    data.update(dict(zip(
        ('context_ht', 'context_x_ht', 'context_y_ht'),
        (context_ht, context_x_ht, context_y_ht)
    )))
    
    # Get methylation-dependent mutation rate
    mutation_rate_ht = hl.read_table(mutation_rate_path).select('mu_snp')
    mutation_rate_ht.write(mutation_rate_local_path,overwrite=args.overwrite)
    data['mutation_rate_ht'] = mutation_rate_ht

    # Build model for observation rate based on coverage
    coverage_ht = hl.read_table(po_coverage_local_path)
    coverage_model, plateau_models = build_models(coverage_ht, args.trimers, True) 
    with open('data/coverage_models.pkl','wb') as fid:
        pickle.dump((coverage_model,plateau_models),fid)  
    data.update(dict(zip(
        ('coverage_model', 'plateau_models'),
        (coverage_model, plateau_models)
    )))

    return data


def get_exomes(exomes_path, gene_intervals):
    # Path to full exome table
    full_exome_ht = hl.read_table(exomes_path)
    # Select relevant columns to avoid getting too much data 
    full_exome_ht = full_exome_ht.select(
        full_exome_ht.freq,
        full_exome_ht.vep.transcript_consequences,
        full_exome_ht.context,
        full_exome_ht.methylation,
        full_exome_ht.coverage
    )
    # Filter by gene
    exome_ht = hl.filter_intervals(full_exome_ht, gene_intervals)
    # Prepare exomes
    exome_ht = prepare_ht(exome_ht)
    return exome_ht


def get_context(context_ht_path, gene_intervals):
    full_context_ht = hl.read_table(context_ht_path)
    full_context_ht = full_context_ht.select(
        full_context_ht.vep.transcript_consequences,
        full_context_ht.context,
        full_context_ht.methylation,
        full_context_ht.coverage
    )
    # Filter this table in similar way to exomes ht
    selected_context_ht = hl.filter_intervals(full_context_ht,gene_intervals)
    selected_context_ht = prepare_ht(selected_context_ht)
    return selected_context_ht

def filter_exomes(full_exome_ht, exomes_local_path, overwrite):
    # filter into X, Y and autosomal regions
    exome_ht = full_exome_ht.filter(full_exome_ht.locus.in_autosome_or_par())
    exome_x_ht = hl.filter_intervals(full_exome_ht, [hl.parse_locus_interval('X')])
    exome_x_ht = exome_x_ht.filter(exome_x_ht.locus.in_x_nonpar())
    exome_y_ht = hl.filter_intervals(full_exome_ht, [hl.parse_locus_interval('Y')])
    exome_y_ht = exome_y_ht.filter(exome_y_ht.locus.in_y_nonpar())
    exome_ht.write(exomes_local_path,overwrite=overwrite)
    exome_x_ht.write(exomes_local_path.replace('.ht', '_x.ht'),overwrite=overwrite)
    exome_y_ht.write(exomes_local_path.replace('.ht', '_y.ht'),overwrite=overwrite)
    return exome_ht, exome_x_ht, exome_y_ht


def filter_context(full_context_ht, context_local_path, overwrite):
    context_ht = full_context_ht.filter(full_context_ht.locus.in_autosome_or_par())
    context_x_ht = hl.filter_intervals(full_context_ht, [hl.parse_locus_interval('X')])
    context_x_ht = context_x_ht.filter(context_x_ht.locus.in_x_nonpar())
    context_y_ht = hl.filter_intervals(full_context_ht, [hl.parse_locus_interval('Y')])
    context_y_ht = context_y_ht.filter(context_y_ht.locus.in_y_nonpar())
    context_ht.write(context_local_path,overwrite=overwrite)
    context_x_ht.write(context_local_path.replace('.ht', '_x.ht'),overwrite=overwrite)
    context_y_ht.write(context_local_path.replace('.ht', '_y.ht'),overwrite=overwrite)
    return context_ht, context_x_ht, context_y_ht


def prepare_exomes_and_context(
        exome_ht: hl.Table,
        context_ht: hl.Table,
        dataset: str = 'gnomad', 
        impose_high_af_cutoff_upfront: bool = True
    ) -> hl.Table:
    '''Prepare exome mutation data for modelling. Untested and currently not used. Why?'''
    # For standard analysis, assume the worst predicted transcript consequence takes place
    exome_ht = add_most_severe_csq_to_tc_within_ht(exome_ht)
    exome_ht = exome_ht.transmute(transcript_consequences=exome_ht.vep.transcript_consequences)
    exome_ht = exome_ht.explode(exome_ht.transcript_consequences)
    
    # Annotate variants with grouping variables. 
    exome_ht, grouping = annotate_constraint_groupings(exome_ht)
    exome_ht = exome_ht.select('context', 'ref', 'alt', 'methylation_level', 'freq', 'pass_filters', *grouping)

    # annotate and filter context table in same way. Exome coverage or genome coverage?
    context_ht, grouping = annotate_constraint_groupings(context_ht)
    context_ht = context_ht.filter(hl.is_defined(context_ht.exome_coverage)).select(
        'context', 'ref', 'alt', 'methylation_level', *grouping)

    # Filter by allele count > 0, allele fraction < cutoff, coverage > 0
    af_cutoff = 0.001
    freq_index = exome_ht.freq_index_dict.collect()[0][dataset]
    def keep_criteria(ht):
        crit = (ht.freq[freq_index].AC > 0) & ht.pass_filters & (ht.coverage > 0)
        if impose_high_af_cutoff_upfront:
            crit &= (ht.freq[freq_index].AF <= af_cutoff)
        return crit
    exome_ht = exome_ht.filter(keep_criteria(exome_ht))

    return exome_ht, context_ht, grouping
