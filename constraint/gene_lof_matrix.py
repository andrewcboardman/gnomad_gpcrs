#!/usr/bin/env python3

__author__ = 'konradk'

from constraint_utils import *

gene_lof_matrix_path = 'gs://gnomad-resources/lof_paper/individual_level/full_gene_lof_matrix_{}.mt'
homozygous_lof_mt_path_raw = 'gs://gnomad-resources/lof_paper/individual_level/all_homozygous_lofs{data_type}.mt'
lofs_by_gene_ht_path_raw = f'{root}/homozygous_lof_summary{{data_type}}.ht'
coverage_summary_ht_path = f'{root}/coverage_summary.ht'


def load_gtf_data():
    ht = hl.experimental.import_gtf('gs://hail-common/references/gencode/gencode.v19.annotation.gtf.bgz', 'GRCh37', True, min_partitions=12)
    ht = ht.annotate(gene_id=ht.gene_id.split('\\.')[0],
                     transcript_id=ht.transcript_id.split('\\.')[0],
                     length=ht.interval.end.position - ht.interval.start.position + 1)
    genes = ht.filter(ht.feature == 'gene').select(
        'gene_id', 'gene_type', 'gene_name', 'length').rename({'length': 'gene_length'}).key_by('gene_id')
    coding_regions = ht.filter(ht.feature == 'CDS').select('gene_id', 'transcript_id', 'transcript_type', 'length', 'level')
    transcripts = coding_regions.group_by('transcript_id', 'transcript_type', 'gene_id',
                                          transcript_level=coding_regions.level).aggregate(
        cds_length=hl.agg.sum(coding_regions.length),
        num_coding_exons=hl.agg.count()
    ).key_by('transcript_id')
    return transcripts.annotate(**genes[transcripts.gene_id])


def create_coverage_summary(overwrite: bool = False):
    gtf = hl.experimental.import_gtf('gs://hail-common/references/gencode/gencode.v19.annotation.gtf.bgz', 'GRCh37', True, min_partitions=100)
    gtf = gtf.annotate(gene_id=gtf.gene_id.split('\\.')[0],
                       transcript_id=gtf.transcript_id.split('\\.')[0])
    gtf = gtf.filter(gtf.feature == 'CDS').select('gene_id', 'transcript_id')
    gtf = gtf.checkpoint(hl.utils.new_temp_file())
    ht = hl.read_table(coverage_ht_path('exomes'))
    ht = ht.annotate(cds=gtf.index(ht.locus, all_matches=True)).explode('cds')
    ht = ht.group_by(gene_id=ht.cds.gene_id, transcript_id=ht.cds.transcript_id).aggregate(
        mean_coverage=hl.agg.mean(ht.mean),
        median_coverage=hl.agg.approx_quantiles(ht.median, 0.5)
    )
    return ht.checkpoint(coverage_summary_ht_path, overwrite=overwrite, _read_if_exists=not overwrite)


def load_gene_expression_data():
    ht = hl.read_table('gs://gnomad-public/papers/2019-tx-annotation/data/GTEx.V7.tx_medians.110818.ht')
    # Filtering out tissues with < 100 samples
    brain_tissues = ht.row.values.filter(lambda x: x.tissue.contains('Brain') &
                                                   (x.tissue != 'Brain-Spinalcord(cervicalc-1)') &
                                                   (x.tissue != 'Brain-Substantianigra')).map(
        lambda y: y.tx_expr_summary)
    return ht.select(
        brain_expression=hl.mean(hl.filter(lambda k: ~hl.is_nan(k), brain_tissues), filter_missing=True)
    ).key_by('transcript_id')


def load_exac_pli_data():
    exac_pli = hl.import_table(f'{root}/old_exac_data/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt.gz', impute=True, force_bgz=True)
    return exac_pli.key_by(transcript=exac_pli.transcript.split('\\.')[0])


def add_rank(ht, field, ascending=True, total_genes=None, bins=10, defined_only=False):
    if total_genes is None:
        if defined_only:
            total_genes = ht.aggregate(hl.agg.count_where(hl.is_defined(ht[field])))
        else:
            total_genes = ht.count()
    rank_field = ht[field] if ascending else -ht[field]
    ht = ht.order_by(rank_field).add_index(f'{field}_rank')
    ht = ht.annotate(**{f'{field}_rank': hl.or_missing(
        hl.is_defined(ht[field]), ht[f'{field}_rank']
    )})
    return ht.annotate(**{
        f'{field}_bin': hl.int(ht[f'{field}_rank'] * bins / total_genes),
        f'{field}_bin_6': hl.int(ht[f'{field}_rank'] * 6 / total_genes)
    })


def select_primitives_from_ht(ht):
    return ht.select(**{x: v for x, v in ht.row_value.items() if
                        v.dtype in {hl.tstr, hl.tint32, hl.tfloat32, hl.tint64, hl.tfloat64, hl.tbool}})


def generate_gene_lof_matrix(mt: hl.MatrixTable, tx_ht: hl.Table, by_transcript: bool = False,
                             filter_an_adj: bool = False, common_var_filt: bool = False, pre_loftee: bool = False,
                             remove_ultra_common: bool = False
                             ) -> hl.MatrixTable:
    filt_criteria = hl.len(mt.filters) == 0
    if filter_an_adj:
        filt_criteria &= get_an_adj_criteria(mt)
    if remove_ultra_common:
        filt_criteria &= mt.freq[0].AF < 0.95
    if common_var_filt:
        filt_criteria &= mt.freq[0].AF < 0.05
    mt = mt.filter_rows(filt_criteria)
    if by_transcript:
        explode_field = 'transcript_consequences'
    else:
        mt = process_consequences(mt)
        explode_field = 'worst_csq_by_gene'
    if pre_loftee:
        lof_cats = hl.literal({"splice_acceptor_variant", "splice_donor_variant", "stop_gained", "frameshift_variant"})
        criteria = lambda x: lof_cats.contains(add_most_severe_consequence_to_consequence(x).most_severe_consequence)
    else:
        criteria = lambda x: (x.lof == 'HC') & hl.is_missing(x.lof_flags)
    lof_csqs = mt.vep[explode_field].filter(criteria)
    mt = mt.select_rows(mt.freq, lof_csqs=lof_csqs)
    mt = mt.explode_rows(mt.lof_csqs)
    annotation_expr = {
        'gene_id': mt.lof_csqs.gene_id,
        'gene': mt.lof_csqs.gene_symbol,
        'indel': hl.is_indel(mt.alleles[0], mt.alleles[1])
    }
    if by_transcript:
        annotation_expr['transcript_id'] = mt.lof_csqs.transcript_id
        annotation_expr['canonical'] = hl.is_defined(mt.lof_csqs.canonical)
    else:
        tx_annotation = annotate_tx_expression_data(mt, tx_ht, mt.lof_csqs).mean_expression
        annotation_expr['expressed'] = hl.case().when(
            tx_annotation >= 0.9, 'high').when(
            tx_annotation > 0.1, 'medium').when(
            hl.is_defined(tx_annotation), 'low').default('missing')
    mt = mt.annotate_rows(**annotation_expr)
    mt.describe()
    return mt.group_rows_by(*list(annotation_expr.keys())).aggregate_rows(
        n_sites=hl.agg.count(),
        n_sites_array=hl.agg.array_sum(mt.freq.map(lambda x: hl.int(x.AC > 0))),
        classic_caf=hl.agg.sum(mt.freq[0].AF),
        max_af=hl.agg.max(mt.freq[0].AF),
        classic_caf_array=hl.agg.array_sum(mt.freq.map(lambda x: x.AF)),
    ).aggregate_entries(
        num_homs=hl.agg.count_where(mt.GT.is_hom_var()),
        num_hets=hl.agg.count_where(mt.GT.is_het()),
        defined_sites=hl.agg.count_where(hl.is_defined(mt.GT))
    ).result()


def generate_gene_lof_summary(mt, collapse_indels: bool = False, by_transcript: bool = False):
    if collapse_indels:
        grouping = ['gene_id', 'gene']
        if by_transcript:
            grouping.extend(['transcript_id', 'canonical'])
        else:
            grouping.append('expressed')
        mt = mt.group_rows_by(*grouping).aggregate_rows(
            n_sites=hl.agg.sum(mt.n_sites),
            n_sites_array=hl.agg.array_sum(mt.n_sites_array),
            classic_caf=hl.agg.sum(mt.classic_caf),
            max_af=hl.agg.max(mt.max_af),
            classic_caf_array=hl.agg.array_sum(mt.classic_caf_array)
        ).aggregate_entries(
            num_homs=hl.agg.sum(mt.num_homs),
            num_hets=hl.agg.sum(mt.num_hets),
            defined_sites=hl.agg.sum(mt.defined_sites)
        ).result()
    ht = mt.annotate_rows(no_lofs=hl.agg.count_where((mt.defined_sites > 0) & (mt.num_homs + mt.num_hets == 0)),
                          obs_het_lof=hl.agg.count_where((mt.num_hets > 0) & (mt.num_homs == 0)),
                          obs_hom_lof=hl.agg.count_where(mt.num_homs > 0),
                          defined=hl.agg.count_where(mt.defined_sites > 0),
                          pop_no_lofs=hl.agg.group_by(
                              mt.meta.pop, hl.agg.count_where((mt.defined_sites > 0) & (mt.num_homs + mt.num_hets == 0))),
                          pop_obs_het_lof=hl.agg.group_by(
                              mt.meta.pop, hl.agg.count_where((mt.num_hets > 0) & (mt.num_homs == 0))),
                          pop_obs_hom_lof=hl.agg.group_by(
                              mt.meta.pop, hl.agg.count_where(mt.num_homs > 0)),
                          pop_defined=hl.agg.group_by(
                              mt.meta.pop, hl.agg.count_where(mt.defined_sites > 0)),
                          ).rows()
    ht = ht.annotate(p=1 - hl.sqrt(hl.float64(ht.no_lofs) / ht.defined),
                     pop_p=hl.dict(hl.array(ht.pop_defined).map(
                         lambda x: (x[0], 1 - hl.sqrt(hl.float64(ht.pop_no_lofs.get(x[0])) / x[1])))))
    ht = ht.annotate(exp_hom_lof=ht.defined * ht.p * ht.p)
    return ht.annotate(oe=ht.obs_hom_lof / ht.exp_hom_lof)


def combine_lof_metrics(gene_lof_matrix, by_transcript: bool = False, pop_specific: bool = False):
    keys = ['gene']
    caf_keys = ['gene']
    caf_drop = ['gene_id', 'oe']
    constraint_path = raw_constraint_ht_path.format(subdir="standard" if by_transcript else "tx_annotation")
    if pop_specific: constraint_path = constraint_path.replace('standard/', 'pop_specific/standard/')
    ht = hl.read_table(constraint_path)
    if by_transcript:
        keys.append('transcript')
        caf_keys.append('transcript_id')
        caf_drop.append('canonical')
    else:
        keys.append('expressed')
        caf_keys.append('expressed')
    constraint_ht = ht.drop(*[x for x in ht.row if x.endswith('_classic') or x.endswith('_with_os')])
    constraint_ht = add_rank(constraint_ht, 'oe_lof_upper', defined_only=True).key_by(*keys)
    caf_ht = hl.read_table(gene_lof_matrix.replace('.mt', '.summary.ht'))
    mapping = caf_ht.freq_index_dict.collect()[0]
    caf_dict = {
        f'classic_caf_{pop}': caf_ht.classic_caf_array[mapping[f'gnomad_{pop}']] for pop in
        map(lambda x: x.lower(), EXOME_POPS)
    }
    caf_dict.update({
        f'p_{pop}': caf_ht.pop_p[pop] for pop in map(lambda x: x.lower(), EXOME_POPS)
    })
    caf_ht = caf_ht.annotate(**caf_dict)
    caf_ht = caf_ht.key_by(*caf_keys).drop(*caf_drop)
    exac_pli = load_exac_pli_data()
    gene_ht = load_gtf_data()
    # expr_ht = load_gene_expression_data()
    coverage_ht = hl.read_table(coverage_summary_ht_path)
    if by_transcript:
        exac = exac_pli[constraint_ht.transcript]
        gene = gene_ht.drop('gene_name')[constraint_ht.transcript]
        # expr_ht = expr_ht.drop('gene_id')
        coverage = coverage_ht[constraint_ht.key]
    else:
        exac = exac_pli.key_by('gene')[constraint_ht.gene]
        gene = gene_ht.key_by('gene_name').drop('transcript_id')[constraint_ht.gene]
        # expr_ht = expr_ht.key_by('gene_id').drop('transcript_id')
        coverage = coverage_ht.key_by('transcript_id').drop('gene_id')[constraint_ht.transcript]
    constraint_ht = constraint_ht.annotate(
        **caf_ht[constraint_ht.key], **gene, **coverage,
        exac_pLI=exac.pLI, exac_obs_lof=exac.n_lof, exac_exp_lof=exac.exp_lof, exac_oe_lof=exac.n_lof / exac.exp_lof
    )
    # If CAF data is missing and LoFs are possible in the gene, set all CAF metrics to zero
    constraint_ht = constraint_ht.annotate(
        **{x: hl.cond(hl.is_missing(constraint_ht[x]) & (constraint_ht.possible_lof > 0),
                      0,
                      constraint_ht[x])
           for x in list(caf_ht.row)
           if x in constraint_ht.row and constraint_ht[x].dtype in {hl.tint32, hl.tfloat32, hl.tint64, hl.tfloat64}}
    )
    return constraint_ht.annotate(
        # **expr_ht[constraint_ht.gene_id]
    ).annotate_globals(freq_meta=caf_ht.index_globals().freq_meta, freq_index_dict=caf_ht.index_globals().freq_index_dict)


def find_in_meta(downsampling, pop, freq_meta):
    for i, x in enumerate(freq_meta):
        if x.get('group') == 'adj' and \
                x.get('pop') == pop and \
                x.get('downsampling') == str(downsampling):
            return i
    print(f'failed on {downsampling} and {pop}')
    return None


def explode_downsamplings(ht, select_canonical=True):
    pop_lengths = get_all_pop_lengths(ht, prefix='exp_lof_')
    downsamplings = list(map(lambda x: x[1], get_downsamplings(ht)))
    freq_meta = ht.freq_meta.collect()[0]
    meta_locs = [(i, pop, downsamplings[i], find_in_meta(downsamplings[i], pop, freq_meta)) for count, pop in pop_lengths for i in range(count)]
    fields = ['data']
    if select_canonical: fields.insert(0, 'canonical')
    ht = ht.annotate(
        data=[
            hl.struct(
                pop=pop, downsampling=downsampling,
                exp_syn=ht[f'exp_syn_{pop}'][i], obs_syn=hl.cond(hl.is_missing(ht[f'obs_syn_{pop}'][i]) & hl.is_defined(ht.obs_syn), 0, ht[f'obs_syn_{pop}'][i]),
                exp_mis=ht[f'exp_mis_{pop}'][i], obs_mis=hl.cond(hl.is_missing(ht[f'obs_mis_{pop}'][i]) & hl.is_defined(ht.obs_mis), 0, ht[f'obs_mis_{pop}'][i]),
                exp_lof=ht[f'exp_lof_{pop}'][i], obs_lof=hl.cond(hl.is_missing(ht[f'obs_lof_{pop}'][i]) & hl.is_defined(ht.obs_lof), 0, ht[f'obs_lof_{pop}'][i]),
                n_sites=ht.n_sites_array[meta_loc],
                caf=ht.classic_caf_array[meta_loc]
            )
            for i, pop, downsampling, meta_loc in meta_locs] + [
            hl.struct(pop='global', downsampling=125748,
                      exp_syn=ht.exp_syn, obs_syn=ht.obs_syn,
                      exp_mis=ht.exp_mis, obs_mis=ht.obs_mis,
                      exp_lof=ht.exp_lof, obs_lof=ht.obs_lof,
                      n_sites=ht.n_sites, caf=ht.classic_caf)
        ]
    ).select(*fields)
    ht = ht.annotate()
    ht = ht.explode('data')
    ht = ht.transmute(**ht.data)
    return oe_confidence_interval(ht, ht.obs_lof, ht.exp_lof, prefix='oe_lof', select_only_ci_metrics=False)


def compute_homozygous_lof(mt, tx_ht):
    mt = mt.filter_rows((hl.len(mt.filters) == 0) & (mt.freq[0].homozygote_count >= 1))
    mt = process_consequences(mt)
    mt = mt.explode_rows(mt.vep.worst_csq_by_gene)
    mt = mt.filter_rows((mt.vep.worst_csq_by_gene.lof == 'HC') & hl.is_missing(mt.vep.worst_csq_by_gene.lof_flags))
    mt = mt.filter_rows(hl.agg.any(mt.GT.is_hom_var() &
                                   (mt.AD[1] / hl.sum(mt.AD) >= 0.8)
                                   ))
    mt = mt.annotate_rows(indel=hl.is_indel(mt.alleles[0], mt.alleles[1]),
                          tx_annotation=annotate_tx_expression_data(mt, tx_ht, mt.vep.worst_csq_by_gene))
    return mt.repartition(100)


def get_gene_lof_path(args):
    extension = 'by_transcript' if args.by_transcript else 'worst'
    if args.filter_an_adj:
        extension += '_an_adj'
    if args.no_loftee:
        extension += '_no_loftee'
    if args.remove_ultra_common:
        extension += '_remove_ultra_common'
    return gene_lof_matrix_path.format(extension)


def get_release_subdir(args):
    subdir = []
    if not args.by_transcript:
        subdir.append('worst_csq')
    if not args.filter_an_adj:
        subdir.append('no_an_adj')
    if args.no_loftee:
        subdir.append('no_loftee')
    if args.pop_specific:
        subdir.append('pop_specific')
    if args.remove_ultra_common:
        subdir.append('remove_ultra_common')
    subdir = '.'.join(subdir)
    if subdir: subdir = f'other_cuts/{subdir}/'
    return subdir


def main(args):
    hl.init(log='/constraint.log')

    gene_lof_matrix = get_gene_lof_path(args)
    subdir = get_release_subdir(args)

    root = 'gs://gnomad-public/papers/2019-flagship-lof/v1.0'
    all_lof_metrics_path = f'{root}/{subdir}gnomad.v2.1.1.lof_metrics.{{extension}}'

    tx_ht = load_tx_expression_data(context=args.genomes)
    if args.calculate_gene_lof_matrix:
        mt = get_gnomad_data('exomes', adj=True, release_samples=True, release_annotations=True)
        mt = generate_gene_lof_matrix(mt, tx_ht, by_transcript=args.by_transcript,
                                      filter_an_adj=args.filter_an_adj, pre_loftee=args.no_loftee,
                                      remove_ultra_common=args.remove_ultra_common)
        mt.write(gene_lof_matrix, args.overwrite)
        send_message(args.slack_channel, 'Gene LoF matrix computed!')

    if args.export_gene_lof_matrix:
        mt = hl.read_matrix_table(gene_lof_matrix)
        ht = generate_gene_lof_summary(mt, not args.dont_collapse_indels, args.by_transcript)
        if args.dont_collapse_indels:
            gene_lof_matrix = gene_lof_matrix.replace('.mt', '_by_indel.mt')
        ht.write(gene_lof_matrix.replace('.mt', '.summary.ht'), args.overwrite)
        hl.read_table(gene_lof_matrix.replace('.mt', '.summary.ht')).drop('classic_caf_array').export(
            gene_lof_matrix.replace('.mt', '.summary.txt.bgz'))

    if args.create_coverage_summary:
        create_coverage_summary(args.overwrite)

    # hail gene_lof_matrix.py --combine_lof_metrics --by_transcript
    if args.combine_lof_metrics:
        ht = combine_lof_metrics(gene_lof_matrix, args.by_transcript, pop_specific=args.pop_specific)
        ht.write(all_lof_metrics_path.format(extension='by_transcript.ht'), args.overwrite)
        # This file has been spot-checked. Of the canonical transcripts,
        # # 10.5% are missing CAF data due to having 0 HC no-flag LoFs
        # # # Many are actually 0, some e.g. single exon flag remove all LoFs
        # # # Some additional ones (SETD1B, TCF4, etc) are DQ'ed bc of AN_Adj
        # # 7.5% are missing ExAC pLI data (AC genes, low coverage, etc)
        # # 2.7% are missing LoF data (only zero LoF possible genes), confirmed by:
        # # # cht.filter(hl.is_missing(cht.obs_lof) & ~cht.gene_issues.contains('no_exp_lof')).show(10)
        # # 1 is missing oe_syn CIs: tiny Y chromosome gene with negative expected
        # # Checked that number of exons is similar to ExAC (minor differences = UTRs)

    if args.export_combined_metrics:
        ht = hl.read_table(all_lof_metrics_path.format(extension='by_transcript.ht'))
        ht = ht.annotate(constraint_flag=hl.delimit(ht.constraint_flag, '|'),
                         chromosome=ht.interval.start.contig,
                         start_position=ht.interval.start.position,
                         end_position=ht.interval.end.position)
        ht = select_primitives_from_ht(ht)
        ht.export(all_lof_metrics_path.format(extension='by_transcript.txt.bgz'))
        if args.by_transcript:
            ht = hl.read_table(all_lof_metrics_path.format(extension='by_transcript.ht'))
            ht = ht.filter(ht.canonical)
            ht = add_rank(ht, 'oe_lof_upper', defined_only=True).drop('canonical').key_by('gene_id', 'gene')
            ht.write(all_lof_metrics_path.format(extension='by_gene.ht'), args.overwrite)
            ht = select_primitives_from_ht(ht)
            ht.export(all_lof_metrics_path.format(extension='by_gene.txt.bgz'))

        ht = hl.read_table(all_lof_metrics_path.format(extension='by_transcript.ht'))
        ht = explode_downsamplings(ht, args.by_transcript)
        ht.export(all_lof_metrics_path.format(extension='downsamplings.txt.bgz').replace('v1.0', 'v1.1'))

    data_type = 'genomes' if args.genomes else 'exomes'
    homozygous_lof_mt_path = homozygous_lof_mt_path_raw.format(data_type='_genomes' if args.genomes else '')
    if args.compute_homozygous_lof:
        mt = get_gnomad_data(data_type, adj=True, non_refs_only=True,
                             release_samples=True, release_annotations=True)
        ht = compute_homozygous_lof(mt, tx_ht)
        ht.write(homozygous_lof_mt_path, args.overwrite)
        send_message(args.slack_channel, 'Homozygous LoFs computed!')

    lofs_by_gene_ht_path = lofs_by_gene_ht_path_raw.format(data_type='_genomes' if args.genomes else '')
    if args.export_homozygous_lof:
        mt = hl.read_matrix_table(homozygous_lof_mt_path)
        if args.include_low_pext:
            homozygous_lof_mt_path = homozygous_lof_mt_path.replace('.mt', '.low_pext.mt')
            lofs_by_gene_ht_path = lofs_by_gene_ht_path.replace('.ht', '.low_pext.ht')
        vep_ht = hl.read_table(annotations_ht_path(data_type, 'vep_csq'))
        if not args.include_low_pext: mt = mt.filter_rows(mt.tx_annotation.mean_expression > 0.1)
        mt = mt.annotate_rows(info=hl.struct(
            AC=mt.freq[0].AC, AN=mt.freq[0].AN, AF=mt.freq[0].AF, n_hom=mt.freq[0].homozygote_count,
            CSQ=vep_ht[mt.row_key].vep
        )).drop('is_missing')
        if args.genomes:
            exome_ht = hl.read_matrix_table(homozygous_lof_mt_path_raw.format(data_type='')).rows()
            mt = mt.filter_rows(hl.is_missing(exome_ht[mt.row_key]))

        hl.export_vcf(mt, homozygous_lof_mt_path.replace('.mt', '.vcf.bgz'), metadata={
            'info': {'CSQ': {'Description': vep_ht.vep_csq_header.collect()[0]}}
        })
        ht = mt.rows()
        # ht.count()  # 3385 variants (110 in genomes)
        ht.select(
            'indel', AC=ht.info.AC, AN=ht.info.AN, AF=ht.info.AF, n_hom=ht.info.n_hom,
            variant_id=hl.delimit([ht.locus.contig, hl.str(ht.locus.position), ht.alleles[0], ht.alleles[1]], '-'),
            csq=ht.vep.worst_csq_by_gene
        ).flatten().export(homozygous_lof_mt_path.replace('.mt', '.txt.bgz'))

        ht = ht.group_by(ht.vep.worst_csq_by_gene.gene_symbol, ht.indel).aggregate(
            total=hl.agg.count(),
            total_twohit=hl.agg.count_where(ht.freq[0].homozygote_count >= 2)
        )
        ht.write(lofs_by_gene_ht_path, args.overwrite)
        hl.read_table(lofs_by_gene_ht_path).export(lofs_by_gene_ht_path.replace('.ht', '.txt.bgz'))

    if args.count_lofs:
        ht = get_gnomad_public_data(data_type)
        result = count_lofs(ht, data_type)
        ht = ht.filter(hl.len(ht.filters) == 0)
        filtered_result = count_lofs(ht, data_type)
        print('\n')
        pprint(dict(result))
        print('\n')
        pprint(dict(filtered_result))


def count_lofs(ht, data_type):
    lof_cats = hl.literal({"splice_acceptor_variant", "splice_donor_variant", "stop_gained", "frameshift_variant"})
    basic_lof_csqs = ht.vep.transcript_consequences.filter(
        lambda x: lof_cats.contains(add_most_severe_consequence_to_consequence(x).most_severe_consequence) &
                  (x.biotype == 'protein_coding')
    )
    basic_canonical_lof_csqs = ht.vep.transcript_consequences.filter(
        lambda x: lof_cats.contains(add_most_severe_consequence_to_consequence(x).most_severe_consequence) &
                  (x.biotype == 'protein_coding') & (x.canonical == 1)
    )
    hc_lof_csqs = ht.vep.transcript_consequences.filter(lambda x: (x.lof == 'HC'))
    hc_noflag_lof_csqs = ht.vep.transcript_consequences.filter(lambda x: (x.lof == 'HC') & hl.is_missing(x.lof_flags))
    hc_canonical_lof_csqs = ht.vep.transcript_consequences.filter(lambda x: (x.lof == 'HC') & (x.canonical == 1))
    hc_noflag_canonical_lof_csqs = ht.vep.transcript_consequences.filter(
        lambda x: (x.lof == 'HC') & hl.is_missing(x.lof_flags) & (x.canonical == 1))
    if data_type == 'genomes':
        sample_sizes = {'female': 6967, 'male': 8741}
    else:
        sample_sizes = {'female': 57787, 'male': 67961}
    result = ht.aggregate(
        hl.struct(
            basic_lofs=hl.agg.count_where(hl.len(basic_lof_csqs) > 0),
            basic_canonical_lofs=hl.agg.count_where(hl.len(basic_canonical_lof_csqs) > 0),

            hc_lofs=hl.agg.count_where(hl.len(hc_lof_csqs) > 0),
            hc_noflag_lofs=hl.agg.count_where(hl.len(hc_noflag_lof_csqs) > 0),
            an_adj_hc_lofs=hl.agg.count_where((hl.len(hc_lof_csqs) > 0) & get_an_adj_criteria(ht, sample_sizes)),
            an_adj_hc_noflag_lofs=hl.agg.count_where(
                (hl.len(hc_noflag_lof_csqs) > 0) & get_an_adj_criteria(ht, sample_sizes)),

            hc_canonical_lofs=hl.agg.count_where(hl.len(hc_canonical_lof_csqs) > 0),
            hc_canonical_noflag_lofs=hl.agg.count_where(hl.len(hc_noflag_canonical_lof_csqs) > 0),
            an_adj_hc_canonical_lofs=hl.agg.count_where(
                (hl.len(hc_canonical_lof_csqs) > 0) & get_an_adj_criteria(ht, sample_sizes)),
            an_adj_hc_canonical_noflag_lofs=hl.agg.count_where(
                (hl.len(hc_noflag_canonical_lof_csqs) > 0) & get_an_adj_criteria(ht, sample_sizes)))
    )
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--by_transcript', help='Use all transcripts instead of worst', action='store_true')
    parser.add_argument('--filter_an_adj', help='Filter variants with less than 80% callrate', action='store_true')
    parser.add_argument('--no_loftee', help='Do not filter with LOFTEE (only for comparison)', action='store_true')
    parser.add_argument('--remove_ultra_common', help='Remove ultra-common (AF > 95%) variants', action='store_true')
    parser.add_argument('--dont_collapse_indels', help='Do not collapse indels when exporting Gene LoF matrix', action='store_true')
    parser.add_argument('--pop_specific', help='Use pop-specific matrix instead', action='store_true')
    parser.add_argument('--include_low_pext', help='Include low pext variants in export', action='store_true')
    parser.add_argument('--calculate_gene_lof_matrix', help='Calculate Gene LoF Matrix', action='store_true')
    parser.add_argument('--export_gene_lof_matrix', help='Export Gene LoF Matrix', action='store_true')
    parser.add_argument('--create_coverage_summary', help='Create coverage summary file', action='store_true')
    parser.add_argument('--combine_lof_metrics', help='Combine all LoF metrics', action='store_true')
    parser.add_argument('--export_combined_metrics', help='Export combined LoF metrics', action='store_true')
    parser.add_argument('--compute_homozygous_lof', help='Compute Homozygous LoFs', action='store_true')
    parser.add_argument('--export_homozygous_lof', help='Export Homozygous LoF summary', action='store_true')
    parser.add_argument('--count_lofs', help='Generate counts of LoFs for paper', action='store_true')
    parser.add_argument('--genomes', help='Use genomes instead of exomes', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
