import hail as hl
from .utils import utils

def summarise_prop_observed(po_ht, summary_path):
    """ Function for drawing final inferences from observed and expected variant counts"""

    # Finish annotation groups    
    classic_lof_annotations = hl.literal({'stop_gained', 'splice_donor_variant', 'splice_acceptor_variant'})
    variant_class = (hl.case()
        .when(classic_lof_annotations.contains(po_ht.annotation) & (po_ht.modifier == 'HC'), 'lof_hc')
        .when(classic_lof_annotations.contains(po_ht.annotation) & (po_ht.modifier == 'LC'), 'lof_lc')
        .when((po_ht.annotation == 'missense_variant') & (po_ht.modifier == 'probably_damaging'), 'mis_pphen')
        .when((po_ht.annotation == 'missense_variant') & (po_ht.modifier != 'probably_damaging'), 'mis_non_pphen')
        .when(po_ht.annotation == 'synonymous_variant', 'syn')
        .default('non-coding variants')
        )
    po_ht = po_ht.annotate(variant_class = variant_class)
    
    # Final aggregation by group 
    groups = ('gene','transcript','canonical','variant_class')    
    agg_expr = {
        'obs': hl.agg.sum(po_ht.observed_variants),
        'exp': hl.agg.sum(po_ht.expected_variants),
        'oe': hl.agg.sum(po_ht.observed_variants) / hl.agg.sum(po_ht.expected_variants),
        'adj_mu': hl.agg.sum(po_ht.adjusted_mutation_rate),
        'raw_mu': hl.agg.sum(po_ht.raw_mutation_rate),
        'poss': hl.agg.sum(po_ht.possible_variants)
    }
    constraint_ht = po_ht.group_by(*groups).aggregate(**agg_expr)
    
    # calculate confidence intervals, join tables and label
    constraint_ht = utils.oe_confidence_interval(constraint_ht, constraint_ht.obs, constraint_ht.exp, select_only_ci_metrics=False)
    constraint_df = constraint_ht.select_globals().to_pandas()
    constraint_df.to_csv(summary_path.replace('.ht','.csv.gz'),compression='gzip')
    return constraint_df

def summarise(paths, data, model):
    if not data:
        data['prop_observed_ht'] = hl.read_table(paths['po_output_path'])
    data['summary'] = summarise_prop_observed(data['prop_observed_ht'], paths['summary_output_path'])