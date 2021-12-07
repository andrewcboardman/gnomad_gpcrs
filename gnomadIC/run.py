import pandas as pd
from .data import *
from .model import *
from .summarise import *

def setup_paths(run_ID):
    root = './data'
    # Google storage paths 
    gs_paths = dict(
        exomes_path = 'gs://gcp-public-data--gnomad/papers/2019-flagship-lof/v1.0/model/exomes_processed.ht/',
        context_path = 'gs://gcp-public-data--gnomad/resources/context/grch37_context_vep_annotated.ht/',
        mutation_rate_path = 'gs://gcp-public-data--gnomad/papers/2019-flagship-lof/v1.0/model/mutation_rate_methylation_bins.ht',
        po_coverage_path = 'gs://gcp-public-data--gnomad/papers/2019-flagship-lof/v1.0/model/prop_observed_by_coverage_no_common_pass_filtered_bins.ht'
    )
    # Local paths
    output_subdir = f'{root}/{run_ID}'
    if not os.path.isdir(output_subdir):
        os.mkdir(output_subdir)
    local_paths = dict(
        # models - shared between runs
        mutation_rate_local_path = f'{root}/models/mutation_rate_methylation_bins.ht',
        po_coverage_local_path = f'{root}/models/prop_observed_by_coverage_no_common_pass_filtered_bins.ht',
        coverage_models_local_path = f'{root}/models/coverage_models.pkl',
        # outputs - specific to run
        exomes_local_path = f'{output_subdir}/exomes.ht',
        context_local_path = f'{output_subdir}/context.ht',        
        possible_variants_ht_path = f'{output_subdir}/possible_transcript_pop.ht',
        po_output_path = f'{output_subdir}/prop_observed.ht',
        finalized_output_path = f'{output_subdir}/constraint.ht',
        summary_output_path = f'{output_subdir}/constraint_final.ht'
    )
    paths = {**gs_paths, **local_paths}

    return paths


def get_gene_intervals(test=False,controls=False):
    # Get Ensembl gene intervals from file

    gpcr_gene_intervals = pd.read_csv('data/Ensembl_Grch37_gpcr_genome_locations.csv')

    gpcr_gene_intervals = gpcr_gene_intervals[['HGNC symbol','Grch37 symbol','Grch37 chromosome','Grch37 start bp','Grch37 end bp']]
    gpcr_gene_intervals['locus_interval_txt'] = gpcr_gene_intervals['Grch37 chromosome'] + ':'  + \
        + gpcr_gene_intervals['Grch37 start bp'].map(str) + '-' + gpcr_gene_intervals['Grch37 end bp'].map(str)
    gpcr_gene_intervals['interval'] = gpcr_gene_intervals['locus_interval_txt'].map(hl.parse_locus_interval)

    if controls:     
        gene_intervals = pd.read_csv('data/ensembl_gene_annotations.txt',sep='\t')
        gene_intervals.columns = ['ensembl_gene_id','Grch37 chromosome','Grch37 start bp','Grch37 end bp','Grch37 symbol']
        gene_intervals[['ensembl_gene_id','Grch37 symbol','Grch37 chromosome','Grch37 start bp','Grch37 end bp']]
        control_gene_intervals = gene_intervals[~gene_intervals['Grch37 symbol'].isin(gpcr_gene_intervals['Grch37 symbol'])].sample(n=500,random_state=0)
        control_gene_intervals.to_csv('data/control_gene_intervals.csv')
        control_gene_intervals['locus_interval_txt'] = control_gene_intervals['Grch37 chromosome'] + ':'  + \
        + control_gene_intervals['Grch37 start bp'].map(str) + '-' + control_gene_intervals['Grch37 end bp'].map(str)
        control_gene_intervals['interval'] = control_gene_intervals['locus_interval_txt'].map(hl.parse_locus_interval)


    if test:
        gpcr_gene_intervals = gpcr_gene_intervals.sample(n=1,random_state=0)
        print(f"{str(gpcr_gene_intervals['HGNC symbol'].values[0])} chosen as test gene")
        return gpcr_gene_intervals['interval'].tolist()
    elif not controls:
        return gpcr_gene_intervals['interval'].tolist() 
    else:
        return gpcr_gene_intervals['interval'].tolist() + control_gene_intervals['interval'].tolist()


def run_tasks(tasks, paths, model, test = False, controls=False):
    '''Runs all requested tasks in specified path'''
    data = {}
    
    if 'download' in tasks:
        # Load gene intervals
        gene_intervals = get_gene_intervals(test,controls)
        # If in test mode only load 1 gene
        print('Getting data from Google Cloud...')
        data = get_data(paths, gene_intervals, model)
        print('Data loaded successfully!')
  
    if 'model' in tasks:
        print('Modelling expected number of variants')
        data = model(paths, data, model)
        print()
    
    if 'summarise' in tasks:
        print('Running aggregation by variant classes')
        data = summarise(paths, data, model)
        print('Aggregated variants successfully!')
