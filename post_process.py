# %%
import pandas as pd

# %%
constraint_by_gene_raw = pd.read_csv('data/gnomad_standard/constraint_final.csv.gz', compression='gzip', index_col=0)
constraint_by_gene_raw

# %%
gpcr_genes = pd.read_csv('data/Ensembl_Grch37_gpcr_genome_locations.csv')
gpcr_genes
gpcr_genes = gpcr_genes[['HGNC symbol','HGNC name','Grch37 symbol']]
gpcr_genes.columns = ['symbol','gene_name_long','gene']

# %%
gpcr_constraint = gpcr_genes.merge(constraint_by_gene_raw, on='gene',how='left')
gpcr_constraint = gpcr_constraint[gpcr_constraint.canonical].drop(columns = ['canonical'])
gpcr_constraint = gpcr_constraint[gpcr_constraint.metric.isin(('lof_hc','mis_pphen'))]
gpcr_constraint = pd.pivot(gpcr_constraint,
    index=['symbol','gene_name_long','gene','transcript'], 
    columns = 'metric', 
    values=['obs','oe'])
gpcr_constraint.columns = ['_'.join(col).strip() for col in gpcr_constraint.columns.values]
gpcr_constraint.to_csv('data/gpcr_constraint_wide.csv')


