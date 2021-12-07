import hail as hl
import pandas as pd

hl.init()

# Select variants in canonical transcript in exomes hail table
exomes_ht = hl.read_table('data/gnomad_standard/exomes.ht')
exomes_ht = exomes_ht.select(*('gene','transcript','canonical','hgvsp','annotation','modifier','freq'))
exomes_ht = exomes_ht.filter(exomes_ht.canonical)
exomes_ht = exomes_ht.transmute(freq_gnomad = exomes_ht.freq[0])
exomes_ht = exomes_ht.transmute(AC = exomes_ht.freq_gnomad.AC, homozygote_count=exomes_ht.freq_gnomad.homozygote_count)
exomes_ht = exomes_ht.drop('canonical','transcript')

# classify variants with standard pipeline
classic_lof_annotations = hl.literal({'stop_gained','frameshift_variant','splice_donor_variant', 'splice_acceptor_variant'})
exomes_ht = exomes_ht.annotate(variant_class = (hl.case()
        .when(classic_lof_annotations.contains(exomes_ht.annotation) & (exomes_ht.modifier == 'HC'), 'lof_hc')
        .when(classic_lof_annotations.contains(exomes_ht.annotation) & (exomes_ht.modifier == 'LC'), 'lof_lc')
        .when((exomes_ht.annotation == 'missense_variant') & (exomes_ht.modifier == 'probably_damaging'), 'mis_pphen')
        .when((exomes_ht.annotation == 'missense_variant') & (exomes_ht.modifier != 'probably_damaging'), 'mis_non_pphen')
        .when(exomes_ht.annotation == 'synonymous_variant', 'syn')
        .default('non-coding')
        ))


gpcr_genes = pd.read_csv('data/Ensembl_Grch37_gpcr_genome_locations.csv',index_col=0)
gpcr_genes = gpcr_genes[['Grch37 symbol','HGNC name']]
gpcr_genes.columns = ['gene','gene_name_long']

lof_variants = exomes_ht.filter(classic_lof_annotations.contains(exomes_ht.annotation)).to_pandas()
print(lof_variants.columns)

lof_variants.columns = ['chromosome','location','alleles','gene','variant_HGVSP_id','variant_type','confidence','n_observed_total','n_observed_homozygous','variant_class']
#lof_variants = gpcr_genes.merge(lof_variants, on='gene',how='left')
lof_variants.to_csv('data/gnomad_standard/lof_variants.csv')

mis_variants = exomes_ht.filter(exomes_ht.annotation == 'missense_variant').to_pandas()
mis_variants.columns = ['chromosome','location','alleles','gene','variant_HGVSP_id','variant_type','polyphen','n_observed_total','n_observed_homozygous','variant_class']
#mis_variants = gpcr_genes.merge(mis_variants, on='gene',how='left')
mis_variants.to_csv('data/gnomad_standard/mis_variants.csv')

