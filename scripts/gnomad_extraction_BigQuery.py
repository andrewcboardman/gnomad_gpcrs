#!/usr/bin/env python
# coding: utf-8

# # BigQuery extraction
# This notebook contains code to download genetic variants in GPCR gene regions from gnomAD v2.1.1 hosted on Google BigQuery

# In[ ]:


import numpy as np
import pandas as pd
import json
from google.cloud import bigquery


# In[2]:


# Setup BigQuery client and convenience function
client = bigquery.Client()

def run_query(query):
    query_job = client.query(query)
    result = query_job.to_dataframe()
    gb_processed = (query_job.total_bytes_billed / 1024 ** 3)
    print(
        'This query processed {} GB of data which is {}% of your 1 TB monthly free quota.'.format(
            gb_processed, round(gb_processed / 1024 * 100, 4)
        )
    )
    return result


# In[3]:


# Query to generate table of column info for the gnomAD variant table
query_template = """
SELECT column_name, field_path, description
FROM `bigquery-public-data`.gnomAD.INFORMATION_SCHEMA.COLUMN_FIELD_PATHS
WHERE table_name = "{GNOMAD_VER}__{CHROM}"
      AND column_name IN (
          SELECT COLUMN_NAME
          FROM `bigquery-public-data`.gnomAD.INFORMATION_SCHEMA.COLUMNS
          WHERE table_name = "{GNOMAD_VER}__{CHROM}")
"""
query = query_template.format(
    GNOMAD_VER='v2_1_1_exomes',
    CHROM='chr21'
)

column_info = run_query(query)
print(
    'There are {} columns in `bigquery-public-data.gnomAD.{}__{}` table'.format(
        len(column_info.index),
        'v2_1_1_exomes',
        'chr21'
    )
)
column_info.to_csv('../data/gnomad_v2_1_1_exomes_column_info.csv')


# In[4]:


# Load table of GPCR gene regions, extracted from Ensembl
gpcr_gene_regions = pd.read_csv('../data/Ensembl_Grch37_gpcr_genome_locations.csv')
gpcr_gene_regions


# In[60]:


# Query to extract genetic variants in region
# NOTE: For v2_1_1 the "variant_type" column must be replaced with "alternate_bases.allele_type AS variant_type"
def query_region(gene_symbol, chromosome, x, y, gnomad_version='v2_1_1_exomes'):
    # Adjust query based on version
    if gnomad_version.startswith('v3'):
      # Variant type (snv, indel, multi-snv, multi-indel, or mixed) is stored under difference columns in V2 and V3
        variant_type_col = 'variant_type'
        extra_columns = ''
    else:
        variant_type_col = 'alternate_bases. allele_type'
        # These vep columns only exist in V2
        extra_columns = 'vep.STRAND AS STRAND, vep.Protein_position AS Protein_pos,'
    # Set query
    query_template = """
    SELECT reference_name AS chromosome, 
           start_position AS genome_pos,
           names AS variant_id,
           reference_bases,
           alternate_bases.alt AS alternate_bases,
           AN AS allele_number,
           alternate_bases.AC AS allele_count,
           alternate_bases.nhomalt AS num_alternate_homozygous,
           vep.Consequence AS vep_consequence,
           vep.IMPACT AS vep_impact,
           vep.SYMBOL AS vep_gene_symbol,
           vep.Gene AS vep_ensembl_gene,
           vep.Feature AS vep_ensembl_transcript,
           vep.ENSP AS vep_ensembl_protein,
           vep.Protein_position as vep_protein_pos,
           vep.Amino_acids as vep_amino_acids,
           vep.DISTANCE as vep_distance_to_transcript,
           vep.SWISSPROT as vep_swissprot_match,
           vep.SIFT as vep_SIFT,
           vep.PolyPhen as vep_PolyPhen
    FROM `bigquery-public-data.gnomAD.{GNOMAD_VER}__{CHROM}` AS main_table,
         main_table.alternate_bases AS alternate_bases,
         alternate_bases.vep AS vep
    WHERE start_position >= {X} AND start_position <= {Y} AND vep.SYMBOL = "{GENE}"
    ORDER BY 1,2
    """
    query = query_template.format(GNOMAD_VER=gnomad_version,
                                  CHROM=chromosome,
                                  GENE=gene_symbol,
                                  X=x, 
                                  Y=y
                                 )

    print(
        'Running Region (Type 1) queries on gnomAD version: {}, chromosome: {}, gene symbol: {}'.format(
            gnomad_version,
            chromosome,
            gene_symbol
        )
    )
    output = run_query(query)
    return output


# In[61]:


# Test region query 
test_DARC = query_region('DARC','chr1','159173097','159176290')


# In[62]:


# Check output of test query
test_DARC


# In[63]:


# Run query for all GPCRs in table
variants_by_gene = []
for i, row in gpcr_gene_regions.iterrows():
    variants = query_region(
        row['Grch37 symbol'],
        'chr'+row['Grch37 chromosome'],
        row['Grch37 start bp'],
        row['Grch37 end bp']
    )
    variants_by_gene.append(variants)


# In[78]:


# Filter non-existent variants and annotate with up-to-date gene symbols
all_variants = pd.concat(variants_by_gene)
all_variants = all_variants[all_variants.allele_count > 0]
all_variants = gpcr_gene_regions[['HGNC symbol','HGNC name','Grch37 symbol']].merge(all_variants,left_on='Grch37 symbol',right_on='vep_gene_symbol')
all_variants


# In[80]:


all_variants.to_csv('../data/gnomad_v2.1.1_gpcr_variants_unfiltered.csv')

