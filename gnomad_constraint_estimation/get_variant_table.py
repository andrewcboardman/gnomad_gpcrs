import argparse
import hail as hl
import numpy as np
import pandas as pd
import json
from google.cloud import bigquery

def query_region(gene_symbol, chromosome, x, y, gnomad_version='v2_1_1_exomes',client):
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

def test(client):
    test_DARC = query_region('DARC','chr1','159173097','159176290',client)
    print(test_DARC)

def iterate_genes(client):
    gpcr_gene_regions = pd.read_csv('../data/Ensembl_Grch37_gpcr_genome_locations.csv')
    gpcr_gene_regions
    variants_by_gene = []
    for i, row in gpcr_gene_regions.iterrows():
        variants = query_region(
            row['Grch37 symbol'],
            'chr'+row['Grch37 chromosome'],
            row['Grch37 start bp'],
            row['Grch37 end bp'],
            client
        )
    variants_by_gene.append(variants)
    all_variants = pd.concat(variants_by_gene)
    all_variants = all_variants[all_variants.allele_count > 0]
    all_variants = gpcr_gene_regions[['HGNC symbol','HGNC name','Grch37 symbol']].merge(all_variants,left_on='Grch37 symbol',right_on='vep_gene_symbol')
    all_variants.to_csv('../data/gnomad_v2.1.1_gpcr_variants_unfiltered.csv')
    

def main(args):
    client = bigquery.Client()
    if args.test:
        test(client)
    else:
        iterate_genes(client)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--test', help='Run tests without actually requesting data',action='store_true')
    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    args = parser.parse_args()
    main(args)