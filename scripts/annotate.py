import requests
import io
from itertools import product
import pandas as pd


alphabet = ['Ala', 'Cys', 'Asp', 'Glu', 'Phe', 'Gly', 'His', 'Ile', 'Lys', 'Leu', 'Met', 'Asn', 'Pro', 'Gln', 'Arg', 'Ser', 'Thr', 'Sec', 'Val', 'Trp', 'Tyr']
gpcr_genes = pd.read_csv('../data/Ensembl_Grch37_gpcr_genome_locations.csv')['Ensembl id Grch37'].to_list()
server = 'http://grch37.rest.ensembl.org'
params={'type':'protein','multiple_sequences':'true'}
headers={ "Content-Type" : "application/json"} 

def get_variants(seq, id, alphabet=alphabet):
    wt_aa = [seq3(aa) for aa in seq]
    vars = []
    for (i, wt), mut in product(enumerate(wt_aa), alphabet):
        vars.append({'hgvsp':id + ':p.'+ wt + str(i + 1) + mut, 'pos':i + 1, 'wt':wt,'mut':mut}) 
    return pd.DataFrame(vars)

variants_all = []

for ensembl_gene in genes:
    ext = f'/sequence/id/{ensembl_gene}'
    
    response = requests.get(server + ext, params=params, headers=headers)
    if response.status_code != 200:
        print(f'Bad request for {ensembl_gene}')
        continue
    else:
        proteins = response.json()
        variants = []
        for protein in proteins:
            variant_df = get_variants(protein['seq'], protein['id'])
            variants.append(variant_df)
        variants = pd.concat(variants)
    
    variants_all.append(variants)

if variants_all:
    variants_all = pd.concat(variants_all)
    variants_all.to_csv('../data/variants.csv',compression='gzip')