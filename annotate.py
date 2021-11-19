import requests
import io
from Bio import SeqIO
from Bio.SeqUtils import seq3
from itertools import product
import pandas as pd

def get_sequence(fasta):
    sequence = ''.join(fasta.split('\n')[1:-1])
    return [seq3(aa) for aa in sequence]

genes = pd.read_csv('data/Ensembl_Grch37_gpcr_genome_locations.csv')['Ensembl id Grch37'].to_list()
alphabet = ['Ala', 'Cys', 'Asp', 'Glu', 'Phe', 'Gly', 'His', 'Ile', 'Lys', 'Leu', 'Met', 'Asn', 'Pro', 'Gln', 'Arg', 'Ser', 'Thr', 'Sec', 'Val', 'Trp', 'Tyr']


for ensembl_gene in genes:
    response = requests.get(f'http://grch37.rest.ensembl.org/sequence/id/{ensembl_gene}?;type=protein;multiple_sequences=true').text
    with io.StringIO(response) as input:
        protein_sequence_records = SeqIO.parse(input)
    for sequence_record in protein_sequence_records:
        protein_id = sequence_record.id
        
        output = []
        for (i, wt), mut in product(enumerate(get_sequence(response)), alphabet):
            output.append({'hgvsp':protein_id + ':p.'+ wt + str(i + 1) + mut, 'pos':i + 1, 'wt':wt,'mut':mut}) 
            variants.append(pd.DataFrame(output))

    variants = pd.concat(variants)

variants.to_csv('../data/variants.csv')