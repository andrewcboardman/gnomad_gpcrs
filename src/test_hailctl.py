import hail as hl
import pandas as pd

hl.init()
ht_exomes_path = 'gs://gcp-public-data--gnomad/release/2.1/ht/exomes/gnomad.exomes.r2.1.sites.ht/'
ht_exomes = hl.read_table(ht_exomes_path)
ht_exomes.describe()
