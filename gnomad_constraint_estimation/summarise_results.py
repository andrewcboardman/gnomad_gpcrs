import argparse
import pickle
import random
import hail as hl
import pandas as pd
from itertools import product
from typing import Dict, List, Optional, Set, Tuple, Any
from gnomad_pipeline.utils import *
from constraint_analysis import *


def load_data_to_summarise(paths):
    final_input_path = paths['final_output_path']

    data = {'finalised_ht':hl.read_table(final_input_path)}
    return data


def summarise(paths, data):
    mut_types = ('lof', 'mis', 'syn','mis_pphen','mis_non_pphen')
    output_var_types = zip(('obs', 'exp', 'oe', 'oe', 'oe'),
                            ('', '', '', '_lower', '_upper'))
    output_vars = product(mut_types,output_var_types)
    data['summary'] = (data['finalised_ht']
        .select(
            'gene','transcript','canonical',
            *[f'{t}_{m}{ci}' for m, (t, ci) in output_vars],
            gene_issues=data['finalised_ht'].constraint_flag
        )
        .select_globals()
    )
    data['summary'].write(
        paths['summary_output_path'], 
        overwrite=args.overwrite
    )
    data['summary'].export(paths['summary_output_path'].replace('.ht', '.txt.bgz'))

    

