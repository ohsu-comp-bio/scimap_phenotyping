import argparse
import ast
import json
import warnings

import pandas as pd
import squidpy as sq
from anndata import read_h5ad


def main(inputs, anndata, output):
    """
    Parameter
    ---------
    inputs : str
        File path to galaxy tool parameter.
    anndata : str
        File path to anndata containing phenotyping info.
    output : str
        File path to output.
    """
    warnings.simplefilter('ignore')

    with open(inputs, 'r') as param_handler:
        params = json.load(param_handler)

    adata = read_h5ad(anndata)

    tool = params['analyses']['selected_tool']
    tool_func = getattr(sq.gr, tool)

    options = params['analyses']['options']

    for k, v in options.items():
        if k == 'genes':    # for spatial_autocorr and sepal
            if v:
                options[k] = [e.strip() for e in v.split(',')]
        elif k == 'radius':    # for spatial_neighbors
            if v:
                options[k] = ast.literal_eval(v)
        elif k == 'numba_parallel':    # for nhood_enrichment and ligrec
            if v == 'false':
                options[k] = False
            elif v == 'true':
                options[k] = True
        elif k == 'interactions':    # for ligrec
            options[k] = pd.read_csv(v, sep="\t")
        elif k == 'max_neighs':
            options[k] = int(v)      # for sepal

        if v in ['', 'none']:
            options[k] = None

    cluster_key = params['analyses'].get('cluster_key')
    if cluster_key:
        tool_func(adata, cluster_key, **options)
    else:
        tool_func(adata, **options)

    adata.write(output)


if __name__ == '__main__':
    aparser = argparse.ArgumentParser()
    aparser.add_argument("-i", "--inputs", dest="inputs", required=True)
    aparser.add_argument("-e", "--output", dest="output", required=True)
    aparser.add_argument("-a", "--anndata", dest="anndata", required=True)

    args = aparser.parse_args()

    main(args.inputs, args.anndata, args.output)
