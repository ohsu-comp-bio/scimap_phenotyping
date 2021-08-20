import argparse
import json
import warnings

import scimap as sm
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
    

    


if __name__ == '__main__':
    aparser = argparse.ArgumentParser()
    aparser.add_argument("-i", "--inputs", dest="inputs", required=True)
    aparser.add_argument("-e", "--output", dest="output", required=True)
    aparser.add_argument("-a", "--anndata", dest="anndata", required=True)

    args = aparser.parse_args()

    main(args.inputs, args.anndata, args.output)
