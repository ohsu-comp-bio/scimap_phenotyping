import argparse
import json
import warnings

import numpy as np
import pandas as pd
from pathlib import Path
from anndata import read_h5ad
from vitessce import (
    VitessceConfig,
    Component as cm,
    AnnDataWrapper,
    OmeTiffWrapper,
    MultiImageWrapper,
)


def main(inputs, output, image, anndata, masks=None):
    """
    Parameter
    ---------
    inputs : str
        File path to galaxy tool parameter.
    output : str
        Output folder for saving web content.
    image : str
        File path to the OME Tiff image.
    anndata : str
        File path to anndata containing phenotyping info.
    masks : str
        File path to the image masks.
    """
    warnings.simplefilter('ignore')

    with open(inputs, 'r') as param_handler:
        params = json.load(param_handler)

    marker = params['marker'].strip()
    from_gate = params['from_gate']
    to_gate = params['to_gate']
    increment = params['increment']
    x_coordinate = params['x_coordinate'].strip() or 'X_centroid'
    y_coordinate = params['y_coordinate'].strip() or 'Y_centroid'

    adata = read_h5ad(anndata)

    # If no raw data is available make a copy
    if adata.raw is None:
        adata.raw = adata

    # Copy of the raw data if it exisits
    if adata.raw is not None:
        adata.X = adata.raw.X

    data = pd.DataFrame(
        np.log1p(adata.X),
        columns = adata.var.index,
        index= adata.obs.index)
    marker_values = data[[marker]].values
    
    # Generate a dataframe with various gates
    def gate (g, d):
        dd = d.copy()
        dd = np.where(dd < g, 0, dd)
        np.warnings.filterwarnings('ignore')
        dd = np.where(dd >= g, 1, dd)
        return dd

    # Identify the list of increments
    gate_names = []
    for num in np.arange(from_gate, to_gate, increment):
        num = round(num, 3)
        key = marker + '--' +str(num)
        adata.obs[key] = gate(num, marker_values)
        gate_names.append(key)
    
    adata.obsm['XY_coordinate'] = adata.obs[[x_coordinate, y_coordinate]].values

    vc = VitessceConfig(name=None,description=None)
    dataset = vc.add_dataset()
    image_wrappers=[OmeTiffWrapper(img_path=image, name='OMETIFF')]
    if masks:
        image_wrappers.append(
            OmeTiffWrapper(img_path=masks, name='MASKS', is_bitmask=True)
        )
    dataset.add_object(MultiImageWrapper(image_wrappers))

    dataset.add_object(
        AnnDataWrapper(
            adata,
            spatial_centroid_obsm='XY_coordinate',
            cell_set_obs=gate_names,
            cell_set_obs_names=[obj[0].upper() + obj[1:] for obj in gate_names],
            expression_matrix="X"
        )
    )
    spatial = vc.add_view(dataset, cm.SPATIAL)
    cellsets = vc.add_view(dataset, cm.CELL_SETS)
    status = vc.add_view(dataset, cm.STATUS)
    lc = vc.add_view(dataset, cm.LAYER_CONTROLLER)

    vc.layout((status / cellsets / lc ) | (spatial) )
    config_dict = vc.export(to='files', base_url='http://localhost', out_dir=output)

    with open(Path(output).joinpath('config.json'), 'w') as f:
        json.dump(config_dict, f, indent=4)


if __name__ == '__main__':
    aparser = argparse.ArgumentParser()
    aparser.add_argument("-i", "--inputs", dest="inputs", required=True)
    aparser.add_argument("-e", "--output", dest="output", required=True)
    aparser.add_argument("-g", "--image", dest="image", required=True)
    aparser.add_argument("-a", "--anndata", dest="anndata", required=True)
    aparser.add_argument("-m", "--masks", dest="masks", required=False)

    args = aparser.parse_args()

    main(args.inputs, args.output, args.image, args.anndata, args.masks)
