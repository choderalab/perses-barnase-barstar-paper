import os
import copy
import pickle
import numpy as np

from simtk.openmm import unit
from tqdm import tqdm_notebook

from collections import OrderedDict
from scipy.stats import pearsonr

import argparse

# Load arguments
parser = argparse.ArgumentParser(description='compute correlations for heatmap')
parser.add_argument('main_dir', type=int, help='main dir')
parser.add_argument('sub_dir', type=int, help='sub dir')
parser.add_argument('phase', type=str, help='phase')
args = parser.parse_args()

def get_correlation_with_du_dlambda_all_replicas(n_replicas, du_dlambda, out_dir, angle=None, res=None, residue_pair=None, is_water=False, position_type='old', is_checkpoint_interval_10=False):
    """
    Get correlation of a dof with du_dlambda across all replicas

    Parameters
    ----------
    n_replicas : int
        number of replicas used in the simulation
    du_dlambda : [n_replicas, n_iterations] np.ndarray
        du/dlambda time series over n_iterations for each replica
    out_dir : str
        path to directory containing data to analyze
    angle : str, default None
        dihedral angle name (i.e., 'phi', 'psi', 'chi1', 'chi2', 'chi3', 'chi4')
    res : str, default None
        chain name and residue id (format: {barnase|barstar}-{residue id})
        example: 'barstar-42'
    residue_pair : str, default None
        residue pair name (format: {barnase|barstar}-{residue id} | {barnase|barstar}-{residue id})
        example: 'barstar-42 | barstar-76 new'
    is_water : boolean, default False
        whether the dof to analyze is neighboring waters
    position_type : str, default 'old'
        the type of positions (i.e., 'old' or 'new')
    checkpoint_interval_10 : boolean, default False
        whether the simulation was run with checkpoint interval of 10 

    Returns
    -------
    pcc : float
        pearson correlation coefficient
    transform_type : str or None
        the transformation applied to the degree of freedom (i.e. 'sin', 'cos' or None)

    """
    
    if angle: # Get angle values
        y = []
        for replica in range(n_replicas):
            with open(os.path.join(out_dir, f"dataframe_dihedral_angles_per_replica{replica}.pickle"), "rb") as f:
                df = pickle.load(f)
            data = df[df["residue"] == res][f"{angle}_{position_type}"].values
            if is_checkpoint_interval_10: # Truncate, if necessary
                data = data[::10]
            y.append(data[:501])
        y = np.concatenate(y)
            
    elif residue_pair: # Get residue distances
        y = []
        for replica in range(n_replicas):
            with open(os.path.join(out_dir, f"dataframe_residue_distances_replica{replica}.pickle"), "rb") as f:
                df = pickle.load(f)
            data = df[residue_pair].values
            if is_checkpoint_interval_10: # Truncate, if necessary
                data = data[::10]
            y.append(data[:501])   
        y = np.concatenate(y)
    
    elif is_water: # Get water values
        with open(os.path.join(out_dir, f"dataframe_waters_per_replica.pickle"), "rb") as f:
            df = pickle.load(f)
        y = []
        for replica in range(n_replicas):
            data = df[f"waters nearby mutating res (radius = 0.5 nm) {position_type}"][df["replica"] == replica].values
            if is_checkpoint_interval_10: # Truncate, if necessary
                data = data[::10]
            y.append(data[:501])
        y = np.concatenate(y)

    if not np.isnan(y).all():
        if angle:
            y_sin = np.sin(y)
            result_sin = pearsonr(du_dlambda, y_sin).statistic
            y_cos = np.cos(y)
            result_cos = pearsonr(du_dlambda, y_cos).statistic
            if abs(result_sin) > abs(result_cos):
                return result_sin, 'sin'
            else:
                return result_cos, 'cos'
        else:
            return pearsonr(du_dlambda, y).statistic, None
    else:
        return None, None

def get_correlations_all_replicas(out_dir, sub_dir, phase, n_replicas, is_checkpoint_interval_10=False):
    """
    Get correlation (over all replicas) for all degrees of freedom

    Parameters
    ----------
    out_dir : str
        path to directory containing data to analyze
    sub_dir : str
        path to sub-directory containing data to analyze
    phase : str
        phase (e.g., 'complex')
    n_replicas : int
        number of replicas used in the simulation
    checkpoint_interval_10 : boolean, default False
        whether the simulation was run with checkpoint interval of 10 

    Returns
    -------
    correlation data : OrderedDict()
        key: name of degree of freedom with the maximum magnitude pearson correlation coefficient in its category
        value: dict containing 'pearsonr' and 'transform_type' as keys

    """

    # Load du_dlambda
    du_dlambda_all = []
    for replica in range(n_replicas):
        with open(os.path.join(out_dir, f"{sub_dir}_{phase}_du_dlambda_replica{replica}.npy"), "rb") as f:
            du_dlambda = np.load(f)
        if is_checkpoint_interval_10:
            du_dlambda = du_dlambda[::10]
        du_dlambda_all.append(du_dlambda[:501]) 
    du_dlambda_all = np.concatenate(du_dlambda_all)
    du_dlambda_all = [val * unit.kilojoules_per_mole.conversion_factor_to(unit.kilocalories_per_mole) for val in du_dlambda_all]

    ## backbone and sidechain dihedral angles ##
    print("Getting dihedral angles")
    
    # Load a dataframe with angles to get the keys
    replica = 0
    with open(os.path.join(out_dir, f"dataframe_dihedral_angles_per_replica{replica}.pickle"), "rb") as f:
        df = pickle.load(f)

    # Get correlations
    d_backbone = {}
    d_sidechain = {}
    for angle in tqdm_notebook(['phi', 'psi', 'chi1', 'chi2', 'chi3', 'chi4']):
        for res in df.residue.unique():
            for position_type in ['old', 'new']:
                
                pearson_r_result, transform_type = get_correlation_with_du_dlambda_all_replicas(n_replicas, 
                                                du_dlambda_all, 
                                                out_dir, 
                                                angle=angle, 
                                                res=res, 
                                                position_type=position_type,
                                                is_checkpoint_interval_10=is_checkpoint_interval_10)

                print(f"{res} {angle}: {pearson_r_result}")
                
                if pearson_r_result is not None:
                    if angle in ['phi', 'psi']:
                        d_backbone[f"{res} {angle} {position_type}"] = {"pearsonr": pearson_r_result, "transform_type": transform_type}
                    else:
                        d_sidechain[f"{res} {angle} {position_type}"] = {"pearsonr": pearson_r_result, "transform_type": transform_type}
 
    ## intra- and inter- interface residue contacts ##
    print("Getting residue disances")
    
    # Load a dataframe with distances to get the keys
    replica = 0
    with open(os.path.join(out_dir, f"dataframe_residue_distances_replica{replica}.pickle"), "rb") as f:
        df = pickle.load(f)
    
    # Get correlations
    d_intra = {}
    d_inter = {}
    for residue_pair in tqdm_notebook(df.keys()[2:]): # skip 'replica' and 'iteration' keys

        pearson_r_result, transform_type = get_correlation_with_du_dlambda_all_replicas(n_replicas, 
                                            du_dlambda_all, 
                                            out_dir, 
                                            residue_pair=residue_pair, 
                                            is_checkpoint_interval_10=is_checkpoint_interval_10)

        print(f"{residue_pair}: {pearson_r_result}")

        res_A, res_B = residue_pair.split(" | ")
        if ("barnase" in res_A and "barnase" in res_B) or ("barstar" in res_A and "barstar" in res_B): 
            d_intra[residue_pair] = {"pearsonr": pearson_r_result, "transform_type": transform_type}
        else:
            d_inter[residue_pair] = {"pearsonr": pearson_r_result, "transform_type": transform_type}

    ## Num waters within 5 angstroms of mutating residue ##
    print("Getting neighboring waters")

    # Get correlations
    d_water = {}
    for position_type in tqdm_notebook(['old', 'new']):

        pearson_r_result, transform_type = get_correlation_with_du_dlambda_all_replicas(n_replicas, 
                                            du_dlambda_all, 
                                            out_dir, 
                                            is_water=True, 
                                            position_type=position_type,
                                            is_checkpoint_interval_10=is_checkpoint_interval_10)

        print(f"{position_type}: {pearson_r_result}")

        d_water[position_type] = {"pearsonr": pearson_r_result, "transform_type": transform_type}

    # Get max for each category
    def get_max_pearsonr(d):
        """
        Get the degree of freedom (and its pcc and transformation type) with the maximum
        magnitude PCC for its category

        Parameters
        ----------
        d : dict
            correlation values (and transformation types) for each degree of freedom in a given category
            key : degree of freedom name
            value : dict containing 'pearsonr' and 'transform_type' as keys

        Returns
        -------
        max_pearsonr_value
            maximum magnitude PCC for its category
        max_pearsonr_key
            name of the degree of freedom with the maximum magnitude PCC for its category
        transform_type : str or None
            transformation type applied to the dihedral angle time series

        """
        d_pearsonr = {k: abs(v['pearsonr']) for k, v in d.items()}
        max_pearsonr_key = max(d_pearsonr, key=d_pearsonr.get) # the name of the degree of freedom with the max pearsonr value
        max_pearsonr_value = d[max_pearsonr_key] # the pearsonr value (without taking absolute value) 
        transform_type = d[max_pearsonr_key]['transform_type'] # the transformation type ('sin', 'cos', or None)
        return max_pearsonr_value, max_pearsonr_key, transform_type

    correlation_data = OrderedDict() 
    for d in [d_backbone, d_sidechain, d_intra, d_inter, d_water]:
        max_pearsonr_value, max_pearsonr_key, transform_type = get_max_pearsonr(d)
        print(f"{max_pearsonr_key} ({transform_type})")
        correlation_data[max_pearsonr_key] = {'pearsonr': max_pearsonr_value, 'transform_type': transform_type}

    return correlation_data

main_dir = args.main_dir
sub_dir = args.sub_dir
phase = args.phase

if main_dir == 45:
    if sub_dir == 11:
        replicate = 0
        n_replicas = 24
        is_checkpoint_interval_10 = False

    elif sub_dir == 10:
        replicate = 1
        n_replicas = 36
        is_checkpoint_interval_10 = True

    elif sub_dir == 9:
        replicate = 1
        n_replicas = 24
        is_checkpoint_interval_10 = True

else:
    replicate = 1 if main_dir == 47 else 0
    is_checkpoint_interval_10 = True if sub_dir in [12, 17] else False
    n_replicas = 24 if sub_dir in [0, 1, 2, 3, 4, 14, 15, 16, 17, 18] else 36

out_dir = f"/data/chodera/zhangi/perses_benchmark/repex/perses-bnbs-paper-fourth-attempt/{main_dir}/{sub_dir}/replicate_{replicate}/"

correlation_data = get_correlations_all_replicas(out_dir, 
                                                 sub_dir, 
                                                 phase, 
                                                 n_replicas, 
                                                 is_checkpoint_interval_10=is_checkpoint_interval_10)

with open(os.path.join(out_dir, f"correlation_data_replicas_all_50ns.pickle"), "wb") as f:
    pickle.dump(correlation_data, f)


