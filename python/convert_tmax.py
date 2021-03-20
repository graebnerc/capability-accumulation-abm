import os
if os.getcwd()[-6:] == "python":
    os.chdir("..")
if os.getcwd()[-23:] != "capability-accumulation":
    print("WARNING: Current working directory: ", os.getcwd())
    exit(1)
import numpy as np
import pandas as pd


def aggregate_tmax(base_path, col_variable, rename_dict=None,
                   relevant_cols=None, save_data=False):
    """Aggregates the result files for cases to a single file
        
    Parameters
    ----------
    base_path : str
        base path the the result files
    col_variable :
    rename_dict : list of dict, optional
        Dictionary to rename the relevant dolumn
    relevant_cols : list
        List of the columns to be kept in the aggregated file
    save_data : bool

    Example
    -------
    path_used = "output/topology/topology_"
    cases = ["BA", "full", "PWC", "random", "regular", "ring"]
    rename_dict_topology = {"regular_4": "Regular", 
                            "ring_2": "Ring",
                            "random_0.25": "Random", 
                            "full": "Complete", 
                            "BA_4": "Barabasi-Albert", 
                            "powerlaw_cluster_4_0.7": "PLC"}
    """
    tmax_data_list = []
    for case in cases:
        # file_name = base_path + case + "_tmax.feather"
        file_name = base_path + case + "_tmax.csv"
        # tmax_data_list.append(pd.read_feather(file_name))
        tmax_data_list.append(pd.read_csv(file_name))
    tmax_data = pd.concat(tmax_data_list, ignore_index=True, sort=False)
    if rename_dict is not None:
        for i in col_variable:
            tmax_data = tmax_data.replace({i: rename_dict})
    if relevant_cols is not None:
        relevant_cols += col_variable
        relevant_cols = list(set(relevant_cols))
        # if col_variable not in relevant_cols:
        # relevant_cols.append(col_variable)
        tmax_data = tmax_data[relevant_cols]
    if save_data is True:
        # tmax_data.to_feather(base_path +"tmax.feather")
        tmax_data.to_csv(base_path + "tmax.csv")
        print("Saved DataFrame in {}".format(base_path + "tmax.feather"))

    return tmax_data


cols_relevant = ["output_all_prdcts", 
                 "share_produced_products", "comp_produced_prod_mean", 
                 "price_produced_prod_mean", "price_produced_prod_std"]

# Topology
col_var = ["prod_space_struc"]
file_path = "output/topology/topology_"
cases = ["BA", "full", "PWC", "random", "regular", "ring"]

rename_dict = {"regular_4": "Regular", 
               "ring_2": "Ring",
               "random_0.25": "Random", 
               "full": "Complete", 
               "BA_4": "Barabasi-Albert", 
               "powerlaw_cluster_4_0.7": "PLC"}

topology_tmax = aggregate_tmax(
    base_path=file_path, col_variable=col_var,
    rename_dict=rename_dict, relevant_cols=cols_relevant, save_data=True)

# Complexity distribution
col_var = ["prod_sp_val_dist"]
file_path = "output/comp_dist/comp_dist_"
cases = ["exp1", "exp01", "exp05", "exp15", "normal01", "uni1", "uni05"] 

rename_dict = {"exponential_[1.0]": "Exp(1)", 
               "exponential_[0.1]": "Exp(0.1)", 
               "exponential_[0.5]": "Exp(0.5)", 
               "exponential_[1.5]": "Exp(1.5)", 
               "normal_[0, 1]": "N(0, 1)", 
               "uniform_[1.0]": "U(1)", 
               "uniform_[0.5]": "U(0.5)"}

distribution_tmax = aggregate_tmax(
    base_path=file_path, col_variable=col_var,
    rename_dict=rename_dict, relevant_cols=cols_relevant, save_data=True)

distribution_tmax_normed = distribution_tmax[["prod_sp_val_dist"]].merge(
    distribution_tmax.groupby(['prod_sp_val_dist']).transform(
        lambda x: (x - x.min()) / (x.max()-x.min())), 
    left_index=True, right_index=True)

distribution_tmax_normed = distribution_tmax_normed.groupby(
    ['prod_sp_val_dist']).agg(
        mean_price=("price_produced_prod_mean", np.mean), 
        sd_price=("price_produced_prod_mean", np.std), 
        mean_price_std=("price_produced_prod_std", np.mean), 
        sd_price_std=("price_produced_prod_std", np.std))
distribution_tmax_normed.reset_index(inplace=True, drop=False)
distribution_tmax_normed.head(2)

# Complexity allocation
col_var = ["prod_sp_val_alct"]
file_path = "output/allocation/allocation_"
cases = ["closeness", "degree", "eigenvector", "random"]
rename_dict = {'closeness_centrality': "Closeness", 
               "degree_centrality": 'Degree', 
               "eigenvector_centrality": 'Eigenvector', 
               "random": 'Random'}

allocation_tmax = aggregate_tmax(
    base_path=file_path, col_variable=col_var,
    rename_dict=rename_dict, relevant_cols=cols_relevant, save_data=True)

# Delta coefficient - reduced values
col_var = ["delta_coefficient"]
file_path = "output/delta_small/delta_small_"
cases = [str(i) for i in [1, 5, 10, 20, 50]]

delta__small_tmax = aggregate_tmax(
    base_path=file_path, col_variable=col_var,
    rename_dict=None, relevant_cols=cols_relevant, save_data=True)

# Maximum range of vision
col_var = ["max_info"]
file_path = "output/range_vis/range_vis_"
cases = ["0" + str(i) for i in range(0, 10, 2)]
cases.append("10")

range_vis = aggregate_tmax(
    base_path=file_path, col_variable=col_var,
    rename_dict=None, relevant_cols=cols_relevant, save_data=True)

# Maximum saturation level
col_var = ["nominal_demand"]
file_path = "output/saturation/saturation_"
cases = [str(i) for i in [25, 50, 75, 100, 125]]

saturation_threshold = aggregate_tmax(
    base_path=file_path, col_variable=col_var,
    rename_dict=None, relevant_cols=cols_relevant, save_data=True)

# Topology and allocation

col_var = ["nominal_demand", "prod_space_struc"]
file_path = "output/topdist/topdist_"
cases = [[str(i)+str(j) for j in ["_degree", "_random"]]
         for i in ["BA", "complete", "random"]]
cases = [item for sublist in cases for item in sublist]
topdist = aggregate_tmax(
    base_path=file_path, col_variable=col_var,
    rename_dict=None, relevant_cols=cols_relevant, save_data=True)

# Topology and vision

col_var = ["nominal_demand", "prod_space_struc"]
file_path = "output/vistop/vistop_"

cases_nb = ["_0" + str(i) for i in range(0, 10, 2)]
cases_nb.append("_10")
cases = [[str(i)+str(j) for j in cases_nb]
         for i in ["BA", "complete", "random"]]
cases = [item for sublist in cases for item in sublist]
topvis = aggregate_tmax(
    base_path=file_path, col_variable=col_var,
    rename_dict=None, relevant_cols=cols_relevant, save_data=True)
