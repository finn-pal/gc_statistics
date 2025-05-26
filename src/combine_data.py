import argparse

import h5py
import numpy as np

from tools.get_clustering_trends import iteration_clustering
from tools.get_galaxy_trends import get_galaxy_data
from tools.get_group_trends import get_group_data

# def save_dict_to_h5py(h5file, path, d):
#     for key, value in d.items():
#         key_path = f"{path}/{key}"
#         if isinstance(value, dict):
#             save_dict_to_h5py(h5file, key_path, value)  # recursive
#         else:
#             # Convert lists to numpy arrays
#             if isinstance(value, list):
#                 value = np.array(value)
#             # Handle strings separately
#             if isinstance(value, str):
#                 h5file.create_dataset(key_path, data=np.string_(value))
#             elif np.isscalar(value):
#                 h5file.create_dataset(key_path, data=value)
#             else:
#                 h5file.create_dataset(key_path, data=value)


def save_dict_to_h5py(h5file, key_path, data):
    for key, value in data.items():
        new_key_path = f"{key_path}/{key}" if key_path else key
        if isinstance(value, dict):
            save_dict_to_h5py(h5file, new_key_path, value)  # recursive
        else:
            try:
                value = np.array(value)
                h5file.create_dataset(new_key_path, data=value)
            except ValueError:
                # Save each sub-element of the list of lists as a separate dataset
                grp = h5file.create_group(new_key_path)
                for i, subvalue in enumerate(value):
                    grp.create_dataset(str(i), data=np.array(subvalue))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--simulation", required=True, type=str, help="simulation name (e.g. m12i)")
    parser.add_argument("-a", "--iteration_low_limit", required=True, type=int, help="lower bound iteration")
    parser.add_argument("-b", "--iteration_up_limit", required=True, type=int, help="upper bound iteration")
    parser.add_argument("-l", "--snap_lim", required=False, default=46, type=int, help="lower snapshot limit")
    parser.add_argument("-c", "--cluster_lim", required=False, default=100, type=int, help="cluster limit")

    args = parser.parse_args()

    it_min = args.iteration_low_limit
    it_max = args.iteration_up_limit
    it_lst = np.arange(it_min, it_max + 1)

    sim = args.simulation

    snap_lim = args.snap_lim
    cluster_lim = args.cluster_lim

    # file locations
    sim_dir = "/Users/z5114326/Documents/simulations/"

    cluster_types = ["birch", "agglom", "kmeans"]
    cluster_variables = ["et", "jr", "jp", "jz"]
    n_clusters = 2

    grp_dict = get_group_data(
        it_lst, sim, sim_dir, snap_lim, cluster_types, cluster_lim, n_clusters, cluster_variables
    )

    cluster_dict = {}
    for cluster_type in cluster_types:
        cluster_dict[cluster_type] = iteration_clustering(
            it_lst, cluster_type, n_clusters, cluster_variables, sim, sim_dir, cluster_lim
        )

    galaxy_info = get_galaxy_data(sim, sim_dir, snap_lim)

    comb_dict = {"clustering": cluster_dict, "groups": grp_dict, "galaxy": galaxy_info}

    save_file = sim_dir + sim + "/" + sim + "_combined.hdf5"  # save location
    with h5py.File(save_file, "w") as f:  # 'w' mode creates or overwrites the file
        save_dict_to_h5py(f, "", comb_dict)
