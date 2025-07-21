import json

import gc_utils
import h5py
import numpy as np
import pandas as pd
import utilities as ut

from tools.get_clustering_trends import group_iteration_clustering


def get_kl(group_id, snap_lst, sim, sim_dir):
    grp = str(group_id)

    data_file = sim_dir + sim + "/kl_data.json"
    with open(data_file, "r") as file:
        kl_dict = json.load(file)

    kl_avg = []
    kl_std = []
    cont_avg = []
    cont_std = []
    for snap in snap_lst:
        snap_msk = np.array(kl_dict[grp]["snaps"]) == snap
        kl_avg_val = np.array(kl_dict[grp]["kl_snap_avg"])[snap_msk]
        kl_std_val = np.array(kl_dict[grp]["kl_snap_std"])[snap_msk]

        if kl_avg_val != "NA":
            kl_avg.append(float(kl_avg_val[0]))
        else:
            kl_avg.append(np.nan)

        if kl_std_val != "NA":
            kl_std.append(float(kl_std_val[0]))
        else:
            kl_std.append(np.nan)

        cont_avg_val = np.array(kl_dict[grp]["cont_snap_avg"])[snap_msk]
        cont_std_val = np.array(kl_dict[grp]["cont_snap_std"])[snap_msk]

        if cont_avg_val != "NA":
            cont_avg.append(float(cont_avg_val[0]))
        else:
            cont_avg.append(np.nan)

        if cont_std_val != "NA":
            cont_std.append(float(cont_std_val[0]))
        else:
            cont_std.append(np.nan)

    return np.array(kl_avg), np.array(kl_std), np.array(cont_avg), np.array(cont_std)


def get_gc_num(group_id, it_lst, snap_lst, proc_data):
    gc_num_avg = []
    gc_num_std = []

    for snap in snap_lst:
        snap_id = gc_utils.snapshot_name(snap)
        gc_num_lst = []

        for it in it_lst:
            it_id = gc_utils.iteration_name(it)
            snp_dat = proc_data[it_id]["snapshots"][snap_id]

            grp_msk = np.abs(snp_dat["group_id"][()]) == group_id
            acc_msk = snp_dat["now_accreted"][()] == 1
            gc_num_lst.append(len(snp_dat["group_id"][grp_msk]))

        gc_num_avg.append(np.nanmean(gc_num_lst))
        gc_num_std.append(np.nanstd(gc_num_lst))

    return np.array(gc_num_avg), np.array(gc_num_std)


def get_dispersion(group_id, it_lst, snap_lst, proc_data):
    gc_dis_rho = []
    gc_dis_the = []
    gc_dis_phi = []
    gc_dis = []

    for snap in snap_lst:
        snap_id = gc_utils.snapshot_name(snap)

        sig_rho = []
        sig_the = []
        sig_phi = []
        sig = []
        for it in it_lst:
            it_id = gc_utils.iteration_name(it)
            snp_dat = proc_data[it_id]["snapshots"][snap_id]

            grp_msk = np.abs(snp_dat["group_id"][()]) == group_id

            vel_sph = snp_dat["vel.sph"][grp_msk]

            sig_rho.append(np.nanstd(vel_sph[:, 0]))
            sig_the.append(np.nanstd(vel_sph[:, 1]))
            sig_phi.append(np.nanstd(vel_sph[:, 2]))

            sig.append(np.linalg.norm([sig_rho[-1], sig_the[-1], sig_phi[-1]]))

        gc_dis_rho.append(np.nanmean(sig_rho))
        gc_dis_the.append(np.nanmean(sig_the))
        gc_dis_phi.append(np.nanmean(sig_phi))
        gc_dis.append(np.nanmean(sig))

    gc_dis_sph = np.column_stack((gc_dis_rho, gc_dis_the, gc_dis_phi))

    return np.array(gc_dis), gc_dis_sph


def get_com(group_id, it_lst, snap_lst, proc_data):
    com_avg = np.empty((0, 3))
    com_err = np.empty((0, 3))

    com_dist_avg = []
    com_dist_err = []

    for snap in snap_lst:
        snap_id = gc_utils.snapshot_name(snap)

        com_std_hold = []
        com_hold = np.empty((0, 3))
        for it in it_lst:
            it_id = gc_utils.iteration_name(it)
            snp_dat = proc_data[it_id]["snapshots"][snap_id]
            grp_msk = np.abs(snp_dat["group_id"][()]) == int(group_id)

            xyz = snp_dat["pos.xyz"][grp_msk]
            mass = 10 ** snp_dat["mass"][grp_msk]

            c_x = np.sum(xyz[:, 0] * mass) / np.sum(mass)
            c_y = np.sum(xyz[:, 1] * mass) / np.sum(mass)
            c_z = np.sum(xyz[:, 2] * mass) / np.sum(mass)

            com = np.array([c_x, c_y, c_z])

            d_c = []
            for pos in xyz:
                d_c.append(ut.coordinate.get_distances(com, pos, total_distance=True))
            d_c = np.array(d_c)

            com_std_hold.append(np.nanstd(d_c))
            com_hold = np.vstack((com_hold, com))

        com_dist_avg.append(np.nanmean(com_std_hold))
        com_dist_err.append(np.nanstd(com_std_hold))

        com_avg = np.vstack((com_avg, np.nanmean(com_hold, axis=0)))
        com_err = np.vstack((com_err, np.nanstd(com_hold, axis=0)))

    return com_avg, com_err, np.array(com_dist_avg), np.array(com_dist_err)


def get_group_data(it_lst, sim, sim_dir, snap_lim, cluster_types, cluster_lim, n_clusters, cluster_variables):
    proc_file = sim_dir + sim + "/" + sim + "_processed.hdf5"
    proc_data = h5py.File(proc_file, "r")  # open processed data file

    data_file = sim_dir + sim + "/gc_groups.json"
    with open(data_file, "r") as file:
        gc_dict = json.load(file)
    group_lst = [int(grp) for grp in gc_dict[sim].keys()]

    pub_data = sim_dir + "snapshot_times_public.txt"
    pub_snaps = pd.read_table(pub_data, comment="#", header=None, sep=r"\s+")
    pub_snaps.columns = [
        "index",
        "scale_factor",
        "redshift",
        "time_Gyr",
        "lookback_time_Gyr",
        "time_width_Myr",
    ]

    snap_lst = pub_snaps[pub_snaps["index"] >= snap_lim]["index"].values
    time_lst = pub_snaps[pub_snaps["index"] >= snap_lim]["time_Gyr"].values

    # main loop
    grp_dict = {}
    grp_dict["other"] = {}
    for group_id in group_lst:
        print("Processing:", group_id)
        kl_avg, kl_std, cont_avg, cont_std = get_kl(group_id, snap_lst, sim, sim_dir)
        gc_num_avg, gc_num_std = get_gc_num(group_id, it_lst, snap_lst, proc_data)
        gc_dis, gc_dis_sph = get_dispersion(group_id, it_lst, snap_lst, proc_data)
        com_avg, com_err, com_dist_avg, com_dist_err = get_com(group_id, it_lst, snap_lst, proc_data)

        grp_dict["other"][str(group_id)] = {}
        grp_dict["other"][str(group_id)]["time"] = time_lst
        grp_dict["other"][str(group_id)]["snaps"] = snap_lst

        grp_dict["other"][str(group_id)]["kl_avg"] = kl_avg
        grp_dict["other"][str(group_id)]["kl_std"] = kl_std
        grp_dict["other"][str(group_id)]["cont_avg"] = cont_avg
        grp_dict["other"][str(group_id)]["cont_std"] = cont_std

        grp_dict["other"][str(group_id)]["gc_num_avg"] = gc_num_avg
        grp_dict["other"][str(group_id)]["gc_num_std"] = gc_num_std

        grp_dict["other"][str(group_id)]["gc_dis"] = gc_dis
        grp_dict["other"][str(group_id)]["gc_dis_sph"] = gc_dis_sph

        grp_dict["other"][str(group_id)]["com_avg"] = com_avg
        grp_dict["other"][str(group_id)]["com_std"] = com_err
        grp_dict["other"][str(group_id)]["com_dist_avg"] = com_dist_avg
        grp_dict["other"][str(group_id)]["com_dist_std"] = com_dist_err

        print("")

    for cluster_type in cluster_types:
        print("Clustering:", cluster_type)
        cluster_dict = group_iteration_clustering(
            it_lst, cluster_type, n_clusters, cluster_variables, sim, sim_dir, cluster_lim
        )

        grp_dict[cluster_type] = {}
        grp_dict[cluster_type]["time"] = cluster_dict["time"]
        grp_dict[cluster_type]["snaps"] = cluster_dict["snap"]
        grp_dict[cluster_type]["in_situ"] = cluster_dict["in_situ"]
        grp_dict[cluster_type]["ex_situ"] = cluster_dict["ex_situ"]

        grp_dict[cluster_type]["groups"] = {}
        for group_id in group_lst:
            grp_dict[cluster_type]["groups"][str(group_id)] = cluster_dict["groups"][str(group_id)]
        # grp_dict[cluster_type]["groups"] = cluster_dict["groups"]

        print("")

    return grp_dict
