import json
import multiprocessing as mp
import os
import sys

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from gc_utils import iteration_name, snapshot_name  # type: ignore
from scipy.interpolate import interp1d


def get_type_flag(gc, grp, gc_survive_lst):
    if gc in gc_survive_lst:
        if grp == 0:
            # formed in-situ and survived to z = 0
            type_flag = 0
        elif grp > 0:
            # formed ex-situ, accreted and survived to z = 0
            type_flag = 2

    elif grp < -2:
        # formed ex-situ and died before accretion
        type_flag = 4

    elif grp > 0:
        # formed ex-situ, accreted but died before z = 0
        type_flag = 3

    elif grp == 0:
        # formed in-situ, but died before z = 0
        type_flag = 1

    else:
        sys.exit("Some GC Missing Type Flag")

    return type_flag


def get_mass_data(it, evolve_mass_loss, sim, sim_dir, mass_dict: dict = {}):
    # get required files

    proc_file = sim_dir + sim + "/" + sim + "_processed.hdf5"
    proc_data = h5py.File(proc_file, "r")  # open processed data file

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

    it_id = iteration_name(it)

    print(it_id)

    it_mass_dict = {}

    src_dat = proc_data[it_id]["source"]
    ana_mask = np.array(src_dat["analyse_flag"]) == 1

    gc_lst = np.array(src_dat["gc_id"])[ana_mask]
    gc_survive_lst = np.array(proc_data[it_id]["snapshots"]["snap600"]["gc_id"])

    for gc in gc_lst:
        gc_id = str(gc)
        it_mass_dict[gc_id] = {}

        idx = np.where(np.array(src_dat["gc_id"])[ana_mask] == gc)[0][0]
        grp = np.array(src_dat["group_id"])[ana_mask][idx]

        it_mass_dict[gc_id]["group_id"] = int(grp)
        it_mass_dict[gc_id]["type_flag"] = get_type_flag(gc, grp, gc_survive_lst)

        time_lst = []
        time_for_lst = []

        log_mass_lst = []  # mass at each time
        mass_loss_lst = []  # mass loss
        mass_loss_det_lst = []  # mass loss but taking t_form as 55% of t_form

        # get formation information
        t_form = np.array(src_dat["form_time"])[ana_mask][idx]

        log_mass_form = np.array(src_dat["logm_tform"])[ana_mask][idx]
        log_mass_form_det = (
            np.log10(1 - evolve_mass_loss) + log_mass_form
        )  # corrected for mass loss by evolution

        # get other details
        t_dis = np.array(src_dat["t_dis"])[ana_mask][idx]

        # update relevant lists
        time_lst.append(t_form)
        time_for_lst.append(t_form - t_form)

        # only for accreted (types 2, 3, 4) GCs
        if it_mass_dict[gc_id]["type_flag"] > 1:
            time_acc_lst = []

            t_acc = np.array(src_dat["t_acc"])[ana_mask][idx]
            time_acc_lst.append(t_form - t_acc)

        log_mass_lst.append(log_mass_form)

        # add mass loss as 0 for first time step
        mass_loss_lst.append(0)
        mass_loss_det_lst.append(0)

        if t_dis == -1:
            snap_lst = pub_snaps[(pub_snaps["time_Gyr"] >= t_form)]["index"]
        else:
            snap_lst = pub_snaps[(pub_snaps["time_Gyr"] >= t_form) & (pub_snaps["time_Gyr"] <= t_dis)][
                "index"
            ]

        # this has been added to ensure that GCs that form and die between snaps are considered in the for loop
        # that gets the details of their death (death loop)
        time = t_form

        for snap in snap_lst:
            snap_id = snapshot_name(snap)
            snap_dat = proc_data[it_id]["snapshots"][snap_id]

            # check, sometimes timing issues with GC formation model
            gc_snap_lst = np.array(snap_dat["gc_id"])
            if gc not in gc_snap_lst:
                continue

            time = pub_snaps[pub_snaps["index"] == snap]["time_Gyr"].values[0]

            snap_idx = np.where(np.array(snap_dat["gc_id"]) == gc)[0][0]
            cur_log_mass = np.array(snap_dat["mass"])[snap_idx]

            pst_log_mass = log_mass_lst[-1]
            pst_mass = 10**pst_log_mass

            cur_mass = 10**cur_log_mass

            # get mass loss (not logged)
            mass_los = pst_mass - cur_mass

            # make initial correction to detectable mass loss list
            if len(log_mass_lst) == 1:
                mass_form_det = 10**log_mass_form_det
                mass_los_det = mass_form_det - cur_mass
            else:
                mass_los_det = mass_los

            # update relevant lists
            log_mass_lst.append(cur_log_mass)

            mass_loss_lst.append(mass_los)
            mass_loss_det_lst.append(mass_los_det)

            time_lst.append(time)
            time_for_lst.append(time - t_form)

            # get accretion time
            if it_mass_dict[gc_id]["type_flag"] > 1:
                time_acc_lst.append(time - t_acc)

        # need to consider full disruption stuff too
        # like if in past list but not in current list mass loss is 100%

        # (death loop)
        # update full disruption details fro GCs that don't survive to z = 0

        # if (mass_dict[it_id][gc]["type_flag"] != 0) and (mass_dict[it_id][gc]["type_flag"] != 2):
        #     last_time = pub_snaps[pub_snaps["time_Gyr"] > time]["time_Gyr"].values[0]

        if t_dis != -1:
            # no longer exists
            # final_log_mass = np.nan
            final_log_mass = np.nan

            # past mass
            pst_log_mass = log_mass_lst[-1]
            pst_mass = 10**pst_log_mass

            # as fully disrupted
            mass_los = pst_mass

            time_lst.append(t_dis)
            time_for_lst.append(t_dis - t_form)

            # get accretion time
            if it_mass_dict[gc_id]["type_flag"] > 1:
                time_acc_lst.append(t_dis - t_acc)

            log_mass_lst.append(final_log_mass)
            mass_loss_lst.append(mass_los)
            mass_loss_det_lst.append(mass_los)

        # update dictionary
        it_mass_dict[gc_id]["time"] = time_lst
        it_mass_dict[gc_id]["form_time"] = time_for_lst

        if it_mass_dict[gc_id]["type_flag"] > 1:
            it_mass_dict[gc_id]["acc_time"] = time_acc_lst

        it_mass_dict[gc_id]["log_mass"] = log_mass_lst
        it_mass_dict[gc_id]["mass_loss"] = mass_loss_lst
        it_mass_dict[gc_id]["mass_loss_detectable"] = mass_loss_det_lst

    # close file
    proc_data.close()

    mass_dict[it_id] = it_mass_dict

    return mass_dict


def create_hdf5(sim: str, data_dict: dict, save_dir: str):
    save_file = save_dir + sim + "_mass_data.hdf5"  # save location

    if not os.path.exists(save_file):
        h5py.File(save_file, "w")

    with h5py.File(save_file, "a") as hdf:
        for it_id in data_dict.keys():
            # print(it_id)
            if it_id in hdf.keys():
                it_group = hdf[it_id]
            else:
                it_group = hdf.create_group(it_id)

            for gc_id in data_dict[it_id].keys():
                # print(gc_id)
                if gc_id in it_group.keys():
                    gc_id_group = it_group[gc_id]
                else:
                    gc_id_group = it_group.create_group(gc_id)

                for key in data_dict[it_id][gc_id].keys():
                    # print(key)
                    if key in gc_id_group.keys():
                        del gc_id_group[key]
                    gc_id_group.create_dataset(key, data=np.array(data_dict[it_id][gc_id][key]))


if __name__ == "__main__":
    # Save to a JSON file
    evolve_mass_loss = 0.45

    it_min = 0
    it_max = 3
    # it_max = 100

    sim = "m12i"

    sim_dir = "/Users/z5114326/Documents/simulations/"
    save_dir = "data/gc_mass_data/"

    it_lst = np.linspace(it_min, it_max, it_max + 1, dtype=int)

    cores = 4

    with mp.Manager() as manager:
        shared_dict = manager.dict()  # Shared dictionary across processes
        args = [(it, evolve_mass_loss, sim, sim_dir, shared_dict) for it in it_lst]

        with mp.Pool(processes=cores, maxtasksperchild=1) as pool:
            pool.starmap(get_mass_data, args, chunksize=1)

        final_dict = dict(shared_dict)

    # for key in final_dict["it000"]["119998556"].keys():
    #     print(key)
    #     print(final_dict["it000"]["119998556"][key])
    create_hdf5(sim, final_dict, save_dir)
