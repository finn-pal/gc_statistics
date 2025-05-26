import json

import gc_utils
import gizmo_analysis as gizmo
import halo_analysis as halo
import numpy as np
import pandas as pd
import utilities as ut


def get_kappa_co(kappa_dat, snap_data):
    kappa_dict = {"snap": [], "time": [], "kappa_co": []}

    for snap_id in kappa_dat.keys():
        snap = int(snap_id[4:])
        time = snap_data[snap_data["index"] == snap]["time_Gyr"].values[0]

        kappa_dict["snap"].append(snap)
        kappa_dict["time"].append(time)

        kappa_dict["kappa_co"].append(kappa_dat[snap_id]["kappa_co_sg"])

    return kappa_dict


def get_halo_details(halt, main_halo, snap_data, snap_lim, fire_dir):
    snap_lst = snap_data[snap_data["index"] >= snap_lim]["index"].values

    gc_utils.block_print()
    hals = halo.io.IO.read_catalogs("index", snap_lst, fire_dir)
    gc_utils.enable_print()

    gal_dict = {"snap": [], "time": [], "num_mergers": [], "prog_idx": [], "prog_mass": [], "acc_rate": []}

    for snap in range(snap_lim, 600 + 1):
        halo_tid = gc_utils.get_halo_prog_at_snap(halt, main_halo, snap)
        idx_init = np.where(halt["tid"] == halo_tid)[0][0]
        idx = idx_init
        idx_lst = []
        mass_lst = []
        for i in range(0, halt["progenitor.number"][idx_init]):
            if i == 0:
                idx = halt["progenitor.main.index"][idx]
            else:
                idx = halt["progenitor.co.index"][idx]

            mass = halt["mass"][idx]

            if halt["am.phantom"][idx] == 0:
                idx_lst.append(idx)
                mass_lst.append(mass)

        gal_dict["snap"].append(snap)
        gal_dict["num_mergers"].append(len(idx_lst) - 1)

        gal_dict["prog_idx"].append(idx_lst)
        gal_dict["prog_mass"].append(mass_lst)

        # accretion rate stuff

        halo_tidx = idx_init
        halo_cidx = halt["catalog.index"][halo_tidx]

        hal_snap = hals[snap]
        acc_rate = hal_snap["accrete.rate.tdyn"][halo_cidx]

        gal_dict["acc_rate"].append(acc_rate)

        # time

        time = snap_data[snap_data["index"] == snap]["time_Gyr"].values[0]
        gal_dict["time"].append(time)

    return gal_dict


def get_gas_frac(halt, main_halo, snap_pub, snap_lim, fire_dir, log_temp_lim=4.5):
    snap_lst = snap_pub[snap_pub["index"] >= snap_lim]["index"].values
    time_lst = snap_pub[snap_pub["index"] >= snap_lim]["time_Gyr"].values

    gas_frac_dict = {"snap": [], "time": [], "cold_gas_frac": []}
    for snap, time in zip(snap_lst, time_lst):
        print("Gas Frac For ", snap)
        print("")
        halo_tid = gc_utils.get_halo_prog_at_snap(halt, main_halo, snap)
        halo_idx = np.where(halt["tid"] == halo_tid)[0][0]
        rad = halt["star.radius.90"][halo_idx]

        gc_utils.block_print()
        part = gizmo.io.Read.read_snapshots(["gas", "star"], "index", snap, fire_dir)
        gc_utils.enable_print()

        gas_dist = ut.particle.get_distances_wrt_center(
            part,
            species=["gas"],
            center_position=halt["position"][halo_idx],
            total_distance=True,
        )

        star_dist = ut.particle.get_distances_wrt_center(
            part,
            species=["star"],
            center_position=halt["position"][halo_idx],
            total_distance=True,
        )

        gas_dist_msk = gas_dist < rad
        star_dist_msk = star_dist < rad

        gas_temp_msk = np.log10(part["gas"]["temperature"]) < log_temp_lim

        cold_gas_mass = np.sum(part["gas"]["mass"][gas_dist_msk & gas_temp_msk])
        baryon_mass = np.sum(part["gas"]["mass"][gas_dist_msk]) + np.sum(part["star"]["mass"][star_dist_msk])

        cold_gas_frac = cold_gas_mass / baryon_mass

        gas_frac_dict["cold_gas_frac"].append(cold_gas_frac)
        gas_frac_dict["snap"].append(snap)
        gas_frac_dict["time"].append(time)

    return gas_frac_dict


def get_galaxy_data(sim, sim_dir, snap_lim):
    fire_dir = sim_dir + sim + "/" + sim + "_res7100"

    snap_file = fire_dir + "/snapshot_times.txt"
    snap_data = pd.read_table(snap_file, comment="#", header=None, sep=r"\s+")
    snap_data.columns = [
        "index",
        "scale_factor",
        "redshift",
        "time_Gyr",
        "lookback_time_Gyr",
        "time_width_Myr",
    ]

    snap_pub_file = sim_dir + "/snapshot_times_public.txt"
    snap_pub = pd.read_table(snap_pub_file, comment="#", header=None, sep=r"\s+")
    snap_pub.columns = [
        "index",
        "scale_factor",
        "redshift",
        "time_Gyr",
        "lookback_time_Gyr",
        "time_width_Myr",
    ]

    kap_dat_file = sim_dir + sim + "/galaxy_kinematics.json"
    with open(kap_dat_file, "r") as file:
        kappa_dat = json.load(file)

    sim_code_file = sim_dir + "simulation_codes.json"
    with open(sim_code_file, "r") as file:
        sim_codes = json.load(file)
    main_halo = sim_codes[sim]["halo"]

    gc_utils.block_print()
    halt = halo.io.IO.read_tree(simulation_directory=fire_dir, species="star")
    gc_utils.enable_print()

    kappa_dict = get_kappa_co(kappa_dat, snap_data)
    halo_dict = get_halo_details(halt, main_halo, snap_data, snap_lim, fire_dir)
    gas_frac_dict = get_gas_frac(halt, main_halo, snap_pub, snap_lim, fire_dir, log_temp_lim=4.5)

    galaxy_info = {"kappa_co": kappa_dict, "halo_data": halo_dict, "gas_frac": gas_frac_dict}

    return galaxy_info
