from collections import defaultdict

import gc_utils
import h5py
import numpy as np
import pandas as pd
from sklearn import cluster as cl
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.preprocessing import StandardScaler


def clustering(X, rnd_grps, n_clusters, cluster_type):
    scaler = StandardScaler()
    X_scal = scaler.fit_transform(X)

    if cluster_type == "kmeans":
        kmeans = cl.KMeans(n_clusters=n_clusters, random_state=0).fit(X_scal)
        labels = kmeans.labels_

    if cluster_type == "agglom":
        agglom = cl.AgglomerativeClustering(n_clusters=n_clusters).fit(X_scal)
        labels = agglom.labels_

    if cluster_type == "birch":
        birch = cl.Birch(n_clusters=n_clusters, threshold=0.2).fit(X_scal)
        labels = birch.labels_

    count = 0
    for label in range(n_clusters):
        label_msk = labels == label

        # Get group ID counts for accreted particles in this cluster
        grp_lst_acc = np.abs(rnd_grps[label_msk])
        unique_acc, counts_acc = np.unique(grp_lst_acc, return_counts=True)
        count_dict_acc = dict(zip(unique_acc, counts_acc))

        # Identify the most common group
        grp_identified = unique_acc[np.argmax(counts_acc)]

        if count_dict_acc[grp_identified] > count:
            count = count_dict_acc[grp_identified]
            grp_label = label
            pred_grp = grp_identified

    return pred_grp, grp_label, labels


def clustering_envelope(it, snap, proc_data, n_clusters, variables, cluster_type):
    it_id = gc_utils.iteration_name(it)
    snap_id = gc_utils.snapshot_name(snap)

    snp_dat = proc_data[it_id]["snapshots"][snap_id]
    acc_msk = snp_dat["now_accreted"][()] == 1

    full_grps = np.abs(snp_dat["group_id"][acc_msk])
    full_gcid = snp_dat["gc_id"][acc_msk]

    rnd_grps = full_grps.copy()
    rnd_gcid = full_gcid.copy()

    unq_grp, cnt_grp = np.unique(full_grps, return_counts=True)
    grp_cnt_dict = dict(zip(unq_grp, cnt_grp))

    res_dict = {
        gc_id: {"True": grp, "Pred": None} if grp_cnt_dict[grp] >= 5 else {"True": -1, "Pred": None}
        for gc_id, grp in zip(full_gcid, full_grps)
    }

    X_init = np.column_stack(
        [
            snp_dat["j.cyl"][:, 2][acc_msk]
            if var == "jz"
            else snp_dat["j.cyl"][:, 1][acc_msk]
            if var == "jp"
            else snp_dat["j.cyl"][:, 0][acc_msk]
            if var == "jr"
            else snp_dat[var][acc_msk]
            for var in variables
        ]
    )

    X = X_init

    while len(X) >= 5:
        pred_grp, grp_label, labels = clustering(X, rnd_grps, n_clusters, cluster_type)

        lab_msk = labels == grp_label
        gc_filt = rnd_gcid[lab_msk]

        for gc_id in gc_filt:
            res_dict[gc_id]["Pred"] = pred_grp

        rnd_grps = rnd_grps[~lab_msk]
        rnd_gcid = rnd_gcid[~lab_msk]
        X = X[~lab_msk]

    # group remaining gc's as other (-1)
    for gc_id in rnd_gcid:
        res_dict[gc_id]["Pred"] = -1

    return res_dict


def calc_entropy(p):
    """Compute Shannon entropy of a probability distribution p."""
    p = p[p > 0]  # filter out zeros to avoid log(0)
    return -np.sum(p * np.log2(p))


def clustering_entropy(true_labels, predicted_labels):
    labels = np.unique(np.concatenate((np.array(true_labels), np.array(predicted_labels))))
    cm = confusion_matrix(true_labels, predicted_labels, labels=labels)
    total = cm.sum()

    weighted_entropies = [] = []
    entropy_by_cluster = {}

    for i, col in enumerate(cm.T):  # loop over predicted clusters
        cluster_label = labels[i]
        cluster_total = np.sum(col)

        if cluster_total == 0:
            entropy_by_cluster[cluster_label] = 0.0
            weighted_entropies.append(0.0)
            continue

        p = col / cluster_total
        e = calc_entropy(p)
        weight = cluster_total / total

        entropy_by_cluster[cluster_label] = e / np.log2(cm.shape[1])
        weighted_entropies.append(weight * e)

    total_entropy = np.sum(weighted_entropies)

    # print(entropy_by_cluster)
    return total_entropy


def snapshot_clustering(it, proc_data, cluster_type, n_clusters, variables, snap_lst):
    snap_dict = {}
    snap_dict["entropy"] = []
    snap_dict["accuracy"] = []
    snap_dict["weight_recall"] = []
    snap_dict["macro_recall"] = []
    snap_dict["weight_precision"] = []
    snap_dict["macro_precision"] = []
    snap_dict["weight_f1"] = []
    snap_dict["macro_f1"] = []

    for snap in snap_lst:
        res_dict = clustering_envelope(it, snap, proc_data, n_clusters, variables, cluster_type)

        true_grp = [res_dict[gc_id]["True"] for gc_id in res_dict.keys()]
        pred_grp = [res_dict[gc_id]["Pred"] for gc_id in res_dict.keys()]

        entropy = clustering_entropy(true_grp, pred_grp)
        num_grps = len(np.unique(np.concatenate((true_grp, pred_grp))))
        entropy_norm = entropy / np.log2(num_grps)
        snap_dict["entropy"].append(entropy_norm)

        report = classification_report(true_grp, pred_grp, output_dict=True, zero_division=0)
        snap_dict["accuracy"].append(report["accuracy"])
        snap_dict["weight_recall"].append(report["weighted avg"]["recall"])
        snap_dict["macro_recall"].append(report["macro avg"]["recall"])
        snap_dict["weight_precision"].append(report["weighted avg"]["precision"])
        snap_dict["macro_precision"].append(report["macro avg"]["precision"])
        snap_dict["weight_f1"].append(report["weighted avg"]["f1-score"])
        snap_dict["macro_f1"].append(report["macro avg"]["f1-score"])

    return snap_dict


def iteration_clustering(it_lst, cluster_type, n_clusters, variables, sim, sim_dir, snap_lim):
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

    snap_lst = pub_snaps["index"].values
    snap_lst = snap_lst[snap_lst >= snap_lim]
    time_lst = [pub_snaps[pub_snaps["index"] == snap]["time_Gyr"].values[0] for snap in snap_lst]

    # ent_it_lst = []
    it_dict = defaultdict(list)

    for it in it_lst:
        snap_dict = snapshot_clustering(it, proc_data, cluster_type, n_clusters, variables, snap_lst)
        for metric in snap_dict.keys():
            it_dict[metric].append(snap_dict[metric])

    cluster_dict = {
        "snap": snap_lst,
        "time": time_lst,
    }

    # Add mean and std of each metric across iterations
    for metric in it_dict.keys():
        cluster_dict[f"{metric}_avg"] = np.nanmean(it_dict[metric], axis=0)
        cluster_dict[f"{metric}_std"] = np.nanstd(it_dict[metric], axis=0)

    return cluster_dict


def averaged_metric(report, group_keys, metric="f1-score", exclude_group="0"):
    total_w_sum = 0.0
    total_m_sum = 0.0
    total_counts = 0

    num_grps = len(group_keys) - 1

    for group in group_keys:
        if group == exclude_group:
            continue

        value = report[group][metric]
        count = report[group]["support"]

        total_w_sum += value * count
        total_m_sum += value
        total_counts += count

    w_avg = total_w_sum / total_counts
    m_avg = total_m_sum / num_grps

    return w_avg, m_avg


def group_snap_clustering(it, proc_data, cluster_type, n_clusters, variables, snap_lst):
    group_dict = {}
    for snap in snap_lst:
        res_dict = clustering_envelope(it, snap, proc_data, n_clusters, variables, cluster_type)

        true_grp = [res_dict[gc_id]["True"] for gc_id in res_dict.keys()]
        pred_grp = [res_dict[gc_id]["Pred"] for gc_id in res_dict.keys()]

        report = classification_report(true_grp, pred_grp, output_dict=True, zero_division=0)

        group_keys = [k for k in report.keys() if not any(c.isalpha() for c in str(k))]
        for group in group_keys:
            if group not in group_dict.keys():
                group_dict[group] = {"snap": [], "recall": [], "precision": [], "f1": []}

            group_dict[group]["snap"].append(snap)
            group_dict[group]["recall"].append(report[group]["recall"])
            group_dict[group]["precision"].append(report[group]["precision"])
            group_dict[group]["f1"].append(report[group]["f1-score"])

        if "ex_situ" not in group_dict.keys():
            group_dict["ex_situ"] = {
                "snap": [],
                "weight_recall": [],
                "macro_recall": [],
                "weight_precision": [],
                "macro_precision": [],
                "weight_f1": [],
                "macro_f1": [],
            }

        w_r, m_r = averaged_metric(report, group_keys, metric="recall")
        w_p, m_p = averaged_metric(report, group_keys, metric="precision")
        w_f, m_f = averaged_metric(report, group_keys, metric="f1-score")

        group_dict["ex_situ"]["snap"].append(snap)
        group_dict["ex_situ"]["weight_recall"].append(w_r)
        group_dict["ex_situ"]["macro_recall"].append(m_r)
        group_dict["ex_situ"]["weight_precision"].append(w_p)
        group_dict["ex_situ"]["macro_precision"].append(m_p)
        group_dict["ex_situ"]["weight_f1"].append(w_f)
        group_dict["ex_situ"]["macro_f1"].append(m_f)

    return group_dict


def group_iteration_clustering(it_lst, cluster_type, n_clusters, variables, sim, sim_dir, snap_lim):
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

    snap_lst = pub_snaps["index"].values
    snap_lst = snap_lst[snap_lst >= snap_lim]
    time_lst = [pub_snaps[pub_snaps["index"] == snap]["time_Gyr"].values[0] for snap in snap_lst]

    it_dict = {}

    for it in it_lst:
        grp_snap_dict = group_snap_clustering(it, proc_data, cluster_type, n_clusters, variables, snap_lst)
        for key in grp_snap_dict.keys():
            if key not in it_dict:
                it_dict[key] = defaultdict(list)

            for metric in grp_snap_dict[key]:
                metric_lst = []
                for snap in snap_lst:
                    grp_snp_lst = np.array(grp_snap_dict[key]["snap"])
                    if snap not in grp_snp_lst:
                        metric_lst.append(np.nan)
                    else:
                        snp_msk = grp_snp_lst == snap
                        val = np.array(grp_snap_dict[key][metric])[snp_msk][0]
                        metric_lst.append(val)

                it_dict[key][metric].append(metric_lst)

    cluster_dict = {"snap": snap_lst, "time": time_lst, "groups": {}}

    # Add mean and std of each metric across iterations
    for key in it_dict.keys():
        if key == "ex_situ":
            cluster_dict[key] = {}
        else:
            cluster_dict["groups"][key] = {}
        for metric in it_dict[key].keys():
            if metric != "snap":
                if key == "ex_situ":
                    cluster_dict[key][f"{metric}_avg"] = np.nanmean(it_dict[key][metric], axis=0)
                    cluster_dict[key][f"{metric}_std"] = np.nanstd(it_dict[key][metric], axis=0)
                else:
                    cluster_dict["groups"][key][f"{metric}_avg"] = np.nanmean(it_dict[key][metric], axis=0)
                    cluster_dict["groups"][key][f"{metric}_std"] = np.nanstd(it_dict[key][metric], axis=0)

    cluster_dict["in_situ"] = {}
    for metric in cluster_dict["groups"]["0"].keys():
        cluster_dict["in_situ"][metric] = cluster_dict["groups"]["0"][metric]

    return cluster_dict
