import argparse
import json
import os
from collections import defaultdict

import agama
import cmasher as cmr
import gc_utils
import h5py
import halo_analysis as halo
import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
import seaborn as sns
import utilities as ut
from astropy.stats import biweight_location
from matplotlib import cm
from matplotlib.animation import PillowWriter
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LogNorm, Normalize
from matplotlib.ticker import LogLocator
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d, zoom
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
from scipy.stats import binned_statistic_2d
from sklearn import cluster as cl
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.preprocessing import StandardScaler


def make_correlation_gif(sim_lst, ):

###########################################################################################

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-s", "--simulations", required=True, type=str, nargs="+", help="simulation name (e.g. m12i)"
    )
    parser.add_argument("-l", "--location", required=True, type=str, help="sim_dir location")
    parser.add_argument("-a", "--iteration_low_limit", required=True, type=int, help="lower bound iteration")
    parser.add_argument("-b", "--iteration_up_limit", required=True, type=int, help="upper bound iteration")
    parser.add_argument("-p", "--snap_lim", required=False, type=int, help="lower snap limit", default=60)
    args = parser.parse_args()

    location = args.location
    if location == "local":
        sim_dir = "../../simulations/"

    elif location == "katana":
        sim_dir = "/srv/scratch/astro/z5114326/simulations/"

    elif location == "one_touch":
        sim_dir = "/Volumes/One_Touch/simulations/"

    elif location == "expansion":
        sim_dir = "/Volumes/Expansion/simulations/"

    elif location == "my_passport":
        sim_dir = "/Volumes/My_Passport_for_Mac/simulations/"

    else:
        raise RuntimeError("Incorrect location provided. Must be local, katana or one_touch.")

    sim_lst = args.simulations
    snap_lim = args.snap_lim

    it_min = args.iteration_low_limit
    it_max = args.iteration_up_limit
    it_rng = it_max - it_min
    it_lst = np.linspace(it_min, it_max, it_rng + 1, dtype=int)

    cluster_type = "birch"

    # variables = ["et", "jr", "jp", "jz"]
    variables = ["et_norm", "jr", "jp", "jz"]
    n_clusters = 2

    agama.setUnits(mass=1, length=1, velocity=1)

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
    snap_lst = pub_snaps["index"].to_numpy()
    time_lst = pub_snaps["time_Gyr"].to_numpy()

    snap_lst = pub_snaps[pub_snaps["index"] >= snap_lim]["index"].to_numpy()
    time_lst = pub_snaps[pub_snaps["index"] >= snap_lim]["time_Gyr"].to_numpy()
