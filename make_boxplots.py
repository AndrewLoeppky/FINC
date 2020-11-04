"""
===============================================================
DATA PROCESSING LIBRARY - FINC OUTPUTS

Author: Andrew Loeppky
UBC/ETH Atmospheric Science, NBD Group
Last Update Nov 4/2020

***Currently under development, dont trust this code yet***

functions for processing raw output from Freezing Ice Nuclei 
Counter (FINC) and generating plots.

TODO:
replicate Killian's heatmap of freeze order, frozen fraction 
curve, create a processing function to produce diff and cumulative
freezing temperatures
===============================================================
"""
# %%
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio


# %%
def get_mat(the_file):
    """
    retrieves freezing temperature vector from matlab workspace
    and returns it as a 288x1 numpy array

    'the_file' (str) is the name of the experiment folder
    """
    the_path = "C:/Users/Owner/UBC_F2020/FINC/outputs/"
    load_this = the_path + the_file + "/workspace.mat"

    mat = sio.loadmat(load_this)
    frzTall = mat["FrzTall"]
    frzTall = np.sort(frzTall, axis=0)
    return frzTall


# %%
def make_hist(frzdata):
    """
    Returns a histogram of freezing events at temp T

    WARNING: NOT NORMALIZED TO N(T)!!! Strongly dependent on
    droplet size
    """
    fig, ax = plt.subplots()
    ax = plt.hist(frzdata, bins=100)
    plt.xlabel("Freezing Temperature ($^oC$)")
    plt.ylabel("Frozen Well Count")

    # return fig


# %%
def make_boxplot(*data):
    """
    given one or more input freezing vectors, make boxplots of
    freezing temperature
    """
    # splat datasets and create list of 1d lists to for plotting
    all_dat = []
    for dat in data:
        dat = dat.reshape(dat.shape[0])
        all_dat.append(dat)

    fig, ax = plt.subplots()
    ax = plt.boxplot(all_dat)
    plt.xlabel("Sample")
    plt.ylabel("Freezing Temp ($^oC$)")


# %%
def main():
    # name of expt runs
    SA1 = "finc_20201007_SA1"
    SA2 = "finc_20201007_SA2"
    SA3 = "finc_20201007_SA3"

    # get the data from matlab
    SA1 = get_mat(SA1)
    SA2 = get_mat(SA2)
    SA3 = get_mat(SA3)

    # make some plots
    make_hist(SA1)
    make_boxplot(SA1, SA2, SA3)


if __name__ == "__main__":
    main()
