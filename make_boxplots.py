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
replicate Killian's heatmap of freeze order

modify functions to return plot objects instead of just 
printing the plots

generate legends with sample names
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
    # frzTall = np.sort(frzTall, axis=0)

    return frzTall


# %%
def make_hist(*frzdata):
    """
    Returns a histogram of freezing events at temp T

    WARNING: NOT NORMALIZED TO N(T)!!! Strongly dependent on
    droplet size
    """
    for data in frzdata:
        plt.hist(data, bins=100)

    plt.xlabel("Freezing Temperature ($^oC$)")
    plt.ylabel("$N_{wells}$")
    plt.title("Frozen Well Count as a Function of Temperature")


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

    plt.boxplot(all_dat)
    plt.xlabel("Sample")
    plt.ylabel("Freezing Temp ($^oC$)")
    plt.title("Freezing Temp Distrubutions")


# %%
def make_ff_curve(data):
    """
    creates a frozen fraction (FF) curve as a function of
    temperature
    """
    data = np.sort(data, axis=0)
    ind = np.linspace(len(data), 1, len(data)) / len(data)

    plt.plot(data, ind)
    plt.xlabel("Temp ($^oC$")
    plt.ylabel("FF")
    plt.title("Frozen Fraction")


# %%
def make_big_K(*data, norm=1, Vwell=10):
    """
    creates a cumulative concentration curve from Vali eq.

    default normalization constant (eg surface area) = 1
    default well volume = 10 microlitre
    """
    # make frozen fraction curve and calculate K from it
    length = len(data[1])
    ind = np.linspace(length, 1, length) / length

    K = -np.log(1 - ind) / (norm * Vwell * 10 ** -6)  # Vali eq.

    for dat in data:
        dat = np.sort(dat, axis=0)
        plt.plot(dat, K)

    plt.yscale("log")
    plt.xlabel("Temp ($^oC$)")
    plt.ylabel("K(T) ($L^{-1}$)")
    plt.title(f"Cumulative Freezing Spectra ({Vwell}$\mu L$ droplet vol)")


# %%
def main():
    # name of expt runs
    SA1 = "finc_20201007_SA1"
    SA2 = "finc_20201007_SA2"  # this one is bunk!
    SA3 = "finc_20201007_SA3"
    lig = "finc_20201028_lignin2"

    # get the data from matlab
    SA1 = get_mat(SA1)
    SA2 = get_mat(SA2)
    SA3 = get_mat(SA3)
    lig = get_mat(lig)

    the_data = SA1, SA3, lig

    # make some plots
    # make_hist(*the_data)
    # make_big_K(*the_data)
    make_boxplot(*the_data)


if __name__ == "__main__":
    main()

# %%
